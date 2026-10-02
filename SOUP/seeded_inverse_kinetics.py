# Bootstrap Uncertainty Estimation - Zymonic Acid Kinetics

# This script estimates confidence intervals for the 10 rate constants of the zymonic acid interconversion network, using bootstrap resampling with convergence monitoring.
# It is meant to be run after Script 1 (inverse_kinetics_optimisation.py) has produced a best-fit set of rate constants via differential evolution.

# WORKFLOW
# Step 1: Run Script 1 to find the best-fit rate constants via differential evolution.
# Step 2: Copy the best-fit parameter values from Script 1 into the best_fit_params array inside run_fast_bootstrap_analysis(), in this file.
# Step 3: Run this script. It will:
#   a) Generate n_bootstrap resampled datasets, by sampling the experimental data with replacement (bootstrap resampling).
#   b) Fit rate constants to each resampled dataset using L-BFGS-B local optimisation, starting from the best-fit parameters. If local optimisation fails, differential evolution is used as a fallback.
#   c) Monitor convergence: every 100 samples, check whether the 95% confidence interval bounds for each parameter have changed by less than 5%. This is the convergence criterion from Perkins et al.
#   d) Report 95% confidence intervals, coefficient of variation, and convergence status for all 10 parameters.
# Step 4: Results and plots are saved to ./Outputs/.


# REPRODUCIBILITY
# A random seed is generated from the current date and time at the start of each run. This seed is used for all bootstrap resampling, so results will differ between runs, but can be reproduced exactly using the saved seed.
# The seed is printed to the terminal and saved to the output JSON file. To reproduce a specific run, set RANDOM_SEED manually at the top of the seed control section below.

# OUTPUTS
# - ./Outputs/bootstrap_distributions.png      : Parameter distribution histograms
# - ./Outputs/bootstrap_results_detailed.json  : Full statistics including seed
# - ./Outputs/publication_ready_results.csv    : Table of CI results


# This script imports from the existing simulation framework: rate_laws.py, differentials.py, data.py, utilities.py.
# The input files (zymonic-atmosphere.dat, zymonic-aqueous.dat) are created automatically by create_zymonic_input_files() before the bootstrap runs.

import os
import sys
import numpy as np
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution
from scipy.integrate import solve_ivp
import time
from datetime import datetime
import json
import warnings

# scipy optimisation warnings from intermediate evaluations are suppressed here.
warnings.filterwarnings('ignore', category=RuntimeWarning)


# SECTION 1: SEED CONTROL FOR REPRODUCIBILITY
# A unique seed is generated from the current timestamp at the start of each run. This ensures:
#   - Results differ between runs (different bootstrap resamples each time).
#   - Any run can be reproduced exactly using the saved seed.
#   - The seed is traceable to a specific date and time.
# The seed is saved in the output JSON file. To reproduce a previous run, replace the generate_seed_from_datetime() call below with:
#     RANDOM_SEED = <saved_seed>
#     SEED_INFO = {'seed': RANDOM_SEED, 'timestamp': '...', 'readable': '...'}

def generate_seed_from_datetime():
    # Purpose: produces a fresh, traceable random seed from the current date and time, so each run's bootstrap resampling is reproducible from its own saved seed.
    # Inputs are: none.
    # Outputs are: seed, an integer seed derived from the current timestamp (converted to microseconds since the Unix epoch, and kept within numpy's valid range of < 2^31); seed_info, a dictionary holding that seed value, an ISO timestamp, and a human-readable datetime string, for logging and reproducibility. A different seed is returned on every call.
    # Called by: directly from module-level code (not inside a function), in this file, to set RANDOM_SEED and SEED_INFO when the module is first imported.
    now = datetime.now()
    seed = int(now.timestamp() * 1000000) % (2**31)

    seed_info = {
        'seed': seed,
        'timestamp': now.isoformat(),
        'readable': now.strftime('%Y-%m-%d %H:%M:%S.%f')
    }

    return seed, seed_info


# The seed is generated at import time, so it stays consistent across the entire run.
RANDOM_SEED, SEED_INFO = generate_seed_from_datetime()


def set_all_seeds(seed=RANDOM_SEED):
    # Purpose: seeds every source of randomness this script uses, so the sequence of bootstrap resamples is deterministic for a given seed.
    # Inputs are: seed, the random seed to use (defaults to RANDOM_SEED, the datetime-derived seed generated above).
    # Outputs are: nothing is returned; every source of randomness this script uses (numpy.random, used for bootstrap resampling via np.random.choice; and Python's built-in random module, set for completeness) is seeded, so all bootstrap resampling after this call produces a deterministic sequence of samples for the given seed.
    # Called by: directly from module-level code (not inside a function), in this file, once at startup.
    np.random.seed(seed)
    import random
    random.seed(seed)

    print(f"Random seed set to {seed}")
    print(f"   Generated from timestamp: {SEED_INFO['readable']}")
    print(f"   All bootstrap resampling will use this seed")
    print(f"   To reproduce these results, use seed={seed}\n")


# The seeds are set immediately at startup.
set_all_seeds(RANDOM_SEED)


# SECTION 2: EXPERIMENTAL DATA
# This is the concentration-time data for each observed species, converted from a percentage of total zymonic acid to a molar concentration (M).
# ZK represents the combined ketone measurement (ZCK and ZOK cannot be resolved experimentally).
# Data source is from Perkins et al.

EXPERIMENTAL_DATA = {
    'ZCD': {
        'time': np.array([0.0, 2.28, 6.09, 7.89, 13.79, 17.76, 21.76, 27.87, 31.94, 38.06,
                         42.15, 56.51, 72.94, 87.32, 118.12, 150.99, 181.80, 245.47, 307.08,
                         430.27, 551.39, 674.53, 857.20, 1039.83, 1138.32]),
        'concentration': np.array([0.0, 27.17, 36.95, 45.65, 54.13, 58.91, 62.17, 63.91, 65.0,
                                 65.86, 66.30, 66.52, 65.86, 65.0, 63.91, 62.17, 60.86, 58.04,
                                 55.65, 51.73, 48.47, 46.52, 42.82, 40.65, 39.78]) / 1000.0
    },
    'ZCE': {
        'time': np.array([0.0, 0.98, 3.32, 7.70, 14.11, 18.36, 24.62, 28.80, 35.00, 39.13,
                         43.27, 57.67, 74.10, 88.48, 152.11, 184.95, 246.52, 308.10, 431.22,
                         552.29, 858.00, 1040.59, 119.27, 675.39, 1139.07]),
        'concentration': np.array([100.0, 71.30, 61.52, 52.17, 43.26, 38.26, 34.56, 31.95, 30.43,
                                 29.34, 28.26, 26.73, 26.08, 25.43, 24.13, 23.47, 22.17, 21.08,
                                 19.34, 18.04, 15.65, 14.56, 24.78, 17.17, 14.13]) / 1000.0
    },
    'ZK': {
        # ZK = ZCK + ZOK (combined ketone, experimentally unresolved)
        'time': np.array([0.0, 1138.44, 1040.01, 857.52, 675.08, 554.11, 431.13, 308.19, 246.73,
                         183.23, 152.52, 119.76, 89.03, 74.70, 60.37, 43.99, 37.90, 25.60, 15.36, 5.12]),
        'concentration': np.array([0.0, 35.65, 34.57, 31.96, 27.83, 25.87, 22.39, 17.83, 15.22,
                                 11.96, 10.22, 8.26, 6.96, 5.65, 4.78, 3.70, 1.52, 1.30, 0.87, 0.22]) / 1000.0
    },
    'ZOD': {
        'time': np.array([0.0, 1139.18, 1040.73, 858.18, 675.64, 552.60, 431.60, 310.61, 247.05,
                         183.48, 152.72, 119.92, 91.21, 74.81, 58.40, 48.15, 37.90, 25.60, 15.36, 5.12]),
        'concentration': np.array([0.0, 10.65, 10.00, 9.57, 8.70, 7.61, 6.52, 5.43, 4.57, 3.70,
                                 3.26, 2.83, 2.39, 2.17, 1.96, 1.74, 1.52, 1.30, 0.87, 0.22]) / 1000.0
    }
}


# SECTION 3: RATE CONSTANT NAMES AND BOUNDS
# There are 10 rate constants, corresponding to the 5 reversible reactions:
#   ZCE <-> ZCD   (k1, k_1)
#   ZCE <-> ZCK   (k2, k_2)
#   ZCD <-> ZCK   (k3, k_3)
#   ZCD <-> ZOD   (k4, k_4)
#   ZOD <-> ZOK   (k5, k_5)
# These bounds constrain the local optimisation during bootstrap fitting, stopping the optimiser wandering into physically unreasonable regions when working with resampled (noisier) datasets.

RATE_CONSTANTS = ['k1', 'k_1', 'k2', 'k_2', 'k3', 'k_3', 'k4', 'k_4', 'k5', 'k_5']

PARAMETER_SPECIFIC_BOUNDS = {
    'k1':  (1e-2, 1e0),
    'k_1': (1e-3, 1e-1),
    'k2':  (1e-3, 1e0),
    'k_2': (1e-1, 1e1),
    'k3':  (1e-10, 1e-6),
    'k_3': (1e-10, 1e-6),
    'k4':  (1e-4, 1e-2),
    'k_4': (1e-3, 1e-1),
    'k5':  (1e-10, 1e-6),
    'k_5': (1e-5, 1e-3),
}

BOUNDS = [PARAMETER_SPECIFIC_BOUNDS[param] for param in RATE_CONSTANTS]

# Every species is weighted equally here; adjust this if a subset of species should dominate the bootstrap fits.
SPECIES_WEIGHTS = {
    'ZCE': 1.0,
    'ZCD': 1.0,
    'ZK':  1.0,
    'ZOD': 1.0
}

print("Inverse Kinetics Parameter Estimation for Zymonic Acid System")
print("=" * 70)
print("Experimental data loaded:")
for species, data in EXPERIMENTAL_DATA.items():
    print(f"   {species}: {len(data['time'])} time points, "
          f"max conc = {max(data['concentration'])*1000:.1f} mM")
print(f"Fitting {len(RATE_CONSTANTS)} rate constants with literature-informed bounds")
print()


# SECTION 4: DATABASE AND INPUT FILE CREATION

def create_zymonic_database(rate_constants_dict, db_filename='zymonic_system.db'):
    # Purpose: writes a temporary reaction database encoding the zymonic network with a given set of candidate rate constants.
    # Inputs are: rate_constants_dict, a mapping of rate constant name to its value in min^-1; db_filename, the filename within ./Inputs/Databases/ (default 'zymonic_system.db').
    # Outputs are: the database name without its .db extension. A temporary database encoding the zymonic reaction network is written as a side effect; each reversible reaction A <-> B is stored as two separate forward-only rows, since the simulation framework does not use the Reverse Rate Constant column directly. The database is meant to be deleted by the caller after each simulation.
    # Called by: run_zymonic_simulation(), in this file.

    reactions = [
        # --- Reaction 1: ZCE <-> ZCD ---
        {'Reactants': '(1)ZCE', 'Products': '(1)ZCD',
         'Forward Rate Constant': rate_constants_dict['k1'],
         'Reverse Rate Constant': 'none',
         'Equilibrium Constant': 'none'},
        {'Reactants': '(1)ZCD', 'Products': '(1)ZCE',
         'Forward Rate Constant': rate_constants_dict['k_1'],
         'Reverse Rate Constant': 'none',
         'Equilibrium Constant': 'none'},

        # --- Reaction 2: ZCE <-> ZCK ---
        {'Reactants': '(1)ZCE', 'Products': '(1)ZCK',
         'Forward Rate Constant': rate_constants_dict['k2'],
         'Reverse Rate Constant': 'none',
         'Equilibrium Constant': 'none'},
        {'Reactants': '(1)ZCK', 'Products': '(1)ZCE',
         'Forward Rate Constant': rate_constants_dict['k_2'],
         'Reverse Rate Constant': 'none',
         'Equilibrium Constant': 'none'},

        # --- Reaction 3: ZCD <-> ZCK ---
        {'Reactants': '(1)ZCD', 'Products': '(1)ZCK',
         'Forward Rate Constant': rate_constants_dict['k3'],
         'Reverse Rate Constant': 'none',
         'Equilibrium Constant': 'none'},
        {'Reactants': '(1)ZCK', 'Products': '(1)ZCD',
         'Forward Rate Constant': rate_constants_dict['k_3'],
         'Reverse Rate Constant': 'none',
         'Equilibrium Constant': 'none'},

        # --- Reaction 4: ZCD <-> ZOD ---
        {'Reactants': '(1)ZCD', 'Products': '(1)ZOD',
         'Forward Rate Constant': rate_constants_dict['k4'],
         'Reverse Rate Constant': 'none',
         'Equilibrium Constant': 'none'},
        {'Reactants': '(1)ZOD', 'Products': '(1)ZCD',
         'Forward Rate Constant': rate_constants_dict['k_4'],
         'Reverse Rate Constant': 'none',
         'Equilibrium Constant': 'none'},

        # --- Reaction 5: ZOD <-> ZOK ---
        {'Reactants': '(1)ZOD', 'Products': '(1)ZOK',
         'Forward Rate Constant': rate_constants_dict['k5'],
         'Reverse Rate Constant': 'none',
         'Equilibrium Constant': 'none'},
        {'Reactants': '(1)ZOK', 'Products': '(1)ZOD',
         'Forward Rate Constant': rate_constants_dict['k_5'],
         'Reverse Rate Constant': 'none',
         'Equilibrium Constant': 'none'},
    ]

    df = pd.DataFrame(reactions)

    db_dir = './Inputs/Databases'
    os.makedirs(db_dir, exist_ok=True)
    db_path = os.path.join(db_dir, db_filename)

    if os.path.exists(db_path):
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    df.to_sql('reactions', conn, index=False)
    conn.close()

    return db_filename.replace('.db', '')


def create_zymonic_input_files():
    # Purpose: writes the atmosphere and aqueous starting-condition files the simulation framework needs, matching the Perkins et al. experimental conditions.
    # Inputs are: none.
    # Outputs are: nothing is returned; the atmosphere and aqueous input files required by the simulation framework are written to ./Inputs/, overwriting any existing files there. ZCE starts at 0.1 M (the sole starting species), every other species starts at 0.0 M, temperature is set to 298.15 K, and pH to 3.0.
    # Called by: the script's entry point, at the bottom of this file (if __name__ == "__main__").

    inputs_dir = './Inputs'
    os.makedirs(inputs_dir, exist_ok=True)

    # There is no gas-phase species in this system, but the simulation framework still requires this file to exist.
    atmosphere_content = """# Atmosphere Data File
# This file contains information about species, their partial pressures, total pressure, total temperature, and the solar spectrum file in use.

# Species and their corresponding partial pressures (in atm) along with Henry's Law Constants
# Format: species, partial_pressure (atm), Henry's Law Constant (M/atm)
species, partial_pressure, constant (M/atm)


# Total pressure (in atm)
Total pressure
1.0

# Temperature (in Kelvin)
Total Temperature
298

# Solar spectrum file (relative or absolute path)
Solar Spectrum File
spectrum.dat
"""

    with open(os.path.join(inputs_dir, 'zymonic-atmosphere.dat'), 'w') as f:
        f.write(atmosphere_content)

    # This sets the aqueous species and their initial concentrations.
    aqueous_content = """# Species and their corresponding concentrations (in M)
# Format: species, concentration (M)
species, concentration
ZCE, 0.1
ZCD, 0.0
ZCK, 0.0
ZOD, 0.0
ZOK, 0.0

# Total Volume (in L)
Total Volume
1

# Temperature (in Kelvin)
Total Temperature
298.15

# pH
pH
3.0
"""

    with open(os.path.join(inputs_dir, 'zymonic-aqueous.dat'), 'w') as f:
        f.write(aqueous_content)

    print("Created zymonic input files:")
    print("   ./Inputs/zymonic-atmosphere.dat")
    print("   ./Inputs/zymonic-aqueous.dat")
    print("   ZCE: 0.1 M (initial concentration)")
    print("   All others: 0.0 M\n")



# SECTION 5: SIMULATION INTERFACE

def run_zymonic_simulation(rate_constants_array, end_time=1200, num_points=500, verbose=False):
    # Purpose: runs one full simulation of the zymonic network for a given set of candidate rate constants, and returns the resulting concentration-time curves.
    # Inputs are: rate_constants_array, the 10 rate constants in RATE_CONSTANTS order; end_time, the simulation end time in minutes (default 1200); num_points, the number of evenly spaced time points to evaluate (default 500); verbose, whether to print diagnostic output at each step.
    # Outputs are: a dictionary with keys 'time', 'ZCE', 'ZCD', 'ZCK', 'ZOD', 'ZOK', 'ZK', or None if the simulation fails for any reason. A temporary database is created with the candidate rate constants, parsed into rate law expressions and differential equations by the existing simulation framework, the starting conditions are read from the .dat input files, the ODE system is integrated with solve_ivp (LSODA), and the species concentrations are returned; the temporary database is deleted either way.
    # Called by: calculate_residuals(), in this file.

    rate_constants_dict = dict(zip(RATE_CONSTANTS, rate_constants_array))
    db_name = create_zymonic_database(rate_constants_dict)

    try:
        from rate_laws import format_species_name, parse_reactants_products
        from differentials import construct_differential_equations, create_bulletproof_ode_system
        from data import read_dat_file, read_aqueous_dat_file
        from utilities import initialise_concentrations

        if verbose:
            print(f"   Running simulation with database: {db_name}")

        # The reactions are read from the database here, and used to build the symbolic rate laws.
        db_path = f'./Inputs/Databases/{db_name}.db'
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query("SELECT * FROM reactions", conn)
        conn.close()

        rate_laws = []
        database_species = set()

        for index, row in df.iterrows():
            forward_rate_constant = row.get('Forward Rate Constant', 'none')
            reverse_rate_constant = row.get('Reverse Rate Constant', 'none')

            reactants_str = str(row.get('Reactants', ''))
            products_str = str(row.get('Products', ''))

            reactants, reactant_species = parse_reactants_products(reactants_str)
            products, product_species = parse_reactants_products(products_str)

            # These build the rate law strings, e.g. "0.1135*ZCE".
            forward_rate_law = None
            if forward_rate_constant != 'none' and reactants:
                forward_rate_law = f"{forward_rate_constant}*{reactants}"  # e.g. "0.1135*ZCE" - rate constant multiplied by the reactant concentration term

            reverse_rate_law = None
            if reverse_rate_constant != 'none' and products:
                reverse_rate_law = f"{reverse_rate_constant}*{products}"  # same pattern, for the reverse direction of this reaction

            rate_laws.append({
                'forward': forward_rate_law,
                'reverse': reverse_rate_law,
                'consumed': reactant_species,
                'produced': product_species
            })

            for s, _ in reactant_species + product_species:
                database_species.add(s)

        database_species = list(database_species)

        # The ODE system is assembled here.
        differential_equations = construct_differential_equations(rate_laws)

        # The physical parameters and initial concentrations are read here.
        species_atm, partial_pressures, henrys_constants = read_dat_file(
            './Inputs/zymonic-atmosphere.dat')
        species_aq, concentrations, total_volume, total_temperature, pH, _ = read_aqueous_dat_file(
            './Inputs/zymonic-aqueous.dat')

        initial_conditions = initialise_concentrations(
            database_species, species_aq, concentrations,
            species_atm, partial_pressures, henrys_constants
        )

        # This guards that ZCE starts at 0.1 M, regardless of the parsing outcome above.
        if 'ZCE' in database_species:
            zce_idx = database_species.index('ZCE')
            if initial_conditions[zce_idx] < 0.09:
                initial_conditions[zce_idx] = 0.1

        ode_system = create_bulletproof_ode_system(
            database_species, differential_equations,
            species_atm, partial_pressures, henrys_constants, pH
        )

        solution = solve_ivp(
            ode_system,
            (0, end_time),
            initial_conditions,
            method='LSODA',
            rtol=1e-6,
            atol=1e-9,
            t_eval=np.linspace(0, end_time, num_points),
            dense_output=True
        )

        if not solution.success:
            if verbose:
                print(f"   Simulation failed: {solution.message}")
            return None

        # The results are packaged here; ZK = ZCK + ZOK is computed as the experimental observable.
        results = {'time': solution.t}

        for species_name in ['ZCE', 'ZCD', 'ZCK', 'ZOD', 'ZOK']:
            if species_name in database_species:
                idx = database_species.index(species_name)
                results[species_name] = solution.y[idx, :]
            else:
                results[species_name] = np.zeros_like(solution.t)

        results['ZK'] = results['ZCK'] + results['ZOK']

        if verbose:
            print(f"   Simulation completed: {len(solution.t)} time points")

        return results

    except Exception as e:
        if verbose:
            print(f"   Simulation error: {e}")
        return None

    finally:
        # The temporary database is deleted here, whether or not the simulation succeeded.
        db_path = f'./Inputs/Databases/{db_name}.db'
        if os.path.exists(db_path):
            os.remove(db_path)


# SECTION 6: OBJECTIVE FUNCTION

def interpolate_simulation_to_experimental_times(sim_results, exp_data):
    # Purpose: resamples the simulated concentration curves onto each species' own experimental time points, so simulation and experiment can be compared directly.
    # Inputs are: sim_results, the output from run_zymonic_simulation(); exp_data, experimental data in the EXPERIMENTAL_DATA format.
    # Outputs are: a dictionary of the simulated concentrations, interpolated (via linear interpolation) onto each species' own experimental time points. For ZK, the combined simulated ketone (ZCK + ZOK) is interpolated onto the experimental ZK time points.
    # Called by: calculate_residuals(), in this file.

    interpolated = {}

    for species in ['ZCE', 'ZCD', 'ZOD']:
        if species in exp_data:
            exp_times = exp_data[species]['time']
            interpolated[species] = np.interp(exp_times, sim_results['time'], sim_results[species])

    if 'ZK' in exp_data:
        exp_times = exp_data['ZK']['time']
        interpolated['ZK'] = np.interp(exp_times, sim_results['time'], sim_results['ZK'])

    return interpolated


def calculate_residuals(rate_constants_array, experimental_data, weights=None):
    # Purpose: scores how well a candidate set of rate constants fits a experimental dataset, as a single weighted residual number.
    # Inputs are: rate_constants_array, the candidate rate constants in RATE_CONSTANTS order; experimental_data, experimental data in the EXPERIMENTAL_DATA format; weights, optional per-species weights (defaults to SPECIES_WEIGHTS).
    # Outputs are: the total weighted sum of squared relative residuals. Each species' residual is normalised by its own mean experimental concentration, so every species has roughly equal influence on the fit regardless of its absolute concentration; a large penalty (1e10) is returned instead if the simulation fails.
    # Called by: bootstrap_with_convergence_monitoring(), in this file, via a lambda (for both the L-BFGS-B and differential evolution fallback fits).

    sim_results = run_zymonic_simulation(rate_constants_array, end_time=1200, verbose=False)

    if sim_results is None:
        return 1e10

    sim_interpolated = interpolate_simulation_to_experimental_times(sim_results, experimental_data)

    if weights is None:
        weights = SPECIES_WEIGHTS

    total_ssr = 0.0

    for species in experimental_data.keys():
        if species in sim_interpolated:
            exp_conc = experimental_data[species]['concentration']
            sim_conc = sim_interpolated[species]

            mean_exp_conc = np.mean(exp_conc)
            if mean_exp_conc > 1e-10:
                ssr = np.sum(((exp_conc - sim_conc) / mean_exp_conc) ** 2)
            else:
                ssr = np.sum((exp_conc - sim_conc) ** 2)

            total_ssr += weights.get(species, 1.0) * ssr

    return total_ssr


# SECTION 7: BOOTSTRAP RESAMPLING AND FITTING

def bootstrap_resample_data(experimental_data, n_bootstrap=100, seed=None):
    # Purpose: generates resampled copies of the experimental dataset, for bootstrap uncertainty estimation.
    # Inputs are: experimental_data, the original experimental data in the EXPERIMENTAL_DATA format; n_bootstrap, the number of bootstrap datasets to generate; seed, an optional seed that, if given, resets numpy's random seed before generating samples (for reproducibility of a specific batch).
    # Outputs are: a list of resampled datasets, each in the EXPERIMENTAL_DATA format. For each bootstrap sample, every species' n time-concentration pairs are resampled with replacement, producing n new pairs, so the resampled dataset has the same number of points as the original but with some points duplicated and others omitted. This variation between datasets is what produces variation in the fitted parameters, from which confidence intervals can be estimated empirically.
    # Called by: bootstrap_with_convergence_monitoring(), in this file.

    if seed is not None:
        np.random.seed(seed)
        print(f"   Bootstrap resampling seed set to {seed}")

    bootstrap_samples = []

    for i in range(n_bootstrap):
        resampled_data = {}

        for species, data in experimental_data.items():
            n_points = len(data['time'])
            indices = np.random.choice(n_points, n_points, replace=True)

            resampled_data[species] = {
                'time': data['time'][indices],
                'concentration': data['concentration'][indices]
            }

        bootstrap_samples.append(resampled_data)

    return bootstrap_samples


def bootstrap_with_convergence_monitoring(best_fit_params, n_bootstrap=2245,
                                          convergence_check_interval=100):
    # Purpose: fits the model to many resampled datasets to build up confidence intervals for every rate constant, periodically checking whether those intervals have stopped changing.
    # Inputs are: best_fit_params, the best-fit parameters from Script 1, used as the starting point for every bootstrap optimisation; n_bootstrap, the total number of bootstrap samples to fit (default 2245, the paper value); convergence_check_interval, the number of samples between convergence checks (default 100).
    # Outputs are: a dictionary of per-parameter statistics (best_fit, mean, std, ci_lower, ci_upper, cv_percent, convergence_status, and the raw bootstrap values), or None if every sample fails. For each resampled dataset, rate constants are first fitted with L-BFGS-B local optimisation starting from best_fit_params (fast, since the best-fit parameters are a good starting point for any resampled dataset); if that fails, differential evolution is retried as a slower fallback. Every convergence_check_interval samples, convergence is assessed by checking whether the 95% confidence interval width for each parameter has changed by less than 5% since the previous check - this is the criterion used in Perkins et al.
    # Called by: run_fast_bootstrap_analysis(), in this file.

    print(f"\nBootstrap with convergence monitoring ({n_bootstrap} samples)...")
    print(f"   Seed: {RANDOM_SEED} (from {SEED_INFO['readable']})")
    print(f"   Convergence criterion: CI width change < 5% between checks")
    print(f"   Check interval: every {convergence_check_interval} samples")

    # Every bootstrap dataset is generated upfront here, using the datetime-based seed.
    bootstrap_data_samples = bootstrap_resample_data(
        EXPERIMENTAL_DATA, n_bootstrap, seed=RANDOM_SEED)

    bootstrap_params = []    # Fitted parameters for each successful sample
    convergence_flags = []   # True = L-BFGS-B converged; False = DE fallback

    # The confidence-interval bounds at each convergence check are stored here, for convergence monitoring.
    convergence_history = {param: [] for param in RATE_CONSTANTS}

    start_time = time.time()

    for i, bootstrap_data in enumerate(bootstrap_data_samples):

        if (i + 1) % 50 == 0:
            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed
            remaining = (n_bootstrap - i - 1) / rate
            print(f"   Sample {i+1}/{n_bootstrap} — "
                  f"rate: {rate:.1f}/min, remaining: {remaining/60:.1f} min")

        try:
            # L-BFGS-B, starting from best_fit_params, is the primary method.
            result = minimize(
                lambda x: calculate_residuals(x, bootstrap_data),
                best_fit_params,
                method='L-BFGS-B',
                bounds=BOUNDS,
                options={
                    'maxiter': 200,
                    'ftol': 1e-8,
                    'gtol': 1e-6
                }
            )

            if result.success:
                bootstrap_params.append(result.x)
                convergence_flags.append(True)

            else:
                # Differential evolution with reduced settings is used here as a fallback, for a heavily resampled dataset with a different landscape where the local method fails to converge.
                result_de = differential_evolution(
                    lambda x: calculate_residuals(x, bootstrap_data),
                    bounds=BOUNDS,
                    maxiter=100,
                    popsize=15,
                    polish=True,
                    seed=RANDOM_SEED + i   # Offset seed per sample for variety
                )

                if result_de.success:
                    bootstrap_params.append(result_de.x)
                    convergence_flags.append(False)
                else:
                    # Both methods failed, so this sample is skipped.
                    continue

        except Exception:
            continue

        # Every convergence_check_interval samples, the current 95% confidence interval bounds are computed and compared with the previous check; a parameter whose interval width has changed by less than 5% is considered converged, confirming that more samples would not substantially change the reported uncertainty.
        if (i + 1) % convergence_check_interval == 0 and len(bootstrap_params) >= 50:
            current_params = np.array(bootstrap_params)

            for j, param_name in enumerate(RATE_CONSTANTS):
                values = current_params[:, j]
                ci_lower = np.percentile(values, 2.5)
                ci_upper = np.percentile(values, 97.5)
                convergence_history[param_name].append((ci_lower, ci_upper))

            # Convergence is only checked once there are at least 2 history points to compare.
            if len(convergence_history[RATE_CONSTANTS[0]]) >= 2:
                converged_params = []

                for param_name in RATE_CONSTANTS:
                    history = convergence_history[param_name]
                    old_range = history[-2][1] - history[-2][0]
                    new_range = history[-1][1] - history[-1][0]

                    if old_range > 0:
                        percent_change = abs(new_range - old_range) / old_range * 100  # how much the 95% CI width has shifted since the last check, as a percentage
                        if percent_change < 5.0:
                            converged_params.append(param_name)

                print(f"   Convergence check at sample {i+1}: "
                      f"{len(converged_params)}/{len(RATE_CONSTANTS)} parameters "
                      f"converged (<5% CI width change)")

    total_time = time.time() - start_time

    if len(bootstrap_params) == 0:
        print("   All bootstrap samples failed.")
        return None

    bootstrap_params = np.array(bootstrap_params)
    local_converged = sum(convergence_flags)

    # The final statistics are computed here.
    param_stats = {}

    print(f"\nBootstrap Results ({len(bootstrap_params)}/{n_bootstrap} successful):")
    print(f"   Total time:        {total_time/60:.1f} minutes")
    print(f"   Local convergence: {local_converged}/{len(bootstrap_params)} "
          f"samples (remainder used DE fallback)")
    print()
    print("   Parameter      Best Fit      Mean +/- Std      95% CI           CV%    Status")
    print("   " + "-" * 90)

    for i, param_name in enumerate(RATE_CONSTANTS):
        values = bootstrap_params[:, i]

        mean_val  = np.mean(values)
        std_val   = np.std(values)
        ci_lower  = np.percentile(values, 2.5)
        ci_upper  = np.percentile(values, 97.5)
        cv_percent = (std_val / mean_val) * 100 if mean_val > 0 else float('inf')  # coefficient of variation: spread relative to the mean, as a percentage

        # The final convergence status is taken from the last two history entries here.
        if len(convergence_history[param_name]) >= 2:
            history = convergence_history[param_name]
            old_range = history[-2][1] - history[-2][0]
            new_range = history[-1][1] - history[-1][0]
            final_change = abs(new_range - old_range) / old_range * 100 if old_range > 0 else float('inf')  # same CI-width-change calculation as above, for the final reported status
            status = "CONVERGED" if final_change < 5.0 else "NOT_CONV"
        else:
            status = "UNKNOWN"

        param_stats[param_name] = {
            'best_fit': best_fit_params[i],
            'mean': mean_val,
            'std': std_val,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'cv_percent': cv_percent,
            'convergence_status': status,
            'values': values
        }

        print(f"   {param_name:12} {best_fit_params[i]:11.2e} "
              f"{mean_val:11.2e}+/-{std_val:.1e} "
              f"[{ci_lower:.1e}, {ci_upper:.1e}] {cv_percent:6.1f}% {status}")

    return param_stats


# SECTION 8: RESULTS VISUALISATION AND SAVING

def plot_parameter_distributions(param_stats, save_plots=True):
    # Purpose: draws and saves a histogram of the bootstrap-fitted values for every rate constant, showing the best fit and its confidence interval.
    # Inputs are: param_stats, the output from bootstrap_with_convergence_monitoring(); save_plots, whether to save the figure to ./Outputs/bootstrap_distributions.png (default True).
    # Outputs are: nothing is returned; a histogram of the bootstrap-fitted values is drawn for each of the 10 rate constants, with the best-fit value from Script 1's primary optimisation marked as a red dashed line, and the 2.5th/97.5th percentile confidence-interval bounds marked as orange dotted lines. The figure title includes the random seed and timestamp, so a plot can always be traced back to the run that produced it. A log x-axis is used for any parameter whose 95% confidence interval spans more than two orders of magnitude.
    # Called by: run_fast_bootstrap_analysis(), in this file.

    fig, axes = plt.subplots(2, 5, figsize=(15, 8))
    fig.suptitle(
        f'Bootstrap Parameter Distributions\n'
        f'Seed={RANDOM_SEED}  Generated: {SEED_INFO["readable"]}',
        fontsize=14, fontweight='bold'
    )

    for i, (param_name, stats) in enumerate(param_stats.items()):
        row, col = i // 5, i % 5
        ax = axes[row, col]

        ax.hist(stats['values'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')

        ax.axvline(stats['best_fit'], color='red', linestyle='--',
                  linewidth=2, label='Best Fit')
        ax.axvline(stats['ci_lower'], color='orange', linestyle=':', alpha=0.7)
        ax.axvline(stats['ci_upper'], color='orange', linestyle=':', alpha=0.7)

        ax.set_xlabel(f'{param_name} (min^-1)')
        ax.set_ylabel('Frequency')
        ax.set_title(param_name)
        ax.legend()

        if stats['ci_upper'] / stats['ci_lower'] > 100:
            ax.set_xscale('log')

    plt.tight_layout()

    if save_plots:
        os.makedirs('./Outputs', exist_ok=True)
        plt.savefig('./Outputs/bootstrap_distributions.png', dpi=300, bbox_inches='tight')
        print("   Plot saved: ./Outputs/bootstrap_distributions.png")

    plt.show()


def save_bootstrap_results_detailed(best_fit_params, param_stats, problematic_params):
    # Purpose: saves the full bootstrap results (statistics, seed, and reproducibility info) to disk as JSON and a CSV table.
    # Inputs are: best_fit_params, the best-fit parameters from Script 1; param_stats, the output from bootstrap_with_convergence_monitoring(); problematic_params, the names of parameters whose best-fit value falls outside its own bootstrap confidence interval (these correspond to the red-shaded parameters in Perkins et al. Figure 7).
    # Outputs are: nothing is returned; two files are written - bootstrap_results_detailed.json (the full statistics, including the seed, timestamp, convergence status, and whether each best-fit parameter falls within its confidence interval, so this specific set of results can be reproduced by setting RANDOM_SEED to the saved value) and publication_ready_results.csv (a tabular summary suitable for a manuscript or report).
    # Called by: run_fast_bootstrap_analysis(), in this file.

    os.makedirs('./Outputs', exist_ok=True)

    results_dict = {
        'timestamp': datetime.now().isoformat(),
        'methodology': 'Bootstrap resampling with convergence monitoring',
        'reproducibility': {
            'random_seed': int(RANDOM_SEED),
            'seed_generated_from': SEED_INFO['readable'],
            'seed_timestamp': SEED_INFO['timestamp'],
            'note': 'To reproduce, set RANDOM_SEED to the value above at the '
                    'top of the seed control section'
        },
        'best_fit_parameters': {
            param: float(val)
            for param, val in zip(RATE_CONSTANTS, best_fit_params)
        },
        'problematic_parameters': problematic_params,
        'bootstrap_statistics': {}
    }

    for param_name, stats in param_stats.items():
        results_dict['bootstrap_statistics'][param_name] = {
            'best_fit': float(stats['best_fit']),
            'mean': float(stats['mean']),
            'std': float(stats['std']),
            'ci_lower_95': float(stats['ci_lower']),
            'ci_upper_95': float(stats['ci_upper']),
            'coefficient_of_variation_percent': float(stats['cv_percent']),
            'convergence_status': stats['convergence_status'],
            'best_fit_within_ci': (
                stats['ci_lower'] <= stats['best_fit'] <= stats['ci_upper']
            )
        }

    with open('./Outputs/bootstrap_results_detailed.json', 'w') as f:
        json.dump(results_dict, f, indent=2)

    table_data = []
    for param in RATE_CONSTANTS:
        stats = param_stats[param]
        within_ci = stats['ci_lower'] <= stats['best_fit'] <= stats['ci_upper']
        table_data.append({
            'Parameter':  param,
            'Best_Fit':   f"{stats['best_fit']:.3e}",
            'CI_95_Lower': f"{stats['ci_lower']:.3e}",
            'CI_95_Upper': f"{stats['ci_upper']:.3e}",
            'CI_Range':   f"[{stats['ci_lower']:.2e}, {stats['ci_upper']:.2e}]",
            'CV_Percent': f"{stats['cv_percent']:.1f}%",
            'Convergence': stats['convergence_status'],
            'Within_CI':  'yes' if within_ci else 'no'
        })

    pd.DataFrame(table_data).to_csv('./Outputs/publication_ready_results.csv', index=False)

    print("\nDetailed results saved:")
    print("   ./Outputs/bootstrap_results_detailed.json")
    print("   ./Outputs/publication_ready_results.csv")


# SECTION 9: MAIN BOOTSTRAP ANALYSIS FUNCTION

def run_fast_bootstrap_analysis():
    # Purpose: the main entry point that runs the full bootstrap uncertainty analysis on the best-fit parameters from Script 1, and saves and plots the results.
    # Inputs are: none.
    # Outputs are: param_stats from bootstrap_with_convergence_monitoring(), or None if the bootstrap fails. This is the primary entry point for this script: it takes the hard-coded best-fit parameters from Script 1 below, runs bootstrap_with_convergence_monitoring() for n_bootstrap samples, checks whether each best-fit parameter falls within its own bootstrap confidence interval (flagging any that don't as 'problematic' - these correspond to the red-shaded parameters in Perkins et al. Figure 7), generates the parameter distribution plots, and saves all results to ./Outputs/.
    # To update the best-fit parameters: after running Script 1, copy the printed best-fit parameter values into the best_fit_params array below, keeping the RATE_CONSTANTS order: [k1, k_1, k2, k_2, k3, k_3, k4, k_4, k5, k_5].
    # Called by: the script's entry point, at the bottom of this file (if __name__ == "__main__").

    print("Bootstrap Analysis (datetime-seeded)")
    print("=" * 60)

    # -------------------------------------------------------------------------
    # THE BEST-FIT PARAMETERS FROM SCRIPT 1 ARE PASTED HERE.
    # These values were obtained by running main_optimisation_only() or
    # fit_with_multiple_starts() in inverse_kinetics_optimisation.py.
    # -------------------------------------------------------------------------
    best_fit_params = np.array([
        1.1135e-01,  # k1  — ZCE -> ZCD
        4.6810e-02,  # k_1 — ZCD -> ZCE
        1.0000e-08,  # k2  — ZCE -> ZCK
        9.9999e-09,  # k_2 — ZCK -> ZCE
        1.0000e-04,  # k3  — ZCD -> ZCK
        1.0000e-03,  # k_3 — ZCK -> ZCD
        1.4257e-03,  # k4  — ZCD -> ZOD
        4.4387e-03,  # k_4 — ZOD -> ZCD
        9.6230e-02,  # k5  — ZOD -> ZOK
        3.0495e-02   # k_5 — ZOK -> ZOD
    ])

    print("Best-fit parameters (from Script 1):")
    for param, value in zip(RATE_CONSTANTS, best_fit_params):
        print(f"   {param} = {value:.4e} min^-1")

    print(f"\nBootstrap settings:")
    print(f"   Random seed:          {RANDOM_SEED}")
    print(f"   Seed generated from:  {SEED_INFO['readable']}")
    print(f"   Each run produces different results (different resamples)")
    print(f"   Seed is saved to output files for reproducibility")

    # This runs the full Perkins et al. paper-style bootstrap (2245 samples); change n_bootstrap here to use fewer samples for a faster test.
    n_bootstrap = 2245
    print(f"\nRunning {n_bootstrap} bootstrap samples...")

    param_stats = bootstrap_with_convergence_monitoring(
        best_fit_params, n_bootstrap=n_bootstrap)

    if param_stats is None:
        print("Bootstrap analysis failed.")
        return None

    # Every best-fit parameter is checked here against its own bootstrap confidence interval: one that falls outside it suggests the optimisation found a different region of parameter space than the bootstrap distribution, which may indicate a local minimum or an identifiability issue.
    print("\nParameter Validation (best-fit value vs bootstrap 95% CI):")
    problematic_params = []

    for param_name, stats in param_stats.items():
        best_fit  = stats['best_fit']
        ci_lower  = stats['ci_lower']
        ci_upper  = stats['ci_upper']
        within_ci = ci_lower <= best_fit <= ci_upper

        status = "GOOD" if within_ci else "OUTSIDE CI"
        print(f"   {param_name}: {best_fit:.3e} in [{ci_lower:.3e}, {ci_upper:.3e}]  {status}")

        if not within_ci:
            problematic_params.append(param_name)

    if problematic_params:
        print(f"\nWARNING: {len(problematic_params)} parameter(s) outside their bootstrap CIs:")
        print(f"   {', '.join(problematic_params)}")
        print("   This suggests potential optimisation or identifiability issues.")
        print("   (Equivalent to the red-shaded parameters in Perkins et al. Fig. 7)")
    else:
        print("\nAll parameters fall within their bootstrap confidence intervals.")
        print("   This indicates good parameter identifiability.")

    plot_parameter_distributions(param_stats, save_plots=True)
    save_bootstrap_results_detailed(best_fit_params, param_stats, problematic_params)

    print(f"\nBootstrap analysis complete.")
    print(f"   Seed {RANDOM_SEED} saved to output files.")
    print(f"   To reproduce these results: set RANDOM_SEED = {RANDOM_SEED} "
          f"in the seed control section.")

    return param_stats


# SECTION 10: ENTRY POINT

if __name__ == "__main__":

    print("\nBootstrap Uncertainty Estimation — Zymonic Acid Kinetics")
    print("=" * 60)
    print("Prerequisites:")
    print("   1. Script 1 (inverse_kinetics_optimisation.py) has been run")
    print("   2. Best-fit parameters have been pasted into best_fit_params")
    print("      in run_fast_bootstrap_analysis() in this file")
    print()
    print(f"Reproducibility:")
    print(f"   Seed:      {RANDOM_SEED}")
    print(f"   Timestamp: {SEED_INFO['readable']}")
    print(f"   Results will differ on each run (different random seed)")
    print(f"   Seed is saved to ./Outputs/ for reproducibility\n")

    create_zymonic_input_files()
    param_stats = run_fast_bootstrap_analysis()
