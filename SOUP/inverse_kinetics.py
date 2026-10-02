# Inverse Kinetics Parameter Estimation for the Zymonic Acid System (can be adapted to more systems)

# PURPOSE
# This script estimates the rate constants for the 5-species zymonic acid interconversion network, by fitting a kinetic simulation to experimental concentration-time data.

# WORKFLOW
# Step 1: The experimental data and parameter bounds are defined at the top of the file. 
# Step 2: A temporary database is created on each function evaluation, encoding the reaction network with the candidate rate constants. This database is read by the existing simulation framework (rate_laws.py, differentials.py, etc.) and deleted after use.
# Step 3: The simulation is run using scipy's solve_ivp (LSODA method). Simulated concentrations are interpolated onto the experimental time points and compared via a weighted sum of squared relative residuals.
# Step 4: An optimisation algorithm (differential evolution, or multi-start differential evolution) minimises the residuals to find the best-fit rate constants.
# Step 5: Results are plotted and saved to ./Outputs/.

# OUTPUTS
# - ./Outputs/zymonic_kinetics_fit.png         : Fit vs experimental data plot
# - ./Outputs/bootstrap_distributions.png      : Bootstrap parameter distributions
# - ./Outputs/zymonic_kinetics_results.json    : Best-fit parameters as JSON
# - ./Outputs/zymonic_kinetics_parameters.csv  : Best-fit parameters as CSV

# RUNNING
# At the bottom of this file, set choice = "2" for a standard single optimisation run, or choice = "3" for multi-start optimisation. choice = "1" runs a quick sanity-check simulation only.

# DEPENDENCIES
# This script imports from the existing simulation framework: rate_laws.py, differentials.py, data.py, utilities.py. These must be importable from the working directory.

# The best-fit parameters produced by this script are intended to be copied manually into Script 2 (seeded_inverse_kinetics.py), for uncertainty estimation via bootstrap resampling.

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

# scipy optimisation warnings that arise from intermediate parameter evaluations during differential evolution (e.g. stiff ODE solver warnings) are suppressed here, since they don't indicate a problem with the final result.
warnings.filterwarnings('ignore', category=RuntimeWarning)


# =============================================================================
# SECTION 1: EXPERIMENTAL DATA
# =============================================================================
# This is the concentration-time data for each observed species, converted from a percentage of total zymonic acid to a molar concentration (M).
# ZK represents the combined ketone measurement (ZCK + ZOK), since these two species cannot be resolved experimentally.

# Data source: Perkins et al.
# Units: time in minutes, concentration in mol/L (M).

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


# =============================================================================
# SECTION 2: RATE CONSTANT NAMES AND BOUNDS
# =============================================================================
# There are 10 rate constants, corresponding to the 5 reversible reactions below.
# Naming convention: k_N is the reverse of kN.

# Reaction network:
#   ZCE <-> ZCD   (k1, k_1)
#   ZCE <-> ZCK   (k2, k_2)
#   ZCD <-> ZCK   (k3, k_3)
#   ZCD <-> ZOD   (k4, k_4)
#   ZOD <-> ZOK   (k5, k_5)

# The bounds are informed by the expected order of magnitude from Perkins et al. Figure 7, and constrain the search to physically realistic regions of parameter space. Units are min^-1 throughout.

RATE_CONSTANTS = ['k1', 'k_1', 'k2', 'k_2', 'k3', 'k_3', 'k4', 'k_4', 'k5', 'k_5']

PARAMETER_SPECIFIC_BOUNDS = {
    'k1':  (1e-2, 1e0),     # ZCE -> ZCD: major forward pathway, moderate-fast
    'k_1': (1e-3, 1e-1),    # ZCD -> ZCE: reverse of major pathway, slower
    'k2':  (1e-3, 1e0),     # ZCE -> ZCK: secondary forward pathway
    'k_2': (1e-1, 1e1),     # ZCK -> ZCE: fast reverse (near equilibrium)
    'k3':  (1e-10, 1e-6),   # ZCD -> ZCK: very slow (literature ~1e-8)
    'k_3': (1e-10, 1e-6),   # ZCK -> ZCD: very slow reverse
    'k4':  (1e-4, 1e-2),    # ZCD -> ZOD: ring opening, slow
    'k_4': (1e-3, 1e-1),    # ZOD -> ZCD: ring closing, faster than opening
    'k5':  (1e-10, 1e-6),   # ZOD -> ZOK: very slow (literature ~5e-9)
    'k_5': (1e-5, 1e-3),    # ZOK -> ZOD: faster reverse (literature ~9e-5)
}

# This converts PARAMETER_SPECIFIC_BOUNDS to the list-of-tuples format required by scipy.optimize.
BOUNDS = [PARAMETER_SPECIFIC_BOUNDS[param] for param in RATE_CONSTANTS]


# =============================================================================
# SECTION 3: SPECIES WEIGHTS
# =============================================================================
# Each species contributes to the total residual with a multiplicative weight.
# A weight of 1.0 treats all species equally; increasing a weight penalises a poor fit to that species more heavily during optimisation.

SPECIES_WEIGHTS = {
    'ZCE': 1.0,
    'ZCD': 1.0,
    'ZK':  1.0,
    'ZOD': 1.0
}

# A summary of the loaded data is printed on import, for a quick sanity check.
print("Inverse Kinetics Parameter Estimation for Zymonic Acid System")
print("=" * 70)
print("Experimental data loaded:")
for species, data in EXPERIMENTAL_DATA.items():
    print(f"   {species}: {len(data['time'])} time points, "
          f"max conc = {max(data['concentration'])*1000:.1f} mM")
print(f"Fitting {len(RATE_CONSTANTS)} rate constants with literature-informed bounds")
print()


# =============================================================================
# SECTION 4: DATABASE AND INPUT FILE CREATION
# =============================================================================

def create_zymonic_database(rate_constants_dict, db_filename='zymonic_system.db'):
    # Purpose: writes a temporary SQLite reaction database encoding the zymonic network with a given set of candidate rate constants.
    # Inputs are: rate_constants_dict, a mapping of rate constant name (e.g. 'k1') to its value in min^-1; db_filename, the filename for the database within ./Inputs/Databases/.
    # Outputs are: the database name without its .db extension, as expected by the simulation framework. A fresh SQLite database encoding the zymonic reaction network is written as a side effect (any existing file at that path is deleted first); it should be deleted by the caller once the simulation has finished with it. Each reversible reaction A <-> B is stored as two separate rows (A -> B and B -> A, each with its own forward rate constant), since the 'Reverse Rate Constant' column is always 'none' - reversibility is handled explicitly, not via that column.
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

    # A clean database is always started from here, to avoid stale data from a previous call.
    if os.path.exists(db_path):
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    df.to_sql('reactions', conn, index=False)
    conn.close()

    return db_filename.replace('.db', '')


def create_zymonic_input_files():
    # Purpose: writes the atmosphere and aqueous starting-condition files the simulation framework needs, matching the Perkins et al. experimental conditions.
    # Inputs are: none.
    # Outputs are: nothing is returned; the atmosphere and aqueous input files required by the simulation framework are written to ./Inputs/, overwriting any existing files there. Temperature is set to 298.15 K and pH to 3.0, to match the experimental conditions, and ZCE is given a starting concentration of 0.1 M (the sole starting species; every other species starts at 0.0 M).
    # Called by: main_optimisation_only(), main_optimisation_with_bootstrap() and test_simulation(), all in this file.

    inputs_dir = './Inputs'
    os.makedirs(inputs_dir, exist_ok=True)

    # There is no gas-phase species in this purely aqueous system, but the simulation framework still requires an atmosphere file to exist in this format.
    atmosphere_content = """# Atmosphere Data File
# This file contains information about species, their partial pressures, total pressure, total temperature, and the solar spectrum file in use.
# Please follow the format below to update or add new data.

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

    # The aqueous file sets the initial concentration of each species: ZCE starts at 0.1 M, and all others start at zero.
    aqueous_content = """# Species and their corresponding concentrations (in M)
# Format: species, concentration
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
    print("   All others: 0.0 M")


# =============================================================================
# SECTION 5: SIMULATION INTERFACE
# =============================================================================

def run_zymonic_simulation(rate_constants_array, end_time=1200, num_points=500, verbose=False):
    # Purpose: runs one full simulation of the zymonic network for a given set of candidate rate constants, and returns the resulting concentration-time curves.
    # Inputs are: rate_constants_array, the 10 rate constants in the order defined by RATE_CONSTANTS; end_time, the simulation end time in minutes (default 1200, i.e. 20 hours); num_points, the number of evenly spaced time points to evaluate (default 500); verbose, whether to print detailed diagnostic output while the simulation runs.
    # Outputs are: a dictionary with keys 'time', 'ZCE', 'ZCD', 'ZCK', 'ZOD', 'ZOK', 'ZK', or None if the simulation fails for any reason. This function is the bridge between the optimiser (which works with plain arrays of numbers) and the existing simulation framework (which reads from databases and .dat files): it converts the rate constants to a temporary SQLite database, parses that database into rate laws and differential equations, reads the starting conditions written by create_zymonic_input_files(), solves the ODE system with scipy's solve_ivp (LSODA), and returns the concentration of all 5 species plus the combined ZK observable, deleting the temporary database afterwards either way. The ZCE starting concentration is checked after being read from file and corrected to 0.1 M if the file parser returns a lower value, as a guard against parsing edge cases.
    # Called by: calculate_residuals(), fit_rate_constants_differential_evolution(), fit_with_multiple_starts(), plot_fit_results() and test_simulation(), all in this file.

    rate_constants_dict = dict(zip(RATE_CONSTANTS, rate_constants_array))
    db_name = create_zymonic_database(rate_constants_dict)

    try:
        from rate_laws import format_species_name, parse_reactants_products
        from differentials import construct_differential_equations, create_bulletproof_ode_system
        from data import read_dat_file, read_aqueous_dat_file
        from utilities import initialise_concentrations

        if verbose:
            print(f"   Running simulation with database: {db_name}")

        # The reaction database is read here, and used to construct rate law expressions directly from the numeric constants currently stored in it, so they reflect the current candidate parameters.
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

            # These build the symbolic rate law strings, e.g. "0.1135*ZCE".
            forward_rate_law = None
            if forward_rate_constant != 'none' and reactants:
                forward_rate_law = f"{forward_rate_constant}*{reactants}"  # e.g. "0.1135*ZCE" i.e., rate constant multiplied by the reactant concentration term

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

        if verbose:
            print(f"   Generated {len(rate_laws)} rate laws for {len(database_species)} species")
            print(f"   Species found: {database_species}")
            non_zero_laws = [law for law in rate_laws if law['forward'] is not None]
            if non_zero_laws:
                print(f"   Sample rate law: {non_zero_laws[0]['forward']}")

        differential_equations = construct_differential_equations(rate_laws)

        if verbose:
            print(f"   Constructed {len(differential_equations)} differential equations")
            sample_species = list(differential_equations.keys())[0] if differential_equations else None
            if sample_species:
                print(f"   Sample equation for {sample_species}: {differential_equations[sample_species]}")

        species_atm, partial_pressures, henrys_constants = read_dat_file(
            './Inputs/zymonic-atmosphere.dat')
        species_aq, concentrations, total_volume, total_temperature, pH, _ = read_aqueous_dat_file(
            './Inputs/zymonic-aqueous.dat')

        if verbose:
            print(f"   Aqueous species from file: {species_aq}")
            print(f"   Concentrations from file:  {concentrations}")

        # The aqueous starting concentrations are mapped here onto the species order used internally by the ODE solver.
        initial_conditions = initialise_concentrations(
            database_species, species_aq, concentrations,
            species_atm, partial_pressures, henrys_constants
        )

        # This guards that ZCE starts at 0.1 M, as specified in the .dat file: if initialise_concentrations() returns a lower value (e.g. due to a parsing mismatch), it is overridden to the correct value here.
        if 'ZCE' in database_species:
            zce_idx = database_species.index('ZCE')
            if initial_conditions[zce_idx] < 0.09:
                initial_conditions[zce_idx] = 0.1
                if verbose:
                    print("   ZCE initial concentration corrected to 0.1 M")

        if verbose:
            print("   Initial conditions (M):")
            for i, species in enumerate(database_species):
                print(f"     {species}: {initial_conditions[i]:.6f}")

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

        # The results are packaged here into a dictionary keyed by species name; ZK (the combined ketone) is computed as ZCK + ZOK, to match the experimental observable.
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
        # The temporary database is always cleaned up here, even if an error occurred above.
        db_path = f'./Inputs/Databases/{db_name}.db'
        if os.path.exists(db_path):
            os.remove(db_path)


# =============================================================================
# SECTION 6: OBJECTIVE FUNCTION
# =============================================================================

def interpolate_simulation_to_experimental_times(sim_results, exp_data):
    # Purpose: resamples the simulated concentration curves onto each species' own experimental time points, so simulation and experiment can be compared directly.
    # Inputs are: sim_results, the output from run_zymonic_simulation(); exp_data, experimental data in the EXPERIMENTAL_DATA format.
    # Outputs are: a dictionary of the simulated concentrations, interpolated (via linear interpolation) onto each species' own experimental time points, since the experimental time points are irregular and differ between species - this makes the simulation directly comparable to the experimental data. For ZK, the combined ketone (ZCK + ZOK) from the simulation is interpolated onto the ZK experimental times.
    # Called by: calculate_residuals(), fit_rate_constants_differential_evolution() and fit_with_multiple_starts(), all in this file.

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
    # Purpose: scores how well a candidate set of rate constants fits the experimental data, as a single weighted residual number.
    # Inputs are: rate_constants_array, the candidate rate constants in RATE_CONSTANTS order; experimental_data, experimental data in the EXPERIMENTAL_DATA format; weights, optional per-species weights (defaults to SPECIES_WEIGHTS if not given).
    # Outputs are: the total weighted sum of squared relative residuals which is the core metric the optimiser minimises. For each species, the residual is normalised by that species' mean experimental concentration (SSR_i = sum(((exp_j - sim_j) / mean(exp))^2)), which stops a species with a high absolute concentration (e.g. ZCD) from dominating the fit at the expense of one present at lower concentrations (e.g. ZOD); the total is the weighted sum of these across every species. If the simulation fails, a large penalty (1e10) is returned instead, so the optimiser moves away from that region of parameter space.
    # Called by: objective_function(), in this file; and fast_bootstrap_from_best_fit(), in this file, via a lambda.

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
                # These are relative (normalised) residuals, so they are scale-independent.
                ssr = np.sum(((exp_conc - sim_conc) / mean_exp_conc) ** 2)
            else:
                # Absolute residuals are used instead as a fallback, for a concentration that is close to zero.
                ssr = np.sum((exp_conc - sim_conc) ** 2)

            total_ssr += weights.get(species, 1.0) * ssr

    return total_ssr


def objective_function(rate_constants_array):
    # Purpose: the function the optimiser actually minimises.
    # Inputs are: rate_constants_array, the candidate rate constants passed in by the optimiser.
    # Outputs are: the total weighted residual from calculate_residuals(), or 1e10 if that call raises any exception - this is a thin wrapper that exists so optimisation can continue even if an individual evaluation fails outright.
    # Called by: fit_rate_constants_differential_evolution(), passed as the function differential_evolution() minimises; and main_optimisation_only(), for the same reason in its backup optimisation attempt; both in this file.
    try:
        return calculate_residuals(rate_constants_array, EXPERIMENTAL_DATA)
    except Exception:
        return 1e10


# =============================================================================
# SECTION 7: OPTIMISATION METHODS
# =============================================================================

def fit_rate_constants_differential_evolution(experimental_data, max_iterations=2000):
    # Purpose: runs the primary global optimisation, searching the full bounded parameter space for the best-fit rate constants.
    # Inputs are: experimental_data, experimental data in the EXPERIMENTAL_DATA format; max_iterations, the maximum number of differential evolution generations (default 2000).
    # Outputs are: the scipy OptimizeResult from differential evolution; result.x holds the best-fit parameters. Differential evolution is a population-based global optimisation algorithm, with many local minima, since it explores the full bounded parameter space before converging. The settings used are: popsize=30 (a large population, for thorough exploration), mutation=(0.3, 1.9) (a wide mutation range - aggressive exploration early on, fine-tuning later), recombination=0.9 (high recombination, mixing candidate solutions freely), polish=True (an L-BFGS-B local refinement step is applied after differential evolution converges), init='sobol' (a Sobol sequence spreads the initial population evenly across the bounded space, which works better than random placement for log-scale parameters), and seed=42 (a fixed seed, so this optimisation step is reproducible). A progress callback reports the residual and convergence at every iteration.
    # Called by: main_optimisation_only() and main_optimisation_with_bootstrap(), both in this file.

    print("Starting differential evolution optimisation...")
    print(f"   Max iterations: {max_iterations}")
    print(f"   Population size: 30 (large, for thorough exploration)")
    print(f"   Species weights: {SPECIES_WEIGHTS}")
    print(f"   Residuals: relative (normalised by mean experimental concentration)")

    start_time = time.time()

    def callback(xk, convergence):
        # Purpose: reports optimisation progress after every generation, and periodically prints the current best parameter estimates.
        # Inputs are: xk, the current best candidate parameter vector; convergence, differential evolution's own convergence measure for this generation.
        # Outputs are: nothing is returned; progress is printed for this generation, and every 3 minutes the current parameter estimates are also printed.
        # Called by: scipy's differential_evolution(), once per generation, since this is passed to it as its callback argument.
        current_residual = objective_function(xk)
        elapsed = time.time() - start_time
        print(f"   Progress: residual={current_residual:.4e}, "
              f"convergence={convergence:.4f}, elapsed={elapsed:.1f}s")

        if elapsed > 0 and int(elapsed) % 180 == 0:
            print("   Current parameter estimates:")
            for param, value in zip(RATE_CONSTANTS, xk):
                magnitude = int(np.log10(value)) if value > 0 else -99  # order of magnitude, e.g. 3 for a value near 1e3; -99 flags a non-positive value
                print(f"     {param} = {value:.2e}  (order 10^{magnitude})")

    result = differential_evolution(
        objective_function,
        bounds=BOUNDS,
        maxiter=max_iterations,
        popsize=30,
        mutation=(0.3, 1.9),
        recombination=0.9,
        atol=1e-12,
        tol=1e-10,
        seed=42,
        callback=callback if max_iterations > 500 else None,
        disp=True,
        workers=1,
        polish=True,
        init='sobol'
    )

    elapsed_time = time.time() - start_time

    print(f"\nOptimisation completed in {elapsed_time:.1f} seconds")
    print(f"   Success: {result.success}")
    print(f"   Final residual: {result.fun:.4e}")
    print(f"   Function evaluations: {result.nfev}")
    print(f"   Iterations: {result.nit}")

    if result.success:
        fitted_constants = dict(zip(RATE_CONSTANTS, result.x))
        print("\nFitted Rate Constants (grouped by order of magnitude):")

        # Parameters are grouped by order of magnitude here, for easier interpretation.
        by_magnitude = {}
        for name, value in fitted_constants.items():
            magnitude = int(np.log10(value)) if value > 0 else -99  # order of magnitude, used as the grouping key below
            by_magnitude.setdefault(magnitude, []).append((name, value))

        for magnitude in sorted(by_magnitude.keys()):
            print(f"   Order 10^{magnitude}:")
            for name, value in by_magnitude[magnitude]:
                print(f"     {name} = {value:.4e} min^-1")

        # The goodness of fit (R-squared) for each species is evaluated here.
        print("\nModel Performance (R² per species):")
        sim_results = run_zymonic_simulation(result.x, end_time=1200, verbose=False)
        if sim_results is not None:
            sim_interpolated = interpolate_simulation_to_experimental_times(
                sim_results, experimental_data)

            for species in experimental_data.keys():
                if species in sim_interpolated:
                    exp_conc = experimental_data[species]['concentration']
                    sim_conc = sim_interpolated[species]

                    ss_res = np.sum((exp_conc - sim_conc) ** 2)
                    ss_tot = np.sum((exp_conc - np.mean(exp_conc)) ** 2)
                    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
                    print(f"   {species}: R² = {r_squared:.4f}")

    return result


def fit_with_multiple_starts(experimental_data, n_starts=3):
    # Purpose: runs several independent optimisations from different random seeds and keeps the best result, to reduce the risk of settling on a poor local minimum.
    # Inputs are: experimental_data, experimental data in the EXPERIMENTAL_DATA format; n_starts, the number of independent optimisation runs (default 3).
    # Outputs are: the scipy OptimizeResult with the lowest residual across all the starts, or None if every run fails. Because differential evolution uses a stochastic population, different random seeds can converge to different local minima; running several independent optimisations and keeping the best result reduces the risk of settling on a poor one. Each run uses a different seed (42, 142, 242, ...) and a reduced iteration count compared to fit_rate_constants_differential_evolution(), trading per-run thoroughness for breadth across starts.
    # Called by: the script's entry point, at the bottom of this file (if __name__ == "__main__"), under choice == "3".

    print(f"Multi-start optimisation with {n_starts} independent starts...")
    print("   Using literature-constrained bounds")

    best_result = None
    best_residual = float('inf')
    all_results = []

    for i in range(n_starts):
        print(f"\n--- Start {i+1}/{n_starts} (seed={42 + i*100}) ---")

        result = differential_evolution(
            objective_function,
            bounds=BOUNDS,
            maxiter=1000,
            popsize=20,
            mutation=(0.5, 1.7),
            recombination=0.9,
            atol=1e-12,
            tol=1e-10,
            seed=42 + i * 100,
            polish=True,
            init='sobol',
            disp=False
        )

        all_results.append(result)

        if result.success and result.fun < best_residual:
            best_residual = result.fun
            best_result = result
            print(f"   New best result. Residual: {result.fun:.4e}")
        elif result.success:
            print(f"   Converged but not best. Residual: {result.fun:.4e}")
        else:
            print(f"   Did not converge. Message: {result.message}")

    print(f"\nMulti-start Summary:")
    print(f"   Successful runs: {sum(1 for r in all_results if r.success)}/{n_starts}")
    print(f"   Best residual:   {best_residual:.4e}")

    if best_result and best_result.success:
        fitted_constants = dict(zip(RATE_CONSTANTS, best_result.x))
        print("\nBest Fitted Rate Constants:")
        for i, (name, value) in enumerate(fitted_constants.items()):
            bounds_range = BOUNDS[i]
            pct = ((value - bounds_range[0]) / (bounds_range[1] - bounds_range[0])) * 100  # where this value sits within its search range, as a percentage (0% = lower bound, 100% = upper bound)
            print(f"   {name} = {value:.4e} min^-1  ({pct:.1f}% of search range)")

        # R-squared and MAPE are evaluated here for each species.
        print("\nFinal Model Performance:")
        sim_results = run_zymonic_simulation(best_result.x, end_time=1200, verbose=False)
        if sim_results is not None:
            sim_interpolated = interpolate_simulation_to_experimental_times(
                sim_results, experimental_data)

            for species in experimental_data.keys():
                if species in sim_interpolated:
                    exp_conc = experimental_data[species]['concentration']
                    sim_conc = sim_interpolated[species]

                    ss_res = np.sum((exp_conc - sim_conc) ** 2)
                    ss_tot = np.sum((exp_conc - np.mean(exp_conc)) ** 2)
                    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
                    mape = np.mean(np.abs((exp_conc - sim_conc) / (exp_conc + 1e-10))) * 100  # mean absolute percentage error; 1e-10 avoids a divide-by-zero at points where exp_conc is ~0
                    print(f"   {species}: R² = {r_squared:.4f}, MAPE = {mape:.1f}%")

    return best_result


# =============================================================================
# SECTION 8: BOOTSTRAP UNCERTAINTY ESTIMATION
# =============================================================================
# The bootstrap functions below are included for completeness, but are not the primary bootstrap tool. For a full bootstrap analysis with convergence monitoring, use Script 2 (inverse_kinetics_bootstrap.py), after obtaining best-fit parameters from this script.

def bootstrap_resample_data(experimental_data, n_bootstrap=100):
    # Purpose: generates resampled copies of the experimental dataset, for bootstrap uncertainty estimation.
    # Inputs are: experimental_data, experimental data in the EXPERIMENTAL_DATA format; n_bootstrap, the number of bootstrap datasets to generate.
    # Outputs are: a list of resampled datasets, each in the EXPERIMENTAL_DATA format. For each bootstrap sample, every species' time and concentration arrays are resampled independently (n points drawn from n, with replacement). This introduces variability into the fitted parameters, which is what allows confidence intervals to be estimated empirically.
    # Called by: fast_bootstrap_from_best_fit(), in this file.

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


def fast_bootstrap_from_best_fit(best_fit_params, n_bootstrap=100, verbose=True):
    # Purpose: estimates the uncertainty on each rate constant, by refitting the model to many resampled versions of the experimental data.
    # Inputs are: best_fit_params, the best-fit parameters from the primary optimisation, used as the starting point for every bootstrap fit; n_bootstrap, the number of bootstrap samples (default 100); verbose, whether to print progress every 20 samples.
    # Outputs are: a dictionary of per-parameter statistics (mean, standard deviation, 95% confidence interval, coefficient of variation, and the raw bootstrap values), or None if every bootstrap sample fails. Each bootstrap dataset is fitted with L-BFGS-B local optimisation, starting from the best-fit parameters; this is faster than a full global optimisation for every sample, since the best-fit parameters are already a good starting point for any resampled dataset. If local optimisation fails, differential evolution with reduced settings is used as a backup. A parameter with a coefficient of variation above 100% may be poorly identifiable from the available data, and may need additional experimental constraints.
    # Called by: main_optimisation_with_bootstrap(), in this file.

    print(f"\nFast bootstrap from best-fit parameters ({n_bootstrap} samples)...")
    print("   Using best-fit parameters as starting points for local optimisation")

    bootstrap_data_samples = bootstrap_resample_data(EXPERIMENTAL_DATA, n_bootstrap)
    bootstrap_params = []
    convergence_flags = []

    start_time = time.time()

    for i, bootstrap_data in enumerate(bootstrap_data_samples):
        if verbose and (i + 1) % 20 == 0:
            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed
            remaining = (n_bootstrap - i - 1) / rate
            print(f"   Sample {i+1}/{n_bootstrap} — "
                  f"rate: {rate:.1f}/min, remaining: {remaining/60:.1f} min")

        try:
            result = minimize(
                lambda x: calculate_residuals(x, bootstrap_data),
                best_fit_params,
                method='L-BFGS-B',
                bounds=BOUNDS,
                options={'maxiter': 100, 'ftol': 1e-6, 'gtol': 1e-5}
            )

            if result.success:
                bootstrap_params.append(result.x)
                convergence_flags.append(True)
            else:
                # Differential evolution with minimal settings is used here as a backup, since local optimisation did not succeed.
                result_de = differential_evolution(
                    lambda x: calculate_residuals(x, bootstrap_data),
                    bounds=BOUNDS,
                    maxiter=50,
                    popsize=10,
                    polish=False,
                    seed=42 + i
                )

                if result_de.success:
                    bootstrap_params.append(result_de.x)
                    convergence_flags.append(False)

        except Exception:
            continue

    total_time = time.time() - start_time

    if len(bootstrap_params) == 0:
        print("   All bootstrap samples failed.")
        return None

    bootstrap_params = np.array(bootstrap_params)
    local_converged = sum(convergence_flags)

    param_stats = {}

    print(f"\nFast Bootstrap Results ({len(bootstrap_params)}/{n_bootstrap} successful):")
    print(f"   Total time: {total_time/60:.1f} min  "
          f"({total_time/len(bootstrap_params):.1f}s per sample)")
    print(f"   Local convergence: {local_converged}/{len(bootstrap_params)}")
    print("   Parameter      Best Fit      Mean +/- Std      95% CI           CV%")
    print("   " + "-" * 80)

    for i, param_name in enumerate(RATE_CONSTANTS):
        values = bootstrap_params[:, i]
        mean_val = np.mean(values)
        std_val = np.std(values)
        ci_lower = np.percentile(values, 2.5)
        ci_upper = np.percentile(values, 97.5)
        cv_percent = (std_val / mean_val) * 100 if mean_val > 0 else float('inf')  # coefficient of variation: spread relative to the mean, as a percentage

        param_stats[param_name] = {
            'best_fit': best_fit_params[i],
            'mean': mean_val,
            'std': std_val,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'cv_percent': cv_percent,
            'values': values
        }

        print(f"   {param_name:12} {best_fit_params[i]:11.2e} "
              f"{mean_val:11.2e}+/-{std_val:.1e} "
              f"[{ci_lower:.1e}, {ci_upper:.1e}] {cv_percent:6.1f}%")

    # Parameters with high or moderate uncertainty are flagged here.
    print("\nParameter Confidence Analysis:")
    well_determined = []
    poorly_determined = []

    for param_name, stats in param_stats.items():
        if stats['cv_percent'] < 20:
            well_determined.append(param_name)
        elif stats['cv_percent'] < 100:
            print(f"   Moderate uncertainty: {param_name} ({stats['cv_percent']:.1f}% CV)")
        else:
            poorly_determined.append(param_name)
            print(f"   High uncertainty:     {param_name} ({stats['cv_percent']:.1f}% CV)")

    if well_determined:
        print(f"   Well-determined: {', '.join(well_determined)}")
    if poorly_determined:
        print("   Note: Parameters with >100% CV may need additional data or constraints")

    return param_stats


def analyse_parameter_correlations(param_stats):
    # Purpose: checks whether any pair of rate constants are so strongly correlated in the bootstrap fits that they cannot be independently determined from this data.
    # Inputs are: param_stats, the output from fast_bootstrap_from_best_fit() (or bootstrap_with_convergence_monitoring(), in Script 2).
    # Outputs are: the full correlation matrix between every pair of parameters (n_params x n_params). A strong correlation between two parameters (|r| > 0.7) means the data cannot independently constrain them - increasing one can be compensated by decreasing the other. This is a symptom of parameters being unidentifiable, and may suggest the reaction network is over-parameterised relative to the available data.
    # Called by: not called anywhere in the repo.

    print("\nParameter Correlation Analysis:")

    param_values = np.array([param_stats[param]['values'] for param in RATE_CONSTANTS])
    correlation_matrix = np.corrcoef(param_values)

    high_correlations = []
    n_params = len(RATE_CONSTANTS)
    for i in range(n_params):
        for j in range(i + 1, n_params):
            corr_val = correlation_matrix[i, j]
            if abs(corr_val) > 0.7:
                high_correlations.append((RATE_CONSTANTS[i], RATE_CONSTANTS[j], corr_val))

    if high_correlations:
        print("   Strong correlations found (|r| > 0.7):")
        for param1, param2, corr in high_correlations:
            direction = "positive" if corr > 0 else "negative"
            print(f"   {param1} <-> {param2}: r = {corr:+.3f} ({direction})")
    else:
        print("   No strong parameter correlations found (good identifiability)")

    return correlation_matrix


# =============================================================================
# SECTION 9: RESULTS VISUALISATION AND REPORTING
# =============================================================================

def plot_fit_results(best_fit_params, experimental_data, param_stats=None, save_plots=True):
    # Purpose: draws and saves the figure comparing the best-fit model against the experimental data for every observed species.
    # Inputs are: best_fit_params, the best-fit parameters to simulate; experimental_data, experimental data in the EXPERIMENTAL_DATA format; param_stats, optional bootstrap statistics from fast_bootstrap_from_best_fit() - if given, the bootstrap parameter distribution plots are also generated; save_plots, whether to save the plots to ./Outputs/ (default True).
    # Outputs are: nothing is returned; a 2x2 panel figure is drawn, comparing the model fit against the experimental data for all four observed species, each panel showing the experimental scatter points, the best-fit model line, and an R-squared annotation.
    # Called by: main_optimisation_only(), main_optimisation_with_bootstrap(), and directly from the entry point (module-level code, not inside a function) at the bottom of this file under choice == "3"; all in this file.

    print("\nGenerating fit comparison plots...")

    sim_results = run_zymonic_simulation(best_fit_params, end_time=1200, num_points=1000)

    if sim_results is None:
        print("   Could not generate simulation for plotting")
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Zymonic Acid Kinetics: Model Fit vs Experimental Data',
                 fontsize=14, fontweight='bold')

    species_list = ['ZCE', 'ZCD', 'ZK', 'ZOD']
    colours = ['#173F3F', '#F37B3A', '#6C9A8B', '#5C4033']

    for i, (species, colour) in enumerate(zip(species_list, colours)):
        row, col = i // 2, i % 2
        ax = axes[row, col]

        if species in experimental_data:
            exp_data = experimental_data[species]

            ax.scatter(exp_data['time'], exp_data['concentration'] * 1000,
                      color=colour, alpha=0.7, s=50,
                      label=f'{species} (Experimental)', zorder=3)

            ax.plot(sim_results['time'], sim_results[species] * 1000,
                   color=colour, linewidth=2, label=f'{species} (Model)', zorder=2)

            ax.set_xlabel('Time (min)')
            ax.set_ylabel('% total zymonic acid')
            ax.set_title(species)
            ax.legend()
            ax.grid(True, alpha=0.3)

            # R-squared is computed and annotated on the panel here.
            exp_conc = experimental_data[species]['concentration']
            sim_interp = np.interp(experimental_data[species]['time'],
                                  sim_results['time'], sim_results[species])
            ss_res = np.sum((exp_conc - sim_interp) ** 2)
            ss_tot = np.sum((exp_conc - np.mean(exp_conc)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

            ax.text(0.05, 0.95, f'R² = {r_squared:.3f}',
                   transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()

    if save_plots:
        os.makedirs('./Outputs', exist_ok=True)
        plt.savefig('./Outputs/zymonic_kinetics_fit.png', dpi=300, bbox_inches='tight')
        print("   Plot saved: ./Outputs/zymonic_kinetics_fit.png")

    plt.show()

    # If bootstrap results are available, the parameter distribution plots are also drawn here.
    if param_stats is not None:
        plot_parameter_distributions(param_stats, save_plots)


def plot_parameter_distributions(param_stats, save_plots=True):
    # Purpose: draws and saves a histogram of the bootstrap-fitted values for every rate constant, showing the best fit and its confidence interval.
    # Inputs are: param_stats, the output from fast_bootstrap_from_best_fit() (or bootstrap_with_convergence_monitoring(), in Script 2); save_plots, whether to save the figure to ./Outputs/bootstrap_distributions.png (default True).
    # Outputs are: nothing is returned; a histogram of the bootstrap-fitted values is drawn for each of the 10 rate constants, with the best-fit value from the primary optimisation marked as a red dashed line, and the 2.5th/97.5th percentiles (the 95% confidence interval bounds) marked as orange dotted lines. A log x-axis is used for any parameter whose 95% confidence interval spans more than two orders of magnitude, so the distribution can be seen clearly.
    # Called by: plot_fit_results(), in this file, when bootstrap statistics are available.

    fig, axes = plt.subplots(2, 5, figsize=(15, 8))
    fig.suptitle('Bootstrap Parameter Distributions', fontsize=14, fontweight='bold')

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


def save_results_to_file(best_fit_params, param_stats=None, optimisation_result=None):
    # Purpose: saves the best-fit rate constants (and bootstrap statistics, if available) to disk as JSON and CSV files.
    # Inputs are: best_fit_params, the best-fit rate constants in RATE_CONSTANTS order; param_stats, optional bootstrap statistics to include in the JSON output; optimisation_result, the optional scipy optimisation result object (used here for its residual and success flag).
    # Outputs are: nothing is returned; two files are written - zymonic_kinetics_results.json (the full results, including a timestamp and bootstrap statistics if given) and zymonic_kinetics_parameters.csv (a tabular summary for easy inspection).
    # Called by: main_optimisation_only(), main_optimisation_with_bootstrap(), and directly from the entry point (module-level code, not inside a function) at the bottom of this file under choice == "3"; all in this file.

    os.makedirs('./Outputs', exist_ok=True)

    results_dict = {
        'timestamp': datetime.now().isoformat(),
        'best_fit_parameters': dict(zip(RATE_CONSTANTS, best_fit_params.tolist())),
        'final_residual': float(optimisation_result.fun) if optimisation_result else None,
        'optimisation_success': bool(optimisation_result.success) if optimisation_result else None
    }

    if param_stats is not None:
        results_dict['bootstrap_statistics'] = {
            param_name: {
                'mean': float(stats['mean']),
                'std': float(stats['std']),
                'ci_lower': float(stats['ci_lower']),
                'ci_upper': float(stats['ci_upper'])
            }
            for param_name, stats in param_stats.items()
        }

    with open('./Outputs/zymonic_kinetics_results.json', 'w') as f:
        json.dump(results_dict, f, indent=2)

    results_df = pd.DataFrame({
        'Parameter': RATE_CONSTANTS,
        'Best_Fit_Value': best_fit_params,
        'Units': ['min^-1'] * len(RATE_CONSTANTS)
    })

    if param_stats is not None:
        results_df['Bootstrap_Mean'] = [param_stats[p]['mean'] for p in RATE_CONSTANTS]
        results_df['Bootstrap_Std']  = [param_stats[p]['std']  for p in RATE_CONSTANTS]
        results_df['CI_Lower']       = [param_stats[p]['ci_lower'] for p in RATE_CONSTANTS]
        results_df['CI_Upper']       = [param_stats[p]['ci_upper'] for p in RATE_CONSTANTS]

    results_df.to_csv('./Outputs/zymonic_kinetics_parameters.csv', index=False)

    print("\nResults saved:")
    print("   ./Outputs/zymonic_kinetics_results.json")
    print("   ./Outputs/zymonic_kinetics_parameters.csv")


# =============================================================================
# SECTION 10: MAIN EXECUTION FUNCTIONS
# =============================================================================
# These functions orchestrate the full analysis pipeline for each choice.

def main_optimisation_only():
    # Purpose: runs the recommended standard workflow: a single thorough optimisation, plotting, and saving the results.
    # Inputs are: none.
    # Outputs are: the best-fit rate constants, or None if every atempt fails. This runs a single, thorough differential evolution optimisation (the recommended starting point), using the settings in fit_rate_constants_differential_evolution(), followed by an automatic local polish step. If the primary optimisation fails to converge, a backup with more aggressive settings (maxiter=5000, popsize=30) is attempted automatically. 
    # Called by: the script's entry point, at the bottom of this file (if __name__ == "__main__"), under choice == "2".

    print("Running standard differential evolution optimisation...")

    create_zymonic_input_files()
    result = fit_rate_constants_differential_evolution(EXPERIMENTAL_DATA, max_iterations=2000)

    if result.success:
        print(f"\nOptimisation succeeded. Final residual: {result.fun:.4e}")
        plot_fit_results(result.x, EXPERIMENTAL_DATA, save_plots=True)
        save_results_to_file(result.x, optimisation_result=result)
        return result.x

    else:
        print("Primary optimisation did not converge. Trying backup settings...")
        result_backup = differential_evolution(
            objective_function,
            bounds=BOUNDS,
            maxiter=5000,
            popsize=30,
            mutation=(0.3, 1.7),
            recombination=0.95,
            atol=1e-12,
            tol=1e-10,
            seed=123,
            polish=True,
            disp=True
        )

        if result_backup.success:
            print(f"Backup optimisation succeeded. Residual: {result_backup.fun:.4e}")
            plot_fit_results(result_backup.x, EXPERIMENTAL_DATA, save_plots=True)
            save_results_to_file(result_backup.x, optimisation_result=result_backup)
            return result_backup.x
        else:
            print("All optimisation attempts failed.")
            return None


def main_optimisation_with_bootstrap():
    # Purpose: runs optimisation followed by a quick bootstrap uncertainty estimate, plotting and saving all the results.
    # Inputs are: none.
    # Outputs are: a tuple (best_fit_params, param_stats), or None if optimisation fails. Step 1 optimises the rate constants via differential evolution (500 iterations); step 2 estimates uncertainty via fast_bootstrap_from_best_fit() (50 samples); step 3 plots and saves all the results. For a more rigorous bootstrap with convergence monitoring (2245 samples), copy the best-fit parameters into Script 2 and run that script instead.
    # Called by: the script's entry point, at the bottom of this file (if __name__ == "__main__"), under choice == "4".

    print("Running optimisation with bootstrap uncertainty estimation...")

    create_zymonic_input_files()

    print("\n" + "="*50)
    print("STEP 1: PARAMETER OPTIMISATION")
    print("="*50)
    result = fit_rate_constants_differential_evolution(EXPERIMENTAL_DATA, max_iterations=500)

    if not result.success:
        print("Optimisation failed — cannot proceed to bootstrap")
        return None

    best_fit_params = result.x

    print("\n" + "="*50)
    print("STEP 2: BOOTSTRAP UNCERTAINTY ESTIMATION (50 samples)")
    print("="*50)
    param_stats = fast_bootstrap_from_best_fit(best_fit_params, n_bootstrap=50)

    print("\n" + "="*50)
    print("STEP 3: RESULTS AND VISUALISATION")
    print("="*50)
    plot_fit_results(best_fit_params, EXPERIMENTAL_DATA, param_stats, save_plots=True)
    save_results_to_file(best_fit_params, param_stats, result)

    return best_fit_params, param_stats


def test_simulation():
    # Purpose: runs a quick sanity-check simulation with arbitrary example parameters, to confirm the simulation framework is working before a lengthy optimisation run.
    # Inputs are: none.
    # Outputs are: True if the test simulation succeeds, False otherwise. This runs a quick test with hard-coded example parameters, to verify that the simulation framework is installed correctly, and that the input files, database creation, and ODE solver all work as expected, before committing to a lengthy optimisation run. The final concentration of all 5 species is printed, and a simple concentration-time plot is shown.
    # Called by: the script's entry point, at the bottom of this file (if __name__ == "__main__"), under choice == "1".

    print("Running test simulation with example parameters...")

    create_zymonic_input_files()

    # These test parameters are arbitrary and just for testing.
    test_params = [0.01, 0.005, 0.02, 0.001, 0.015, 0.008, 0.001, 0.0005, 0.003, 0.002]

    print("Test parameters:")
    for name, val in zip(RATE_CONSTANTS, test_params):
        print(f"   {name} = {val}")

    results = run_zymonic_simulation(test_params, end_time=1200, verbose=True)

    if results is not None:
        print("\nTest simulation successful. Final concentrations (mM):")
        for species in ['ZCE', 'ZCD', 'ZCK', 'ZOD', 'ZOK']:
            print(f"   {species}: {results[species][-1] * 1000:.3f} mM")

        # Every species with a non-negligible concentration is plotted here.
        plt.figure(figsize=(10, 6))
        for species in ['ZCE', 'ZCD', 'ZCK', 'ZOD', 'ZOK']:
            if np.max(results[species]) > 1e-10:
                plt.plot(results['time'], results[species] * 1000,
                        label=species, linewidth=2)

        plt.xlabel('Time (min)')
        plt.ylabel('Concentration (mM)')
        plt.title('Test Simulation')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

        return True

    else:
        print("Test simulation failed.")
        return False


# =============================================================================
# SECTION 11: ENTRY POINT
# =============================================================================

if __name__ == "__main__":

    print("Choose an analysis option:")
    print("   1 — Test simulation only (quick sanity check, ~1 minute)")
    print("   2 — Standard optimisation (single DE run, ~10-30 minutes)")
    print("   3 — Multi-start optimisation (best-of-5 DE runs, ~8-12 hours)")
    print("   4 — Optimisation + quick bootstrap (50 samples, ~2 hours total)")
    print()
    print("   NOTE: For a full 2245-sample bootstrap with convergence monitoring,")
    print("   copy the best-fit parameters from option 2 or 3 into Script 2")
    print("   (inverse_kinetics_bootstrap.py) and run that script separately.")
    print()

    # choice is set here, so the script runs non-interactively.
    choice = "2"

    if choice == "1":
        test_simulation()

    elif choice == "2":
        best_params = main_optimisation_only()
        if best_params is not None:
            print("\nStandard optimisation complete.")
            print("To proceed to full bootstrap analysis:")
            print("   1. Copy the best-fit parameters above into Script 2")
            print("   2. Run: python inverse_kinetics_bootstrap.py")

    elif choice == "3":
        print("Running multi-start optimisation (5 independent starts)...")
        create_zymonic_input_files()

        result = fit_with_multiple_starts(EXPERIMENTAL_DATA, n_starts=5)

        if result and result.success:
            plot_fit_results(result.x, EXPERIMENTAL_DATA, save_plots=True)
            save_results_to_file(result.x, optimisation_result=result)
            print(f"\nMulti-start optimisation complete.")
            print(f"   Final residual: {result.fun:.4e}")
            print("To proceed to full bootstrap analysis:")
            print("   1. Copy the best-fit parameters above into Script 2")
            print("   2. Run: python inverse_kinetics_bootstrap.py")
        else:
            print("Multi-start optimisation failed.")

    elif choice == "4":
        result = main_optimisation_with_bootstrap()
        if result is not None:
            print("\nOptimisation + bootstrap complete.")

    else:
        print("Invalid choice. Set choice to '1', '2', '3', or '4'.")
