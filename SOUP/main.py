# This program runs a chemical kinetics simulation. This is the OG script from which everything runs :)
# The chemistry type (equilibrium, Fenton, or prebiotic) is detected from the database name.
# A solver is then chosen to match the system's stiffness (how differently the reaction rates vary in speed): equilibrium uses a high-precision Radau method, Fenton uses a simpler BDF method, and prebiotic uses a "bulletproof" solver built for extreme stiffness.

import os
import sys
import numpy as np
import re
import time

import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# These modules are always required.
from rate_laws import generate_rate_laws_with_database_type
from differentials import construct_differential_equations
from data import read_dat_file, read_aqueous_dat_file
from utilities import initialise_concentrations, update_concentrations_from_files
from plotting import plot_results

# These modules add extra features (the main solvers, long-term integration, and reaction-rate logging); if one is missing, the program still runs, just without that feature.
try:
    from differentials import (
        ode_system, update_concentrations,
        create_bulletproof_ode_system, solve_with_bulletproof_method
    )
    CORE_SYSTEMS_AVAILABLE = True
    print("All core ODE systems available")
except ImportError as e:
    CORE_SYSTEMS_AVAILABLE = False
    print(f"Core systems import error: {e}")

try:
    from long_term_solver import solve_long_term_chemistry, create_high_performance_ode_system
    LONG_TERM_SOLVER_AVAILABLE = True
    print("Long-term solver available")
except ImportError:
    LONG_TERM_SOLVER_AVAILABLE = False
    print("Long-term solver not available")

try:
    from reaction_rates_logger import ReactionRatesLogger
    REACTION_LOGGING_AVAILABLE = True
    print("Reaction logging available")
except ImportError:
    REACTION_LOGGING_AVAILABLE = False
    print("Reaction logging not available")


def create_simple_fenton_ode_system(database_species, differential_equations,
                                   species_atm, partial_pressures, henrys_constants, pH, rate_laws=None):
    # Purpose: builds the ODE function used to run Fenton-chemistry simulations, with reaction-rate logging attached to it.
    # Inputs are: the list of species in the database, the differential equation for each species, the atmospheric species/partial pressures/Henry's constants used for gas exchange, the solution pH, and (optionally) the rate laws, used only to switch on reaction-rate logging.
    # Outputs are: a function (fenton_ode) that calculates the rate of change of every species at a given time, for use with scipy's solve_ivp.
    # Called by: main(), when the Fenton solver is selected.
    print("   Creating simplified Fenton ODE system...")

    # Reaction-rate logging is set up here, if the logging module is available and rate laws were supplied.
    reaction_logger = None
    integration_start_time = None
    total_time = None

    if rate_laws is not None and REACTION_LOGGING_AVAILABLE:
        try:
            reaction_logger = ReactionRatesLogger(
                database_species, differential_equations, rate_laws,
                output_file='./Outputs/reaction_rates.html'
            )
            print("   - Reaction rate logging enabled")
        except Exception as e:
            print(f"   - Reaction logging setup failed: {e}")
    else:
        print("   - Reaction rate logging disabled")

    def simple_update_concentrations_fenton(y, database_species, species_atm, partial_pressures, henrys_constants, pH):
        # Purpose: applies Henry's law, iron mass balance, and negative-concentration clean-up before each rate calculation.
        # Inputs are: y, the current concentration of every species; database_species, the list of species names; species_atm, partial_pressures and henrys_constants, the atmospheric data used for Henry's law; pH, the solution pH (not used directly here, kept only for a consistent function signature).
        # Outputs are: y, the same array with atmospheric species topped up via Henry's law, the FeIII concentration recalculated from iron mass balance, and any negative values reset to zero.
        # Called by: fenton_ode(), in this file.
        y = np.array(y, copy=True)

        # A species that also exists in the atmosphere is topped up to at least the concentration Henry's law predicts for it.
        for i, species in enumerate(database_species):
            if species in species_atm:
                index = species_atm.index(species)
                if henrys_constants[index] is not None:
                    try:
                        equilibrium_concentration = partial_pressures[index] * float(henrys_constants[index])
                        if y[i] < equilibrium_concentration:
                            y[i] = equilibrium_concentration
                    except (ValueError, TypeError):
                        pass

        # Iron(III) is not integrated directly; instead its concentration is recalculated here from mass balance, as the total iron added (0.0002 M) minus whatever has already become intermediates I1/I2 or iron(II).
        if 'FeIII' in database_species and 'FeII' in database_species and 'I1' in database_species and 'I2' in database_species:
            try:
                fe_iii_idx = database_species.index('FeIII')
                i1_conc = y[database_species.index('I1')]
                i2_conc = y[database_species.index('I2')]
                fe_ii_conc = y[database_species.index('FeII')]

                y[fe_iii_idx] = max(0.0, 0.0002 - (i1_conc + i2_conc) - fe_ii_conc)  # total iron (0.0002 M) minus what has already become I1, I2, or Fe(II)

            except (ValueError, IndexError):
                pass

        # A concentration cannot physically be negative, so any negative value here, however small, is reset to exactly zero; this deliberately allows a radical to reach true zero rather than a tiny negative number. Don't love this. Working on fixing it. If worried run mass balance checks.
        for i in range(len(y)):
            if y[i] < 0:
                if y[i] > -1e-15:  # Very small negative
                    y[i] = 0.0     # Allow true zero for radicals
                else:
                    y[i] = 0.0     # Larger negative, set to zero

        return y

    def safe_pow(x, p):
        # Purpose: raises a number to a power without crashing on a negative or zero base.
        # Inputs are: x, a number; p, the power to raise it to.
        # Outputs are: x raised to the power p, or 0.0 if x is zero, negative, or the calculation fails - this avoids errors from raising a negative or zero concentration to a fractional power.
        # Called by: called from inside the rate-law expressions themselves, once replace_powers_with_safe_pow() has rewritten every "**" in them as a safe_pow(...) call; those expressions are then run by fenton_ode() via eval(), with safe_pow supplied as part of eval()'s namespace.
        try:
            if x <= 0:
                return 0.0
            return x ** p
        except Exception:
            return 0.0

    def replace_powers_with_safe_pow(expr):
        # Purpose: rewrites every "**" power operation in a rate-law expression to go through safe_pow() instead.
        # Inputs are: expr, a rate-law expression as a string, e.g. containing a term like "FeII**2".
        # Outputs are: the same expression with every "**" power operation rewritten as a call to safe_pow(...), e.g. "safe_pow(FeII, 2)".
        # Called by: fenton_ode(), in this file.
        pattern = re.compile(r'([a-zA-Z_][a-zA-Z0-9_]*)\s*\*\*\s*([-+]?[0-9]*\.?[0-9]+)')  # matches a species name followed by **exponent, e.g. "FeII**2"
        def replacer(match):
            # Purpose: builds the replacement text for one matched "**" power operation.
            # Inputs are: match, a regex match object for one "base**exponent" occurrence found by the pattern above.
            # Outputs are: the matched term rewritten as a safe_pow(base, exponent) call.
            # Called by: pattern.sub(), just below, once per match found in expr.
            base = match.group(1)
            exponent = match.group(2)
            return f"safe_pow({base}, {exponent})"
        return pattern.sub(replacer, expr)

    def create_eval_context(database_species, y):
        # Purpose: builds the variable namespace a rate-law expression needs to be evaluated.
        # Inputs are: database_species, the list of species names; y, their current concentrations.
        # Outputs are: a dictionary mapping each species name to its concentration in y, used as the variable namespace when a rate-law expression is evaluated.
        # Called by: fenton_ode(), in this file.
        return {species: y[i] for i, species in enumerate(database_species)}

    def fenton_ode(t, y):
        # Purpose: this is the ODE function itself and it calculates every species' rate of change at time t, for the solver to call repeatedly.
        # Inputs are: t, the current simulation time; y, the current concentration of every species.
        # Outputs are: dy_dt, the rate of change of every species at this time and concentration, calculated from the differential equations built earlier by construct_differential_equations().
        # Called by: scipy's solve_ivp(), in main(), once create_simple_fenton_ode_system() has returned this function.
        nonlocal integration_start_time, reaction_logger, total_time

        # The start time is recorded the first time this runs, so reaction-rate logging can report elapsed time rather than absolute time.
        if integration_start_time is None:
            integration_start_time = t

        # Concentrations are adjusted (Henry's law, iron mass balance, negative clean-up) before the rates are calculated.
        y_updated = simple_update_concentrations_fenton(y, database_species, species_atm,
                                                       partial_pressures, henrys_constants, pH)

        dy_dt = [0] * len(database_species)

        for i, species in enumerate(database_species):
            if species in differential_equations:
                right = differential_equations[species].split('=')[1]  # keeps only the right-hand side of "dspecies_dt = ..."
                safe_right = replace_powers_with_safe_pow(right)
                context = create_eval_context(database_species, y_updated)
                try:
                    dy_dt[i] = eval(safe_right, {"safe_pow": safe_pow}, context)
                except Exception as e:
                    print(f"Error in species '{species}' at t={t}: {right}")
                    raise e

        # If logging is switched on and the total run time is known, the reaction rates at this timepoint are recorded (only at a limited number of points, decided by reaction_logger).
        if reaction_logger and total_time:
            current_time = t - integration_start_time
            if reaction_logger.should_log_time(current_time, total_time):
                reaction_logger.log_timepoint(current_time, y_updated, dy_dt)

        return dy_dt

    # set_total_time and finalise_logging are attached to fenton_ode below, so the caller can control logging without needing direct access to reaction_logger.
    def set_total_time(time):
        # Purpose: tells the logger how long the whole run will be, so it can space out its logging points.
        # Inputs are: time, the total length of the integration in seconds, used so the logger knows how to space out its logging points.
        # Outputs are: nothing is returned; total_time is updated, and a message is printed if logging is enabled.
        # Called by: main(), via simple_ode_system.set_total_time(duration).
        nonlocal total_time
        total_time = time
        if reaction_logger:
            print(f"   - Fenton logging: {time/3600:.1f} hours, 50 timepoints")

    def finalise_logging():
        # Purpose: writes out the reaction-rate log, if logging was switched on.
        # Inputs are: none.
        # Outputs are: nothing is returned; if logging was enabled, the reaction-rate log is written out to an HTML file.
        # Called by: main(), via simple_ode_system.finalise_logging().
        if reaction_logger:
            reaction_logger.finalise_html()
            print("Fenton reaction rates analysis saved to: ./Outputs/reaction_rates.html")

    # Attach methods to function object
    fenton_ode.set_total_time = set_total_time
    fenton_ode.finalise_logging = finalise_logging

    return fenton_ode


def get_solver_config(database_name):
    # Purpose: picks the right solver settings for the detected chemistry type.
    # Inputs are: database_name, the name of the database file (without .db) chosen by the user.
    # Outputs are: a dictionary describing which solver to use (method, tolerances, and a short description), chosen by matching keywords in the dataase name.
    # Called by: main().
    database_name_lower = database_name.lower()

    if 'equilibrium' in database_name_lower:
        return {
            'method': 'Radau',
            'atol': 1e-8,
            'rtol': 1e-14,
            'solver_type': 'legacy',
            'description': 'High precision Radau for equilibrium studies'
        }
    elif 'fenton' in database_name_lower:
        return {
            'method': 'BDF',
            'atol': 1e-8,
            'rtol': 1e-10,
            'solver_type': 'fenton',
            'description': 'Simple BDF for Fenton chemistry'
        }
    else:
        # Any database name that doesn't match "equilibrium" or "fenton" is assumed to be prebiotic chemistry, and gets the universal robust ("bulletproof") solver.
        return {
            'method': 'bulletproof',
            'solver_type': 'bulletproof',
            'description': 'Universal robust system'
        }


def main():
    # Purpose: runs the whole simulation end to end - from asking the user for a database, through solving the chemistry, to saving and plotting the results.
    # Inputs are: none directly; the database name, end time, and other run parameters are collected from the user via input() prompts while this runs.
    # Outputs are: nothing is returned; the simulation is run end-to-end, results are saved to ./Outputs/species_concentrations.md, and a plot is shown and optionally saved.
    # Called by: the script's entry point, at the bottom of this file (if __name__ == "__main__").
    print("Chemical Kinetics Simulation")
    print("="*60)

    # The database name is used to work out which chemistry is being modelled.
    database_name = input("Please enter the database name (without .db): ")

    # Determine chemistry type from database name
    database_name_lower = database_name.lower()
    if 'equilibrium' in database_name_lower:
        database_type = 'equilibrium'
    elif 'fenton' in database_name_lower:
        database_type = 'fenton'
    else:
        database_type = 'prebiotic'

    print(f"Detected chemistry type: {database_type}")

    # The atmosphere and aqueous starting-condition files are read here, before the rate laws are generated, since a Fenton run's starting H2O2 concentration (read from the aqueous file) is needed below to choose the correct stiffness floor. A database-specific file is used if one exists, otherwise the generic default file is used instead.
    print("\nReading input files...")
    atmosphere_input = f'./Inputs/{database_name.lower()}-atmosphere.dat'
    aqueous_input = f'./Inputs/{database_name.lower()}-aqueous.dat'

    if not os.path.exists(atmosphere_input):
        atmosphere_input = './Inputs/atmosphere.dat'
    if not os.path.exists(aqueous_input):
        aqueous_input = './Inputs/aqueous.dat'

    species_atm, partial_pressures, henrys_constants = read_dat_file(atmosphere_input)
    species_aq, concentrations, total_volume, total_temperature, pH, solar_spectrum_file = read_aqueous_dat_file(aqueous_input)

    print("Atmospheric data:")
    print(f"  Species: {len(species_atm)}")
    print("Aqueous data:")
    print(f"  Species: {len(species_aq)}")
    print(f"  pH: {pH}")

    # Set this to True to use the concentration-dependent Fenton stiffness floor described in the original methods write-up (k_min = 60 s^-1 below 1 mM starting H2O2, 200 s^-1 at or above), which requires the starting H2O2 concentration to be read before the rate laws are generated. Set to False to use the original fixed floor (200 s^-1 regardless of concentration) exactly as before this option existed. I would like to rectify this in the future...
    USE_ADAPTIVE_FENTON_KMIN = True
    # For a Fenton run, the starting H2O2 concentration is picked out here (species 'HHOO' in the aqueous file), so generate_rate_laws_with_database_type() below can use the concentration-dependent stiffness floor. This stays None (giving the original, unchanged fixed-floor behaviour) for every other chemistry type, or if USE_ADAPTIVE_FENTON_KMIN is set to False above.
    initial_h2o2 = None
    if database_type == 'fenton' and USE_ADAPTIVE_FENTON_KMIN and 'HHOO' in species_aq:
        initial_h2o2 = concentrations[species_aq.index('HHOO')]
        print(f"  Starting H2O2: {initial_h2o2:.3e} M (adaptive Fenton stiffness floor enabled)")

    # The rate laws (forward/reverse rate constants and stoichiometry for every reaction) are generated from the reaction database.
    print(f"Generating rate laws for {database_type} chemistry...")
    rate_laws, database_species, equilibrium_expressions = generate_rate_laws_with_database_type(database_type, database_name, initial_h2o2)

    def format_rate_law(rate_law_string):
        # Purpose: makes a rate-law string readable for printing, by swapping the PLUS/MINUS placeholders back to +/-.
        # Inputs are: rate_law_string, a rate-law expression as a string, or None.
        # Outputs are: the same string with the placeholder text "PLUS" and "MINUS" swapped back to "+" and "-", or None if the input was None.
        # Called by: main(), when printing the sample rate laws below.
        if rate_law_string is None:
            return None
        return re.sub(r'MINUS', '-', re.sub(r'PLUS', '+', rate_law_string))

    def format_species(species_list):
        # Purpose: makes a list of species names readable for printing, by swapping the PLUS/MINUS placeholders back to +/-.
        # Inputs are: species_list, a list of (species_name, stoichiometric_coefficient) pairs.
        # Outputs are: the same list, with each species name's "PLUS"/"MINUS" placeholders swapped back to "+"/"-".
        # Called by: main(), when printing the sample rate laws below.
        return [(re.sub(r'MINUS', '-', re.sub(r'PLUS', '+', species[0])), species[1]) for species in species_list]

    print("\nSample rate laws:")
    for index, law in enumerate(rate_laws[2:5]):  # Only 3 reactions are shown here, to keep the output readable.
        forward_print = format_rate_law(law['forward'])
        reverse_print = format_rate_law(law['reverse'])
        consumed_print = format_species(law['consumed'])
        produced_print = format_species(law['produced'])

        print(f"{index+3}: Forward: {forward_print}")
        print(f"    Reverse: {reverse_print}")
        print(f"    {consumed_print} → {produced_print}")

    # The differential equation for each species is built next, combining every reaction that produces or consumes it.
    print("\n🔬 Constructing differential equations...")
    differential_equations = construct_differential_equations(rate_laws)
    print(differential_equations)

    print(f"Generated equations for {len(differential_equations)} species")

    # The starting concentration of every species is worked out from the input files (already read above).
    print("\nInitialising concentrations...")
    concs = initialise_concentrations(database_species, species_aq, concentrations, species_atm, partial_pressures, henrys_constants)

    start_time = 0
    end_time = float(input("\nEnter the end time (seconds): "))
    time_span = (start_time, end_time)
    initial_conditions = concs

    # A default number of data points is suggested (at least 500, or one per 10 seconds of simulated time), but the user can override it.
    default_points = max(500, int(end_time / 10))
    num_points_input = input(f"Enter number of data points (default {default_points}): ").strip()
    if num_points_input:
        try:
            num_points = int(num_points_input)
        except ValueError:
            print(f"Invalid input, using default: {default_points}")
            num_points = default_points
    else:
        num_points = default_points

    t_eval = np.linspace(start_time, end_time, num_points)

    # The solver (and its tolerance) is chosen based on the chemistry type detected from the database name.
    solver_config = get_solver_config(database_name)
    print(f"\nSolver configuration:")
    print(f"   {solver_config['description']}")

    # Every starting concentration is checked: a negative, NaN, or infinite value is not physically possible, so it is reset to zero. Don't love this.
    print("\nValidating initial conditions...")
    fixed_count = 0
    for i, (species, conc) in enumerate(zip(database_species, initial_conditions)):
        if conc < 0:
            print(f"Fixed negative concentration: {species.replace('PLUS', '+').replace('MINUS', '-')} = {conc:.2e} → 0.0")
            initial_conditions[i] = 0.0
            fixed_count += 1
        elif np.isnan(conc) or np.isinf(conc):
            print(f"Fixed invalid concentration: {species.replace('PLUS', '+').replace('MINUS', '-')} = {conc} → 0.0")
            initial_conditions[i] = 0.0
            fixed_count += 1

    if fixed_count > 0:
        print(f"Fixed {fixed_count} problematic initial concentrations")
    else:
        print("All initial concentrations are valid")

    start_clock = time.time()
    print(f"\n Starting integration...")
    print(f"   Chemistry: {database_type}")
    print(f"   Time span: {start_time} to {end_time} seconds ({(end_time-start_time)/3600:.2f} hours)")
    print(f"   Solver: {solver_config['method']} ({solver_config['solver_type']})")

    # The correct solver is dispatched below, depending on the chemistry type and how long the simulated run is.
    if not CORE_SYSTEMS_AVAILABLE:
        print("ERROR: Core ODE systems not available!")
        sys.exit(1)

    # A run longer than a week (604800 seconds) is treated as long-term, and is only handled specially for the bulletproof solver, since long-term chunking has only been built for that one.
    duration = end_time - start_time
    if duration > 604800 and LONG_TERM_SOLVER_AVAILABLE and solver_config['solver_type'] != 'bulletproof':  # More than 1 week, non-bulletproof only
        print(f"Long-term integration detected ({duration/86400:.1f} days)")

        if solver_config['solver_type'] == 'bulletproof':
            print("   Using high-performance universal robust system for long-term")
            bulletproof_ode_system = create_bulletproof_ode_system(
                database_species, differential_equations,
                species_atm, partial_pressures, henrys_constants, pH, rate_laws
            )
            bulletproof_ode_system.set_total_time(duration)

            hp_ode_system = create_high_performance_ode_system(bulletproof_ode_system)
            solution = solve_long_term_chemistry(hp_ode_system, time_span, initial_conditions)

            if hasattr(bulletproof_ode_system, 'finalise_logging'):
                bulletproof_ode_system.finalise_logging()
        else:
            print("   Long-term integration only available for universal robust systems")
            print("   Proceeding with standard integration")

    if duration <= 604800 or not LONG_TERM_SOLVER_AVAILABLE or solver_config['solver_type'] == 'bulletproof':
        # Standard integration based on solver type
        if solver_config['solver_type'] == 'legacy':
            # Equilibrium chemistry uses the Radau method, which has proven reliable for this chemistry type.
            print("Using legacy ode_system with Radau...")

            solution = solve_ivp(ode_system, time_span, initial_conditions,
                               method='Radau', atol=1e-8, rtol=1e-14,
                               args=(database_species, differential_equations,
                                    species_atm, partial_pressures, henrys_constants, pH))

        elif solver_config['solver_type'] == 'fenton':
            # Fenton chemistry uses the simplified BDF approach built above.
            print("Using simplified Fenton approach with BDF...")

            simple_ode_system = create_simple_fenton_ode_system(
                database_species, differential_equations,
                species_atm, partial_pressures, henrys_constants, pH, rate_laws
            )

            simple_ode_system.set_total_time(duration)

            solution = solve_ivp(
                simple_ode_system,
                time_span,
                initial_conditions,
                method='BDF',
                atol=1e-8,
                rtol=1e-10,
                t_eval=t_eval
            )

            simple_ode_system.finalise_logging()

        elif solver_config['solver_type'] == 'bulletproof':
            # The bulletproof solver handles any database not recognised as equilibrium or Fenton, and escalates to a stiffer method automatically if needed.
            print("Using universal robust system...")

            bulletproof_ode_system = create_bulletproof_ode_system(
                database_species, differential_equations,
                species_atm, partial_pressures, henrys_constants, pH, rate_laws
            )

            bulletproof_ode_system.set_total_time(duration)

            # The user is asked whether the system should behave as open (species lost to, and replenished from, the atmosphere via Henry's law) or closed (no atmospheric exchange).
            try:
                _ans = input("Open system (Henry's law atmospheric replenishment)? [Y/n]: ").strip().lower()
                _open = _ans not in ('n', 'no')
            except EOFError:
                _open = True
            if not _open:
                _n = len(database_species)
                bulletproof_ode_system.min_concentrations = {i: 0.0 for i in range(_n)}
                print("   Running as CLOSED system (no atmospheric replenishment)")
            else:
                print("   Running as OPEN system (Henry's law floors active)")

            try:
                solution = solve_with_bulletproof_method(bulletproof_ode_system, time_span, initial_conditions)
            except Exception as e:
                print(f"Bulletproof method failed: {e}")
                # If the main bulletproof method fails, the extreme-stiffness handler is tried as a last resort, since it is slower and is only needed occasionally.
                try:
                    from extreme_stiffness_solver import handle_extreme_stiffness_efficiently
                    print("   Escalating to extreme stiffness handler...")
                    solution = handle_extreme_stiffness_efficiently(
                        bulletproof_ode_system, time_span, initial_conditions,
                        database_species, rate_laws
                    )
                except Exception as e2:
                    print(f"Extreme stiffness handler also failed: {e2}")
                    sys.exit(1)

            if hasattr(bulletproof_ode_system, 'finalise_logging'):
                bulletproof_ode_system.finalise_logging()

        else:
            print(f"Unknown solver type: {solver_config['solver_type']}")
            sys.exit(1)

    if not solution.success:
        print(f"Integration failed: {solution.message}")
        sys.exit(1)

    end_clock = time.time()
    elapsed_time = end_clock - start_clock
    print(f"Integration completed in {elapsed_time:.2f} seconds")

    # If the solver produced a continuous (dense) solution, that is used to evaluate results at evenly spaced points; otherwise, the solver's own time points are used as they are.
    if hasattr(solution, 'sol') and solution.sol is not None:
        t_eval_final = np.linspace(time_span[0], time_span[1], num_points)
        results = solution.sol(t_eval_final).T
        time_points = t_eval_final
        print(f"Using dense output with {len(time_points)} points")
    else:
        results = solution.y.T
        time_points = solution.t
        print(f"Using solver time points: {len(time_points)} points")

    # NaN/Inf values are checked for here; negative concentrations are no longer clamped at this stage, since that is now handled inside the solver itself via event-based restarts.
    if np.any(np.isnan(results)) or np.any(np.isinf(results)):
        print("WARNING: NaN or Inf values detected in results!")
        results = np.where(np.isnan(results) | np.isinf(results), 1e-15, results)
        print("Replaced with minimum concentrations")
    else:
        print("All concentrations are valid!")

    print(f"\nIntegration Summary:")
    print(f"   Chemistry type: {database_type}")
    print(f"   Integration time: {duration/3600:.2f} hours")
    print(f"   Solver used: {solver_config['method']} ({solver_config['solver_type']})")
    print(f"   Data points: {len(time_points)}")
    print(f"   Species tracked: {len(database_species)}")
    print(f"   Final concentrations: min={np.min(results[-1]):.2e}, max={np.max(results[-1]):.2e}")

    # Any species whose concentration changed by more than 10% over the run is reported as a significant change.
    significant_changes = []
    for i, species in enumerate(database_species):
        initial = initial_conditions[i]
        final = results[-1][i]
        if abs(final - initial) / max(initial, 1e-12) > 0.1:  # Changed by more than 10%
            significant_changes.append((species, initial, final))

    if significant_changes:
        print(f"Species with significant changes (>{len(significant_changes)} total):")
        for species, initial, final in significant_changes[:5]:  # Only the first 5 are shown, to keep the output readable.
            formatted_name = species.replace('PLUS', '+').replace('MINUS', '-')
            percent_change = ((final - initial) / max(initial, 1e-15)) * 100
            print(f"  {formatted_name}: {initial:.2e} → {final:.2e} ({percent_change:+.1f}%)")

    def format_species_name(name):
        # Purpose: converts a species name from the internal database format to the display format, for writing results and printing.
        # Inputs are: name, a species name in the internal database format, e.g. 'HPLUS1'.
        # Outputs are: the same name converted to the display format, e.g. 'H+1'.
        # Called by: main(), when writing the results table and printing the list of available species below.
        return name.replace("MINUS", "-").replace("PLUS", "+")

    output_file = './Outputs/species_concentrations.md'
    os.makedirs('./Outputs', exist_ok=True)

    # Results are written out as a Markdown table, one row per timepoint.
    with open(output_file, 'w') as f:
        formatted_headers = [format_species_name(s) for s in database_species]
        header = '| time | ' + ' | '.join(formatted_headers) + ' |\n'
        separator = '|---' * (len(database_species) + 1) + '|\n'
        f.write(header)
        f.write(separator)

        for t, concs in zip(time_points, results):
            line = f"| {t} | " + ' | '.join([f"{c:.6e}" for c in concs]) + ' |\n'
            f.write(line)

    print(f"Results saved to {output_file}")

    # Automated sanitty checks are run on the saved results, if the checks module is available.
    try:
        _checks_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tests')
        if _checks_dir not in sys.path:
            sys.path.insert(0, _checks_dir)
        from checks import run_all as _run_checks
        _run_checks(output_file)
    except Exception as _e:
        print(f"  [checks skipped: {_e}]")

    formatted_species = [species.replace("MINUS", "-").replace("PLUS", "+") for species in database_species]
    print("Available species:", formatted_species)

    selected_species = input("Enter species to plot (comma-separated) or 'all': ").split(',')
    plot_results(time_points, results, database_species, selected_species)

    print(f"\nSimulation completed successfully!")


if __name__ == "__main__":
    main()
