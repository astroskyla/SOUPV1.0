# This is a diagnostic tool for the chemical kinetics simulation.
# It is meant to be run before the main simulation, to check the reaction database, initial conditions, and solver settings for anything likely to cause trouble.

import numpy as np
import sqlite3
import pandas as pd
import os
import sys


def diagnose_rate_constants(rate_laws):
    # Purpose: checks the rate constants for anything unusualy fast or slow.
    # Inputs are: rate_laws, the list of rate-law dictionaries.
    # Outputs are: nothing is returned; the range and spread of the forward and reverse rate constants are printed, along with a warning for any that are extremely high (>1e8) or extremely low (<1e-10).
    # Called by: run_full_diagnostic(), in this file.
    print("=== RATE CONSTANT ANALYSIS ===")

    forward_rates = []
    reverse_rates = []

    for i, law in enumerate(rate_laws):
        if law['forward']:
            try:
                # The rate constant is the leading number before the first "*" in the expression.
                rate_str = law['forward'].split('*')[0].strip()
                rate = float(rate_str)
                forward_rates.append((i, rate, law['forward']))
            except:
                print(f"Warning: Could not parse forward rate: {law['forward']}")

        if law['reverse'] and law['reverse'] != 'None':
            try:
                rate_str = law['reverse'].split('*')[0].strip()
                rate = float(rate_str)
                reverse_rates.append((i, rate, law['reverse']))
            except:
                print(f"Warning: Could not parse reverse rate: {law['reverse']}")

    if forward_rates:
        rates = [r[1] for r in forward_rates]
        print(f"\nForward Rate Constants:")
        print(f"  Count: {len(rates)}")
        print(f"  Range: {min(rates):.2e} to {max(rates):.2e}")
        if len(rates) > 1:
            print(f"  Span: {max(rates)/min(rates):.1e} orders of magnitude")

        extreme_high = [(i, r, expr) for i, r, expr in forward_rates if r > 1e8]
        extreme_low = [(i, r, expr) for i, r, expr in forward_rates if r < 1e-10]

        if extreme_high:
            print(f"  WARNING: {len(extreme_high)} extremely high rates (>1e8):")
            for i, r, expr in extreme_high[:3]:  # Only the first 3 are shown, to keep the output readable.
                print(f"    Reaction {i}: {r:.2e}")

        if extreme_low:
            print(f"  WARNING: {len(extreme_low)} extremely low rates (<1e-10):")
            for i, r, expr in extreme_low[:3]:
                print(f"    Reaction {i}: {r:.2e}")

    if reverse_rates:
        rates = [r[1] for r in reverse_rates]
        print(f"\nReverse Rate Constants:")
        print(f"  Count: {len(rates)}")
        print(f"  Range: {min(rates):.2e} to {max(rates):.2e}")
        if len(rates) > 1:
            print(f"  Span: {max(rates)/min(rates):.1e} orders of magnitude")


def diagnose_initial_conditions(database_species, initial_conditions):
    # Purpose: checks the starting concentrations for anything zero, negative, or unusually high.
    # Inputs are: database_species, the list of species names; initial_conditions, their starting concentrations.
    # Outputs are: nothing is returned; the number of zero, negative, and unusually high (>1 M) starting concentrations is printed, along with the overall concentration range.
    # Called by: run_full_diagnostic(), in this file.
    print("\n=== INITIAL CONDITIONS ANALYSIS ===")

    print(f"Total species: {len(database_species)}")

    zero_species = []
    negative_species = []
    high_species = []

    for i, (species, conc) in enumerate(zip(database_species, initial_conditions)):
        if conc < 0:
            negative_species.append((species, conc))
        elif conc == 0:
            zero_species.append((species, conc))
        elif conc > 1.0:  # Very high concentration
            high_species.append((species, conc))

    print(f"Zero concentrations: {len(zero_species)}")
    print(f"Negative concentrations: {len(negative_species)}")
    if negative_species and len(negative_species) <= 5:
        for species, conc in negative_species:
            formatted_name = species.replace('PLUS', '+').replace('MINUS', '-')
            print(f"  {formatted_name}: {conc}")

    print(f"High concentrations (>1M): {len(high_species)}")
    if high_species:
        for species, conc in high_species[:5]:
            formatted_name = species.replace('PLUS', '+').replace('MINUS', '-')
            print(f"  {formatted_name}: {conc:.3f} M")

    nonzero_concs = [c for c in initial_conditions if c > 0]
    if nonzero_concs:
        min_conc = min(nonzero_concs)
        max_conc = max(nonzero_concs)
        print(f"Concentration range: {min_conc:.2e} to {max_conc:.2e} M")
        if min_conc > 0:
            print(f"Dynamic range: {max_conc/min_conc:.1e} orders of magnitude")


def diagnose_henry_law(species_atm, partial_pressures, henrys_constants):
    # Purpose: checks the Henry's law setup for any species with an unusually high equilibrium concentration.
    # Inputs are: species_atm, partial_pressures and henrys_constants, the atmospheric species and their Henry's law data.
    # Outputs are: nothing is returned; each species' Henry's law equilibrium concentration is printed, along with a warning for any that come out unusually high (>0.1 M), since a high one is likely to dominate the system.
    # Called by: run_full_diagnostic(), in this file.
    print("\n=== HENRY'S LAW ANALYSIS ===")

    print(f"Atmospheric species: {len(species_atm)}")

    valid_henry = 0
    high_equilibrium = []

    for i, (species, pp, hc) in enumerate(zip(species_atm, partial_pressures, henrys_constants)):
        if hc is not None and hc != 'N/A':
            try:
                equilibrium_conc = float(pp) * float(hc)
                print(f"  {species}: P={pp} atm, H={hc} M/atm → C_eq={equilibrium_conc:.2e} M")
                valid_henry += 1

                if equilibrium_conc > 0.1:  # Very high equilibrium concentration
                    high_equilibrium.append((species, equilibrium_conc))

            except (ValueError, TypeError):
                print(f"  {species}: Invalid Henry constant '{hc}'")
        else:
            print(f"  {species}: No Henry constant")

    print(f"Valid Henry's Law species: {valid_henry}")

    if high_equilibrium:
        print(f"WARNING: High equilibrium concentrations (may dominate system):")
        for species, conc in high_equilibrium:
            print(f"  {species}: {conc:.2e} M")


def diagnose_stiffness(rate_laws):
    # Purpose: estimates how numerically stiff the system is, and suggests a suitable solver.
    # Inputs are: rate_laws, the list of rate-law dictionaries.
    # Outputs are: nothing is returned; the stiffness ratio (fastest rate divided by slowest rate) is printed, along with a rough recommendation for which solver is likely to cope with it.
    # Called by: run_full_diagnostic(), in this file.
    print("\n=== STIFFNESS ANALYSIS ===")

    all_rates = []

    for law in rate_laws:
        for rate_expr in [law['forward'], law['reverse']]:
            if rate_expr and rate_expr != 'None':
                try:
                    rate_str = rate_expr.split('*')[0].strip()
                    rate = float(rate_str)
                    if rate > 0:
                        all_rates.append(rate)
                except:
                    pass

    if all_rates and len(all_rates) > 1:
        min_rate = min(all_rates)
        max_rate = max(all_rates)
        stiffness_ratio = max_rate / min_rate

        print(f"Rate constant range: {min_rate:.2e} to {max_rate:.2e}")
        print(f"Stiffness ratio: {stiffness_ratio:.1e}")

        if stiffness_ratio > 1e12:
            print("  EXTREMELY STIFF - expect solver difficulties")
            print("      Recommend: Use bulletproof system or extreme stiffness handler")
        elif stiffness_ratio > 1e8:
            print("  VERY STIFF - may need careful solver selection")
            print("      Recommend: Use BDF or Radau methods")
        elif stiffness_ratio > 1e6:
            print("  MODERATELY STIFF - should be manageable")
            print("      Recommend: Standard stiff solvers work well")
        else:
            print("  OW STIFFNESS - should integrate easily")
    else:
        print("  Could not determine stiffness (insufficient rate data)")


def diagnose_species_balance(rate_laws, database_species):
    # Purpose: flags species that are only ever produced, only ever consumed, or never involved in any reaction.
    # Inputs are: rate_laws, the list of rate-law dictionaries; database_species, the list of species names.
    # Outputs are: nothing is returned; species that are only ever produced (and will accumulate), only ever consumed (and will deplete), or that appear in no reaction at all (inactive), are printed.
    # Called by: run_full_diagnostic(), in this file.
    print("\n=== SPECIES BALANCE ANALYSIS ===")

    species_production = {}
    species_consumption = {}

    for species in database_species:
        species_production[species] = []
        species_consumption[species] = []

    for i, law in enumerate(rate_laws):
        for species, coeff in law.get('consumed', []):
            if species in species_consumption:
                species_consumption[species].append(i)

        for species, coeff in law.get('produced', []):
            if species in species_production:
                species_production[species].append(i)

    only_produced = []
    only_consumed = []
    inactive = []

    for species in database_species:
        produced = len(species_production[species])
        consumed = len(species_consumption[species])

        if produced > 0 and consumed == 0:
            only_produced.append(species)
        elif consumed > 0 and produced == 0:
            only_consumed.append(species)
        elif produced == 0 and consumed == 0:
            inactive.append(species)

    print(f"Species only produced (will accumulate): {len(only_produced)}")
    if only_produced:
        for species in only_produced[:5]:  # Only the first 5 are shown, to keep the output readable.
            formatted_name = species.replace('PLUS', '+').replace('MINUS', '-')
            print(f"  {formatted_name}")

    print(f"Species only consumed (will deplete): {len(only_consumed)}")
    if only_consumed:
        for species in only_consumed[:5]:  # Only the first 5 are shown, to keep the output readable.
            formatted_name = species.replace('PLUS', '+').replace('MINUS', '-')
            print(f"  {formatted_name}")

    print(f"Inactive species (no reactions): {len(inactive)}")
    if inactive:
        for species in inactive[:5]:  # Only the first 5 are shown, to keep the output readable.
            formatted_name = species.replace('PLUS', '+').replace('MINUS', '-')
            print(f"  {formatted_name}")


def test_ode_evaluation(database_species, differential_equations, initial_conditions,
                       species_atm, partial_pressures, henrys_constants, pH):
    # Purpose: runs the ODE system once, as a smoke test, and checks the result for NaN/infinite/oversized derivatives.
    # Inputs are: database_species, differential_equations, initial_conditions, species_atm, partial_pressures, henrys_constants and pH, the same simulation setup used to run the actual simulation.
    # Outputs are: nothing is returned; the legacy ODE system is evaluated once at the starting concentrations, and the result is checked for NaN, infinite, or unusually large (>1e6) derivatives.
    # Called by: run_full_diagnostic(), in this file.
    print("\n=== ODE EVALUATION TEST ===")

    try:
        from differentials import ode_system

        print("Testing legacy ODE system evaluation at t=0...")
        derivatives = ode_system(0.0, initial_conditions, database_species, differential_equations,
                               species_atm, partial_pressures, henrys_constants, pH)

        nan_count = np.sum(np.isnan(derivatives))
        inf_count = np.sum(np.isinf(derivatives))
        large_count = np.sum(np.abs(derivatives) > 1e6)

        print(f"  NaN derivatives: {nan_count}")
        print(f"  Infinite derivatives: {inf_count}")
        print(f"  Very large derivatives (>1e6): {large_count}")

        if nan_count == 0 and inf_count == 0:
            print("  ✅ ODE evaluation successful")

            nonzero_derivs = [d for d in derivatives if abs(d) > 1e-15]
            if nonzero_derivs:
                print(f"  Active species: {len(nonzero_derivs)}/{len(derivatives)}")
                print(f"  Derivative range: {min(nonzero_derivs):.2e} to {max(nonzero_derivs):.2e}")
        else:
            print("  ODE evaluation has numerical issues")

            # Only the first 5 problem species are shown, to keep the output readable.
            problem_count = 0
            for i, (species, deriv) in enumerate(zip(database_species, derivatives)):
                if (np.isnan(deriv) or np.isinf(deriv)) and problem_count < 5:
                    formatted_name = species.replace('PLUS', '+').replace('MINUS', '-')
                    print(f"    Problem with {formatted_name}: derivative = {deriv}")
                    problem_count += 1

    except ImportError:
        print("  Could not import differentials module - skipping ODE test")
    except Exception as e:
        print(f"  ODE evaluation failed: {e}")


def diagnose_database_file(database_name):
    # Purpose: checks that the reaction database file exists and can be read.
    # Inputs are: database_name, the database filename without its .db extension.
    # Outputs are: True if the database file exists and could be read successfully (a short summary of it is printed), False otherwise.
    # Called by: run_full_diagnostic(), in this file.
    print("=== DATABASE FILE ANALYSIS ===")

    db_path = f'./Inputs/Databases/{database_name}.db'

    if not os.path.exists(db_path):
        print(f"Database file not found: {db_path}")
        return False

    try:
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query("SELECT * FROM reactions", conn)
        conn.close()

        print(f"Database file loaded successfully")
        print(f"   Total reactions: {len(df)}")
        print(f"   Columns: {list(df.columns)}")

        if 'Reactants' in df.columns and 'Products' in df.columns:
            empty_reactants = df['Reactants'].isna().sum()
            empty_products = df['Products'].isna().sum()
            print(f"   Empty reactants: {empty_reactants}")
            print(f"   Empty products: {empty_products}")

        return True

    except Exception as e:
        print(f"Database file error: {e}")
        return False


def run_full_diagnostic():
    # Purpose: runs every diagnostic check in this file, in turn, and prints a summary at the end.
    # Inputs are: none directly; the database name is collected from the user via input() while this runs.
    # Outputs are: nothing is returned; every diagnostic check in this file is run in turn, and a summary with general recommendations is printed at the end.
    # Called by: the script's entry point, at the bottom of this file (if __name__ == "__main__").
    print("CHEMICAL KINETICS DIAGNOSTIC TOOL")
    print("=" * 60)

    database_name = input("Enter database name (without .db): ")

    if not diagnose_database_file(database_name):
        print("Cannot proceed without valid database file.")
        return

    try:
        from rate_laws import generate_rate_laws_with_database_type
        from differentials import construct_differential_equations
        from data import read_dat_file, read_aqueous_dat_file
        from utilities import initialise_concentrations

        database_name_lower = database_name.lower()
        if 'equilibrium' in database_name_lower:
            database_type = 'equilibrium'
        elif 'prebiotic' in database_name_lower:
            database_type = 'prebiotic'
        elif 'fenton' in database_name_lower:
            database_type = 'fenton'
        else:
            database_type = 'equilibrium'

        print(f"\nDetected chemistry type: {database_type}")

        print("\nLoading rate laws...")
        rate_laws, database_species, equilibrium_expressions = generate_rate_laws_with_database_type(
            database_type, database_name
        )

        print("Loading input files...")
        atmosphere_input = f'./Inputs/{database_type}-atmosphere.dat'
        aqueous_input = f'./Inputs/{database_type}-aqueous.dat'

        if not os.path.exists(atmosphere_input):
            atmosphere_input = './Inputs/atmosphere.dat'
        if not os.path.exists(aqueous_input):
            aqueous_input = './Inputs/aqueous.dat'

        species_atm, partial_pressures, henrys_constants = read_dat_file(atmosphere_input)
        species_aq, concentrations, total_volume, total_temperature, pH, solar_spectrum_file = read_aqueous_dat_file(aqueous_input)

        initial_conditions = initialise_concentrations(
            database_species, species_aq, concentrations,
            species_atm, partial_pressures, henrys_constants
        )

        print("Constructing differential equations...")
        differential_equations = construct_differential_equations(rate_laws)

        print("\n" + "="*60)
        print("RUNNING DIAGNOSTIC TESTS")
        print("="*60)

        diagnose_rate_constants(rate_laws)
        diagnose_initial_conditions(database_species, initial_conditions)
        diagnose_henry_law(species_atm, partial_pressures, henrys_constants)
        diagnose_stiffness(rate_laws)
        diagnose_species_balance(rate_laws, database_species)
        test_ode_evaluation(database_species, differential_equations, initial_conditions,
                           species_atm, partial_pressures, henrys_constants, pH)

        print("\n" + "="*60)
        print("DIAGNOSTIC SUMMARY")
        print("="*60)
        print("Diagnostic completed successfully")
        print("\nGeneral Recommendations:")
        print("1. Check any species marked as 'only produced' or 'only consumed'")
        print("2. If stiffness ratio > 1e8, consider shorter integration times or robust solvers")
        print("3. Monitor Henry's Law species")
        print("4. Verify that high/low rate constants are physically reasonable")
        print("5. For Fenton systems: Simple BDF solver usually works well")
        print("6. For equilibrium systems: High-precision Radau solver recommended")
        print("7. For prebiotic systems: Use bulletproof solver with extreme stiffness handling")

        print(f"\nSystem appears ready for {database_type} chemistry simulation!")

    except ImportError as e:
        print(f"\nImport error: {e}")
        print("Please ensure all required modules are available")
    except Exception as e:
        print(f"\nDiagnostic failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    run_full_diagnostic()
