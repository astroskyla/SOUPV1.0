# This module generates rate laws (forward/reverse rate constants and stoichiometry) from a reaction database.
# The rate constants used depend on which chemistry type the database is for: equilibrium chemistry gets high rates for fast equilibration, Fenton chemistry gets high rates for fast radical chemistry, and prebiotic chemistry gets moderate rates suited to complex reaction networks.

import pandas as pd
import sqlite3
import os
import re
import numpy as np

# These species need a non-negativity blend (see regord_pow in differentials.py) even at reaction order 0.
# HCN is a special case: reaction 10 (AMS formation) is zero-order in HCN by experiment, but the reaction must not be allowed to consume HCN below zero. This species must also have a matching entry in REGORD_SPECIES_CSTAR in differentials.py.
ZERO_ORDER_REGORD_SPECIES = {'HCN'}


def format_species_name(species):
    # Purpose: turns a species name into a valid Python variable name.
    # Inputs are: species, a species name as written in the reaction database, e.g. 'H+1'.
    # Outputs are: the same name with '+' replaced by 'PLUS', '-' replaced by 'MINUS', and any spaces removed, e.g. 'HPLUS1' - this is needed because '+', '-' and spaces are not valid characters in a Python variable name, and species names are used as variable names in the generated rate-law expressions.
    # Called by: parse_reactants_products(), in this file. It is also imported by inverse_kinetics.py and seeded_inverse_kinetics.py, but is not actually called there - both files define their own local version instead.
    species = species.replace('+', 'PLUS').replace('-', 'MINUS')
    return species.replace(' ', '')


def parse_reactants_products(reactants_products, apply_regord=False):
    # Purpose: turns a reactants/products string from the database into a rate-law expression plus a list of species and coefficients.
    # Inputs are: reactants_products, a string listing a reaction's reactants or products as written in the database, e.g. "2 FeII + H2O2"; apply_regord, whether to wrap fractional or negative reaction orders (n < 1, n != 0) with regord_pow(species, n, 'species') instead of species**n, which keeps the rate law well-behaved (Lipschitz, i.e. its slope stays finite) as a concentration approaches zero. Integer exponents, and exponents of 1, are left as plain species**n.
    # Outputs are: rate_law, the multiplied-together rate-law expression as a string, e.g. "FeII*H2O2"; species_list, a list of (species_name, stoichiometric_coefficient) pairs.
    # Called by: generate_rate_laws_with_database_type(), in this file; and directly (module-level code, not inside a function) in inverse_kinetics.py and seeded_inverse_kinetics.py.
    if not isinstance(reactants_products, str):
        return "", []

    components = reactants_products.strip().split(' + ')
    rate_law = []
    species_list = []

    # This pattern handles an optional leading sign/coefficient, optional parentheses, and an optional exponent that may itself be negative, e.g. "SO3-2^-0.6".
    pattern = r'\(?\s*([+-]?[\d.]+)?\s*\)?\s*([A-Za-z0-9\+\-]+)(?:\^([+-]?[\d.]+))?'  # captures three groups: an optional leading coefficient, the species name, and an optional exponent after "^"

    for component in components:
        match = re.match(pattern, component.strip())
        if match:
            coefficient_str, species, exponent = match.groups()
            species = format_species_name(species)

            try:
                coefficient = float(coefficient_str) if coefficient_str else 1.0
            except:
                print(f"Warning: Could not convert coefficient '{coefficient_str}' for species '{species}'. Skipping.")
                continue

            if species != 'H2O':  # Water is left out of rate laws, since its concentration is treated as constant (activity ~ 1).
                if exponent is not None:
                    exp_val = float(exponent)
                    if apply_regord and exp_val < 1 and exp_val != 0:
                        # A fractional or negative reaction order uses the concentration-dependent blend instead of a plain power.
                        rate_law.append(f"regord_pow({species}, {exp_val}, '{species}')")
                    elif apply_regord and exp_val == 0 and species in ZERO_ORDER_REGORD_SPECIES:
                        # A zero-order species that is still consumed also uses the blend, so its rate falls to zero as the species itself runs out.
                        rate_law.append(f"regord_pow({species}, {exp_val}, '{species}')")
                    else:
                        rate_law.append(f"({species}**{exponent})")
                elif coefficient == 1.0:
                    rate_law.append(f"{species}")
                else:
                    rate_law.append(f"({species}**{coefficient})")
                # The stoichiometric coefficient is kept alongside the species name, for use later when building the ODE terms.
                species_list.append((species, coefficient))
        else:
            print(f"Regex failed to parse: '{component.strip()}'")

    return '*'.join(rate_law), species_list


def calculate_equilibrium_rate_constants(equilibrium_constant, database_type, initial_h2o2=None):
    # Purpose: turns a single equilibrium constant into a matching pair of forward/reverse rate constants, sized appropriately for the chemistry type.
    # Inputs are: equilibrium_constant, the K_eq value for a reaction; database_type, which chemistry this reaction belongs to ('equilibrium', 'fenton', 'prebiotic', or anything else, which falls back to a conservative default); initial_h2o2, the starting H2O2 concentration (M) for a Fenton run - optional, and only used by the 'fenton' branch below. When None (the default), the 'fenton' branch uses its original fixed stiffness floor, unchanged from before this parameter existed; when given a number, it switches to the concentration-dependent floor described in the original methods write-up for this model.
    # Outputs are: forward_rate_constant and reverse_rate_constant, a pair of rate constants whose ratio equals K_eq (kf/kr = K_eq, as required for thermodynamic consistency), chosen within a range of magnitudes suited to the chemistry type.
    # Called by: generate_rate_laws_with_database_type(), in this file.
    try:
        K_eq = float(equilibrium_constant)
    except (ValueError, TypeError):
        print(f"Warning: Invalid equilibrium constant: {equilibrium_constant}")
        return 0.001, 0.001

    if database_type == 'equilibrium' or database_type == 'sai':
        # Equilibrium chemistry uses high rates, so the system reaches equilibrium quickly relative to the timescales of interest.
        rate_limit = 5e8

        reverse_rate_constant = rate_limit  # Start with maximum
        forward_rate_constant = K_eq * reverse_rate_constant

        if forward_rate_constant > rate_limit:
            forward_rate_constant = rate_limit
            reverse_rate_constant = forward_rate_constant / K_eq

        forward_rate_constant = min(forward_rate_constant, rate_limit)
        reverse_rate_constant = min(reverse_rate_constant, rate_limit)

        return forward_rate_constant, reverse_rate_constant

    elif database_type == 'fenton':
        # Fenton chemistry also uses high rates, since radical reactions are fast, but keeps them within bounds to avoid making the system unnecessarily stiff.
        base_rate = 1e4      # characteristic magnitude for fast chemistry
        max_rate = 5e7       # upper bound

        if initial_h2o2 is None:
            # This is the original fixed stiffness floor and still what every caller gets unless it explicitly opts into the concentration-dependent floor below.
            min_rate = 2e2
        else:
            # Below 1 mM starting H2O2, the system is stiff enough that a higher minimum rate is needed to avoid numerical issues, but above 1 mM the system is less stiff and can tolerate a lower minimum rate. This concentration-dependent floor is described in the original methods write-up for this model.
            min_rate = 60.0 if initial_h2o2 < 1e-3 else 2e2

        # The forward rate is first picked from a fast but bounded range.
        forward_rate = np.clip(base_rate, min_rate, max_rate)

        # The reverse rate is then set so that kf/kr = K_eq (thermodynamic consistency).
        reverse_rate = forward_rate / K_eq

        # If that reverse rate falls outside the allowed bounds, both rates are rescaled together so the ratio is preserved.
        if reverse_rate < min_rate:
            reverse_rate = min_rate
            forward_rate = reverse_rate * K_eq
        elif reverse_rate > max_rate:
            reverse_rate = max_rate
            forward_rate = reverse_rate * K_eq

        # The forward rate is clipped again in case the rescaling above pushed it back outside the bounds. Only relevant for extreme K_eq values, but included for completeness.
        forward_rate = np.clip(forward_rate, min_rate, max_rate)

        return forward_rate, reverse_rate


    elif database_type in ('prebiotic'):
        # Prebiotic networks use moderate rates, since these systems are large and combining many fast rates would make the equations very stiff.
        base_rate = 0.1
        max_rate = 10.0
        min_rate = 1e-15  # This floor only guards against a literal K_eq of zero, and is deliberately far below 1e-6 so it never distorts a genuinely small reverse rate.

        if K_eq > 100:
            forward_rate = min(max_rate, base_rate * 10)
            reverse_rate = forward_rate / K_eq          # No floor is applied here: a small reverse rate is not, by itself, a sign of stiffness.
        elif K_eq < 0.01:
            reverse_rate = min(max_rate, base_rate * 10)
            forward_rate = reverse_rate * K_eq          # No floor is applied here either, for the same reason.
        else:
            scaling_factor = min(5, max(0.2, np.sqrt(abs(K_eq))))  # grows or shrinks the base rate with K_eq, but clamped to a 5x-0.2x range so it never swings too far
            forward_rate = base_rate * scaling_factor
            reverse_rate = max(min_rate, forward_rate / K_eq)

        forward_rate = np.clip(forward_rate, min_rate, max_rate)
        reverse_rate = np.clip(reverse_rate, min_rate, max_rate)

        return forward_rate, reverse_rate

    else:
        # Any database type not recognised above falls back to this conservative method.
        print(f"Warning: Unknown database type '{database_type}' - using conservative rates")
        base_rate = 1e3
        rate_limit = 1e8

        if K_eq > 1e8:
            forward_rate = min(rate_limit, base_rate * 10)
            reverse_rate = forward_rate / K_eq
        elif K_eq < 1e-8:
            reverse_rate = min(rate_limit, base_rate * 10)
            forward_rate = reverse_rate * K_eq
        else:
            forward_rate = base_rate * min(10, max(0.1, np.sqrt(K_eq)))
            reverse_rate = forward_rate / K_eq

        return min(forward_rate, rate_limit), min(reverse_rate, rate_limit)


def generate_rate_laws_with_database_type(database_type, database_name=None, initial_h2o2=None):
    # Purpose: reads a reaction database and turns it into the rate laws and species list the rest of the simulation needs.
    # Inputs are: database_type, which chemistry the database represents ('equilibrium', 'prebiotic', or 'fenton'); database_name, the database filename without its .db extension (if not given, the user is prompted for it); initial_h2o2, the starting H2O2 concentration (M), passed straight through to calculate_equilibrium_rate_constants() for 'fenton' databases - optional, and has no effect for any other database_type (see that function for what it does with it).
    # Outputs are: rate_laws, a list of dictionaries (one per reaction) each holding the forward/reverse rate-law expressions and the species consumed/produced; database_species, the list of every unique species found across all reactions; equilibrium_expressions, a list of the K_eq = products/reactants expressions built for reactions that have an equilibrium constant.
    # Called by: main() in main.py; run_full_diagnostic() in diagnose.py; run_one() in run_sai_miyakawa_networks.py; _do() and _eff_order() in tests/prebiotic_tests.py; run_fe_cn_simulation() in Photoaquation/photoaquation_checker_simplified.py, Photoaquation/photoaquation_checker_simplified-V2.py and Photoaquation/photoaquation_rate.py; run_simulation() in Photoaquation/fe-cn-inverse-checker.py and Photoaquation/fecn_inversekinetics.py; directly from module-level code (not inside a function) in prebiotic_scenario_runner.py; and, within this file, by generate_rate_laws_from_markdown() and the test block at the bottom of this file.
    folder_path = './Inputs/Databases'

    if database_name is None:
        file_name = input("Please enter the database name (without .db): ")
    else:
        file_name = database_name

    db_path = os.path.join(folder_path, file_name + '.db')
    updated_db_path = os.path.join(folder_path, file_name + '_updated.db')

    print(f"Processing database: {file_name}.db")
    print(f"Chemistry type: {database_type}")

    # The original database is read in first.
    try:
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query("SELECT * FROM reactions", conn)
        conn.close()
        print(f"Successfully read {len(df)} reactions from database")
    except Exception as e:
        print(f"Error reading database: {e}")
        raise

    # Column names are cleaned up (whitespace stripped, lower-cased), and any leftover id column is dropped.
    df.rename(columns=lambda x: x.strip().lower(), inplace=True)
    if 'id' in df.columns:
        df.drop(columns=['id'], inplace=True)

    # Column names are then standardised to a consistent set of labels used throughout the rest of this function.
    column_map = {
        'reactants': 'Reactants',
        'products': 'Products',
        'forward_rate_constant': 'Forward Rate Constant',
        'reverse_rate_constant': 'Reverse Rate Constant',
        'equilibrium_constant': 'Equilibrium Constant'
    }
    df.rename(columns=column_map, inplace=True)

    df['Reactants'] = df['Reactants'].fillna('')
    df['Products'] = df['Products'].fillna('')

    def is_invalid(value):
        # Purpose: checks whether a database cell should be treated as missing.
        # Inputs are: value, a cell value read from the database (a rate constant or equilibrium constant).
        # Outputs are: True if the value is missing or unusable (None, NaN, or one of the text placeholders 'none'/'nan'/''/'n/a'), False otherwise.
        # Called by: generate_rate_laws_with_database_type(), in this file (several places, both above and below where this is defined).
        if value is None:
            return True
        if isinstance(value, str):
            val_lower = value.strip().lower()
            return val_lower in ('none', 'nan', '', 'n/a')
        if isinstance(value, (int, float)):
            return np.isnan(value) if not isinstance(value, str) else False
        return True

    # Any reaction missing both a forward and a reverse rate constant, but with an equilibrium constant given instead, has its rate constants calculated from that equilibrium constant.
    equilibrium_processed = 0

    print(f"Calculating rate constants for {database_type} chemistry...")

    for index, row in df.iterrows():
        forward_val = row.get('Forward Rate Constant')
        reverse_val = row.get('Reverse Rate Constant')

        # Only process if both forward and reverse rates are missing/invalid
        if is_invalid(forward_val) and is_invalid(reverse_val):
            equilibrium_val = row.get('Equilibrium Constant')

            if not is_invalid(equilibrium_val):
                try:
                    equilibrium_constant = float(str(equilibrium_val).strip())

                    forward_rate_constant, reverse_rate_constant = calculate_equilibrium_rate_constants(
                        equilibrium_constant, database_type, initial_h2o2
                    )

                    df.at[index, 'Forward Rate Constant'] = forward_rate_constant
                    df.at[index, 'Reverse Rate Constant'] = reverse_rate_constant

                    equilibrium_processed += 1

                    # Only the first few conversions are printed, as a sanity check, to keep the output readable.
                    if equilibrium_processed <= 3:
                        print(f"   K_eq={equilibrium_constant:.2e} → k_f={forward_rate_constant:.2e}, k_r={reverse_rate_constant:.2e}")

                except Exception as e:
                    print(f"Error processing equilibrium constant at row {index}: {e}")
                    continue

    print(f"Processed {equilibrium_processed} equilibrium constants using {database_type} method")

    # The database, now with rate constants filled in, is saved to a temporary "_updated" copy, which is re-read below and deleted again at the end of this function.
    if os.path.exists(updated_db_path):
        os.remove(updated_db_path)
    try:
        conn_updated = sqlite3.connect(updated_db_path)
        df.to_sql("reactions", conn_updated, index=False)
        conn_updated.close()
        print(f"Updated database saved to: {updated_db_path}")
    except Exception as e:
        print(f"Error saving updated database: {e}")
        raise

    conn_updated = sqlite3.connect(updated_db_path)
    df = pd.read_sql_query("SELECT * FROM reactions", conn_updated)
    conn_updated.close()

    rate_laws = []
    equilibrium_expressions = []
    database_species = set()
    parsing_errors = 0

    print("Converting database to rate laws...")
    # For prebiotic networks, fractional/negative reaction orders are wrapped with regord_pow (see differentials.py) so the rate law stays Lipschitz (its slope stays finite) as a concentration approaches zero.
    apply_regord = (database_type == 'prebiotic')

    for index, row in df.iterrows():
        try:
            forward_rate_constant = str(row.get('Forward Rate Constant', 'none')).strip()
            reverse_rate_constant = str(row.get('Reverse Rate Constant', 'none')).strip()
            equilibrium_constant = str(row.get('Equilibrium Constant', '')).strip()

            if forward_rate_constant.lower() == 'none':
                forward_rate_constant = None
            if reverse_rate_constant.lower() == 'none':
                reverse_rate_constant = None

            reactants, reactant_species = parse_reactants_products(row.get('Reactants', ''), apply_regord=apply_regord)
            products, product_species = parse_reactants_products(row.get('Products', ''), apply_regord=apply_regord)

            forward_rate_law = None
            if forward_rate_constant and reactants:
                forward_rate_law = f"{forward_rate_constant}*{reactants}"

            reverse_rate_law = None
            if not is_invalid(reverse_rate_constant) and products:
                reverse_rate_law = f"{reverse_rate_constant}*{products}"

            # An equilibrium expression (K_eq = products/reactants) is recorded too, where an equilibrium constant is available, for reference/diagnostics rather than for the simulation itself.
            if equilibrium_constant and equilibrium_constant.lower() not in ('none', 'nan', ''):
                try:
                    if product_species and reactant_species:
                        product_expr = ' * '.join([f"{s}**{c}" for s, c in product_species])
                        reactant_expr = ' * '.join([f"{s}**{c}" for s, c in reactant_species])
                        equilibrium_expressions.append(f"{equilibrium_constant} = ({product_expr}) / ({reactant_expr})")
                except Exception as e:
                    print(f"Warning: Error creating equilibrium expression for row {index}: {e}")

            rate_laws.append({
                'forward': forward_rate_law,
                'reverse': reverse_rate_law,
                'consumed': reactant_species,
                'produced': product_species
            })

            for s, _ in reactant_species + product_species:
                database_species.add(s)

        except Exception as e:
            parsing_errors += 1
            print(f"Error parsing row {index}: {e}")
            continue

    print(f"Successfully processed {len(rate_laws)} reactions")
    if parsing_errors > 0:
        print(f"Warning: {parsing_errors} reactions had parsing errors")
    print(f"Found {len(database_species)} unique species")

    # The temporary "_updated" database file is no longer needed, now that its contents have been read into rate_laws.
    try:
        os.remove(updated_db_path)
        print(f"Temporary database cleaned up")
    except Exception as e:
        print(f"Warning: Could not delete temporary file: {e}")

    print(f"Rate laws generated successfully for {database_type} chemistry!")

    return rate_laws, list(database_species), equilibrium_expressions


def generate_rate_laws_from_markdown():
    # Purpose: a legacy fallback that always generates equilibrium-chemistry rate laws.
    # Inputs are: none.
    # Outputs are: the same as generate_rate_laws_with_database_type('equilibrium') - this function only exists to give the old equilibrium-only behaviour a name to fall back on.
    # Called by: not called anywhere in the repo.
    print("Using legacy rate law generation - defaulting to equilibrium chemistry")
    return generate_rate_laws_with_database_type('equilibrium')


def query_rate_laws_by_species(rate_laws, species):
    # Purpose: finds every reaction that produces, or consumes, a given species.
    # Inputs are: rate_laws, the list of rate-law dictionaries produced by generate_rate_laws_with_database_type(); species, the species name to search for.
    # Outputs are: producing_rate_laws and consuming_rate_laws, the reactions (from rate_laws) that produce, and that consume, the given species.
    # Called by: construct_differential_equations(), in differentials.py.
    producing_rate_laws = []
    consuming_rate_laws = []

    for rate_law in rate_laws:
        if any(species == s for s, _ in rate_law['produced']):
            producing_rate_laws.append(rate_law)

        if any(species == s for s, _ in rate_law['consumed']):
            consuming_rate_laws.append(rate_law)

    return producing_rate_laws, consuming_rate_laws


def analyse_rate_law_statistics(rate_laws, database_type):
    # Purpose: prints a summary of the rate constants generated, as a sanity check.
    # Inputs are: rate_laws, the list of rate-law dictionaries; database_type, the chemistry type, used only for the printed heading.
    # Outputs are: nothing is returned; summary statistics (range, median, and an estimated stiffness ratio) for the forward and reverse rate constants are printed.
    # Called by: the test block at the bottom of this file (if __name__ == "__main__").
    print(f"\nRate Law Statistics for {database_type} chemistry:")

    forward_rates = []
    reverse_rates = []

    for law in rate_laws:
        if law['forward']:
            try:
                rate_str = law['forward'].split('*')[0]
                rate = float(rate_str)
                forward_rates.append(rate)
            except:
                pass

        if law['reverse']:
            try:
                rate_str = law['reverse'].split('*')[0]
                rate = float(rate_str)
                reverse_rates.append(rate)
            except:
                pass

    if forward_rates:
        print(f"   Forward rates: {len(forward_rates)} reactions")
        print(f"     Range: {min(forward_rates):.2e} to {max(forward_rates):.2e}")
        print(f"     Median: {np.median(forward_rates):.2e}")

    if reverse_rates:
        print(f"   Reverse rates: {len(reverse_rates)} reactions")
        print(f"     Range: {min(reverse_rates):.2e} to {max(reverse_rates):.2e}")
        print(f"     Median: {np.median(reverse_rates):.2e}")

    if forward_rates and reverse_rates:
        all_rates = forward_rates + reverse_rates
        # The stiffness ratio (fastest rate divided by slowest rate) gives a rough sense of how numerically demanding this reaction network will be to integrate.
        stiffness_ratio = max(all_rates) / min(all_rates) if min(all_rates) > 0 else float('inf')
        print(f"   Estimated stiffness ratio: {stiffness_ratio:.1e}")


if __name__ == "__main__":
    # Running this file directly tests rate-law generation on its own, without running a full simulation.
    print("Rate Laws Module Test")
    print("=" * 40)

    database_name = input("Enter database name (without .db): ")
    database_type = input("Enter chemistry type (equilibrium/prebiotic/fenton): ").lower()

    if database_type not in ['equilibrium', 'prebiotic', 'fenton']:
        print("Invalid chemistry type - defaulting to equilibrium")
        database_type = 'equilibrium'

    try:
        rate_laws, database_species, equilibrium_expressions = generate_rate_laws_with_database_type(
            database_type, database_name
        )

        analyse_rate_law_statistics(rate_laws, database_type)

        print(f"\nTest completed successfully!")
        print(f"   Generated {len(rate_laws)} rate laws")
        print(f"   Found {len(database_species)} species")
        print(f"   Created {len(equilibrium_expressions)} equilibrium expressions")

    except Exception as e:
        print(f"Test failed: {e}")
        import traceback
        traceback.print_exc()
