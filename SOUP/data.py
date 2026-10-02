# This script reads the input data files used by the model (atmosphere and aqueous starting conditions).
# It does not read the SOUP reaction database file itself, which is selected separately by the user.

def normalise_species_name(name):
    # Purpose: converts a species name's '+'/'-' characters into the 'PLUS'/'MINUS' form the reaction database uses.
    # Inputs are: name, a species name as written in an input file, e.g. 'H+1'.
    # Outputs are: the same name with '+' replaced by 'PLUS' and '-' replaced by 'MINUS', e.g. 'HPLUS1', since the reaction database stores species under these names rather than using the '+'/'-' characters directly.
    # Called by: read_dat_file() and read_aqueous_dat_file(), both in this file.
    return name.replace('+', 'PLUS').replace('-', 'MINUS')

def read_dat_file(file_path):
    # Purpose: reads an atmosphere .dat file and returns its species, partial pressures, and Henry's law constants.
    # Inputs are: file_path, the path to an atmosphere .dat file listing species, partial pressures and Henry's law constants.
    # Outputs are: species_atm, partial_pressures and henrys_constants, three lists (one entry per atmospheric species) built by reading the file line by line and stopping once a "Total pressure" line is reached.
    # Called by: main() in main.py; run_zymonic_simulation() in inverse_kinetics.py and seeded_inverse_kinetics.py; run_full_diagnostic() in diagnose.py; _do() in tests/prebiotic_tests.py; run_fe_cn_simulation() in Photoaquation/photoaquation_checker_simplified.py, Photoaquation/photoaquation_checker_simplified-V2.py and Photoaquation/photoaquation_rate.py; run_simulation() in Photoaquation/fe-cn-inverse-checker.py and Photoaquation/fecn_inversekinetics.py; and directly from the scenario loop (module-level code, not inside a function) in prebiotic_scenario_runner.py.
    with open(file_path, 'r') as file:
        lines = file.readlines()

    species_atm = []
    partial_pressures = []
    henrys_constants = []

    for line in lines:
        if line.startswith('#') or line.strip() == '':
            continue
        if line.startswith('species, partial_pressure, constant'):
            continue
        if line.startswith('Total pressure'):
            break
        parts = line.strip().split(', ')
        if len(parts) == 3:
            species_name = normalise_species_name(parts[0])
            species_atm.append(species_name)
            partial_pressures.append(float(parts[1]))
            henrys_constants.append(parts[2] if parts[2] != 'N/A' else None)

    return species_atm, partial_pressures, henrys_constants

def read_aqueous_dat_file(file_path):
    # Purpose: reads an aqueous .dat file and returns its species/concentrations plus the overall solution conditions.
    # Inputs are: file_path, the path to an aqueous .dat file listing species and starting concentrations, plus overall conditions (total volume, temperature, pH, and optionally a solar spectrum file (this file does not do anything yet...)).
    # Outputs are: species_aq and concentrations, one entry per aqueous species; and total_volume, total_temperature, pH and solar_spectrum_file, the overall conditions read from their respective lines (each is left as None if that line is not present in the file).
    # Called by: main() in main.py; run_zymonic_simulation() in inverse_kinetics.py and seeded_inverse_kinetics.py; run_full_diagnostic() in diagnose.py; _do() in tests/prebiotic_tests.py; run_fe_cn_simulation() in Photoaquation/photoaquation_checker_simplified.py, Photoaquation/photoaquation_checker_simplified-V2.py and Photoaquation/photoaquation_rate.py; run_simulation() in Photoaquation/fe-cn-inverse-checker.py and Photoaquation/fecn_inversekinetics.py; and directly from the scenario loop (module-level code, not inside a function) in prebiotic_scenario_runner.py.
    with open(file_path, 'r') as file:
        lines = file.readlines()

    species_aq = []
    concentrations = []
    total_volume = None
    total_temperature = None
    pH = None
    solar_spectrum_file = None

    for i, line in enumerate(lines):
        if line.startswith('#') or line.strip() == '':
            continue
        if line.startswith('species, concentration'):
            continue
        if line.startswith('Total Volume'):
            total_volume = float(lines[i + 1].strip().replace(',', ''))
        elif line.startswith('Total Temperature'):
            total_temperature = float(lines[i + 1].strip())
        elif line.startswith('pH'):
            pH = float(lines[i + 1].strip())
        elif line.startswith('Solar Spectrum File'):
            solar_spectrum_file = lines[i + 1].strip()
        else:
            parts = line.strip().split(', ')
            if len(parts) == 2:
                species_name = normalise_species_name(parts[0])
                species_aq.append(species_name)
                concentrations.append(float(parts[1]))

    return species_aq, concentrations, total_volume, total_temperature, pH, solar_spectrum_file
