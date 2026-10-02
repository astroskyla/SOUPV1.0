# This script holds utility functions that are used across multiple modules. This was written forvever ago and other things should probably be added here.

def initialise_concentrations(database_species, species_aq, concentrations, species_atm, partial_pressures, henrys_constants):
    # Purpose: works out the starting concentration of every species, combining the aqueous input file with any Henry's law top-up from the atmosphere.
    # Inputs are: database_species, the full list of species tracked by the simulation; species_aq and concentrations, the aqueous species names and starting concentrations read from the aqueous input file; species_atm, partial_pressures and henrys_constants, the atmospheric species and their Henry's law data, used to set the starting concentration of any species that also exists in the atmosphere.
    # Outputs are: concs, a list of starting concentrations (one per entry in database_species), in the same order as database_species.
    # Called by: main() in main.py; run_zymonic_simulation() in inverse_kinetics.py and seeded_inverse_kinetics.py; run_full_diagnostic() in diagnose.py; _do() in tests/prebiotic_tests.py; run_fe_cn_simulation() in Photoaquation/photoaquation_checker_simplified.py, Photoaquation/photoaquation_checker_simplified-V2.py and Photoaquation/photoaquation_rate.py; run_simulation() in Photoaquation/fe-cn-inverse-checker.py and Photoaquation/fecn_inversekinetics.py; and directly from the scenario loop (module-level code, not inside a function) in prebiotic_scenario_runner.py.
    concs = [0] * len(database_species)

    # Every species that appears in the aqueous input file gets its starting concentration from there.
    for i, species in enumerate(database_species):
        if species in species_aq:
            index = species_aq.index(species)
            concs[i] += concentrations[index]

    # A species that also exists in the atmosphere gets its concentration topped up to the level Henry's law predicts for that partial pressure, if that level is higher than what was just set above.
    for i, species in enumerate(database_species):
        if species in species_atm:
            index = species_atm.index(species)
            if henrys_constants[index] is not None:
                delta = partial_pressures[index] * float(henrys_constants[index]) - concs[i]
                concs[i] += delta

    return concs

def update_concentrations_from_files(database_species, concs, new_species, new_concentrations):
    # Purpose: merges a new set of species/concentrations into the existing lists, overwriting any species already present and appending any new ones.
    # Inputs are: database_species and concs, the current list of species and their concentrations; new_species and new_concentrations, a list of species names and concentrations to merge in.
    # Outputs are: concs, updated in place so that any species already in database_species has its concentration overwritten, and any new species is appended to both lists.
    # Called by: not called anywhere else in the repo (imported in main.py but currently unused).
    for species, concentration in zip(new_species, new_concentrations):
        if species not in database_species:
            database_species.append(species)
            concs.append(concentration)
        else:
            index = database_species.index(species)
            concs[index] = concentration

    return concs
