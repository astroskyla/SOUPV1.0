import os
import matplotlib.pyplot as plt

def format_to_internal(species_name):
    # Purpose: converts a species name from the display format to the internal database format.
    # Inputs are: species_name, a species name in the user-friendly display format, e.g. 'H+1'.
    # Outputs are: the same name converted to the internal database format used by the reaction database, e.g. 'HPLUS1'.
    # Called by: plot_results(), in this file.
    return species_name.replace("+", "PLUS").replace("-", "MINUS")

def format_to_display(species_name):
    # Purpose: converts a species name from the internal database format back to the display format.
    # Inputs are: species_name, a species name in the internal database format, e.g. 'HPLUS1'.
    # Outputs are: the same name converted back to the user-friendly display format, e.g. 'H+1'.
    # Called by: plot_results(), in this file.
    return species_name.replace("PLUS", "+").replace("MINUS", "-")

def plot_results(time_points, results, database_species, selected_species):
    # Purpose: repeatedly plots the requested species' concentrations over time, letting the user save the plot and choose new species, until they say they're done.
    # Inputs are: time_points, the list of times the simulation was evaluated at; results, the concentration of every species at each of those times; database_species, the full list of species in the simulation; selected_species, the species the user wants plotted (or 'all').
    # Outputs are: nothing is returned; a plot is drawn on screen and, if the user chooses to, saved to the Outputs folder as a .png file.
    # Called by: main(), in main.py.
    while True:
        # 'all' is treated as a request to plot every species in the database.
        if 'all' in [species.strip().lower() for species in selected_species]:
            selected_species = list(database_species)

        for species in selected_species:
            species = species.strip()

            # The species name is converted to the internal database format so it can be matched against database_species, which is stored that way.
            internal_species = format_to_internal(species)

            try:
                i = list(map(format_to_internal, database_species)).index(internal_species)  # finds this species' position by converting every database name to the same format and matching against it
                species_concs = results[:, i]

                # The species name is converted back to the display format for the plot legend.
                formatted_species = format_to_display(species)

                plt.plot(time_points, species_concs, label=formatted_species)
            except ValueError:
                print(f"Species '{species}' not found in the database.")

        plt.xlabel('Time')
        plt.ylabel('Concentration')

        pltyscale = input("Do you want to plot the y-axis on a log scale? (yes/no): ").strip().lower()
        plt.yscale('log' if pltyscale == 'yes' else 'linear')
        plt.legend()

        save_plot = input("Do you want to save the plot? (yes/no): ").strip().lower()
        if save_plot == 'yes':
            file_name = input("Enter the file name to save the plot (without extension): ").strip()
            output_folder = 'Outputs'
            os.makedirs(output_folder, exist_ok=True)
            file_path = os.path.join(output_folder, f"{file_name}.png")
            plt.savefig(file_path)
            print(f"Plot saved as {file_path}")

        plt.show()

        another_plot = input("Do you want to plot another species? (yes/no): ").strip().lower()
        if another_plot != 'yes':
            break

        selected_species = input("Enter the species you want to plot (comma-separated): ").split(',')
