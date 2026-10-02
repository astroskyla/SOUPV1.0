# Plots the zymonic acid case study

import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib import rcParams

# ======= Theme configuration =======
def set_theme(theme='light'):
    # Purpose: switches matplotlib's global style between light and dark.
    # Inputs are: theme, either 'light' or 'dark'.
    # Outputs are: nothing is returned; matplotlib's global plotting style (fonts, sizes, and colours) is updated to suit the chosen theme.
    # Called by: create_zymonic_plots(), in this file.
    if theme == 'dark':
        rcParams.update({
            'font.family': 'Helvetica',
            'font.size': 14,
            'figure.facecolor': '#000000',
            'axes.facecolor': '#000000',
            'axes.edgecolor': 'white',
            'axes.labelcolor': 'white',
            'text.color': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'legend.facecolor': '#1a1a1a',
            'legend.edgecolor': '#444444'
        })
    else:  # light theme
        rcParams.update({
            'font.family': 'Helvetica',
            'font.size': 14,
            'figure.facecolor': 'white',
            'axes.facecolor': 'white',
            'axes.edgecolor': 'black',
            'axes.labelcolor': 'black',
            'text.color': 'black',
            'xtick.color': 'black',
            'ytick.color': 'black'
        })

def get_colour_scheme(theme='light'):
    # Purpose: returns the colour palette to use for the chosen theme.
    # Inputs are: theme, either 'light' or 'dark'.
    # Outputs are: a dictionary of colours (model, scatter) to use for the plot elements, suited to the chosen theme.
    # Called by: create_zymonic_plots(), in this file.
    if theme == 'dark':
        return {
            'model': '#ff9966',      # Brighter orange for dark background
            'scatter': '#7fc4b0'     # Brighter teal for dark background
        }
    else:  # light theme
        return {
            'model': '#f57b3b',      # Original orange
            'scatter': '#6c9a8b'     # Original teal
        }

# ======= Main plotting function =======
def create_zymonic_plots(theme='light'):
    # Purpose: draws the model-vs-literature zymonic acid comparison figure and saves it.
    # Inputs are: theme, either 'light' or 'dark', which selects the plot's colour scheme.
    # Outputs are: nothing is returned; a 2x2 figure comparing the model's zymonic acid species (ZK, ZOD, ZCE, ZCD) against the literature-provided data is drawn and saved as both .png and .pdf files.
    # Called by: the script's entry point, at the bottom of this file (if __name__ == "__main__"), once for each theme.

    set_theme(theme)
    colour_scheme = get_colour_scheme(theme)

    # This is the literature-provided data (time in minutes, scattered concentration as a percentage), used as a reference to compare the model's own output against.
    data = {
        'ZCDtime': [0.0, 2.280617933, 6.096512183, 7.892976589, 13.79837554, 17.76078993, 21.76779742, 27.87068004, 31.94139194, 38.06975633, 42.15957955, 56.51218347, 72.9415512, 87.32600733, 118.1270903, 150.9985667, 181.8060201, 245.4785794, 307.0871158, 430.2787068, 551.3999044, 674.5341615, 857.2065616, 1039.834369, 1138.321389],
        'ZCDscat': [0.0, 27.17391304, 36.95652174, 45.65217391, 54.13043478, 58.91304348, 62.17391304, 63.91304348, 65.0, 65.86956522, 66.30434783, 66.52173913, 65.86956522, 65.0, 63.91304348, 62.17391304, 60.86956522, 58.04347826, 55.65217391, 51.73913043, 48.47826087, 46.52173913, 42.82608696, 40.65217391, 39.7826087],
        'ZCEtime': [0.0, 0.987418379, 3.325370282, 7.701863354, 14.1168976, 18.36598184, 24.6281255, 28.80713489, 35.00557414, 39.13999044, 43.27440675, 57.6779742, 74.10734193, 88.48542762, 152.1133939, 184.953018, 246.5297022, 308.1000159, 431.2279025, 552.2917662, 858.0028667, 1040.598821, 119.2737697, 675.394171, 1139.073101],
        'ZCEscat': [100.0, 71.30434783, 61.52173913, 52.17391304, 43.26086957, 38.26086957, 34.56521739, 31.95652174, 30.43478261, 29.34782609, 28.26086957, 26.73913043, 26.08695652, 25.43478261, 24.13043478, 23.47826087, 22.17391304, 21.08695652, 19.34782609, 18.04347826, 15.65217391, 14.56521739, 24.7826087, 17.17391304, 14.13043478],
        'ZKtime': [0.0, 1138.442427, 1040.012741, 857.5250836, 675.0820194, 554.1137124, 431.1387164, 308.1955725, 246.7335563, 183.2393693, 152.5211021, 119.7579232, 89.02691511, 74.7061634, 60.37267081, 43.9942666, 37.90412486, 25.60280299, 15.35913362, 5.121834687],
        'ZKscat': [0.0, 35.65217391, 34.56521739, 31.95652174, 27.82608696, 25.86956522, 22.39130435, 17.82608696, 15.2173913, 11.95652174, 10.2173913, 8.260869565, 6.956521739, 5.652173913, 4.782608696, 3.695652174, 1.52173913, 1.304347826, 0.869565217, 0.217391304],
        'ZODtime': [0.0, 1139.175028, 1040.732601, 858.1812391, 675.6426183, 552.5975474, 431.6037586, 310.6099697, 247.0457079, 183.4814461, 152.7249562, 119.9171843, 91.21197643, 74.80809046, 58.40420449, 48.15416468, 37.90412486, 25.60280299, 15.35913362, 5.121834687],
        'ZODscat': [0.0, 10.65217391, 10.0, 9.565217391, 8.695652174, 7.608695652, 6.52173913, 5.434782609, 4.565217391, 3.695652174, 3.260869565, 2.826086957, 2.391304348, 2.173913043, 1.956521739, 1.739130435, 1.52173913, 1.304347826, 0.869565217, 0.217391304]
    }

    df_zcd = pd.DataFrame({'time': data['ZCDtime'], 'concentration': data['ZCDscat']})
    df_zce = pd.DataFrame({'time': data['ZCEtime'], 'concentration': data['ZCEscat']})
    df_zk = pd.DataFrame({'time': data['ZKtime'], 'concentration': data['ZKscat']})
    df_zod = pd.DataFrame({'time': data['ZODtime'], 'concentration': data['ZODscat']})

    directory = '../Outputs'

    species_data = {'ZK': [], 'ZOD': [], 'ZCE': [], 'ZCD': []}

    # Every markdown results file in the output directory is read in and folded into species_data below.
    for filename in os.listdir(directory):
        if filename.endswith('.md'):
            file_path = os.path.join(directory, filename)

            df = pd.read_csv(file_path, sep="|", skipinitialspace=True, header=0, skiprows=[1], engine='python')

            # The first and last columns are empty, since the markdown table format has a leading and trailing "|" on every row.
            df = df.drop(df.columns[[0, -1]], axis=1)

            df.columns = df.columns.str.strip()

            # ZOK and ZCK (the open and closed ketone forms) are combined into a single ZK species, since they aren't distinguished in the literature data being compared against.
            df['ZK'] = (df['ZOK'] + df['ZCK'])

            df = df.drop(columns=['ZOK', 'ZCK'])

            # Concentrations are converted to a percentage of total zymonic acid (0.1 M), to match the units the literature data is given in.
            species_columns = df.columns[1:]  # Exclude the 'time' column
            df[species_columns] = df[species_columns].apply(lambda x: (x / 0.1) * 100)

            for species in species_data.keys():
                if species in df.columns:
                    species_data[species].append((df['time'], df[species]))

    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    for time, concentration in species_data['ZK']:
        axs[0, 0].plot(time, concentration, color=colour_scheme['model'], label='ZK')
    axs[0, 0].scatter(df_zk['time'], df_zk['concentration'], color=colour_scheme['scatter'],
                      s=40, label='Provided ZK', edgecolors='none')
    axs[0, 0].set_title('Zymonic Closed & Open Ketone', font='Helvetica')
    axs[0, 0].set_xlabel('time (min)')
    axs[0, 0].set_ylabel('% total zymonic acid')

    for time, concentration in species_data['ZOD']:
        axs[0, 1].plot(time, concentration, color=colour_scheme['model'], label='ZOD')
    axs[0, 1].scatter(df_zod['time'], df_zod['concentration'], color=colour_scheme['scatter'],
                      s=40, label='Provided ZOD', edgecolors='none')
    axs[0, 1].set_title('Zymonic Open Diol', font='Helvetica')
    axs[0, 1].set_xlabel('time (min)')
    axs[0, 1].set_ylabel('% total zymonic acid')

    for time, concentration in species_data['ZCE']:
        axs[1, 0].plot(time, concentration, color=colour_scheme['model'], label='ZCE')
    axs[1, 0].scatter(df_zce['time'], df_zce['concentration'], color=colour_scheme['scatter'],
                      s=40, label='Provided ZCE', edgecolors='none')
    axs[1, 0].set_title('Zymonic Closed Enol', font='Helvetica')
    axs[1, 0].set_xlabel('time (min)')
    axs[1, 0].set_ylabel('% total zymonic acid')

    for time, concentration in species_data['ZCD']:
        axs[1, 1].plot(time, concentration, color=colour_scheme['model'], label='ZCD')
    axs[1, 1].scatter(df_zcd['time'], df_zcd['concentration'], color=colour_scheme['scatter'],
                      s=20, label='Provided ZCD', edgecolors='none')
    axs[1, 1].set_title('Zymonic Closed Diol', font='Helvetica')
    axs[1, 1].set_xlabel('time (min)')
    axs[1, 1].set_ylabel('% total zymonic acid')

    plt.tight_layout()

    suffix = '_dark' if theme == 'dark' else ''
    bg_colour = 'black' if theme == 'dark' else 'white'

    plt.savefig(f'zymonic{suffix}.png', dpi=400, bbox_inches='tight',
                facecolor=bg_colour, edgecolor='none')
    plt.savefig(f'zymonic{suffix}.pdf', bbox_inches='tight',
                facecolor=bg_colour, edgecolor='none')

    print(f"{theme.capitalize()} theme plots saved:")
    print(f"   zymonic{suffix}.png")
    print(f"   zymonic{suffix}.pdf")

    plt.close()
# ======= Main execution =======
if __name__ == "__main__":
    print("Creating light theme plots...")
    create_zymonic_plots(theme='light')

    print("\nCreating dark theme plots...")
    create_zymonic_plots(theme='dark')

    print("\nAll plots created successfully!")
