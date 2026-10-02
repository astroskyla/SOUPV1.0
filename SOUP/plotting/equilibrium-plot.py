import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib import rcParams

# ======= Theme configuration =======
def set_theme(theme='light'):
    # Purpose: switches matplotlib's global style between light and dark.
    # Inputs are: theme, either 'light' or 'dark'.
    # Outputs are: nothing is returned; matplotlib's global plotting style (fonts, sizes, and colours) is updated to suit the chosen theme.
    # Called by: create_comparison_plots(), in this file.
    if theme == 'dark':
        rcParams.update({
            'font.family': 'Helvetica',
            'font.size': 16,
            'axes.labelsize': 18,
            'axes.titlesize': 20,
            'xtick.labelsize': 16,
            'ytick.labelsize': 16,
            'mathtext.fontset': 'custom',
            'mathtext.rm': 'Helvetica',
            'mathtext.it': 'Helvetica:italic',
            'mathtext.bf': 'Helvetica:bold',
            'figure.facecolor': '#000000',
            'axes.facecolor': '#000000',
            'axes.edgecolor': 'white',
            'axes.labelcolor': 'white',
            'text.color': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white'
        })
    else:  # light theme
        rcParams.update({
            'font.family': 'Helvetica',
            'font.size': 16,
            'axes.labelsize': 18,
            'axes.titlesize': 20,
            'xtick.labelsize': 16,
            'ytick.labelsize': 16,
            'mathtext.fontset': 'custom',
            'mathtext.rm': 'Helvetica',
            'mathtext.it': 'Helvetica:italic',
            'mathtext.bf': 'Helvetica:bold',
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
    # Outputs are: a dictionary of colours (GWB, SOUP, residual, edge, line) to use for the plot elements, suited to the chosen theme.
    # Called by: create_comparison_plots(), in this file.
    if theme == 'dark':
        return {
            'GWB': '#6c9a8b',       # Original teal
            'SOUP': '#f57b3b',      # Original orange
            'residual': '#183f3f',  # Dark teal for residuals
            'edge': 'white',
            'line': 'white'
        }
    else:  # light theme
        return {
            'GWB': '#6c9a8b',       # Original teal
            'SOUP': '#f57b3b',      # Original orange
            'residual': '#183f3f',  # Dark teal for residuals
            'edge': 'black',
            'line': 'black'
        }

# --- Input Data ---
# These are the equilibrium concentrations (M) being compared between GWB (Geochemist's Workbench, an established reference model) and SOUP, at three pH values, so the two models' predictions can be checked against each other. Each species maps to a (GWB value, SOUP value) pair.
data = {
    'pH7': {
        'Ca+2': (0.03941, 1.93E-02),
        'HPO4-2': (0.0237, 7.42E-03),
        'H+1': (1.22E-07, 1.00E-07),
        'CO2': (7.12E-08, 6.20E-08),
        'H2PO4-1': (0.01566, 1.19E-02),
        'PO4-3': (5.10E-07, 3.64E-08),
        'H2P2O7-2': (2.19E-06, 2.27E-07),
        'HP2O7-3': (1.92E-05, 1.38E-06),
        'H3PO4': (1.66E-07, 1.67E-07),
        'OH-1': (1.38E-07, 1.00E-07),
        'CaOH+1': (3.93E-08, 3.95E-08),
        'CaHCO3+1': (9.55E-08, 8.58E-08),
        'CaPO4-1': (2.00E-03, 2.02E-03),
        'CaHPO4': (5.86E-02, 7.86E-02),
    },
    'pH9': {
        'Ca+2': (0.01702, 7.02E-03),
        'CaCO3': (1.46E-05, 1.39E-05),
        'HCO3-1': (3.99E-05, 2.67E-05),
        'HPO4-2': (0.01689, 6.76E-03),
        'CO2': (7.12E-08, 6.27E-08),
        'CO3-2': (4.16E-06, 1.21E-06),
        'H2PO4-1': (0.000114, 1.08E-04),
        'PO4-3': (3.51E-05, 3.32E-06),
        'HP2O7-3': (9.68E-08, 1.15E-08),
        'OH-1': (1.38E-05, 1e-5),
        'CaOH+1': (1.72E-06, 1.44E-06),
        'CaHCO3+1': (4.19E-06, 3.12E-06),
        'CaPO4-1': (0.06406, 6.71E-02),
        'CaHPO4': (0.0189, 2.61E-02),
    },
    'pH5': {
        'Ca+2': (0.09365, 7.89E-02),
        'HPO4-2': (0.001718, 4.88E-04),
        'H+1': (1.24E-05, 1.00E-05),
        'CO2': (7.12E-08, 6.20E-08),
        'H2PO4-1': (0.09164, 7.82E-02),
        'H2P2O7-2': (8.79E-05, 9.81E-06),
        'HP2O7-3': (1.10E-05, 5.97E-07),
        'H3PO4': (9.16E-05, 1.10E-04),
        'OH-1': (1.00E-09, 1.00E-09),
        'CaPO4-1': (2.30E-06, 5.43E-06),
        'CaHPO4': (6.35E-03, 2.11E-02),
    }
}

# --- Helper to Format Chemical Labels ---
def format_species_label(species):
    # Purpose: turns a plain species name into a matplotlib maths label with proper subscripts/superscripts.
    # Inputs are: species, a species name in display format, e.g. 'HPO4-2'.
    # Outputs are: a matplotlib-ready label string with the number and charge rendered as subscript/superscript, e.g. "$HPO_{4}^{2-}$", or the species name unchanged if it doesn't match the expected pattern.
    # Called by: create_comparison_plots(), in this file, when labelling the x-axis of each plot.
    pattern = r"^([A-Za-z0-9]+)([\+\-]\d*)?$"  # splits the name into the formula part and an optional trailing charge, e.g. "HPO4" + "-2"
    match = re.match(pattern, species)
    if not match:
        return species
    formula = re.sub(r'(\d+)', r'_{\1}', match.group(1))  # wraps every digit run in "_{...}" so matplotlib renders it as a subscript
    charge = match.group(2)
    if charge:
        sign = charge[0]
        number = charge[1:] if len(charge) > 1 else ''
        if number == '' or number == '1':
            charge_sup = f"^{{{sign}}}"  # a single +/- charge, e.g. "^{-}"
        else:
            charge_sup = f"^{{{number}{sign}}}"  # a multi-charge, e.g. "^{2-}"
    else:
        charge_sup = ''
    return f"${formula}{charge_sup}$"

# ======= Main plotting function =======
def create_comparison_plots(theme='light'):
    # Purpose: draws and saves the GWB-vs-SOUP concentration and residual comparison figures.
    # Inputs are: theme, either 'light' or 'dark', which selects the plot's colour scheme.
    # Outputs are: nothing is returned; two figures are drawn (a bar chart of concentrations, and a bar chart of the log10 ratio between GWB and SOUP) and saved as both .png and .pdf files.
    # Called by: the script's entry point, at the bottom of this file (if __name__ == "__main__"), once for each theme.

    set_theme(theme)
    colour_scheme = get_colour_scheme(theme)

    # pH values are plotted in this fixed order, rather than dictionary order.
    ph_order = ['pH5', 'pH7', 'pH9']

    # --- Plot 1: Bar Plots of Concentrations ---
    fig, axs = plt.subplots(3, 1, figsize=(12, 10), sharey=True)
    width = 0.35
    colours = {'GWB': colour_scheme['GWB'], 'SOUP': colour_scheme['SOUP']}

    for ax, pH in zip(axs, ph_order):
        species_data = data[pH]
        species = list(species_data.keys())
        GWB_vals = [max(species_data[s][0], 1e-10) for s in species]
        SOUP_vals = [max(species_data[s][1], 1e-10) for s in species]
        x = np.arange(len(species))

        ax.bar(x - width/2, GWB_vals, width, label='GWB', color=colours['GWB'],
               edgecolor=colour_scheme['edge'], linewidth=0.8)
        ax.bar(x + width/2, SOUP_vals, width, label='SOUP', color=colours['SOUP'],
               edgecolor=colour_scheme['edge'], linewidth=0.8)

        ax.set_xticks(x)
        ax.set_xticklabels([format_species_label(s) for s in species],
                           rotation=45, ha='right', fontsize=16)
        ax.set_yscale('log')
        ax.set_ylabel('ln(concentration / M)', fontsize=18)
        ax.set_title(f"pH {pH[2:]}", fontsize=20)

    axs[-1].set_xlabel('Species', fontsize=18)
    plt.tight_layout()

    suffix = '_dark' if theme == 'dark' else ''
    bg_colour = 'black' if theme == 'dark' else 'white'

    fig.savefig(f'bar_plot_concentrations{suffix}.png', dpi=400,
                bbox_inches='tight', facecolor=bg_colour, edgecolor='none')
    fig.savefig(f'bar_plot_concentrations{suffix}.pdf',
                bbox_inches='tight', facecolor=bg_colour, edgecolor='none')

    print(f"✅ {theme.capitalize()} theme bar plots saved:")
    print(f"   bar_plot_concentrations{suffix}.png")
    print(f"   bar_plot_concentrations{suffix}.pdf")

    plt.close()

    # --- Plot 2: Residuals (log10 Ratio) ---
    fig, axs = plt.subplots(3, 1, figsize=(12, 9), sharey=True)

    for ax, pH in zip(axs, ph_order):
        species_data = data[pH]
        species = list(species_data.keys())

        GWB_vals = np.array([max(species_data[s][0], 1e-10) for s in species])
        SOUP_vals = np.array([max(species_data[s][1], 1e-10) for s in species])

        residuals = np.log10(GWB_vals / SOUP_vals)
        x = np.arange(len(species))

        ax.bar(x, residuals, color=colour_scheme['residual'],
               edgecolor=colour_scheme['edge'], linewidth=0.8)
        ax.axhline(0, color=colour_scheme['line'], linestyle='--', linewidth=1)

        ax.set_xticks(x)
        ax.set_xticklabels([format_species_label(s) for s in species],
                           rotation=45, ha='right', fontsize=16)
        ax.set_ylabel('log₁₀(GWB / SOUP)', fontsize=18)
        ax.set_title(f"pH {pH[2:]}", fontsize=20)

    axs[-1].set_xlabel('Species', fontsize=18)
    plt.tight_layout()

    fig.savefig(f'residuals_log_ratio{suffix}.png', dpi=400,
                bbox_inches='tight', facecolor=bg_colour, edgecolor='none')
    fig.savefig(f'residuals_log_ratio{suffix}.pdf',
                bbox_inches='tight', facecolor=bg_colour, edgecolor='none')

    print(f"{theme.capitalize()} theme residuals plots saved:")
    print(f"   residuals_log_ratio{suffix}.png")
    print(f"   residuals_log_ratio{suffix}.pdf")

    plt.close()

# ======= Main execution =======
if __name__ == "__main__":
    print("Creating light theme plots...")
    create_comparison_plots(theme='light')

    print("\nCreating dark theme plots...")
    create_comparison_plots(theme='dark')

    print("\nAll plots created successfully!")
