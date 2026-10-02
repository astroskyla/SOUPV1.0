# This plots the Fenton reaction model validation figure, comparing this model's results against the Gallard & De Laat model and experimental data.

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams
from matplotlib.lines import Line2D

def set_theme(theme='light'):
    # Purpose: switches matplotlib's global style between light and dark.
    # Inputs are: theme, either 'light' or 'dark'.
    # Outputs are: nothing is returned; matplotlib's global plotting style (fonts, sizes, and colours) is updated to suit the chosen theme.
    # Called by: create_fenton_plots(), in this file.
    if theme == 'dark':
        rcParams.update({
            'font.family': 'sans-serif',
            'font.sans-serif': ['Helvetica', 'Arial'],
            'font.size': 14,
            'axes.titlesize': 16,
            'axes.labelsize': 14,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'legend.fontsize': 14,
            'legend.title_fontsize': 14,
            'figure.dpi': 400,
            'axes.linewidth': 1.0,
            'axes.spines.top': False,
            'axes.spines.right': False,
            'xtick.major.size': 4,
            'ytick.major.size': 4,
            'axes.grid': False,
            # Dark theme specific
            'figure.facecolor': '#000000',
            'axes.facecolor': '#000000',
            'axes.edgecolor': 'white',
            'axes.labelcolor': 'white',
            'text.color': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'grid.color': '#333333',
            'legend.facecolor': '#1a1a1a',
            'legend.edgecolor': '#444444'
        })
    else:  # light theme
        rcParams.update({
            'font.family': 'sans-serif',
            'font.sans-serif': ['Helvetica', 'Arial'],
            'font.size': 14,
            'axes.titlesize': 16,
            'axes.labelsize': 14,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'legend.fontsize': 14,
            'legend.title_fontsize': 14,
            'figure.dpi': 400,
            'axes.linewidth': 1.0,
            'axes.spines.top': False,
            'axes.spines.right': False,
            'xtick.major.size': 4,
            'ytick.major.size': 4,
            'axes.grid': False,
            # Light theme specific
            'figure.facecolor': 'white',
            'axes.facecolor': 'white',
            'axes.edgecolor': 'black',
            'axes.labelcolor': 'black',
            'text.color': 'black',
            'xtick.color': 'black',
            'ytick.color': 'black'
        })

def get_colour_scheme(theme='light'):
    # Purpose: returns the colour palette (one per H2O2 concentration, plus supporting colours) for the chosen theme.
    # Inputs are: theme, either 'light' or 'dark'.
    # Outputs are: a dictionary of colours to use for the plot elements (one per H2O2 concentration level, plus experimental/label/legend colours), suited to the chosen theme.
    # Called by: create_fenton_plots(), in this file.
    if theme == 'dark':
        return {
            'concentrations': {
                '200uM': '#b8d4a8',   # Light green
                '1mM': '#7fc4b0',     # Light teal
                '10mM': '#4fa3a3',    # Bright teal
                '100mM': '#d4925e',   # Light brown/orange
                '1M': '#ff9966'       # Bright orange (highest conc)
            },
            'experimental': 'white',
            'label_bg': '#1a1a1a',
            'legend_face': '#1a1a1a',
            'legend_edge': '#444444'
        }
    else:  # light theme
        return {
            'concentrations': {
                '200uM': '#dde7d8',   # Light gray-green
                '1mM': '#6c9a8b',     # Medium teal
                '10mM': '#173f3f',    # Dark teal-green
                '100mM': '#5c4033',   # Dark brown
                '1M': '#f37b3a'       # Orange-red
            },
            'experimental': 'black',
            'label_bg': 'white',
            'legend_face': 'white',
            'legend_edge': 'lightgray'
        }

# Styling
WHITE_LINESTYLE = '-'            # Solid line for White model
GALLARD_LINESTYLE = ':'          # Dotted line for Gallard model

def load_data_file(file_path):
    # Purpose: reads a .dat or .md data file into a consistently-formatted DataFrame.
    # Inputs are: file_path, the path to a .dat or .md data file.
    # Outputs are: a pandas DataFrame with whitespace-stripped column names (and, for a .md file, its extra leading/trailing empty columns dropped, and its time column converted from seconds to minutes if it looks like it's in seconds), or None if the file couldn't be read or has an unsupported extension.
    # Called by: create_fenton_plots(), in this file.
    try:
        if file_path.endswith('.dat'):
            data = pd.read_csv(file_path, delimiter=',')
            data.columns = data.columns.str.strip()
        elif file_path.endswith('.md'):
            data = pd.read_csv(file_path, sep='|', skipinitialspace=True,
                             header=0, skiprows=[1], engine='python')
            data = data.drop(data.columns[[0, -1]], axis=1)
            data.columns = data.columns.str.strip()
            # A time column with a maximum value over 200 is assumed to be in seconds, and is converted to minutes for these plots.
            if 'time' in data.columns and data['time'].max() > 200:
                data['time'] = data['time'] / 60
        else:
            print(f"Unsupported file format: {file_path}")
            return None

        return data
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None

def create_fenton_plots(theme='light'):
    # Purpose: draws the 2x2 Fenton chemistry validation figure and saves it.
    # Inputs are: theme, either 'light' or 'dark', which selects the plot's colour scheme.
    # Outputs are: nothing is returned; a 2x2 figure comparing this model's Fenton chemistry results against the Gallard & De Laat model and experimental data (atrazine degradation, radical concentrations, Fe(II), and OH-) is drawn and saved as both .png and .pdf files.
    # Called by: the script's entry point, at the bottom of this file (if __name__ == "__main__"), once for each theme.

    set_theme(theme)
    colours = get_colour_scheme(theme)

    # A figure sized to suit an A4 document.
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 9))

    concentrations = ['200uM', '1mM', '10mM', '100mM', '1M']
    conc_labels = ['200 μM', '1 mM', '10 mM', '100 mM', '1 M']

    markers = ['o', 's', '^', 'D', 'v']

    # These track which data sources actually had data plotted, so the legend only lists sources that are really present.
    has_experimental = False
    has_gallard = False
    has_white = False

    # --- Plot 1: Atrazine Degradation ---
    for i, (conc, label) in enumerate(zip(concentrations, conc_labels)):
        conc_colour = colours['concentrations'][conc]

        exp_file = f"../Outputs/Digitised_Data/atrazine-exp-{conc}H2O2.dat"
        if os.path.exists(exp_file):
            exp_data = load_data_file(exp_file)
            if exp_data is not None and 'time' in exp_data.columns and 'atrazine' in exp_data.columns:
                edge_colour = 'black' if theme == 'dark' else 'white'
                ax1.scatter(exp_data['time'], exp_data['atrazine'],
                           color=colours['experimental'], marker=markers[i], s=30, alpha=0.8,
                           edgecolors=edge_colour, linewidths=0.5, zorder=3)
                has_experimental = True

        # The Gallard & De Laat model's own predicted values, for comparison.
        gallard_file = f"../Outputs/Digitised_Data/atrazine-{conc}H2O2.dat"
        if os.path.exists(gallard_file):
            gallard_data = load_data_file(gallard_file)
            if gallard_data is not None and 'time' in gallard_data.columns and 'atrazine' in gallard_data.columns:
                ax1.plot(gallard_data['time'], gallard_data['atrazine'],
                        color=conc_colour, linestyle=GALLARD_LINESTYLE, linewidth=2.0,
                        alpha=0.9)
                has_gallard = True

        # This model's own predicted values (labelled "White" after the model), normalised to the starting atrazine concentration.
        white_file = f"../Outputs/FentonData/species_concentrations_modelled_{conc}H2O2.md"
        if os.path.exists(white_file):
            white_data = load_data_file(white_file)
            if white_data is not None and 'time' in white_data.columns and 'atrazine' in white_data.columns:
                if white_data['atrazine'].iloc[0] > 0:
                    relative_atrazine = white_data['atrazine'] / white_data['atrazine'].iloc[0]
                    ax1.plot(white_data['time'], relative_atrazine,
                            color=conc_colour, linestyle=WHITE_LINESTYLE, linewidth=2.0,
                            alpha=0.9)
                    has_white = True

    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('[Atrazine]/[Atrazine]₀')
    ax1.set_title('Relative Concentration of Atrazine')
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, 1.05)

    # --- Plot 2: Radical Concentrations ---
    for conc in concentrations:
        conc_colour = colours['concentrations'][conc]

        radical_file = f"../Outputs/Digitised_Data/HO2+O2_{conc}H2O2.dat"
        if os.path.exists(radical_file):
            radical_data = load_data_file(radical_file)
            if radical_data is not None and 'time' in radical_data.columns:
                radical_col = None
                if 'HOOdot + OOdot' in radical_data.columns:
                    radical_col = 'HOOdot + OOdot'
                elif 'HO2+O2' in radical_data.columns:
                    radical_col = 'HO2+O2'
                if radical_col:
                    ax2.semilogy(radical_data['time'], radical_data[radical_col],
                               color=conc_colour, linestyle=GALLARD_LINESTYLE, linewidth=2.0,
                               alpha=0.9)

        white_file = f"../Outputs/FentonData/species_concentrations_modelled_{conc}H2O2.md"
        if os.path.exists(white_file):
            white_data = load_data_file(white_file)
            if white_data is not None and 'time' in white_data.columns:
                if 'HOOdot' in white_data.columns and 'OOdot' in white_data.columns:
                    total_radicals = white_data['HOOdot'] + white_data['OOdot']
                    mask = total_radicals > 1e-15
                    ax2.semilogy(white_data['time'][mask], total_radicals[mask],
                               color=conc_colour, linestyle=WHITE_LINESTYLE, linewidth=2.0,
                               alpha=0.9)

    ax2.set_xlabel('Time (min)')
    ax2.set_ylabel('[HO₂•] + [O₂•⁻] (M)')
    ax2.set_title('Concentration of Radicals')
    ax2.set_xlim(0, 40)
    ax2.set_ylim(1e-12, 1e-6)

    # --- Plot 3: Fe(II) Concentration ---
    for conc in concentrations:
        conc_colour = colours['concentrations'][conc]

        fe2_file = f"../Outputs/Digitised_Data/Fe2_{conc}H2O2.dat"
        if os.path.exists(fe2_file):
            fe2_data = load_data_file(fe2_file)
            if fe2_data is not None and 'time' in fe2_data.columns:
                fe_col = 'Fe2' if 'Fe2' in fe2_data.columns else 'FeII'
                if fe_col in fe2_data.columns:
                    ax3.plot(fe2_data['time'], fe2_data[fe_col],
                           color=conc_colour, linestyle=GALLARD_LINESTYLE, linewidth=2.0,
                           alpha=0.9)

        white_file = f"../Outputs/FentonData/species_concentrations_modelled_{conc}H2O2.md"
        if os.path.exists(white_file):
            white_data = load_data_file(white_file)
            if white_data is not None and 'time' in white_data.columns and 'FeII' in white_data.columns:
                ax3.plot(white_data['time'], white_data['FeII'],
                        color=conc_colour, linestyle=WHITE_LINESTYLE, linewidth=2.0,
                        alpha=0.9)

    ax3.set_xlabel('Time (min)')
    ax3.set_ylabel('Fe(II) (M)')
    ax3.set_title('Concentration of Fe(II)')
    ax3.set_xlim(0, 40)
    ax3.set_ylim(0, 7.5e-7)
    ax3.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

    # --- Plot 4: OH• Concentration ---
    for conc in concentrations:
        conc_colour = colours['concentrations'][conc]

        oh_file = f"../Outputs/Digitised_Data/OH_{conc}H2O2.dat"
        if os.path.exists(oh_file):
            oh_data = load_data_file(oh_file)
            if oh_data is not None and 'time' in oh_data.columns and 'OHdot' in oh_data.columns:
                ax4.plot(oh_data['time'], oh_data['OHdot'],
                        color=conc_colour, linestyle=GALLARD_LINESTYLE, linewidth=2.0,
                        alpha=0.9)

        white_file = f"../Outputs/FentonData/species_concentrations_modelled_{conc}H2O2.md"
        if os.path.exists(white_file):
            white_data = load_data_file(white_file)
            if white_data is not None and 'time' in white_data.columns and 'OHdot' in white_data.columns:
                ax4.plot(white_data['time'], white_data['OHdot'],
                        color=conc_colour, linestyle=WHITE_LINESTYLE, linewidth=2.0,
                        alpha=0.9)

    ax4.set_xlabel('Time (min)')
    ax4.set_ylabel('[•OH] (M)')
    ax4.set_title('Concentration of •OH')
    ax4.set_xlim(0, 40)
    ax4.set_ylim(0, 1.6e-12)
    ax4.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

    # --- Legends ---
    source_elements = []
    legend_colour = 'white' if theme == 'dark' else 'gray'

    if has_experimental:
        source_elements.append(
            Line2D([0], [0], marker='o', color='white' if theme == 'light' else 'black',
                   markerfacecolor=colours['experimental'],
                   markersize=6, markeredgecolor='white' if theme == 'light' else 'black',
                   linestyle='None', label='Experimental Data')
        )
    if has_gallard:
        source_elements.append(
            Line2D([0], [0], color=legend_colour, linestyle=GALLARD_LINESTYLE, linewidth=2.5,
                   label='Gallard & De Laat Model')
        )
    if has_white:
        source_elements.append(
            Line2D([0], [0], color=legend_colour, linestyle=WHITE_LINESTYLE, linewidth=2.5,
                   label='This Work')
        )

    conc_elements = [
        Line2D([0], [0], color=colours['concentrations']['200uM'], linestyle='-', linewidth=3,
               label='200 μM H₂O₂'),
        Line2D([0], [0], color=colours['concentrations']['1mM'], linestyle='-', linewidth=3,
               label='1 mM H₂O₂'),
        Line2D([0], [0], color=colours['concentrations']['10mM'], linestyle='-', linewidth=3,
               label='10 mM H₂O₂'),
        Line2D([0], [0], color=colours['concentrations']['100mM'], linestyle='-', linewidth=3,
               label='100 mM H₂O₂'),
        Line2D([0], [0], color=colours['concentrations']['1M'], linestyle='-', linewidth=3,
               label='1 M H₂O₂')
    ]

    # Data sources legend (first row under plots)
    if source_elements:
        legend1 = fig.legend(handles=source_elements,
                             loc='lower center',
                             bbox_to_anchor=(0.5, -0.05),
                             ncol=len(source_elements),
                             frameon=True,
                             fancybox=False,
                             shadow=False,
                             facecolor=colours['legend_face'],
                             title='Data Sources',
                             title_fontsize=14)
        legend1.get_frame().set_alpha(0.95)
        legend1.get_frame().set_linewidth(0.8)
        legend1.get_frame().set_edgecolor(colours['legend_edge'])

    # Concentration legend (second row under plots)
    legend2 = fig.legend(handles=conc_elements,
                         loc='lower center',
                         bbox_to_anchor=(0.5, -0.15),
                         ncol=len(conc_elements),
                         frameon=True,
                         fancybox=False,
                         shadow=False,
                         facecolor=colours['legend_face'],
                         title='H₂O₂ Concentrations',
                         title_fontsize=14)
    legend2.get_frame().set_alpha(0.95)
    legend2.get_frame().set_linewidth(0.8)
    legend2.get_frame().set_edgecolor(colours['legend_edge'])

    plt.tight_layout()
    plt.subplots_adjust(top=0.94, bottom=0.1, left=0.08, right=0.96,
                       hspace=0.35, wspace=0.3)

    # Panel labels are placed in figure coordinates, flush with the left edge
    # of each panel's y-axis label, so they line up regardless of how wide
    # that panel's tick labels are (e.g. the scientific-notation offset on C/D).
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax, panel_label in zip((ax1, ax2, ax3, ax4), ('(A)', '(B)', '(C)', '(D)')):
        ylabel_bbox = ax.yaxis.label.get_window_extent(renderer=renderer)
        x_fig = fig.transFigure.inverted().transform((ylabel_bbox.x0, 0))[0]
        y_fig = ax.get_position().y1 + 0.02
        fig.text(x_fig, y_fig, panel_label, fontsize=16, fontweight='bold',
                  verticalalignment='bottom', horizontalalignment='left')

    suffix = '_dark' if theme == 'dark' else ''
    bg_colour = 'black' if theme == 'dark' else 'white'

    plt.savefig(f'../Outputs/Fenton_Model_Validation{suffix}.png', dpi=400,
                bbox_inches='tight', facecolor=bg_colour, edgecolor='none')
    plt.savefig(f'../Outputs/Fenton_Model_Validation{suffix}.pdf',
                bbox_inches='tight', facecolor=bg_colour, edgecolor='none')

    print(f"{theme.capitalize()} theme plots saved:")
    print(f"   ../Outputs/Fenton_Model_Validation{suffix}.png")
    print(f"   ../Outputs/Fenton_Model_Validation{suffix}.pdf")

    plt.close()

if __name__ == "__main__":
    os.makedirs('../Outputs', exist_ok=True)

    print("Creating light theme plots...")
    create_fenton_plots(theme='light')

    print("\nCreating dark theme plots...")
    create_fenton_plots(theme='dark')

    print("All plots created successfully!")
