# This is the multi-scenario prebiotic plot.
# It reads the .npz files from Outputs/PrebioticData/ (produced by prebiotic_scenario_runner.py) and generates two figures, one per iron level. npz chosen here as it is faster to load than the original ,md method.
# Each figure has 4 panels (HCN, Formamide & Formate, ADMS & AMS, Fe Complexes).
# A solid line means 1 ppb HCN; a dashed line means 100 ppb HCN.
# The shaded band is +/-1 standard deviation across the 3 runs.
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib import rcParams

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, '..', 'Outputs', 'PrebioticData')
PLOT_DIR   = os.path.join(SCRIPT_DIR, '..', 'Outputs', 'PrebioticData')
os.makedirs(PLOT_DIR, exist_ok=True)

FLOOR = 1e-11   # This matches the solver's BASE_ATOL - concentrations below this tolerance are numerical noise, so they are not plotted.

SPECIES_LABELS = {
    'HCN':             'HCN',
    'HCONH2':          r'HCONH$_2$',
    'HCO2MINUS1':      r'HCO$_2^-$',
    'NH2CHS2O6MINUS2': r'NH$_2$CH(SO$_3^-$)$_2$',
    'NH2CH2SO3MINUS1': r'NH$_2$CH$_2$SO$_3^-$',
    'FeCN6MINUS4':     r'Fe(CN)$_6^{4-}$',
    'FeCN5H2OMINUS3':  r'Fe(CN)$_5$(H$_2$O)$^{3-}$',
    'FePLUS2':         r'Fe$^{2+}$',
}

PANELS = [
    ('HCN',                ['HCN']),
    ('Formamide & Formate', ['HCONH2', 'HCO2MINUS1']),
    ('ADMS & AMS',          ['NH2CHS2O6MINUS2', 'NH2CH2SO3MINUS1']),
    ('Fe Complexes',        ['FeCN6MINUS4', 'FeCN5H2OMINUS3', 'FePLUS2']),
]


def get_theme_cfg(theme):
    # Purpose: switches matplotlib's global style and returns the matching colour palette.
    # Inputs are: theme, either 'light' or 'dark'.
    # Outputs are: a dictionary of colours to use for the plot elements (species, legend, the HCN line style key), suited to the chosen theme. matplotlib's global plotting style is also updated as a side effect.
    # Called by: plot_figure(), in this file.
    if theme == 'dark':
        rcParams.update({
            'font.family': 'Helvetica', 'font.size': 14,
            'figure.facecolor': '#000000', 'axes.facecolor': '#000000',
            'axes.edgecolor': 'white', 'axes.labelcolor': 'white',
            'text.color': 'white', 'xtick.color': 'white',
            'ytick.color': 'white', 'grid.color': '#333333',
            'legend.facecolor': '#1a1a1a', 'legend.edgecolor': '#444444',
        })
        return {
            'species': ["#7fc4b0", "#C4A484", "#ff9966", "#ff8c42", "#40c9a2", "#9b8bda"],
            'legend_face': '#1a1a1a', 'legend_edge': '#444444',
            'hcn_line': 'white',
        }
    else:
        rcParams.update({
            'font.family': 'Helvetica', 'font.size': 14,
            'figure.facecolor': 'white', 'axes.facecolor': 'white',
            'axes.edgecolor': 'black', 'axes.labelcolor': 'black',
            'text.color': 'black', 'xtick.color': 'black', 'ytick.color': 'black',
        })
        return {
            'species': ["#6C9A8B", "#5C4033", "#F37B3A", "#D95F02", "#1B9E77", "#7570B3"],
            'legend_face': 'white', 'legend_edge': 'gray',
            'hcn_line': 'black',
        }


def load_npz(label):
    # Purpose: loads one scenario's saved results from its .npz file.
    # Inputs are: label, the scenario label used as the .npz filename (e.g. 'LowHCN_LowFe').
    # Outputs are: t, the time points in days; mean and std, the mean and standard deviation of every species' concentration across the repeated runs; idx, a dictionary mapping each species name to its row index in mean/std.
    # Called by: plot_figure(), in this file.
    path = os.path.join(DATA_DIR, f'{label}.npz')
    d    = np.load(path, allow_pickle=True)
    t    = d['t'] / 86400 # Need to change if anyone ever wants to not plot days OR used a different time unit in the solver. 
    idx  = {s: i for i, s in enumerate(d['species'])}
    return t, d['mean'], d['std'], idx


def plot_figure(iron_label, low_hcn_label, high_hcn_label, theme='light'):
    # Purpose: draws and saves the 4-panel low-vs-high-HCN comparison figure for one iron level.
    # Inputs are: iron_label, a display label for the iron level ('HighFe' or 'LowFe'), used in the output filename and to choose the Fe-panel y-axis range; low_hcn_label and high_hcn_label, the .npz scenario labels for the 1 ppb and 100 ppb HCN runs at this iron level; theme, either 'light' or 'dark'.
    # Outputs are: nothing is returned; a 4-panel figure comparing the low- and high-HCN scenarios at this iron level is drawn and saved as a .png file.
    # Called by: the script's entry point, at the bottom of this file (if __name__ == "__main__"), once per iron level and theme.
    cs      = get_theme_cfg(theme)
    colours = cs['species']

    t1, mean1, std1, idx1 = load_npz(low_hcn_label)
    t2, mean2, std2, idx2 = load_npz(high_hcn_label)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    # Each panel's y-axis limits are tuned to the actual range of data it shows.
    d_top = 0.1 if iron_label == 'HighFe' else 1e-5
    YLIMS = [
        (8e-12, 5e-6),      # A: 1ppb HCN at ~9e-9, 100ppb at ~9e-7
        (FLOOR, 5e-2),      # B: formamide starts near FLOOR; formate tops ~6e-3
        (FLOOR, 1e-5),      # C: AMS peaks at 1.6e-6 in all scenarios
        (FLOOR, d_top),     # D: Fe Complexes
    ]

    for pi, (ax, letter, (title, sp_list)) in enumerate(
            zip(axes, 'abcd', PANELS)):

        for ci, sp in enumerate(sp_list):
            col = colours[ci % len(colours)]
            m1 = m2 = None

            if sp in idx1:
                m1 = np.maximum(mean1[idx1[sp], :], FLOOR)
                ax.plot(t1, m1, color=col, linestyle='-', linewidth=1.8)
            if sp in idx2:
                m2 = np.maximum(mean2[idx2[sp], :], FLOOR)
                ax.plot(t2, m2, color=col, linestyle='--', linewidth=1.8)

            # The range between the two HCN scenarios is shaded, to show how much this species' trajectory depends on the starting HCN level.
            if m1 is not None and m2 is not None:
                ax.fill_between(t1,
                                np.minimum(m1, m2),
                                np.maximum(m1, m2),
                                color=col, alpha=0.20)

        ax.set_yscale('log')
        ax.set_ylim(*YLIMS[pi])
        ax.set_title(title, fontsize=14, fontweight='bold', pad=6)
        ax.grid(False)
        ax.set_xlabel('Time (days)', fontsize=14)
        ax.set_ylabel('Concentration (M)', fontsize=14)

        ax.text(-0.10, 1.05, f'({letter.upper()})', transform=ax.transAxes,
                fontsize=14, fontweight='bold', va='bottom', ha='left')

        handles = [
            plt.Line2D([0], [0], color=colours[ci % len(colours)], lw=3,
                       label=SPECIES_LABELS.get(sp, sp))
            for ci, sp in enumerate(sp_list)
            if sp in idx1 or sp in idx2
        ]
        if handles:
            leg = ax.legend(handles=handles, fontsize=12,
                            loc='best', frameon=True)
            leg.get_frame().set_facecolor(cs['legend_face'])
            leg.get_frame().set_edgecolor(cs['legend_edge'])
            leg.get_frame().set_alpha(0.9)

        # Only panel D of the LowFe figure gets this inset: a linear-scale, nanomolar zoom on the fast Fe transient at the very start of the run, which is otherwise too small to see on the main log-scale axis.
        if pi == 3 and iron_label == 'LowFe':
            axins = ax.inset_axes([0.13, 0.55, 0.43, 0.38])
            t_sec = t1 * 86400
            INSET_MAX_S = 10_000
            mask = t_sec <= INSET_MAX_S
            for ci2, sp2 in enumerate(sp_list):
                col2 = colours[ci2 % len(colours)]
                if sp2 in idx1:
                    axins.plot(t_sec[mask], mean1[idx1[sp2], mask] * 1e9,
                               color=col2, linestyle='-', linewidth=1.2)
                if sp2 in idx2:
                    axins.plot(t_sec[mask], mean2[idx2[sp2], mask] * 1e9,
                               color=col2, linestyle='--', linewidth=1.2)
            axins.set_xlim(0, INSET_MAX_S)
            axins.set_ylim(0, 1.3)
            axins.set_xlabel('Time (s)', fontsize=9)
            axins.set_ylabel('(nM)', fontsize=9)

            axins.tick_params(labelsize=8)
            axins.xaxis.get_major_formatter().set_useOffset(False)
            ins_face = '#111111' if theme == 'dark' else '#f5f5f5'
            axins.set_facecolor(ins_face)
            for spine in axins.spines.values():
                spine.set_edgecolor('gray')
                spine.set_linewidth(0.8)

    # A single legend for the whole figure encodes the HCN level via line style (solid/dashed), rather than repeating it on every panel.
    hcn_col = cs['hcn_line']
    gl_handles = [
        mlines.Line2D([], [], color=hcn_col, linestyle='-',  linewidth=2, label='1 ppb HCN'),
        mlines.Line2D([], [], color=hcn_col, linestyle='--', linewidth=2, label='100 ppb HCN'),
    ]
    gl = fig.legend(handles=gl_handles, loc='lower center', ncol=2,
                    frameon=True, fontsize=13,
                    title='Initial Atmospheric HCN Concentration')
    gl.get_frame().set_facecolor(cs['legend_face'])
    gl.get_frame().set_edgecolor(cs['legend_edge'])
    gl.get_frame().set_alpha(0.9)

    plt.tight_layout(rect=[0, 0.08, 1, 0.95])

    suffix    = '_dark' if theme == 'dark' else ''
    bg_colour = 'black'  if theme == 'dark' else 'white'
    out_path = os.path.join(PLOT_DIR, f'Prebiotic_{iron_label}{suffix}.png')
    plt.savefig(out_path, dpi=400, bbox_inches='tight',
                facecolor=bg_colour, edgecolor='none')
    print(f"Saved: {out_path}")
    plt.close()


if __name__ == '__main__':
    for theme in ['light', 'dark']:
        plot_figure('HighFe', 'LowHCN_HighFe', 'HighHCN_HighFe', theme=theme)
        plot_figure('LowFe',  'LowHCN_LowFe',  'HighHCN_LowFe',  theme=theme)
    print("\nAll plots saved.")
