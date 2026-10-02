# This script runs the prebiotic chemistry across a grid of scenarios, to see how HCN and iron levels affect the outcome.
# It runs all 4 combinations of low/high HCN and low/high Fe2+, 3 times each (to get a mean and standard deviation), without editing any input files i.e., parameters are changed in memory only.
# Results (mean and standard deviation on a shared time grid) are saved to Outputs/PrebioticData/.
import sys, os, gc, time, subprocess
import numpy as np
from scipy.integrate import solve_ivp

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SOUP_DIR   = SCRIPT_DIR
sys.path.insert(0, SOUP_DIR)
os.chdir(SOUP_DIR)

from rate_laws     import generate_rate_laws_with_database_type
from differentials import construct_differential_equations, create_bulletproof_ode_system, REGORD_SPECIES_CSTAR
from data          import read_dat_file, read_aqueous_dat_file
from utilities     import initialise_concentrations

# -- Configuration ------------------------------------------------------------
DB_NAME   = 'prebiotic'
DB_TYPE   = 'prebiotic'
END_TIME  = 7 * 86400       # 7 days
N_RUNS    = 3
N_GRID    = 500             # timepoints on common grid
BASE_ATOL = 1e-11
BASE_RTOL = 1e-8

# 1 ppb HCN gives a partial pressure of 1e-9 atm, which equilibrates to about 9 nM.
# 100 ppb HCN gives a partial pressure of 1e-7 atm, which equilibrates to about 900 nM.
SCENARIOS = [
    {'hcn_atm': 1e-9,  'fe_M': 1e-9,  'label': 'LowHCN_LowFe',  'hcn_ppb': 1,   'fe_nM': 1},
    {'hcn_atm': 1e-9,  'fe_M': 10e-3, 'label': 'LowHCN_HighFe', 'hcn_ppb': 1,   'fe_mM': 10},
    {'hcn_atm': 1e-7,  'fe_M': 1e-9,  'label': 'HighHCN_LowFe', 'hcn_ppb': 100, 'fe_nM': 1},
    {'hcn_atm': 1e-7,  'fe_M': 10e-3, 'label': 'HighHCN_HighFe','hcn_ppb': 100, 'fe_mM': 10},
]

OUT_DIR = os.path.join(SOUP_DIR, 'Outputs', 'PrebioticData')
os.makedirs(OUT_DIR, exist_ok=True)

COMMON_TIMES = np.linspace(0, END_TIME, N_GRID)

# -- Memory check --------------------------------------------------------------
# This looks for other Python processes already running, since this script can use a lot of memory and running several copies at once risks running out. Can bypass this if you know about hardware...
result = subprocess.run(['pgrep', '-a', 'python3'], capture_output=True, text=True)
other = [l for l in result.stdout.strip().split('\n')
         if l and str(os.getpid()) not in l and l.strip()]
if other:
    print("WARNING: Other Python processes detected:")
    for l in other:
        print(f"  {l}")
    print("Consider killing them to free memory before proceeding.\n")
else:
    print("Memory check: no other Python processes running.\n")

# -- Rate laws are loaded once here, and shared across every scenario below, since they don't change between scenarios. --
print("Loading rate laws...", flush=True)
rate_laws, database_species, _ = generate_rate_laws_with_database_type(DB_TYPE, DB_NAME)
diff_eqs = construct_differential_equations(rate_laws)
n_species = len(database_species)
print(f"  {n_species} species loaded.\n", flush=True)

sp_idx_map = {sp: i for i, sp in enumerate(database_species)}

# -- The base input files are also read once here; each scenario below only overrides the HCN and Fe2+ values, working from copies of these base lists. --
atm_file = os.path.join(SOUP_DIR, 'Inputs', f'{DB_NAME}-atmosphere.dat')
aq_file  = os.path.join(SOUP_DIR, 'Inputs', f'{DB_NAME}-aqueous.dat')

species_atm_base, pp_base, hc_base = read_dat_file(atm_file)
species_aq_base, conc_base, vol, temp, pH, _ = read_aqueous_dat_file(aq_file)

# HCN is located in the atmosphere list.
hcn_atm_idx = species_atm_base.index('HCN')

# Fe2+ is located in the aqueous list (normalised to 'FePLUS2').
fe_aq_idx = species_aq_base.index('FePLUS2')

# -- Every species gets its own absolute tolerance here, tightened for species with a blended (regord) rate law near their cstar, for the same reason described in differentials.py --
atol_arr = np.full(n_species, BASE_ATOL)
for sp, cstar_val in REGORD_SPECIES_CSTAR.items():
    if sp in sp_idx_map:
        tight = max(cstar_val * 1e-8, 1e-13)   # This floor stops the tolerance being tightened too far for a species with a very small cstar, such as HCN.
        if tight < BASE_ATOL:
            atol_arr[sp_idx_map[sp]] = tight

# -- Main loop ------------------------------------------------------------------
total_start = time.time()

for scenario in SCENARIOS:
    label   = scenario['label']
    hcn_ppb = scenario['hcn_ppb']
    fe_str  = f"1 nM" if scenario['fe_M'] < 1e-6 else f"10 mM"

    out_path = os.path.join(OUT_DIR, f'{label}.npz')
    if os.path.exists(out_path):
        print(f"SKIP {label} — {out_path} already exists.", flush=True)
        continue

    print("=" * 60)
    print(f"SCENARIO: {label}  (HCN = {hcn_ppb} ppb, Fe²⁺ = {fe_str})")
    print("=" * 60, flush=True)

    # The base partial-pressure and concentration lists are copied here, so this scenario's changes don't affect the base lists used by the next scenario.
    pp   = list(pp_base)
    conc = list(conc_base)
    pp[hcn_atm_idx]  = scenario['hcn_atm']
    conc[fe_aq_idx]  = scenario['fe_M']

    # The starting concentrations are rebuilt using this scenario's modified HCN and Fe2+ values.
    y0 = np.maximum(
        np.array(initialise_concentrations(
            database_species, species_aq_base, conc,
            species_atm_base, pp, hc_base), dtype=float),
        0.0)

    run_traces = []

    for run_idx in range(1, N_RUNS + 1):
        print(f"\n  Run {run_idx}/{N_RUNS}...", flush=True)
        t_run = time.time()

        # rate_laws is passed in so the analytic Jacobian is compiled; reaction-rate logging is then switched off immediately below, since it isn't needed here and would otherwise add memory overhead across N_RUNS repeats. Again, feel free to turn on if your PC can manage it. You can also edit how much logging is done.
        ode_sys = create_bulletproof_ode_system(
            database_species, diff_eqs,
            species_atm_base, pp, hc_base, pH, rate_laws)
        ode_sys.set_total_time(END_TIME)
        ode_sys.reaction_logger = None

        sol = solve_ivp(
            ode_sys, (0, END_TIME), y0.copy(),
            method='BDF', rtol=BASE_RTOL, atol=atol_arr,
            max_step=3600.0, dense_output=False, jac=ode_sys.jacobian)

        elapsed = time.time() - t_run
        status  = 'SUCCESS' if sol.success else 'FAILED'
        print(f"  {status} — {len(sol.t)} steps, {elapsed:.1f}s", flush=True)

        if not sol.success:
            print(f"  Solver message: {sol.message}")

        # Any CN- excursion below zero is reported here, as a sanity check on the numerical error.
        cn_key = 'CNMINUS1'
        if cn_key in sp_idx_map:
            cn_min = sol.y[sp_idx_map[cn_key], :].min()
            if cn_min < -BASE_ATOL * 10:
                print(f"  CN⁻ min excursion: {cn_min:.3e} M")

        # The result is interpolated onto the shared time grid here, and then deleted straight away, to keep memory use down across the N_RUNS repeats.
        from scipy.interpolate import interp1d
        trace = np.zeros((n_species, N_GRID))
        for si in range(n_species):
            f = interp1d(sol.t, sol.y[si, :], bounds_error=False,
                         fill_value=(sol.y[si, 0], sol.y[si, -1]))
            trace[si, :] = f(COMMON_TIMES)
        run_traces.append(trace)

        del sol, ode_sys
        gc.collect()

    # The mean and sample standard deviation are taken across the N_RUNS repeats of this scenario.
    stack = np.stack(run_traces, axis=0)   # (N_RUNS, n_species, N_GRID)
    mean  = stack.mean(axis=0)
    std   = stack.std(axis=0, ddof=1)

    np.savez(out_path,
             t=COMMON_TIMES,
             mean=mean,
             std=std,
             species=np.array(database_species, dtype=object))
    print(f"\n  Saved: {out_path}", flush=True)
    gc.collect()

total_elapsed = time.time() - total_start
print(f"\nAll scenarios complete. Total wall time: {total_elapsed/60:.1f} min")
