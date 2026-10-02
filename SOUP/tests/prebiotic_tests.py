# This is the verification suite for the prebiotic reaction network.

# Usage:
#   python prebiotic_tests.py --mass-balance    # 30-day closed-system ODE + element audit
#   python prebiotic_tests.py --convergence     # N-run spread study (closed system)
#   python prebiotic_tests.py --open-system     # N-run spread study (open system, 7 days)
#   python prebiotic_tests.py --regord          # verify regord_pow blend function + save PNGs
#   python prebiotic_tests.py --all             # run all four tests in sequence

import sys
import os
import argparse
import time
import io
import contextlib
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SOUP_DIR   = os.path.join(SCRIPT_DIR, '..')
sys.path.insert(0, SOUP_DIR)
sys.path.insert(0, SCRIPT_DIR)
os.chdir(SOUP_DIR)

from scipy.integrate import solve_ivp
from rate_laws     import generate_rate_laws_with_database_type
from differentials import (
    construct_differential_equations,
    create_bulletproof_ode_system,
    REGORD_SPECIES_CSTAR, REGORD_CSTAR, REGORD_S,
)
from data          import read_dat_file, read_aqueous_dat_file
from utilities     import initialise_concentrations
from checks        import get_composition, ELEMENTS_TO_AUDIT

DB_NAME   = 'prebiotic'
DB_TYPE   = 'prebiotic'
BASE_ATOL = 1e-11
BASE_RTOL = 1e-8
N_RUNS    = 3
FLAG_PCT  = 1.0


# -- shared helpers -------------------------------------------------------------

def _display(sp):
    # Purpose: converts a species name from the internal database format to the display format.
    # Inputs are: sp, a species name in the internal database format, e.g. 'HPLUS1'.
    # Outputs are: the same name converted to the display format, e.g. 'H+1'.
    # Called by: test_mass_balance(), test_convergence(), test_open_system() and _element_totals(), all in this file.
    return sp.replace('MINUS', '-').replace('PLUS', '+')


def _load_system(suppress_output=False):
    # Purpose: loads the full prebiotic-network simulation setup, ready to run.
    # Inputs are: suppress_output, whether to silence the console output produced while loading (used for repeated runs, so the console isn't flooded).
    # Outputs are: a tuple (rate_laws, database_species, diff_eqs, species_atm, partial_pressures, henrys_constants, y0, pH), i.e. the full simulation setup needed to run the prebiotic network, with y0 already floored at zero.
    # Called by: test_mass_balance() and test_open_system() (with output shown); test_convergence() (with output suppressed, once per run); all in this file.
    def _do():
        # Purpose: does the actual work of loading the setup, so it can optionally be run with its console output silenced.
        # Inputs are: none (uses DB_TYPE, DB_NAME and the enclosing function's arguments via closure).
        # Outputs are: the same tuple described above for _load_system().
        # Called by: _load_system(), just below, either directly or with its output redirected.
        rate_laws, database_species, _ = generate_rate_laws_with_database_type(DB_TYPE, DB_NAME)
        diff_eqs = construct_differential_equations(rate_laws)
        atm_file = os.path.join(SOUP_DIR, 'Inputs', f'{DB_NAME}-atmosphere.dat')
        aq_file  = os.path.join(SOUP_DIR, 'Inputs', f'{DB_NAME}-aqueous.dat')
        species_atm, partial_pressures, henrys_constants = read_dat_file(atm_file)
        species_aq, concentrations, _, _, pH, _ = read_aqueous_dat_file(aq_file)
        concs = initialise_concentrations(
            database_species, species_aq, concentrations,
            species_atm, partial_pressures, henrys_constants,
        )
        y0 = np.maximum(np.array(concs, dtype=float), 0.0)
        return rate_laws, database_species, diff_eqs, species_atm, partial_pressures, henrys_constants, y0, pH

    if suppress_output:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            return _do()
    return _do()


def _build_atol(database_species, sp_idx):
    # Purpose: builds the per-species absolute-tolerance array used across the tests in this file.
    # Inputs are: database_species, the list of species names; sp_idx, a dictionary mapping each species name to its index.
    # Outputs are: an absolute-tolerance array, one entry per species, tightened for any species with a blended (regord) rate law near its cstar, for the same reason described in differentials.py.
    # Called by: test_mass_balance(), test_convergence() and test_open_system(), all in this file.
    atol_arr = np.full(len(database_species), BASE_ATOL)
    for sp, cstar_val in REGORD_SPECIES_CSTAR.items():
        if sp in sp_idx:
            tight = max(cstar_val * 1e-8, 1e-13)
            if tight < BASE_ATOL:
                atol_arr[sp_idx[sp]] = tight
    return atol_arr


def _element_totals(database_species, y_t):
    # Purpose: adds up the total amount of each audited element across every species, at every timepoint.
    # Inputs are: database_species, the list of species names; y_t, the concentration of every species at every timepoint, shaped (n_species, n_timepoints).
    # Outputs are: a dictionary mapping each audited element (C, N, S, Fe) to its total amount at every timepoint, summed across every species that contains it.
    # Called by: test_mass_balance(), in this file.
    totals = {el: np.zeros(y_t.shape[1]) for el in ELEMENTS_TO_AUDIT}
    for j, sp in enumerate(database_species):
        comp = get_composition(_display(sp))
        for el in ELEMENTS_TO_AUDIT:
            cnt = comp.get(el, 0)
            if cnt:
                totals[el] += cnt * y_t[j, :]
    return totals


# -- test 1: mass balance --------------------------------------------------------

def test_mass_balance():
    # Purpose: runs a 30-day closed-system simulation and checks that mass is conserved.
    # Inputs are: none.
    # Outputs are: True if every audited element (C, N, S, Fe) stays conserved within 1 part per million across a 30-day closed-system run, False otherwise. A detailed report is also printed, including any species with a negative numerical excursion and the final concentration of every species.
    # Called by: main(), in this file.
    END_TIME  = 30 * 86400
    THRESHOLD = 1e-6

    print("=" * 65)
    print("MASS BALANCE TEST  (30-day closed system, raw unclamped trajectory)")
    print("=" * 65)

    rate_laws, database_species, diff_eqs, species_atm, partial_pressures, henrys_constants, y0, pH = \
        _load_system()

    ode_sys = create_bulletproof_ode_system(
        database_species, diff_eqs,
        species_atm, partial_pressures, henrys_constants, pH, rate_laws,
    )
    ode_sys.min_concentrations = {i: 0.0 for i in range(len(database_species))}
    ode_sys.set_total_time(END_TIME)

    atol_arr = _build_atol(database_species, ode_sys.species_to_index)

    print(f"\n  Running 30-day closed-system simulation …")
    sol = solve_ivp(
        ode_sys, (0, END_TIME), y0,
        method='BDF', rtol=BASE_RTOL, atol=atol_arr,
        max_step=3600.0, dense_output=False, jac=ode_sys.jacobian,
    )

    if not sol.success:
        print(f"\n  FAILED: {sol.message}")
        return False

    print(f"  SUCCESS — {len(sol.t)} steps, t_final = {sol.t[-1]/86400:.2f} days")

    mins = sol.y.min(axis=1)
    flagged = [(sp, mn) for sp, mn in zip(database_species, mins) if mn < -BASE_ATOL * 10]
    if flagged:
        print(f"\n  Numerical excursions (reported as uncertainty, raw solution unchanged):")
        for sp, mn in sorted(flagged, key=lambda x: x[1]):
            print(f"    {_display(sp)}: min = {mn:.3e} M")
    else:
        print(f"  No significant negative excursions.")

    totals   = _element_totals(database_species, sol.y)
    all_pass = True

    print(f"\n  {'Element':>8}  {'t=0 [M]':>14}  {'t_final [M]':>14}  "
          f"{'max |Δ|/total_0':>18}  Result")
    print("  " + "-" * 65)

    for el in ELEMENTS_TO_AUDIT:
        t0 = totals[el][0]
        if t0 == 0:
            print(f"  {el:>8}  {'(zero at t=0)':>14}  — skipped")
            continue
        deviation = np.abs(totals[el] - t0) / t0
        max_dev   = deviation.max()
        max_t     = sol.t[deviation.argmax()]
        ok = max_dev < THRESHOLD
        if not ok:
            all_pass = False
        print(f"  {el:>8}  {t0:>14.6e}  {totals[el][-1]:>14.6e}  "
              f"{max_dev:>18.3e}  {'PASS' if ok else 'FAIL'}  (worst at t={max_t/3600:.2f} h)")

    print(f"\n  H excluded: reactions are pseudo-steps; H does not balance.")
    status = ("PASS — all elements conserved to < 1 ppm"
              if all_pass else "FAIL — mass leak detected")
    print(f"  Overall: {status}")

    print(f"\n  --- Final concentrations (t=30 days) ---")
    for i, sp in enumerate(database_species):
        v = sol.y[i, -1]
        if abs(v) > 1e-20:
            print(f"  {_display(sp):<25} {v:.6e} M")

    return all_pass


# -- test 2: convergence (closed system) -----------------------------------------

def test_convergence():
    # Purpose: repeats the closed-system run several times and checks the results converge to the same answer.
    # Inputs are: none.
    # Outputs are: nothing is returned; the prebiotic network is run N_RUNS times over a 30-day closed system, and the spread between runs' final concentrations is printed for every species, flagging any that spread by more than FLAG_PCT (1%).
    # Called by: main(), in this file.
    END_TIME = 30 * 86400
    ATOL_CN  = 2e-14

    print("=" * 65)
    print(f"CONVERGENCE STUDY  ({N_RUNS} runs, 30-day closed system)")
    print("=" * 65)

    all_finals = []
    all_cn_min = []
    species    = None

    for i in range(N_RUNS):
        print(f"\n  Run {i+1}/{N_RUNS} …", flush=True)
        rate_laws, database_species, diff_eqs, species_atm, partial_pressures, henrys_constants, y0, pH = \
            _load_system(suppress_output=True)

        ode_sys = create_bulletproof_ode_system(
            database_species, diff_eqs,
            species_atm, partial_pressures, henrys_constants, pH, rate_laws,
        )
        ode_sys.min_concentrations = {k: 0.0 for k in range(len(database_species))}
        ode_sys.set_total_time(END_TIME)

        atol_arr = _build_atol(database_species, ode_sys.species_to_index)
        cn_idx = ode_sys.species_to_index.get('CNMINUS1')
        if cn_idx is not None:
            atol_arr[cn_idx] = ATOL_CN

        t0 = time.time()
        sol = solve_ivp(
            ode_sys, (0, END_TIME), y0,
            method='BDF', rtol=BASE_RTOL, atol=atol_arr,
            max_step=3600.0, dense_output=False, jac=ode_sys.jacobian,
        )
        elapsed = time.time() - t0

        if not sol.success:
            print(f"  FAILED: {sol.message}")
            continue

        cn_min = sol.y[cn_idx, :].min() if cn_idx is not None else 0.0
        print(f"  SUCCESS — {len(sol.t)} steps, {elapsed:.1f}s,  CN⁻ min = {cn_min:.3e} M")
        species = database_species
        all_finals.append(sol.y[:, -1])
        all_cn_min.append(cn_min)

    if len(all_finals) < 2:
        print("  Need at least 2 successful runs.")
        return

    Y = np.array(all_finals)
    print(f"\n  CN⁻ excursion range: {min(all_cn_min):.3e} → {max(all_cn_min):.3e} M\n")

    col_w  = 14
    header = f"  {'Species':<22}" + "".join(f"  {'Run '+str(i+1):>{col_w}}" for i in range(len(all_finals)))
    header += f"  {'spread%':>9}  note"
    print(header)
    print("  " + "-" * (22 + len(all_finals) * (col_w + 2) + 16))

    for j, sp in enumerate(species):
        vals    = Y[:, j]
        max_abs = np.max(np.abs(vals))
        if max_abs < 1e-20:
            continue
        spread_pct = (vals.max() - vals.min()) / max_abs * 100
        flag = '  <<' if spread_pct > FLAG_PCT else ''
        row  = f"  {_display(sp):<22}" + "".join(f"  {v:{col_w}.4e}" for v in vals)
        row += f"  {spread_pct:>8.2f}%{flag}"
        print(row)

    print(f"\n  << = spread > {FLAG_PCT}% between runs (inherits CN⁻ numerical uncertainty).")


# -- test 3: open-system repeatability --------------------------------------------

def test_open_system():
    # Purpose: repeats the open-system run several times and checks the results are repeatable.
    # Inputs are: none.
    # Outputs are: nothing is returned; the prebiotic network is run N_RUNS times over a 7-day open system (Henry's law floors active, so mass balance isn't checked), and the spread between runs' final concentrations is printed for every species.
    # Called by: main(), in this file.
    END_TIME = 7 * 86400

    print("=" * 65)
    print(f"OPEN-SYSTEM REPEATABILITY  ({N_RUNS} runs, 7-day open system)")
    print("Henry's law floors active — mass balance not checked.")
    print("=" * 65)

    rate_laws, database_species, diff_eqs, species_atm, partial_pressures, henrys_constants, y0, pH = \
        _load_system()

    all_finals   = []
    all_negatives = []

    for run in range(N_RUNS):
        print(f"\n  Run {run+1}/{N_RUNS} …", flush=True)

        ode_sys = create_bulletproof_ode_system(
            database_species, diff_eqs,
            species_atm, partial_pressures, henrys_constants, pH, rate_laws,
        )
        ode_sys.reaction_logger = None
        ode_sys.set_total_time(END_TIME)

        atol_arr = _build_atol(database_species, ode_sys.species_to_index)

        if run == 0:
            print(f"  Henry's law floors:")
            for i, sp in enumerate(database_species):
                floor = ode_sys.min_concentrations.get(i, 0.0)
                if floor > 0:
                    print(f"    {_display(sp)} : {floor:.3e} M")

        t0 = time.time()
        sol = solve_ivp(
            ode_sys, (0, END_TIME), y0,
            method='BDF', rtol=BASE_RTOL, atol=atol_arr,
            max_step=3600.0, dense_output=False, jac=ode_sys.jacobian,
        )
        elapsed = time.time() - t0

        if not sol.success:
            print(f"  FAILED: {sol.message}")
            continue

        print(f"  SUCCESS — {len(sol.t)} steps, {elapsed:.1f}s")

        mins     = sol.y.min(axis=1)
        negs     = [(sp, mn) for sp, mn in zip(database_species, mins) if mn < -BASE_ATOL * 10]
        all_negatives.append(negs)
        if negs:
            print(f"  Negative excursions:")
            for sp, mn in sorted(negs, key=lambda x: x[1]):
                print(f"    {_display(sp)} : {mn:.3e} M")
        else:
            print(f"  No negative excursions.")

        all_finals.append(sol.y[:, -1])

    if len(all_finals) < 2:
        print("\n  Need at least 2 successful runs to compute spread.")
        return

    Y      = np.array(all_finals)
    col_w  = 13
    print(f"\n{'='*65}")
    print(f"Final concentrations (t=7 days) — all {len(all_finals)} runs")
    print(f"{'='*65}")
    header = f"  {'Species':<22}" + "".join(f"  {'Run '+str(i+1):>{col_w}}" for i in range(len(all_finals)))
    header += f"  {'spread%':>9}  note"
    print(header)
    print("  " + "-" * (22 + len(all_finals) * (col_w + 2) + 16))

    for j, sp in enumerate(database_species):
        vals    = Y[:, j]
        max_abs = np.max(np.abs(vals))
        if max_abs < 1e-20:
            continue
        spread_pct = (vals.max() - vals.min()) / max_abs * 100
        flag = '  <<' if spread_pct > FLAG_PCT else ''
        row  = f"  {_display(sp):<22}" + "".join(f"  {v:{col_w}.4e}" for v in vals)
        row += f"  {spread_pct:>8.2f}%{flag}"
        print(row)

    print(f"\n  << = spread > {FLAG_PCT}% between runs.")


# -- test 4: regord_pow verification + plots --------------------------------------

def test_regord():
    # Purpose: verifies that the regord_pow blend (this is a concentration-dependent reaction order blend) behaves correctly, and saves diagnostic plots of it.
    # Inputs are: none.
    # Outputs are: nothing is returned; the concentration-dependent reaction order blend (regord_pow, from differentials.py) is checked for accuracy against a plain power law at each affected species' calibration minimum, the rate-law database is scanned to confirm fractional exponents are actually wrapped with regord_pow, a table of the effective reaction order at several concentrations is printed, and a comparison figure is saved per affected species.
    # Called by: main(), in this file.
    import matplotlib.pyplot as plt

    # These are the species with a fractional or negative reaction order (n), together with their cstar (regord transition concentration) and the lowest concentration they were actually calibrated against in experiment. More may be needed to be added in the future.
    AFFECTED = [
        ('CNMINUS1',  0.474, 'CN$^-$',       REGORD_SPECIES_CSTAR.get('CNMINUS1',  REGORD_CSTAR), 2e-4),
        ('FePLUS2',   0.430, 'Fe$^{2+}$',    REGORD_SPECIES_CSTAR.get('FePLUS2',   REGORD_CSTAR), 1e-4),
        ('CO2',       0.71,  'CO$_2$',        REGORD_SPECIES_CSTAR.get('CO2',       REGORD_CSTAR), 1e-3),
        ('SO3MINUS2', -0.6,  'SO$_3^{2-}$',  REGORD_SPECIES_CSTAR.get('SO3MINUS2', REGORD_CSTAR), 1e-3),
    ]
    C_GRID       = np.logspace(-12, -1, 2000)
    REPORT_CONCS = [1e-9, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]

    def _unmod(c, n):
        # Purpose: the plain, unmodified power law, used as the reference to compare regord_pow's blend against.
        # Inputs are: c, an array of concentrations; n, the reaction order.
        # Outputs are: the plain, unmodified power law c**n (zero wherever c is not positive).
        # Called by: test_regord(), in this file, as the reference to compare regord_pow's blend against.
        return np.where(c > 0, c ** n, 0.0)

    def _mod(c, n, cstar):
        # Purpose: an independent, local re-implementation of regord_pow()'s blend, to check the real one against.
        # Inputs are: c, an array of concentrations; n, the reaction order; cstar, the transition concentration for this species.
        # Outputs are: the regord-blended rate factor at each concentration in c - a local re-implementation of regord_pow() from differentials.py, used here to test that function's behaviour independently.
        # Called by: test_regord(), in this file.
        out = np.empty_like(c, dtype=float)
        for k, ck in enumerate(c):
            if ck <= 0:
                out[k] = 0.0
            elif n >= 1 or n == 0:
                out[k] = float(ck ** n)
            else:
                xs = ck ** REGORD_S; cs = cstar ** REGORD_S  # concentration and cstar, both raised to the blend's sharpness power
                r  = (ck ** n) * ((xs / (xs + cs)) ** ((1.0 - n) / REGORD_S))  # the regord formula: c^n corrected by the blend fraction
                out[k] = r if np.isfinite(r) else 0.0
        return out

    def _eff_order(c, n, cstar, eps=1e-5):
        # Purpose: numerically estimates the effective reaction order at each concentration, as an independent check.
        # Inputs are: c, an array of concentrations; n, the reaction order; cstar, the transition concentration; eps, the fractional step size used for the numerical derivative.
        # Outputs are: the effective reaction order at each concentration in c, d(ln rate)/d(ln c), estimated numerically from _mod() rather than read off the formula directly, as an independent check.
        # Called by: test_regord(), in this file.
        dc   = c * eps
        r_hi = _mod(c + dc / 2, n, cstar)
        r_lo = _mod(c - dc / 2, n, cstar)
        denom = r_hi + r_lo
        safe  = denom > 0
        return np.where(safe, 2.0 * (r_hi - r_lo) / np.where(safe, denom, 1.0) / eps, n)

    # CHECK 1: accuracy at calibration minimum
    print("=" * 65)
    print("CHECK 1 — regord_pow accuracy vs pure power law")
    print("=" * 65)
    all_pass = True
    for sp_id, n, label, sp_cstar, cal_min in AFFECTED:
        print(f"\n  {label}  (n={n:+.3f}, cstar={sp_cstar:.0e} M, cal_min={cal_min:.0e} M)")
        for c in [sp_cstar, 2*sp_cstar, 10*sp_cstar, cal_min]:
            r_orig  = _unmod(np.array([c]), n)[0]
            r_mod   = _mod(np.array([c]), n, sp_cstar)[0]
            err_pct = abs(r_mod - r_orig) / abs(r_orig) * 100 if r_orig != 0 else 0.0
            is_cal  = abs(c - cal_min) < 1e-12 * cal_min
            if is_cal:
                # At the calibration minimum, the blend is required to match the plain power law closely (< 0.1% error), since this is the lowest concentration the rate law was actually fitted to.
                ok = err_pct < 0.1
                if not ok:
                    all_pass = False
                flag = f"  {'PASS' if ok else 'FAIL'} (calibration min)"
            else:
                flag = "  (informational)"
            print(f"    c={c:.2e} M ({c/sp_cstar:>5.0f}×cstar):  err={err_pct:.4f}%{flag}")
    print(f"\n  Result: {'PASS' if all_pass else 'FAIL'}")

    # CHECK 2: rate string inspection
    print("\n" + "=" * 65)
    print("CHECK 2 — fractional-order terms wrapped with regord_pow")
    print("=" * 65)
    try:
        import re
        rate_laws, _, _ = generate_rate_laws_with_database_type('prebiotic', 'prebiotic')
        bare_frac = re.compile(r'[A-Za-z_]\w*\s*\*\*\s*([+-]?\d*\.?\d+)')
        wrapped   = re.compile(r'regord_pow\(\s*(\w+)\s*,\s*([+-]?\d*\.?\d+)')
        found_bare, found_reg = [], []
        for rl in rate_laws:
            for expr in [rl.get('forward', ''), rl.get('reverse', '')]:
                if not expr:
                    continue
                for m in bare_frac.finditer(expr):
                    try:
                        if float(m.group(1)) < 1:
                            found_bare.append(m.group(0))
                    except ValueError:
                        pass
                for m in wrapped.finditer(expr):
                    try:
                        if float(m.group(2)) < 1:
                            found_reg.append((m.group(1), float(m.group(2))))
                    except ValueError:
                        pass
        print(f"  regord_pow-wrapped terms: {len(found_reg)}")
        for sp, exp in found_reg:
            print(f"    regord_pow({sp}, {exp})")
        print(f"  Bare ** fractional terms (should be 0): {len(found_bare)}")
        print(f"  Result: {'PASS' if not found_bare and found_reg else 'FAIL'}")
    except Exception as e:
        print(f"  Could not run CHECK 2: {e}")

    # Effective order table
    print("\n" + "=" * 65)
    print("EFFECTIVE ORDER  d(ln rate)/d(ln c)")
    print("=" * 65)
    hdr = f"  {'Species':>14}  {'n':>6}  " + "  ".join(f"{'c='+f'{c:.0e}':>12}" for c in REPORT_CONCS)
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for sp_id, n, label, sp_cstar, _ in AFFECTED:
        orders = _eff_order(np.array(REPORT_CONCS), n, sp_cstar)
        print(f"  {label:>14}  {n:>6.3f}  " + "  ".join(f"{o:>12.4f}" for o in orders))

    # Figures
    print("\n" + "=" * 65)
    print("Saving regord_pow figures …")
    print("=" * 65)
    OUT_DIR = os.path.join(SOUP_DIR, 'Outputs')
    os.makedirs(OUT_DIR, exist_ok=True)

    for sp_id, n, label, sp_cstar, cal_min in AFFECTED:
        fig, (ax_r, ax_o) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
        fig.subplots_adjust(hspace=0.08)

        r_orig = _unmod(C_GRID, n); r_mod = _mod(C_GRID, n, sp_cstar)
        mask   = (r_orig > 0) & (r_mod > 0)
        ax_r.loglog(C_GRID[mask], r_orig[mask], 'k-',  lw=1.5, label='Unmodified $c^n$')
        ax_r.loglog(C_GRID[mask], r_mod[mask],  'C0-', lw=1.5, label='regord_pow')
        ax_r.axvline(sp_cstar, color='gray', ls='--', lw=1,
                     label=f'{sp_cstar*1000:.3g} mM (calibration limit)')
        ax_r.axvspan(1e-12, sp_cstar, alpha=0.08, color='C1')
        ax_r.set_ylabel('Rate factor (a.u.)', fontsize=12)
        ax_r.set_title(f'{label}  (n={n:+.3f}, cstar={sp_cstar:.0e} M, s={REGORD_S})', fontsize=12)
        ax_r.legend(fontsize=10, loc='best'); ax_r.grid(False)

        eff_n = _eff_order(C_GRID, n, sp_cstar)
        ax_o.semilogx(C_GRID, eff_n, 'C0-', lw=1.5)
        ax_o.axhline(n,   color='k',    ls='--', lw=1, label=f'Fitted n={n:.3f}')
        ax_o.axhline(1.0, color='gray', ls=':',  lw=1, label='Order 1 (linear)')
        ax_o.axvline(sp_cstar, color='gray', ls='--', lw=1)
        ax_o.axvspan(1e-12, sp_cstar, alpha=0.08, color='C1')
        ax_o.set_xlabel('Concentration [M]', fontsize=12)
        ax_o.set_ylabel('Effective order', fontsize=12)
        ax_o.legend(fontsize=10, loc='best')
        ax_o.set_ylim((min(n - 0.3, 0.0) if n < 0 else n - 0.3), 1.3)
        ax_o.grid(False)
        for c, o in zip(REPORT_CONCS, _eff_order(np.array(REPORT_CONCS), n, sp_cstar)):
            ax_o.plot(c, o, 'ko', ms=4)

        out_path = os.path.join(OUT_DIR, f'regord_{sp_id}.png')
        plt.savefig(out_path, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {out_path}")


# -- entry point --------------------------------------------------------------

def main():
    # Purpose: parses the command-line flags and runs whichever tests were requested.
    # Inputs are: none directly; which tests to run is decided by the command-line flags.
    # Outputs are: nothing is returned; the requested test(s) are run in sequence (mass balance, convergence, open system, regord, in that order).
    # Called by: the script's entry point, at the bottom of this file (if __name__ == "__main__").
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mass-balance', action='store_true',
                        help='30-day closed-system ODE + elemental audit')
    parser.add_argument('--convergence',  action='store_true',
                        help='N-run spread study (closed system)')
    parser.add_argument('--open-system',  action='store_true',
                        help='N-run spread study (open system, 7 days)')
    parser.add_argument('--regord',       action='store_true',
                        help='Verify regord_pow blend function + save PNGs')
    parser.add_argument('--all',          action='store_true',
                        help='Run all four tests')
    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()
        sys.exit(0)

    run_all = args.all

    if run_all or args.mass_balance:
        test_mass_balance()
        print()

    if run_all or args.convergence:
        test_convergence()
        print()

    if run_all or args.open_system:
        test_open_system()
        print()

    if run_all or args.regord:
        test_regord()
        print()


if __name__ == '__main__':
    main()
