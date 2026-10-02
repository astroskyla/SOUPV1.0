# This module runs automated diagnostic checks after a simulation finishes.
# It is called automatically by main.py after each run, and can also be run standalone.
# Each species' atom counts are worked out from its formula-style name; the only manual overrides needed are for metal complexes where the number of CN groups is ambiguous from the name alone (e.g. FeCN6-4 parsed naively gives C=1, not C=6).

import os
import re
import sys
import numpy as np


# ---------------------------------------------------------------------------
# Composition derivation
# ---------------------------------------------------------------------------

ELEMENTS_TO_AUDIT = ['C', 'N', 'S', 'Fe']
# H is excluded, since the reactions are pseudo-steps and H does not balance mechanistically.

# These are overrides for species whose formula name is ambiguous for simple parsing.
# FeCN6 parses naively as Fe=1, C=1, N=6, but it should be Fe=1, C=6, N=6, since it has six cyanide ligands.
COMPOSITION_OVERRIDES = {
    'FeCN6-4':    {'C': 6, 'N': 6, 'Fe': 1},
    'FeCN5H2O-3': {'C': 5, 'N': 5, 'Fe': 1},
}


def _parse_formula(name):
    # Purpose: counts each element's atoms in a formula-style species name.
    # Inputs are: name, a formula-style species name, e.g. 'CaHCO3+1'.
    # Outputs are: a dictionary mapping each element symbol to how many atoms of it are in the formula, worked out by stripping the trailing charge notation (e.g. '-1', '+2') and then scanning for element symbols followed by an optional count.
    # Called by: get_composition(), in this file.
    formula = re.sub(r'[+-]\d+$', '', name.strip())
    atoms = {}
    for elem, count in re.findall(r'([A-Z][a-z]?)(\d*)', formula):
        atoms[elem] = atoms.get(elem, 0) + (int(count) if count else 1)
    return atoms


def get_composition(species_name):
    # Purpose: returns a species' elemental composition, using a manual override if one exists.
    # Inputs are: species_name, a species name, either formula-style or one of the manual overrides above.
    # Outputs are: a dictionary mapping each element symbol to its atom count, from COMPOSITION_OVERRIDES if the species is listed there, otherwise parsed directly from its name.
    # Called by: check_mass_balance(), in this file (twice); _element_totals(), in tests/prebiotic_tests.py.
    if species_name in COMPOSITION_OVERRIDES:
        return COMPOSITION_OVERRIDES[species_name]
    return _parse_formula(species_name)


# ---------------------------------------------------------------------------
# Markdown table parser (shared by all checks)
# ---------------------------------------------------------------------------

def _parse_md_table(path):
    # Purpose: reads a saved markdown results table into headers and a numeric array.
    # Inputs are: path, the path to a markdown file containing a pipe-delimited table (as saved by main.py).
    # Outputs are: headers, the list of column names; data, a numpy array of every numeric row. Raises a ValueError if no table could be parsed from the file.
    # Called by: check_mass_balance(), in this file.
    headers = None
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith('|'):
                continue
            cells = [c.strip() for c in line.split('|')[1:-1]]
            if headers is None:
                headers = cells
                continue
            if set(cells[0]) <= {'-'}:
                continue
            try:
                rows.append([float(c) for c in cells])
            except ValueError:
                continue
    if headers is None or not rows:
        raise ValueError(f"Could not parse table from {path}")
    return headers, np.array(rows)


# ---------------------------------------------------------------------------
# CHECK: elemental mass balance
# ---------------------------------------------------------------------------

def check_mass_balance(md_path, threshold=1e-6):
    # Purpose: checks that C, N, S and Fe are all conserved across the run, within a small tolerance.
    # Inputs are: md_path, the path to the saved species_concentrations.md results table; threshold, the maximum fractional deviation from the starting total that still counts as conserved (default one part per million).
    # Outputs are: True if every audited element (C, N, S, Fe) stays conserved within threshold across the whole run, False otherwise. A formatted table of results is also printed.
    # Called by: run_all(), in this file.
    print("=" * 65)
    print("Mass balance check — elemental conservation")
    print(f"  Source: {md_path}")
    print("=" * 65)

    if not os.path.exists(md_path):
        print(f"  ERROR: output file not found — run the simulation first.")
        return False

    headers, data = _parse_md_table(md_path)
    times = data[:, 0]
    species_headers = headers[1:]   # skip 'time' column
    conc_data = data[:, 1:]

    print(f"\n  Time span : {times[0]:.2e} → {times[-1]:.2e} s  "
          f"({times[-1]/86400:.1f} days)")
    print(f"  Data pts  : {len(times)}")
    print(f"  Species   : {len(species_headers)}")

    # Every species' contribution to each audited element's total is added up here, across every timepoint.
    totals = {el: np.zeros(len(times)) for el in ELEMENTS_TO_AUDIT}
    for j, sp in enumerate(species_headers):
        comp = get_composition(sp)
        for el in ELEMENTS_TO_AUDIT:
            cnt = comp.get(el, 0)
            if cnt:
                totals[el] += cnt * conc_data[:, j]

    print(f"\n  {'Element':>7}  {'t=0 [M]':>14}  {'t_final [M]':>14}  "
          f"{'max |Δ|/total₀':>16}  Result")
    print("  " + "-" * 64)

    all_pass = True
    for el in ELEMENTS_TO_AUDIT:
        t0 = totals[el][0]
        if t0 == 0:
            print(f"  {el:>7}  {'(zero at t=0)':>14}  — skipped")
            continue
        deviation = np.abs(totals[el] - t0) / t0
        max_dev = deviation.max()
        worst_t = times[deviation.argmax()]
        ok = max_dev < threshold
        if not ok:
            all_pass = False
        flag = "PASS" if ok else "FAIL"
        print(f"  {el:>7}  {t0:>14.6e}  {totals[el][-1]:>14.6e}  "
              f"{max_dev:>16.3e}  {flag}  (worst t={worst_t/3600:.2f} h)")

    print(f"\n  Note: H excluded — reactions are pseudo-steps; H does not balance.")
    status = "PASS — all elements conserved to < 1 ppm" if all_pass else "FAIL — see above"
    print(f"  Overall: {status}\n")
    return all_pass


# ---------------------------------------------------------------------------
# Entry point: run all checks
# ---------------------------------------------------------------------------

def run_all(md_path):
    # Purpose: runs every check in this file and reports whether they all passed.
    # Inputs are: md_path, the path to the saved species_concentrations.md results table.
    # Outputs are: True only if every check below passes, False if any fail.
    # Called by: main(), in main.py (imported there as _run_checks); the test block at the bottom of this file (if __name__ == "__main__").
    results = []
    results.append(check_mass_balance(md_path))
    # Future checks can be added here, in the same way, e.g.:
    #   results.append(check_charge_balance(md_path))
    return all(results)


if __name__ == '__main__':
    # Standalone usage: python3 plotting/checks.py [path/to/species_concentrations.md]
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(script_dir, '..', 'Outputs', 'species_concentrations.md')

    ok = run_all(path)
    sys.exit(0 if ok else 1)
