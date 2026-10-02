# This module builds and solves the differential equations for the simulation.
# Three ODE systems are provided: a legacy system (high precision, for equilibrium chemistry), a simple system (a clean BDF approach, for Fenton chemistry), and a bulletproof system (robust handling, for prebiotic chemistry and other reaction networks).

import re
import os
import numpy as np
from datetime import datetime
from rate_laws import query_rate_laws_by_species
import sys
import warnings


def format_species_name(name):
    # Purpose: converts a species name from the internal database format to the display format.
    # Inputs are: name, a species name in the internal database format, e.g. 'HPLUS1'.
    # Outputs are: the same name converted to the display format, e.g. 'H+1'.
    # Called by: safe_derivative_calculation() and replace_y_with_species(), both in this file; build_initial_conditions() and save_md(), both in run_sai_miyakawa_networks.py.
    return name.replace("PLUS", "+").replace("MINUS", "-")


def construct_differential_equations(rate_laws):
    # Purpose: builds the differential equation for every species, from the reactions that produce or consume it.
    # Inputs are: rate_laws, the list of rate-law dictionaries produced by generate_rate_laws_with_database_type() in rate_laws.py.
    # Outputs are: differential_equations, a dictionary mapping each species to its differential equation as a string (e.g. "dFeII_dt = ..."), built by combining every reaction that produces or consumes it. A handful of species that are treated as constant (water, H+, OH-, OH, H) are given a fixed derivative of zero.
    # Called by: main() in main.py; run_full_diagnostic() in diagnose.py; run_one() in run_sai_miyakawa_networks.py; _do() in tests/prebiotic_tests.py; run_zymonic_simulation() in inverse_kinetics.py and seeded_inverse_kinetics.py; run_fe_cn_simulation() in Photoaquation/photoaquation_checker_simplified.py, Photoaquation/photoaquation_checker_simplified-V2.py and Photoaquation/photoaquation_rate.py; run_simulation() in Photoaquation/fe-cn-inverse-checker.py and Photoaquation/fecn_inversekinetics.py; and directly from module-level code (not inside a function) in prebiotic_scenario_runner.py.
    unique_species = set()
    for rate_law in rate_laws:
        for species, _ in rate_law['consumed']:
            unique_species.add(species)
        for species, _ in rate_law['produced']:
            unique_species.add(species)

    differential_equations = {}
    for species in unique_species:
        producing_rate_laws, consuming_rate_laws = query_rate_laws_by_species(rate_laws, species)
        equation_terms = []

        for rate_law in consuming_rate_laws:
            for s, coefficient in rate_law['consumed']:
                if s == species:
                    if rate_law['forward'] is not None:
                        equation_terms.append(f"(-{abs(coefficient)}) * ({rate_law['forward']})")
                    if rate_law['reverse']:
                        equation_terms.append(f"(+{abs(coefficient)}) * ({rate_law['reverse']})")

        for rate_law in producing_rate_laws:
            for s, coefficient in rate_law['produced']:
                if s == species:
                    if rate_law['forward'] is not None:
                        equation_terms.append(f"(+{abs(coefficient)}) * ({rate_law['forward']})")
                    if rate_law['reverse']:
                        equation_terms.append(f"(-{abs(coefficient)}) * ({rate_law['reverse']})")

        if equation_terms:
            differential_equations[species] = f"d{species}_dt = " + ' + '.join(equation_terms)

    # These species are held constant, so their derivative is fixed at zero rather than being built from reactions. This can be toggled.
    differential_equations['H2O'] = 'dH2O_dt = 0'
    differential_equations['HPLUS1'] = 'dHPLUS1_dt = 0'
    differential_equations['OHMINUS1'] = 'dOHMINUS1_dt = 0'
    differential_equations['OH'] = 'dOH_dt = 0'
    differential_equations['H'] = 'dH_dt = 0'

    return differential_equations


# =============================================================================
# LEGACY SYSTEM - High precision for equilibrium chemistry
# =============================================================================

def safe_pow(x, p):
    # Purpose: raises a number to a power without crashing on a negative or zero base.
    # Inputs are: x, a number; p, the power to raise it to.
    # Outputs are: x raised to the power p, or 0.0 if x is zero, negative, or the calculation fails - this avoids errors from raising a negative or zero concentration to a fractional power.
    # Called by: from inside the rate-law expressions themselves, once replace_powers_with_safe_pow() has rewritten every "**" in them as a safe_pow(...) call; those expressions are then run by ode_system() via eval(), with safe_pow supplied as part of eval()'s namespace, all in this file.
    try:
        if x <= 0:
            return 0.0
        return x ** p
    except Exception:
        return 0.0


def replace_powers_with_safe_pow(expr):
    # Purpose: rewrites every "**" power operation in a rate-law expression to go through safe_pow() instead.
    # Inputs are: expr, a rate-law expression as a string, e.g. containing a term like "FeII**2".
    # Outputs are: the same expression with every "**" power operation rewritten as a call to safe_pow(...), e.g. "safe_pow(FeII, 2)".
    # Called by: ode_system(), in this file.
    pattern = re.compile(r'([a-zA-Z_][a-zA-Z0-9_]*)\s*\*\*\s*([-+]?[0-9]*\.?[0-9]+)')  # matches a species name followed by **exponent, e.g. "FeII**2"

    def replacer(match):
        # Purpose: turns one regex match (a species and its exponent) into a safe_pow(...) call.
        # Inputs are: match, a regex match object with the species name as group 1 and the exponent as group 2.
        # Outputs are: the replacement text, e.g. "safe_pow(FeII, 2)".
        # Called by: pattern.sub(), just below, once per "**" match found in expr.
        base = match.group(1)
        exponent = match.group(2)
        return f"safe_pow({base}, {exponent})"

    return pattern.sub(replacer, expr)


def update_concentrations(y, database_species, species_atm, partial_pressures, henrys_constants, pH):
    # Purpose: applies Henry's law and iron mass balance (for fenton), and cleans up negative concentrations, before each rate calculation.
    # Inputs are: y, the current concentration of every species; database_species, the list of species names; species_atm, partial_pressures and henrys_constants, the atmospheric data used for Henry's law; pH, the solution pH (not used directly here, kept only for a consistent function signature).
    # Outputs are: y, the same array with atmospheric species topped up via Henry's law, the FeIII concentration recalculated from iron mass balance, and any negative values reset to zero.
    # Called by: ode_system(), in this file. It is also imported by main.py, but is not actually called there. main.py uses its own local Fenton-specific version instead.
    for i, species in enumerate(database_species):
        # A species that also exists in the atmosphere is topped up to at least the concentration Henry's law predicts for it.
        if species in species_atm:
            index = species_atm.index(species)
            if henrys_constants[index] is not None:
                equilibrium_concentration = partial_pressures[index] * float(henrys_constants[index])
                if y[i] < equilibrium_concentration:
                    y[i] = equilibrium_concentration

        # Iron(III) is not integrated directly; instead its concentration is recalculated here from mass balance, as the total iron added (0.0002 M) minus whatever has already become intermediates I1/I2 or iron(II).
        if species == 'FeIII':
            try:
                fe_ii_conc = y[database_species.index('FeII')] if 'FeII' in database_species else 0.0
                i1_conc = y[database_species.index('I1')] if 'I1' in database_species else 0.0
                i2_conc = y[database_species.index('I2')] if 'I2' in database_species else 0.0
                y[i] = max(0.0, 0.0002 - (i1_conc + i2_conc) - fe_ii_conc)  # total iron (0.0002 M) minus what has already become I1, I2, or Fe(II)
            except (ValueError, IndexError):
                pass

    y[y < 0] = 0.0
    return y


def create_eval_context(database_species, y):
    # Purpose: builds the variable namespace a rate-law expression needs to be evaluated.
    # Inputs are: database_species, the list of species names; y, their current concentrations.
    # Outputs are: a dictionary mapping each species name to its concentration in y, used as the variable namespace when a rate-law expression is evaluated.
    # Called by: ode_system(), in this file.
    return {species: y[i] for i, species in enumerate(database_species)}


def ode_system(t, y, database_species, differential_equations, species_atm, partial_pressures,
               henrys_constants, pH):
    # Purpose: this is the ODE function itself for the legacy (equilibrium-chemistry) system. It calculates every species' rate of change at time t.
    # Inputs are: t, the current simulation time; y, the current concentration of every species; database_species, the list of species names; differential_equations, the equations built by construct_differential_equations(); species_atm, partial_pressures, henrys_constants and pH, the atmospheric/solution data needed to update concentrations before the rates are calculated.
    # Outputs are: dy_dt, the rate of change of every species at this time and concentration. This is the legacy, high-precision system used for equilibrium chemistry.
    # Called by: main(), in main.py; test_ode_evaluation(), in diagnose.py.
    y = update_concentrations(y, database_species, species_atm, partial_pressures, henrys_constants, pH)

    dy_dt = [0] * len(database_species)

    for i, species in enumerate(database_species):
        if species in differential_equations:
            right = differential_equations[species].split('=')[1]  # keeps only the right-hand side of "dspecies_dt = ..."
            safe_right = replace_powers_with_safe_pow(right)
            context = create_eval_context(database_species, y)
            try:
                dy_dt[i] = eval(safe_right, {"safe_pow": safe_pow}, context)
            except Exception as e:
                print(f"Error in species '{species}' at t={t}: {right}")
                raise e

    return dy_dt


# =============================================================================
# BULLETPROOF SYSTEM - Robust handling for prebiotic chemistry
# =============================================================================

def safe_pow_bulletproof(x, p):
    # Purpose: raises a number to a power without crashing on a bad base, more defensively than safe_pow().
    # Inputs are: x, a number; p, the power to raise it to.
    # Outputs are: x raised to the power p; 0.0 (or 1.0, if p is not positive) if x is zero or negative, or if the result is NaN/infinite, or if the calculation fails. This is a more defensive version of safe_pow(), used by the bulletproof system.
    # Called by: from inside the rate-law expressions themselves, once replace_powers_with_safe_pow_bulletproof() has rewritten every "**" in them as a safe_pow(...) call; those expressions are then run by UniversalRobustODESystem.safe_derivative_calculation() via eval(), in this file.
    if x <= 0:
        return 0.0 if p > 0 else 1.0
    try:
        result = x ** p
        if np.isnan(result) or np.isinf(result):
            return 0.0 if p > 0 else 1.0
        return result
    except:
        return 0.0 if p > 0 else 1.0


# ---------------------------------------------------------------------------
# Concentration-dependent reaction order (regord)
# ---------------------------------------------------------------------------
# Some reactions have an empirical rate order n below 1, including negative values, which is only valid where it was actually measured.
# Taken literally at low concentration, this makes the rate law non-Lipschitz (its slope goes to infinity), which the solver cannot integrate through reliably.
# The effective order is therefore blended from n, valid down to REGORD_CSTAR (the lowest concentration calibrated against experiment), to 1 (ordinary first-order decay) below that, so the rate falls cleanly to zero as a species is exhausted rather than diverging.

# Formula: c**n * (c**s / (c**s + cstar**s))**((1-n)/s)

# REGORD_CSTAR or REGORD_S can be changed here, and per-species overrides can be added to REGORD_SPECIES_CSTAR, e.g. {'CNMINUS1': 1e-4, 'SO3MINUS2': 1e-6}.

# cstar is set a hundred times below the lowest calibrated concentration, so the blend is essentially complete (<0.01% error) by the time any measured concentration is reached; it must never sit inside the calibrated range, since at c = cstar the blend is only halfway and the rate sits 10-40% below the true power law.

# Calibration minimums (from experiment):
#   CN-: 0.2 mM -> cstar = 2e-6 M (100x below)
#   Fe2+: 0.1 mM -> cstar = 1e-6 M (100x below)
#   SO32-: 1 mM -> cstar = 1e-5 M (100x below)
#   CO2: 1 mM -> cstar = 1e-5 M (100x below)
REGORD_CSTAR = 1e-5          # [M] default transition concentration (fallback)
REGORD_S = 2                 # sharpness of transition (higher = sharper)

REGORD_SPECIES_CSTAR = {
    'CNMINUS1':   2e-6,   # CN-:   calibrated >= 0.2 mM; cstar 100x below
    'FePLUS2':    1e-6,   # Fe2+:  calibrated >= 0.1 mM; cstar 100x below
    'SO3MINUS2':  1e-5,   # SO32-: calibrated >= 1 mM;   cstar 100x below
    'CO2':        1e-5,   # CO2:   calibrated >= 1 mM;   cstar 100x below
    'HCN':        1e-11,  # HCN:   non-negativity guard for reaction 10 in prebiotic.db (AMS, zero-order in HCN); 1e-10 converged, 1e-11 chosen for negligible suppression
}


def regord_pow(x, p, species=None):
    # Purpose: a regularised power function that keeps reaction rates well-behaved near zero concentration.
    # Inputs are: x, a concentration; p, the reaction order to raise it to; species, the species name (optional), used to look up a per-species cstar in REGORD_SPECIES_CSTAR.
    # Outputs are: x**p for p >= 1, unchanged. For p < 1 (and p != 0), the concentration-dependent blend described above is applied instead, so the effective order rises to 1 below cstar, keeping the rate Lipschitz (finite-sloped) at zero. For the special case p == 0 with a species that has a registered cstar, a smooth blend x^s / (x^s + cstar^s) is returned instead of a discontinuous jump from 0 to 1 at x = 0 - this keeps the Jacobian bounded and stops the BDF solver's step size collapsing in supply-limited zero-order reactions.
    # Called by: build_analytic_jacobian()'s jacobian() closure, in this file; and from inside rate-law expressions built by parse_reactants_products() in rate_laws.py, when those expressions are evaluated by UniversalRobustODESystem.safe_derivative_calculation(), also in this file.
    if x <= 0:
        return 0.0
    if p >= 1 or p == 0:
        if p == 0 and species and species in REGORD_SPECIES_CSTAR:
            cstar = REGORD_SPECIES_CSTAR[species]
            s = REGORD_S
            xs = x ** s
            cs = cstar ** s
            return xs / (xs + cs)  # smooth 0-to-1 ramp: near 0 well below cstar, near 1 well above it
        try:
            r = x ** p
            return r if np.isfinite(r) else (0.0 if p > 0 else 1.0)
        except Exception:
            return 0.0 if p > 0 else 1.0
    cstar = REGORD_SPECIES_CSTAR.get(species, REGORD_CSTAR) if species else REGORD_CSTAR
    s = REGORD_S
    xs = x ** s
    cs = cstar ** s
    blend = xs / (xs + cs)  # 0 well below cstar, 1 well above it. This is what smoothly switches the exponent from n to 1
    try:
        r = (x ** p) * (blend ** ((1.0 - p) / s))  # the formula from the block comment above: x^p, corrected by the blend so the effective order rises to 1 as x falls below cstar
        return r if np.isfinite(r) else 0.0
    except Exception:
        return 0.0


def replace_powers_with_safe_pow_bulletproof(expr):
    # Purpose: rewrites every "**" power operation in a rate-law expression to go through safe_pow_bulletproof() instead.
    # Inputs are: expr, a rate-law expression as a string, e.g. containing a term like "FeII**2".
    # Outputs are: the same expression with every "**" power operation rewritten as a call to safe_pow(...), e.g. "safe_pow(FeII, 2)". This is the bulletproof-system equivalent of replace_powers_with_safe_pow().
    # Called by: UniversalRobustODESystem.safe_derivative_calculation(), in this file.
    pattern = re.compile(r'([a-zA-Z_][a-zA-Z0-9_]*)\s*\*\*\s*([-+]?[0-9]*\.?[0-9]+)')  # matches a species name followed by **exponent, e.g. "FeII**2"
    def replacer(match):
        # Purpose: turns one regex match (a species and its exponent) into a safe_pow(...) call.
        # Inputs are: match, a regex match object with the species name as group 1 and the exponent as group 2.
        # Outputs are: the replacement text, e.g. "safe_pow(FeII, 2)".
        # Called by: pattern.sub(), just below, once per "**" match found in expr.
        base = match.group(1)
        exponent = match.group(2)
        return f"safe_pow({base}, {exponent})"
    return pattern.sub(replacer, expr)


_FIXED_SPECIES = frozenset({'H2O', 'HPLUS1', 'OHMINUS1', 'OH', 'H'})


def _rate_law_exponents(rate_law_str, species_to_idx):
    # Purpose: reads the exponent each species appears with in a rate-law expression, whichever of the three ways it's written.
    # Inputs are: rate_law_str, a rate-law expression as a string; species_to_idx, a dictionary mapping every species name to its index.
    # Outputs are: a dictionary mapping each species that appears in the expression to its exponent. Three forms are recognised: regord_pow(SPECIES, EXP, ...), for fractional/negative reaction orders; SPECIES**EXP, an explicit exponent (integer or the older format); and a bare SPECIES, which implies an exponent of 1.
    # Called by: _compile_rate_law_terms(), in this file.
    exps = {}
    for s in species_to_idx:
        # A regord_pow() call is checked for first, since fractional exponents in prebiotic networks are written that way.
        m = re.search(r'regord_pow\(\s*' + re.escape(s) + r'\s*,\s*([+-]?\d+\.?\d*)', rate_law_str)  # looks for "regord_pow(species, exponent" and captures the exponent
        if m:
            exps[s] = float(m.group(1))
            continue
        m = re.search(r'\b' + re.escape(s) + r'\*\*([+-]?\d+\.?\d*)', rate_law_str)  # looks for "species**exponent" and captures the exponent
        if m:
            exps[s] = float(m.group(1))
        elif re.search(r'\b' + re.escape(s) + r'\b', rate_law_str):  # the species appears with no exponent at all, so it's implicitly to the power of 1
            exps[s] = 1.0
    return exps


def _compile_rate_law_terms(rate_laws, database_species):
    # Purpose: pre-parses every reaction's rate law into a fast-to-reuse form, ready for the Jacobian to be built from.
    # Inputs are: rate_laws, the list of rate-law dictionaries; database_species, the list of species names.
    # Outputs are: compiled, a list of (k, rate_reactants, contributions) tuples, one per forward or reverse reaction step, where k is the rate constant, rate_reactants is a list of (species_index, exponent) pairs read from the rate-law expression, and contributions is a list of (affected_species_index, signed_stoichiometric_coefficient) pairs for every non-fixed species that reaction step changes. This pre-compiled form is what build_analytic_jacobian() uses to build the Jacobian without re-parsing the rate-law strings on every call.
    # Called by: UniversalRobustODESystem.__init__(), in this file.
    species_to_idx = {s: i for i, s in enumerate(database_species)}
    compiled = []

    for rl in rate_laws:
        # Forward step: k_fwd * product of reactant^exponent.
        if rl['forward'] is not None and rl['consumed']:
            try:
                k_fwd = float(rl['forward'].split('*')[0])
            except (ValueError, IndexError):
                continue

            fwd_exps = _rate_law_exponents(rl['forward'], species_to_idx)
            reactants = [(species_to_idx[s], exp) for s, exp in fwd_exps.items()]

            contributions = []
            for s, coeff in rl['consumed']:
                if s in species_to_idx and s not in _FIXED_SPECIES:
                    contributions.append((species_to_idx[s], -abs(coeff)))
            for s, coeff in rl['produced']:
                if s in species_to_idx and s not in _FIXED_SPECIES:
                    contributions.append((species_to_idx[s], +abs(coeff)))

            if reactants and contributions:
                compiled.append((k_fwd, reactants, contributions))

        # Reverse step: k_rev * product of product^exponent.
        if rl['reverse'] is not None and rl['produced']:
            try:
                k_rev = float(rl['reverse'].split('*')[0])
            except (ValueError, IndexError):
                continue

            rev_exps = _rate_law_exponents(rl['reverse'], species_to_idx)
            reactants_rev = [(species_to_idx[s], exp) for s, exp in rev_exps.items()]

            contributions_rev = []
            for s, coeff in rl['consumed']:
                if s in species_to_idx and s not in _FIXED_SPECIES:
                    contributions_rev.append((species_to_idx[s], +abs(coeff)))
            for s, coeff in rl['produced']:
                if s in species_to_idx and s not in _FIXED_SPECIES:
                    contributions_rev.append((species_to_idx[s], -abs(coeff)))

            if reactants_rev and contributions_rev:
                compiled.append((k_rev, reactants_rev, contributions_rev))

    return compiled


def build_analytic_jacobian(compiled_terms, n_species, database_species=None):
    # Purpose: builds a function that calculates the Jacobian analytically, so the solver doesn't have to estimate it numerically.
    # Inputs are: compiled_terms, the list produced by _compile_rate_law_terms(); n_species, the total number of species; database_species, the list of species names (optional; used to look up each species' cstar for the regord blend).
    # Outputs are: jacobian, a function of (t, y) that returns the Jacobian matrix (the matrix of partial derivatives of each species' rate of change with respect to every other species' concentration). For an integer or super-linear exponent (p >= 1), the derivative is sign * (p / y_j) * rate_value. For a fractional/negative exponent (p < 1, the regord path), the derivative instead uses an effective order that accounts for the regord blend (derived from d/dy[y^(p+s)/(y^s+cstar^s)] = (p+s-s*blend)*f(y)/y). Supplying an analytic Jacobian like this lets the solver take larger, more confident steps than it could with a numerically estimated one.
    # Called by: UniversalRobustODESystem.__init__(), in this file.
    # Each species' cstar is looked up once here, when the Jacobian is built, so it always stays consistent with regord_pow().
    s = REGORD_S
    if database_species is not None:
        cstar_by_idx = [
            REGORD_SPECIES_CSTAR.get(sp, REGORD_CSTAR)
            for sp in database_species
        ]
    else:
        cstar_by_idx = None

    def jacobian(t, y):
        # Purpose: calculates the actual Jacobian matrix for the current concentrations.
        # Inputs are: t, the current simulation time (accepted for a consistent signature with solve_ivp's jac callable, not used directly here); y, the current concentration of every species.
        # Outputs are: J, the n_species x n_species matrix of partial derivatives, built up term by term from compiled_terms.
        # Called by: automatically, whenever the solver needs the Jacobian - this is the function returned by build_analytic_jacobian(), and is set as the .jacobian attribute used throughout this file and in solve_with_bulletproof_method().
        J = np.zeros((n_species, n_species))
        for k, reactants, contributions in compiled_terms:
            rate_value = k
            for idx, exp in reactants:
                rate_value *= regord_pow(max(y[idx], 0.0), exp)  # builds up the reaction rate, one reactant's concentration term at a time

            for i, sign_coeff in contributions:
                for j, exp_j in reactants:
                    yj = y[j]
                    if yj <= 0 or rate_value == 0:
                        continue
                    if exp_j >= 1 or exp_j == 0:
                        J[i, j] += sign_coeff * (exp_j / yj) * rate_value  # ordinary power-law derivative: d(y^n)/dy = n*y^(n-1), written as (n/y) * y^n
                    else:
                        cstar = cstar_by_idx[j] if cstar_by_idx is not None else REGORD_CSTAR
                        xs = yj ** s  # y raised to the blend's sharpness power s
                        blend = xs / (xs + cstar ** s)  # 0 near y = 0, rising to 1 well above cstar
                        # This is the correct derivative of y^n * y^s/(y^s+cstar^s): d/dy = (n + s*(1-blend)) * regord_pow(y)/y.
                        eff_order = exp_j + s * (1.0 - blend)  # the reaction order effectively rises toward 1 as y falls below cstar
                        J[i, j] += sign_coeff * (eff_order / yj) * rate_value  # same derivative form as above, but using this blended order
        return J
    return jacobian


class UniversalRobustODESystem:
    # This is the "bulletproof" ODE system: it handles any database type not specifically recognised as equilibrium or Fenton, and is built to prevent negative concentrations under a wide range of conditions.
    # Instantiated by: create_bulletproof_ode_system(), in this file.
    def __init__(self, database_species, differential_equations,
                 species_atm, partial_pressures, henrys_constants, pH, rate_laws=None):
        # Purpose: sets up a new bulletproof ODE system, ready to be handed to a solver.
        # Inputs are: database_species, differential_equations, species_atm, partial_pressures, henrys_constants and pH, the same simulation setup used throughout this file; rate_laws, the rate-law list (optional), used to enable reaction-rate logging and to build the analytic Jacobian.
        # Outputs are: none directly; this sets up the object's state (Henry's law floors, reaction-rate logging, and the analytic Jacobian, if available).
        # Called by: automatically, when create_bulletproof_ode_system() creates a new UniversalRobustODESystem, in this file.
        self.database_species = database_species
        self.differential_equations = differential_equations
        self.species_atm = species_atm
        self.partial_pressures = partial_pressures
        self.henrys_constants = henrys_constants
        self.pH = pH

        # Henry's law equilibrium concentrations are pre-computed once here, rather than recalculated on every timestep.
        self.henry_equilibrium = {}
        for i, species in enumerate(species_atm):
            if henrys_constants[i] is not None and henrys_constants[i] != 'N/A':
                try:
                    self.henry_equilibrium[species] = partial_pressures[i] * float(henrys_constants[i])
                    print(f"Henry's Law: {species} equilibrium = {self.henry_equilibrium[species]:.2e} M")
                except (ValueError, TypeError):
                    pass

        self.species_to_index = {species: i for i, species in enumerate(database_species)}

        # Every species gets a floor concentration: the Henry's law equilibrium value for an atmospheric species, or zero otherwise. The derivative is adjusted (see __call__ below) whenever a species is at or below its floor and still falling.
        self.min_concentrations = {}
        for i, species in enumerate(database_species):
            if species in self.henry_equilibrium:
                self.min_concentrations[i] = self.henry_equilibrium[species]
            else:
                self.min_concentrations[i] = 0.0

        # Reaction-rate logging is set up here, if rate laws were supplied and the logging module is available.
        self.reaction_logger = None
        self.total_time = None
        self.integration_start_time = None

        if rate_laws is not None:
            try:
                from reaction_rates_logger import ReactionRatesLogger
                self.reaction_logger = ReactionRatesLogger(
                    database_species, differential_equations, rate_laws
                )
                print("   - Reaction rate logging enabled")
            except ImportError:
                print("   - Reaction rate logging not available")

        # An analytic Jacobian is only built for prebiotic and other unknown database types, since that is where it is most needed and rate_laws is available.
        self.jacobian = None
        self._compiled_terms = []
        if rate_laws is not None:
            try:
                compiled = _compile_rate_law_terms(rate_laws, database_species)
                self._compiled_terms = compiled
                _raw_jac = build_analytic_jacobian(compiled, len(database_species), database_species)
                _K_HENRY = 1e-3
                _floors = self.min_concentrations

                def jacobian(t, y, _rj=_raw_jac, _fl=_floors, _k=_K_HENRY):
                    # Purpose: wraps the raw analytic Jacobian to also account for the Henry's law replenishment term.
                    # Inputs are: t, the current simulation time; y, the current concentration of every species; _rj, _fl and _k, the raw Jacobian function, floor concentrations, and Henry's law rate constant, captured as default arguments so they can't change after this closure is built.
                    # Outputs are: J, the Jacobian from _rj(t, y), with an extra diagonal correction added for every species currently below its Henry's law floor.
                    # Called by: automatically, whenever the solver needs the Jacobian - this becomes self.jacobian, used throughout this file and in solve_with_bulletproof_method().
                    J = _rj(t, y)
                    # The Henry's law replenishment term (see __call__ below) also contributes to the Jacobian's diagonal: d/dy_i [k_H*(floor-y_i)] = -k_H, whenever a species is below its floor.
                    for i, floor in _fl.items():
                        if floor > 0 and y[i] < floor:
                            J[i, i] -= _k
                    return J

                self.jacobian = jacobian
                print(f"   - Analytic Jacobian compiled ({len(compiled)} terms)")
            except Exception as e:
                print(f"   - Analytic Jacobian compilation failed: {e}")

        # This tracks which error messages have already been shown, so the same warning is not printed on every timestep.
        self.warning_count = {}

    def apply_hard_caps(self, y):
        # Purpose: floors every concentration at zero. With new solver system this is not really used anymore... 
        # Inputs are: y, a concentration array.
        # Outputs are: y, with every value floored at 0.0, to prevent a negative concentration.
        # Called by: not called anywhere in the repo (defined but currently unused).
        y_capped = np.array(y, dtype=float)
        np.maximum(y_capped, 0.0, out=y_capped)
        return y_capped

    def safe_derivative_calculation(self, y_safe):
        # Purpose: calculates how fast every species is currently changing.
        # Inputs are: y_safe, the current concentration of every species.
        # Outputs are: dy_dt, the rate of change of every species, calculated by evaluating each species' differential equation; if evaluating a particular species' equation fails or produces NaN/infinity, that species' derivative is set to 0.0 instead, and an error is printed (only the first two times, to avoid flooding the console).
        # Called by: __call__(), in this file.
        dy_dt = np.zeros_like(y_safe)

        context = {species: y_safe[i] for i, species in enumerate(self.database_species)}
        _eval_globals = {"safe_pow": safe_pow_bulletproof, "regord_pow": regord_pow}

        for i, species in enumerate(self.database_species):
            if species in self.differential_equations:
                equation = self.differential_equations[species]

                try:
                    _, right_side = equation.split('=', 1)
                    right_side = right_side.strip()

                    if right_side == '0':
                        dy_dt[i] = 0.0
                        continue

                    # replace_powers_with_safe_pow_bulletproof() wraps any remaining plain "**" (integer exponents); any regord_pow(...) calls already present in the rate law string are left as they are.
                    safe_expr = replace_powers_with_safe_pow_bulletproof(right_side)
                    derivative = eval(safe_expr, _eval_globals, context)

                    if np.isnan(derivative) or np.isinf(derivative):
                        derivative = 0.0

                    dy_dt[i] = derivative

                except Exception as e:
                    dy_dt[i] = 0.0
                    error_key = f"error_{species}"
                    if error_key not in self.warning_count:
                        self.warning_count[error_key] = 0
                    if self.warning_count[error_key] < 2:
                        formatted_name = format_species_name(species)
                        print(f"  Error in {formatted_name}: {str(e)[:50]}...")
                        self.warning_count[error_key] += 1

        return dy_dt

    def apply_physical_constraints(self, y, dy_dt):
        # Purpose: stops a species going any lower once it's already at its floor concentration. One again, not really used anymore due to new solving technique...
        # Inputs are: y, the current concentration of every species; dy_dt, their calculated rate of change.
        # Outputs are: dy_dt, with the derivative set to zero for any species that is at (or within 1% of) its floor concentration and still decreasing.
        # Called by: not called anywhere in the repo (defined but currently unused).
        dy_dt_constrained = np.array(dy_dt)

        for i in range(len(y)):
            current_conc = y[i]
            current_deriv = dy_dt_constrained[i]
            min_conc = self.min_concentrations.get(i, 0.0)

            if current_conc <= min_conc * 1.01 and current_deriv < 0:
                dy_dt_constrained[i] = 0.0

        return dy_dt_constrained

    def __call__(self, t, y):
        # Purpose: this is the ODE function itself for the bulletproof system i.e., it's what the solver calls at every step.
        # Inputs are: t, the current simulation time; y, the current concentration of every species.
        # Outputs are: dy_dt, the rate of change of every species, including a Henry's law replenishment term for any atmospheric species below its floor.
        # Called by: automatically, whenever this object is passed as the right-hand-side function to scipy's solve_ivp() (e.g. inside solve_with_bulletproof_method(), in this file) or to a related solver elsewhere in the repo.
        if self.integration_start_time is None:
            self.integration_start_time = t

        dy_dt = self.safe_derivative_calculation(y)

        # For an atmospheric species below its Henry's law floor, a smooth first-order dissolution flux is added, proportional to how far below the floor it is. This replaces an earlier hard clamp (forcing the derivative to exactly zero), which created a discontinuity that caused the BDF/Radau solvers' step size to collapse; a continuous right-hand side lets the solver extrapolate through the floor without rejecting steps. The rate constant of 1e-3 per second gives equilibration in around 17 minutes - fast relative to the chemistry (hours to days), but not so fast that it makes the system stiffer.
        _K_HENRY = 1e-3
        for i, floor in self.min_concentrations.items():
            if floor > 0:
                deficit = floor - y[i]
                if deficit > 0:
                    dy_dt[i] += _K_HENRY * deficit

        if self.reaction_logger and self.total_time:
            if self.reaction_logger.should_log_time(t - self.integration_start_time, self.total_time):
                self.reaction_logger.log_timepoint(t - self.integration_start_time, y, dy_dt)

        return dy_dt

    def set_total_time(self, total_time):
        # Purpose: tells the logger how long the whole run will be, so it can space out its logging points.
        # Inputs are: total_time, the total length of the integration in seconds, used so the logger knows how to space out its logging points.
        # Outputs are: nothing is returned; self.total_time is updated, and a message is printed if logging is enabled.
        # Called by: callers of create_bulletproof_ode_system() that switch on logging, e.g. main() in main.py, via bulletproof_ode_system.set_total_time(duration).
        self.total_time = total_time
        if self.reaction_logger:
            print(f"   - Universal robust logging: {total_time/3600:.1f} hours, 50 timepoints")

    def finalise_logging(self):
        # Purpose: writes out the reaction-rate log, if logging was switched on.
        # Inputs are: none.
        # Outputs are: nothing is returned; if logging was enabled, the reaction-rate log is written out to an HTML file.
        # Called by: callers of create_bulletproof_ode_system() that switch on logging, e.g. main() in main.py, via bulletproof_ode_system.finalise_logging().
        if self.reaction_logger:
            self.reaction_logger.finalise_html()


def create_bulletproof_ode_system(database_species, differential_equations,
                                species_atm, partial_pressures, henrys_constants, pH, rate_laws=None):
    # Purpose: creates a new bulletproof ODE system for the given simulation setup.
    # Inputs are: the same simulation setup used throughout this file (database_species, differential_equations, species_atm, partial_pressures, henrys_constants, pH), plus rate_laws (optional).
    # Outputs are: a new UniversalRobustODESystem instance, ready to be passed to a solver as the right-hand-side function.
    # Called by: main() in main.py (twice); run_zymonic_simulation() in inverse_kinetics.py and seeded_inverse_kinetics.py; run_one() in run_sai_miyakawa_networks.py; test_mass_balance(), test_convergence() and test_open_system() in tests/prebiotic_tests.py; run_fe_cn_simulation() in Photoaquation/photoaquation_checker_simplified.py, Photoaquation/photoaquation_checker_simplified-V2.py and Photoaquation/photoaquation_rate.py; run_simulation() in Photoaquation/fe-cn-inverse-checker.py and Photoaquation/fecn_inversekinetics.py; and directly from module-level code (not inside a function) in prebiotic_scenario_runner.py.
    return UniversalRobustODESystem(
        database_species, differential_equations,
        species_atm, partial_pressures, henrys_constants, pH, rate_laws
    )



def solve_with_bulletproof_method(ode_system, time_span, initial_conditions):
    # Purpose: runs the integration, trying progressively more conservative solver settings until one succeeds.
    # Inputs are: ode_system, a UniversalRobustODESystem instance; time_span, the (start_time, end_time) tuple; initial_conditions, the starting concentration of every species.
    # Outputs are: solution, the result from scipy's solve_ivp() (or, if every attempt fails, an object with success=False and a message explaining the failure).
    # Called by: main(), in main.py; handle_extreme_stiffness_efficiently(), in extreme_stiffness_solver.py.
    from scipy.integrate import solve_ivp

    start_time, end_time = time_span
    duration = end_time - start_time

    # The solver settings get more conservative (larger tolerance, smaller maximum step) as the run gets longer, since a longer run has more opportunity to accumulate error or run into stiffness.
    if duration <= 3600:  # <= 1 hour
        configs = [
            {'method': 'Radau', 'rtol': 1e-6, 'atol': 1e-9, 'max_step': 360.0},
            {'method': 'BDF', 'rtol': 1e-6, 'atol': 1e-9, 'max_step': 360.0},
        ]
        tier_name = "high_precision"
    elif duration <= 86400:  # <= 1 day
        configs = [
            {'method': 'BDF', 'rtol': 1e-6, 'atol': 1e-9, 'max_step': 3600.0},
        ]
        tier_name = "balanced"
    elif duration <= 604800:  # <= 1 week
        configs = [
            {'method': 'BDF', 'rtol': 1e-6, 'atol': 1e-9, 'max_step': 21600.0},
        ]
        tier_name = "long_run"
    else:  # > 1 week
        configs = [
            {'method': 'BDF', 'rtol': 1e-8, 'atol': 1e-11, 'max_step': 3600.0},
        ]
        tier_name = "very_long_run"

    print(f"Using {tier_name} solver settings for {duration/3600:.1f}-hour simulation")

    # The analytic Jacobian (only available for prebiotic/unknown database types) is used by default, since it is faster than letting the solver estimate one numerically; this can be skipped when running interactively.
    jac = None
    if hasattr(ode_system, 'jacobian') and ode_system.jacobian is not None:
        if sys.stdin.isatty():
            try:
                answer = input("Use analytic Jacobian? (faster, recommended) [Y/n]: ").strip().lower()
            except EOFError:
                answer = ''
        else:
            answer = 'y'  # non-interactive: default to yes
        if answer in ('', 'y', 'yes'):
            jac = ode_system.jacobian
            print("   Using analytic Jacobian")
        else:
            print("   Using finite-difference Jacobian")

    initial_conditions = np.maximum(initial_conditions, 0.0)

    base_atol = configs[0].get('atol', 1e-11)

    # A tighter atol (absolute tolerance i.e., the smallest change in concentration the solver is required to resolve) is set for species with a blended rate law near cstar.
    # Their rate law changes shape sharply in that region, and the BDF solver can overshoot through zero if its steps are too large there.
    # Setting atol to roughly 1e-8 x cstar for these species forces small steps through the transition region, well below where CN-'s quasi-steady-state concentration (~4e-15 M) sits.
    sp_names = getattr(ode_system, 'database_species', [])
    sp_idx_map = getattr(ode_system, 'species_to_index', {})
    if sp_names and sp_idx_map:
        atol_arr = np.full(len(sp_names), base_atol)
        for sp, cstar_val in REGORD_SPECIES_CSTAR.items():
            if sp in sp_idx_map:
                tight = cstar_val * 1e-8  # a hundred-millionth of cstar which is tight enough to force small steps through the blend's transition region
                if tight < base_atol:
                    atol_arr[sp_idx_map[sp]] = tight
        atol_to_use = atol_arr
    else:
        atol_to_use = base_atol

    for i, config in enumerate(configs):
        try:
            print(f"Attempt {i+1}: {config['method']}")

            solution = solve_ivp(
                ode_system,
                time_span,
                initial_conditions,
                dense_output=True,
                jac=jac,
                **{**config, 'atol': atol_to_use},
            )

            if solution.success:
                print(f"Success with {config['method']}!")
                print(f"   Integration steps: {len(solution.t)}")

                # A concentration that dips slightly negative is reported here as numerical uncertainty rather than corrected: the raw solution is self-consistent (atoms still balance), and clamping it would break that consistency. Reporting it lets the numerical error bar on the affected species be judged directly.
                sp_names = getattr(ode_system, 'database_species', [])
                if sp_names:
                    mins = solution.y.min(axis=1)
                    flagged = [
                        (sp, mn) for sp, mn in zip(sp_names, mins)
                        if mn < -base_atol * 10   # Excursions this small are just floating-point noise below atol, so they are ignored.
                    ]
                    if flagged:
                        print(f"\n  Numerical excursions (reported as uncertainty, not corrected):")
                        for sp, mn in sorted(flagged, key=lambda x: x[1]):
                            display = sp.replace('MINUS', '-').replace('PLUS', '+')
                            print(f"    {display}: min = {mn:.3e} M")
                        print()

                return solution

            else:
                print(f"   {config['method']} failed: {solution.message}")

        except Exception as e:
            print(f"   {config['method']} crashed: {str(e)[:100]}")
            continue

    # If every standard method above has failed, the extreme-stiffness handler is tried as a last resort.
    print("All standard methods failed - trying extreme stiffness handling...")
    try:
        from extreme_stiffness_solver import handle_extreme_stiffness_efficiently
        return handle_extreme_stiffness_efficiently(
            ode_system, time_span, initial_conditions,
            ode_system.database_species if hasattr(ode_system, 'database_species') else [],
            None
        )
    except Exception as e:
        print(f"Extreme stiffness handler failed: {e}")

    # If nothing above worked, a failure result is returned rather than raising an exception, so the caller can handle it the same way as any other failed solve_ivp() result.
    class FailedSolution:
        # This is a stand-in result object, shaped like a failed scipy solve_ivp() result, so calling code can handle it the same way.
        def __init__(self):
            # Purpose: marks this result as a failure.
            # Inputs are: none.
            # Outputs are: none directly; self.success and self.message are set.
            # Called by: automatically, when a FailedSolution is created, just below.
            self.success = False
            self.message = "All universal robust solver attempts failed"

    return FailedSolution()


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def split_top_level_terms(expr):
    # Purpose: splits an expression into its top-level +/- separated terms, ignoring +/- inside brackets.
    # Inputs are: expr, an expression string, e.g. "(-1) * (k1*A) + (+1) * (k2*B)".
    # Outputs are: a list of its top-level "+"/"-" separated terms, ignoring any "+"/"-" that appears inside parentheses.
    # Called by: not called anywhere in the repo (defined but currently unused).
    terms = []
    current_term = ''
    depth = 0
    for char in expr:
        if char == '(':
            depth += 1
        elif char == ')':
            depth -= 1
        if depth == 0 and char in '+-' and current_term:
            terms.append(current_term.strip())
            current_term = char
        else:
            current_term += char
    if current_term:
        terms.append(current_term.strip())
    return terms


def replace_y_with_species(expression, database_species):
    # Purpose: makes a "y[i]"-indexed expression readable, by swapping in species names.
    # Inputs are: expression, a string containing references like "y[0]", "y[1]", etc.; database_species, the list of species names in the same order as those indices.
    # Outputs are: the same expression with every "y[i]" replaced by the display-format name of the corresponding species, e.g. "y[0]" -> "H+1".
    # Called by: not called anywhere in the repo (defined but currently unused).
    for i, species in enumerate(database_species):
        formatted_species = format_species_name(species)
        expression = expression.replace(f"y[{i}]", formatted_species)
    return expression
