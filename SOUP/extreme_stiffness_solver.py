# This module handles extremely stiff reaction systems (ones where reaction rates vary hugely in speed), for when the usual bulletproof solver runs into difficulty.
# It is a fallback of last resort: conservative implicit solver methods with relaxed tolerances and small steps, used only when needed.
# Thinking of ditching as just super erroneous. 

import numpy as np
from scipy.integrate import solve_ivp
import warnings
import time


class ExtremeStiffnessSolver:
    # This solver is for extremely stiff systems: it uses conservative implicit solver methods with relaxed tolerances, trading some precision for a better chance of succeeding at all.

    def __init__(self, stiffness_ratio):
        # Purpose: decides which of the three fallback strategies to use, based on how stiff the system is.
        # Inputs are: stiffness_ratio, the ratio of the fastest to the slowest reaction rate in the system, used to decide how conservative the solver needs to be.
        # Outputs are: none directly; self.integration_method is set based on stiffness_ratio.
        # Called by: handle_extreme_stiffness_efficiently(), in this file.
        self.stiffness_ratio = stiffness_ratio
        self.integration_method = self.select_method(stiffness_ratio)

    def select_method(self, stiffness_ratio):
        # Purpose: maps a stiffness ratio onto a named integration strategy.
        # Inputs are: stiffness_ratio, the fastest-to-slowest reaction rate ratio.
        # Outputs are: a string naming the integration strategy to use: "ultra_conservative" above 1e15, "conservative_implicit" above 1e12, and "standard_stiff" otherwise.
        # Called by: __init__(), in this file.
        if stiffness_ratio > 1e15:
            return "ultra_conservative"
        elif stiffness_ratio > 1e12:
            return "conservative_implicit"
        else:
            return "standard_stiff"

    def solve_extreme_stiff(self, ode_system, time_span, initial_conditions, database_species):
        # Purpose: dispatches to the right solve_* method for this system's stiffness level.
        # Inputs are: ode_system, the ODE system to integrate; time_span, the (start_time, end_time) tuple; initial_conditions, the starting concentration of every species; database_species, the list of species names (accepted for a consistent signature, not used directly here).
        # Outputs are: the solution from whichever of the three solve_* methods below matches self.integration_method.
        # Called by: handle_extreme_stiffness_efficiently(), in this file.
        print(f"EXTREME STIFFNESS: {self.stiffness_ratio:.0e}")
        print(f"Using method: {self.integration_method}")

        if self.integration_method == "ultra_conservative":
            return self.solve_ultra_conservative(ode_system, time_span, initial_conditions)
        elif self.integration_method == "conservative_implicit":
            return self.solve_conservative_implicit(ode_system, time_span, initial_conditions)
        else:
            return self.solve_standard_stiff(ode_system, time_span, initial_conditions)

    def solve_ultra_conservative(self, ode_system, time_span, initial_conditions):
        # Purpose: the most cautious solve attempt, for the very stiffest systems.
        # Inputs are: ode_system, time_span and initial_conditions, as above.
        # Outputs are: a successful solution from the first of Radau, BDF or LSODA that succeeds, using very relaxed tolerances and very small steps; raises a RuntimeError if all three fail.
        # Called by: solve_extreme_stiff(), for the most extreme stiffness ratios; and as a last-resort fallback from solve_conservative_implicit(), both in this file.
        print("   Using ultra-conservative integration...")

        duration = time_span[1] - time_span[0]

        # These settings favour stability over precision: a relaxed tolerance and very small maximum step.
        settings = {
            'rtol': 1e-2,  # Relaxed tolerance for stability. At this point accuracy is out the window.
            'atol': 1e-4,
            'max_step': min(0.1, duration / 10000),  # Very small steps
            'first_step': 1e-8
        }

        methods = ['Radau', 'BDF', 'LSODA']

        for method in methods:
            try:
                print(f"   Trying ultra-conservative {method}...")

                solution = solve_ivp(
                    ode_system,
                    time_span,
                    initial_conditions,
                    method=method,
                    **settings,
                    dense_output=False
                )

                if solution.success:
                    print(f"   Ultra-conservative {method} succeeded!")
                    print(f"   Integration steps: {len(solution.t)}")
                    solution.y = np.maximum(solution.y, 0.0)  # Final safety
                    return solution
                else:
                    print(f"   {method} failed: {solution.message}")

            except Exception as e:
                print(f"   {method} crashed: {str(e)[:100]}")
                continue

        raise RuntimeError("All ultra-conservative methods failed")

    def solve_conservative_implicit(self, ode_system, time_span, initial_conditions):
        # Purpose: a moderately cautious solve attempt, for moderately extreme stiffness.
        # Inputs are: ode_system, time_span and initial_conditions, as above.
        # Outputs are: a successful solution from the first configuration that succeeds (Radau, then BDF, then a relaxed LSODA), each less strict than a standard solve but not as relaxed as solve_ultra_conservative(); falls back to solve_ultra_conservative() if every configuration fails.
        # Called by: solve_extreme_stiff(), for moderately extreme stiffness ratios; and as an escalation from solve_standard_stiff(), both in this file.
        print("   Using conservative implicit integration...")

        duration = time_span[1] - time_span[0]

        # These settings are conservative, but not as extreme as solve_ultra_conservative()'s.
        configs = [
            {
                'method': 'Radau',
                'rtol': 1e-3,
                'atol': 1e-5,
                'max_step': min(1.0, duration / 1000),
                'first_step': 1e-6,
                'desc': 'Conservative Radau'
            },
            {
                'method': 'BDF',
                'rtol': 1e-3,
                'atol': 1e-5,
                'max_step': min(2.0, duration / 500),
                'first_step': 1e-6,
                'desc': 'Conservative BDF'
            },
            {
                'method': 'LSODA',
                'rtol': 1e-2,
                'atol': 1e-4,
                'max_step': min(5.0, duration / 200),
                'desc': 'Relaxed LSODA'
            }
        ]

        for config in configs:
            try:
                print(f"   Trying {config['desc']}...")

                solution = solve_ivp(
                    ode_system,
                    time_span,
                    initial_conditions,
                    **{k: v for k, v in config.items() if k not in ['desc']},
                    dense_output=False
                )

                if solution.success:
                    print(f"   {config['desc']} succeeded!")
                    print(f"   Integration steps: {len(solution.t)}")
                    solution.y = np.maximum(solution.y, 0.0)
                    return solution
                else:
                    print(f"   {config['method']} failed: {solution.message}")

            except Exception as e:
                print(f"   {config['method']} crashed: {str(e)[:100]}")
                continue

        # If none of the conservative configurations succeed, the ultra-conservative method is tried as a last resort.
        print("   Conservative methods failed - trying ultra-conservative...")
        return self.solve_ultra_conservative(ode_system, time_span, initial_conditions)

    def solve_standard_stiff(self, ode_system, time_span, initial_conditions):
        # Purpose: a standard stiff-solver attempt, for the least extreme stiffness this fallback handles.
        # Inputs are: ode_system, time_span and initial_conditions, as above.
        # Outputs are: a successful solution from the first of Radau, BDF or LSODA that succeeds, using standard stiff-solver tolerances; escalates to solve_conservative_implicit() if every configuration fails.
        # Called by: solve_extreme_stiff(), for the least extreme stiffness ratios that still reach this solver, in this file.
        print("   Using standard stiff integration...")

        configs = [
            {'method': 'Radau', 'rtol': 1e-6, 'atol': 1e-9, 'max_step': 0.1},
            {'method': 'BDF', 'rtol': 1e-6, 'atol': 1e-9, 'max_step': 1.0},
            {'method': 'LSODA', 'rtol': 1e-5, 'atol': 1e-8, 'max_step': 2.0}
        ]

        for config in configs:
            try:
                print(f"   Trying standard {config['method']}...")

                solution = solve_ivp(
                    ode_system, time_span, initial_conditions,
                    dense_output=False, **config
                )

                if solution.success:
                    print(f"   Standard {config['method']} succeeded!")
                    return solution

            except Exception as e:
                print(f"   {config['method']} crashed: {str(e)[:100]}")
                continue

        # If none of the standard configurations succeed, the more conservative method is tried next.
        print("   Standard methods failed - escalating to conservative...")
        return self.solve_conservative_implicit(ode_system, time_span, initial_conditions)


def handle_extreme_stiffness_efficiently(ode_system, time_span, initial_conditions, database_species, rate_laws=None):
    # Purpose: the entry point for this fallback which works out how stiff the system is, and routes it to the right solver.
    # Inputs are: ode_system, the ODE system to integrate; time_span, the (start_time, end_time) tuple; initial_conditions, the starting concentration of every species; database_species, the list of species names; rate_laws, the rate-law list (optional), used to estimate how stiff the system is.
    # Outputs are: a solution, either from the extreme-stiffness solver above (if the system is judged extremely stiff) or from the standard bulletproof solver (if it is not).
    # Called by: main(), in main.py, as a last resort when the bulletproof method has already failed; and solve_with_bulletproof_method(), in differentials.py, for the same reason.
    # The stiffness ratio (fastest rate divided by slowest rate, from the rate laws) is estimated here, as a rough gauge of how numerically demanding the system is.
    if rate_laws:
        all_rates = []
        for law in rate_laws:
            for rate_expr in [law['forward'], law['reverse']]:
                if rate_expr and rate_expr != 'None':
                    try:
                        rate_str = rate_expr.split('*')[0].strip()
                        rate = float(rate_str)
                        if rate > 0:
                            all_rates.append(rate)
                    except:
                        pass

        if len(all_rates) >= 2:
            stiffness_ratio = max(all_rates) / min(all_rates)
        else:
            stiffness_ratio = 1e16  # Assume very stiff
    else:
        stiffness_ratio = 1e16  # Assume very stiff

    duration = time_span[1] - time_span[0]
    print(f"Extreme stiffness analysis:")
    print(f"   Stiffness ratio: {stiffness_ratio:.1e}")
    print(f"   Integration time: {duration}s ({duration/3600:.2f} hours)")

    if stiffness_ratio > 1e12:
        print("   Extreme stiffness detected - using conservative methods")
        solver = ExtremeStiffnessSolver(stiffness_ratio)
        return solver.solve_extreme_stiff(ode_system, time_span, initial_conditions, database_species)
    else:
        # A system that is not extremely stiff is handled by the standard bulletproof method instead, since it doesn't need this fallback's extra caution.
        print("   System not extremely stiff - using standard universal robust method")
        from differentials import solve_with_bulletproof_method
        return solve_with_bulletproof_method(ode_system, time_span, initial_conditions)


def test_extreme_stiffness_solver():
    # Purpose: prints a summary of what this module does, as a quick sanity check when run directly.
    # Inputs are: none.
    # Outputs are: nothing is returned; a short description of this module's approach and stiffness thresholds is printed.
    # Called by: the test block at the bottom of this file (if __name__ == "__main__").
    print("Extreme Stiffness Solver - Simplified Version")
    print("=" * 50)
    print("This version handles extreme stiffness through:")
    print("1. Conservative implicit methods (Radau, BDF, LSODA)")
    print("2. Relaxed tolerances for stability")
    print("3. Small time steps")
    print("4. No complex radical approximations")
    print()
    print("Stiffness thresholds:")
    print("- > 1e15: Ultra-conservative")
    print("- > 1e12: Conservative implicit")
    print("- < 1e12: Standard stiff")
    print()
    print("Extreme stiffness solver ready")


if __name__ == "__main__":
    test_extreme_stiffness_solver()
