# This module solves chemistry simulations that run for longer than a week.
# Rather than integrating the whole period in one go, it splits the run into adaptive time chunks, which keeps long-duration simulations efficient and manageable.

import numpy as np
from scipy.integrate import solve_ivp
import time


def solve_long_term_chemistry(ode_system, time_span, initial_conditions, target_points=1000):
    # Purpose: runs a long simulation in manageable time chunks and stitches the results back together.
    # Inputs are: ode_system, the ODE system function (or callable object) to integrate; time_span, the (start_time, end_time) tuple in seconds; initial_conditions, the starting concentration of every species; target_points, the desired number of output points (currently unused, kept for a consistent function signature).
    # Outputs are: a combined solution object covering the whole time span, built by joining together the results from each individual time chunk.
    # Called by: main(), in main.py, for integrations longer than a week.
    start_time, end_time = time_span
    total_time = end_time - start_time

    print(f"Long-term solver for {total_time/86400:.1f} days ({total_time/31536000:.2f} years)")

    # The chunk size is chosen based on the total run length: week-long chunks for a run of more than 30 days, day-long chunks for a run of more than a week, and 10 even chunks otherwise.
    if total_time > 86400 * 30:  # More than 30 days
        chunk_size = 86400 * 7  # 1 week chunks
        print(f"   Using 1-week chunks ({total_time/chunk_size:.0f} chunks total)")
    elif total_time > 86400 * 7:  # More than 1 week
        chunk_size = 86400  # 1 day chunks
        print(f"   Using 1-day chunks ({total_time/chunk_size:.0f} chunks total)")
    else:
        chunk_size = total_time / 10  # 10 chunks for shorter periods
        print(f"   Using {chunk_size/3600:.1f}-hour chunks")

    all_times = [start_time]
    all_solutions = [np.array(initial_conditions)]
    current_time = start_time
    current_state = np.array(initial_conditions)

    # LSODA is used here since it copes well with a mixture of stiff and non-stiff behaviour over a long run, and the tolerances are relaxed slightly for speed.
    solver_config = {
        'method': 'LSODA',  # Best for long-term mixed stiff/non-stiff
        'rtol': 1e-4,       # Relaxed for speed. Probably not ideal based on the concentrations you are looking at. Treatwitha pinch of salt...
        'atol': 1e-7,       # Reasonable absolute tolerance
        'max_step': chunk_size / 100,  # Allow larger steps within chunks
        'first_step': None,  # Let solver choose
    }

    chunk_count = 0
    total_chunks = int(np.ceil(total_time / chunk_size))

    print(f"   Solver: {solver_config['method']}, rtol={solver_config['rtol']:.0e}, atol={solver_config['atol']:.0e}")
    print(f"   Max step size: {solver_config['max_step']/3600:.1f} hours")
    print()

    while current_time < end_time:
        chunk_count += 1
        chunk_start = current_time
        chunk_end = min(current_time + chunk_size, end_time)
        chunk_duration = chunk_end - chunk_start

        print(f"Chunk {chunk_count}/{total_chunks}: t={chunk_start/86400:.1f} to {chunk_end/86400:.1f} days", end=" ")

        chunk_start_clock = time.time()

        try:
            chunk_solution = solve_ivp(
                ode_system,
                (chunk_start, chunk_end),
                current_state,
                dense_output=True,
                **solver_config
            )

            if not chunk_solution.success:
                print(f"FAILED: {chunk_solution.message}")
                # If the default settings fail for this chunk, a tighter, more conservative retry is attempted before giving up on it.
                print("   Retrying with conservative settings...")

                conservative_config = solver_config.copy()
                conservative_config['rtol'] = 1e-5
                conservative_config['atol'] = 1e-8
                conservative_config['max_step'] = chunk_duration / 1000

                chunk_solution = solve_ivp(
                    ode_system,
                    (chunk_start, chunk_end),
                    current_state,
                    dense_output=True,
                    **conservative_config
                )

                if not chunk_solution.success:
                    print(f"Conservative retry also failed: {chunk_solution.message}")
                    break

            # The final state and time of this chunk become the starting point for the next one.
            current_state = chunk_solution.y[:, -1]
            current_time = chunk_solution.t[-1]

            # Only up to 100 points per chunk are kept, to avoid the combined result using excessive memory over a long run.
            chunk_points = min(100, len(chunk_solution.t))  # Max 100 points per chunk
            if chunk_points > 2:
                indices = np.linspace(0, len(chunk_solution.t)-1, chunk_points, dtype=int)
                selected_times = chunk_solution.t[indices[1:]]  # Skip first (duplicate)
                selected_solutions = chunk_solution.y[:, indices[1:]].T
            else:
                selected_times = [chunk_solution.t[-1]]
                selected_solutions = [chunk_solution.y[:, -1]]

            all_times.extend(selected_times)
            all_solutions.extend(selected_solutions)

            chunk_wall_time = time.time() - chunk_start_clock
            steps_taken = len(chunk_solution.t)

            print(f"{chunk_wall_time:.1f}s ({steps_taken} steps)")

            # Any concentration that drifted slightly negative between chunks is clamped back to zero here. Not a fan of thisbe aware an maybe change in future?
            negative_count = np.sum(current_state < -1e-10)
            if negative_count > 0:
                print(f"   {negative_count} negative concentrations detected")
                current_state = np.maximum(current_state, 0.0)

            # Progress is reported every 10 chunks, with a rough time-remaining estimate based on how long the last chunk took.
            if chunk_count % 10 == 0:
                avg_time_per_chunk = chunk_wall_time
                remaining_chunks = total_chunks - chunk_count
                estimated_time_remaining = avg_time_per_chunk * remaining_chunks

                print(f"   Progress: {chunk_count/total_chunks*100:.1f}% complete")
                print(f"   Estimated time remaining: {estimated_time_remaining/60:.1f} minutes")
                print(f"   Memory usage: {len(all_times)} time points stored")

        except Exception as e:
            print(f"CRASHED: {str(e)}")
            print("   This chunk failed completely - stopping integration")
            break

    print(f"\n Long-term integration completed: {current_time/86400:.1f} days")
    print(f"   Final time reached: {current_time:.0f} seconds")
    print(f"   Total data points: {len(all_times)}")

    combined_times = np.array(all_times)
    combined_solutions = np.array(all_solutions).T

    class LongTermSolution:
        # This mimics the parts of scipy's usual solve_ivp() result that the rest of the code relies on (t, y, success, message), since the combined long-term result isn't produced by a single solve_ivp() call.
        def __init__(self, t, y, success, message="Long-term integration completed"):
            # Purpose: stores the combined result so it looks like a normal solve_ivp() result.
            # Inputs are: t, the combined time points; y, the combined concentrations; success, whether the run reached (or got close to) end_time; message, a description of the outcome.
            # Outputs are: none directly; self.t, self.y, self.success, self.message and self.sol are all set.
            # Called by: automatically, when a LongTermSolution is created, at the bottom of solve_long_term_chemistry().
            self.t = t
            self.y = y
            self.success = success
            self.message = message
            self.sol = None  # No dense output for memory efficiency

        def __call__(self, t_eval):
            # Purpose: evaluates the solution at arbitrary times, standing in for a normal solve_ivp() result's dense output.
            # Inputs are: t_eval, the times to evaluate the solution at.
            # Outputs are: the concentrations at those times, found by simple linear interpolation between the stored points (this stands in for the dense output that a normal solve_ivp() result would provide).
            # Called by: whichever code asks this solution object for values at specific times, since it behaves like a normal solve_ivp() result.
            return np.interp(t_eval, self.t, self.y.T).T

    success = current_time >= end_time * 0.95  # Accept if we got close
    return LongTermSolution(combined_times, combined_solutions, success)


def create_high_performance_ode_system(base_ode_system):
    # Purpose: wraps an ODE system so its call speed can be monitored during a long run.
    # Inputs are: base_ode_system, the ODE system function (or callable object) to wrap.
    # Outputs are: a HighPerformanceODESystem instance that behaves exactly like base_ode_system, but also times each call and periodically prints a performance summary.
    # Called by: main(), in main.py, before a long-term integration.

    class HighPerformanceODESystem:
        def __init__(self, base_system):
            # Purpose: stores the wrapped ODE system and sets up the performance-tracking state.
            # Inputs are: base_system, the ODE system to wrap.
            # Outputs are: none directly; self.base_system and the tracking attributes are set.
            # Called by: automatically, when a HighPerformanceODESystem is created, in create_high_performance_ode_system().
            self.base_system = base_system
            self.call_count = 0
            self.last_report_time = time.time()
            self.evaluation_times = []

        def __call__(self, t, y):
            # Purpose: calls the wrapped ODE system, timing it and periodically printing a performance summary.
            # Inputs are: t, the current simulation time; y, the current concentration of every species.
            # Outputs are: result, whatever the wrapped ODE system returns for (t, y), unchanged.
            # Called by: automatically, whenever this wrapper is passed as the right-hand-side function to a solver.
            eval_start = time.time()

            result = self.base_system(t, y)

            eval_time = time.time() - eval_start
            self.evaluation_times.append(eval_time)
            self.call_count += 1

            # A performance summary is printed at most every 5 minutes, so it doesn't flood the console during a long run.
            current_time = time.time()
            if current_time - self.last_report_time > 300:  # Every 5 minutes
                avg_eval_time = np.mean(self.evaluation_times[-1000:])  # Last 1000 calls
                print(f"   ODE evaluations: {self.call_count} total, {avg_eval_time*1000:.2f}ms avg")
                print(f"   Current simulation time: t={t/86400:.2f} days")
                self.last_report_time = current_time

            return result

        # Any attribute not defined on this wrapper (e.g. .jacobian, .database_species) is forwarded through to the wrapped system, so this wrapper is a drop-in replacement for it.
        def __getattr__(self, name):
            # Purpose: forwards any attribute this wrapper doesn't itself define through to the wrapped ODE system.
            # Inputs are: name, the attribute name being looked up.
            # Outputs are: the corresponding attribute from self.base_system.
            # Called by: automatically, by Python, whenever code accesses an attribute on this wrapper that isn't defined directly on it (e.g. .jacobian, .database_species).
            return getattr(self.base_system, name)

    return HighPerformanceODESystem(base_ode_system)


def estimate_integration_time(ode_system, initial_conditions, target_duration):
    # Purpose: estimates how long a full run will take, by timing a short test run first.
    # Inputs are: ode_system, the ODE system to test; initial_conditions, the starting concentrations; target_duration, the full simulated duration (in seconds) that the estimate is being made for.
    # Outputs are: an estimate of how long the full integration will take (in seconds), based on timing a short test run of 1 hour or 1% of the target duration (whichever is smaller); returns None if the test integration fails or crashes.
    # Called by: not called anywhere in the repo.
    print("Estimating integration performance...")

    test_duration = min(3600, target_duration / 100)  # Test 1 hour or 1% of target

    start_time = time.time()

    try:
        test_solution = solve_ivp(
            ode_system,
            (0, test_duration),
            initial_conditions,
            method='LSODA',
            rtol=1e-4,
            atol=1e-7
        )

        wall_time = time.time() - start_time

        if test_solution.success:
            steps_per_second = len(test_solution.t) / test_duration
            time_per_simulated_second = wall_time / test_duration

            estimated_total_time = target_duration * time_per_simulated_second

            print(f"   Test integration: {test_duration/3600:.1f} hours simulated in {wall_time:.1f}s")
            print(f"   Solver steps per simulated hour: {steps_per_second * 3600:.0f}")
            print(f"   Estimated time for full integration: {estimated_total_time/3600:.1f} hours")

            if estimated_total_time > 7200:  # More than 2 hours
                print("   Full integration may take a very long time")
                print("   Consider using shorter time periods or coarser tolerances")

            return estimated_total_time
        else:
            print(f"   Test failed: {test_solution.message}")
            return None

    except Exception as e:
        print(f"   Test crashed: {e}")
        return None
