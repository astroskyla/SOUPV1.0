# Running the Model

This guide explains how to set up and run the model.

---

## Input Files

1. **Navigate to the [Database Overview](Database.md#database)**  
    - Browse the existing databases to find data related to your reaction network of interest.
    - If you can’t find relevant data or notice missing information, you can create a new database by following the instructions in the [database section](Database.md#database).

2. **Prepare Input Files**

- After cloning the repository from GitHub, you should have a folder named `SOUP` in the main project directory.
- Within the `SOUP` folder, locate the `Inputs` folder. Here, you will find the `atmosphere.dat` and `aqueous.dat` files.
- Edit these files to define the environment and starting concentrations for your simulation.

<details class=validation>
<summary><strong>aqueous.dat</strong></summary>

```
# Aqueous Solution Data File
# This file contains information about species, their concentrations and their Henry's Law Constants as well as the properties of the aqueous solution.
# Please follow the format below to update or add new data.

# Species and their corresponding partial pressures (in atm)
# Format: species, concentration (M)
species, concentration
A, 0.1
B, 0

# Total Volume (in L)
Total Volume
1

# Temperature (in Kelvin)
# This is currently set as an average but a temperature spectrum could be used in the future
Total Temperature
298.15

# pH
# This is currently set as an average pH but a spectrum could be used in the future 
pH
3.0

```

</details>

<details class=validation>
<summary><strong>atmosphere.dat</strong></summary>

```
# Atmosphere Data File
# This file contains information about species, their partial pressures, total pressure, total temperature, and the solar spectrum file in use.
# Please follow the format below to update or add new data.

# Species and their corresponding partial pressures (in atm) along with Henry's Law Constants
# Format: species, partial_pressure (atm), Henry's Law Constant (M/atm)
species, partial_pressure, constant (M/atm)
H2O, 0.001, N/A
O2, 0.2095, N/A
N2, 0.7808, N/A
Ar, 9.3e-3, N/A
SO2, 0.0, 1.47910838817
CO2, 2.0e-6, 3.10e-2
H2S, 0.0, 0.087


# Total pressure (in atm)
# This is currently set as however a pressure spectrum could be used in the future
Total pressure
1.0

# Temperature (in Kelvin)
# This is currently set as an average but a temperature spectrum could be used in the future
Total Temperature
288

# Solar spectrum file (relative or absolute path)
# There is currently no solar spectrum but it can be added in the future to work with specific photochemistry
Solar Spectrum File
spectrum.dat

```

</details>

**Important**: Maintain the format strictly to avoid errors during simulation.

## Running the Model

3. **Execute the Script**  
    - Once your input files are configured and you have selected a database, run the `main.py` script from the `SOUP` folder in your terminal:

        ```
        python main.py
        ```

    - You’ll be prompted to select a database by entering its title as shown in the database overview.

4. **Complete the Prompts**  
    - Respond to the following prompts to customise your simulation:

        - **Enter the end time (seconds)** you'd like to evolve your system over.
        - **Enter number of data points** (press enter to use the default).
        - **Open system (Henry's law atmospheric replenishment)?** This is only asked for prebiotic-type databases.
        - **Enter species to plot (comma-separated) or 'all'**

    - `species_concentrations.md` is saved to `Outputs/` automatically once the simulation finishes. For each species you plot, you'll then be asked whether to use a log-scale y-axis, whether to save the plot as an image, and if so, what to name it.

5. **Analyse the Results**  
    - Once the simulation is complete, review the output plots and analyse the results as needed.
    - You can plot the results yourself by using the species_concentrations.md output file.
    - You can probe the reaction rates at each time step by looking at the reaction_rates.html output file.
    - After each run, a mass-balance check runs automatically for carbon, nitrogen, sulfur, and iron (hydrogen is excluded, since reactions are treated as pseudo-steps), printing a pass/fail summary to the terminal. It won't stop the run if something fails.
    

---

## Model Architecture

### Overview
The model is organised into a core simulation engine (`SOUP/`), a set of standalone analysis/plotting scripts built on top of it, and a small Flask app for browsing and editing the reaction database.

### Core simulation (`SOUP/`)
- **`main.py`** — entry point; detects the chemistry type from the selected database, picks a matching solver, and runs the simulation.
- **`rate_laws.py`** — builds forward/reverse rate laws and stoichiometry from a reaction database.
- **`differentials.py`** — builds and solves the ODE system (separate legacy, simple, and bulletproof solvers for equilibrium, Fenton, and prebiotic chemistry respectively).
- **`data.py`** — reads the atmosphere/aqueous input data files.
- **`plotting.py`** — plots simulation output concentrations over time.
- **`utilities.py`** — shared helper functions (e.g. initial concentration setup) used across the other modules.
- **`diagnose.py`** — pre-flight diagnostic checks on a database/input setup before running a full simulation.
- **`extreme_stiffness_solver.py`** — fallback solver for reaction systems too stiff for the normal bulletproof solver.
- **`long_term_solver.py`** — splits multi-week simulations into adaptive time chunks.
- **`inverse_kinetics.py`** — fits rate constants for the 5-species zymonic acid network to experimental data via optimisation.
- **`seeded_inverse_kinetics.py`** — bootstrap uncertainty estimation for the zymonic acid rate constants fitted by `inverse_kinetics.py`.
- **`prebiotic_scenario_runner.py`** — runs the prebiotic network across a grid of HCN/Fe²⁺ scenarios and saves the results.
- **`run_sai_miyakawa_networks.py`** — batch runner comparing the Sai and Miyakawa HCN hydrolysis networks across pH/sulfur scenarios.
- **`reaction_rates_logger.py`** — generates an HTML report of reaction fluxes over the course of a simulation.

### Tests (`SOUP/tests/`)
- **`checks.py`** — automated mass-balance/composition checks run after each simulation.
- **`prebiotic_tests.py`** — verification suite for the prebiotic network (mass balance, convergence, open-system, and `regord_pow` checks).

### Plotting scripts (`SOUP/plotting/`)
- **`equilibrium-plot.py`** — equilibrium chemistry comparison plots.
- **`plot-fenton.py`** — Fenton reaction validation plot against the Gallard & De Laat model and experimental data.
- **`prebiotic_scenario_plot.py`** — plots the multi-scenario prebiotic results produced by `prebiotic_scenario_runner.py`.
- **`zymonic-plot.py`** — plots the zymonic acid case study results.

### Database generators (`SOUP/Inputs/Databases/`)
- **`mk-fenton.py`** — generates the corrected Fenton reaction database.
- **`Sai_Miyakawa/make_network_dbs.py`** — generates the four Sai/Miyakawa comparison databases used by `run_sai_miyakawa_networks.py`.

### Other analysis scripts
- **`SOUP/EXPERIMENTAL_RAW_DATA/parsing.py`** — parses raw experimental data filenames/CSVs.
- **`SOUP/Outputs/SaiData/sai-plot.py`** — publication figure comparing the Sai and Miyakawa prebiotic chemistry results.

### Database app (`Flask/`)
- **`flask-app.py`** — the Flask web app for browsing, forking, and editing the reaction database described in [Database](Database.md).
- **`sql-db.py`** — (re)builds `master.db` from `database.sql`.

### Interactions
- `main.py` ties together `rate_laws.py`, `differentials.py`, `data.py`, `utilities.py`, and `plotting.py` for a single simulation run.
- The standalone scripts above build on this core (calling into `data.py`/`utilities.py`/`differentials.py` directly) to run batches of simulations, fit parameters, or generate figures without going through `main.py`.
- The Flask app manages the SQLite database files that the core simulation reads from.

---

## Installation

Dependencies are managed with [`uv`](https://docs.astral.sh/uv/); run `uv sync` to install everything listed in `pyproject.toml`.

---