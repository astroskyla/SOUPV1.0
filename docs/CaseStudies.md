# Case Studies

This page demonstrates how to run the four main case studies that showcase the SOUP model's capabilities across different chemical domains. Each case study is automatically optimised by the model's chemistry detection and specialised solver systems.

---

## Equilibrium Chemistry

**Demonstrates:** A comparison between Geochemist's Workbench (GWB) and the SOUP Model

### Setup and Execution

1. **Prepare your database:** Use or create a database with "equilibrium" in the filename (e.g., `my_equilibrium_system.db`). The model will automatically detect this as equilibrium chemistry. Note that to replicate the exact case study you should use `equilibrium.db`.

2. **Configure input files:** Set up your `atmosphere.dat` and `aqueous.dat` files with appropriate initial conditions for your equilibrium system. Once again note that to replicate the exact case study you should use `equilibrium-aqueous.dat` and `equilibrium-atmosphere.dat`.

3. **Run the simulation:**
    ```bash
    python main.py
    ```
    - When prompted, enter your equilibrium database name
    - **Recommended duration:** 1000+ seconds for full equilibration
    - The model automatically uses the high-precision Radau solver optimised for equilibrium chemistry

4. **What happens automatically:**
    - Chemistry detection triggers equilibrium-specific rate constants (up to 5×10⁸ s⁻¹)
    - Legacy solver system activates with maximum precision settings (rtol=1e-14)

5. **Analyse results:** The model outputs detailed equilibrium species distributions and generates reaction rate reports showing how the system reaches equilibrium. To generate figures as shown in the case study run:
```bash
cd SOUP/plotting
python equilibrium-plot.py
```

**Expected behavior:** Rapid initial changes as the system equilibrates, followed by stable steady-state concentrations that represent the true equilibrium.

---

## Fenton Chemistry

**Demonstrates:** Fast radical chemistry with iron catalysis.

### Setup and Execution

1. **Prepare your database:** Use or create a database with "fenton" in the filename. The model automatically optimises for radical chemistry. Note that to replicate the exact case study you should use `fenton.db`.
 

2. **Configure initial conditions:** 
    - Set up your `fenton-aqueous.dat` file with appropriate iron concentrations
    - **Key parameter:** Adjust H₂O₂ concentration (typically 0.001 M or 0.0002 M) to control reaction intensity

3. **Run the simulation:**
    ```bash
    python main.py
    ```
    - Enter your fenton database name when prompted
    - **Recommended duration:** 6000 minutes (allows full radical chemistry evolution)

4. **What happens automatically:**
    - Fenton chemistry detection triggers optimised radical rate constants
    - Simple BDF solver activates with settings that allow natural radical behavior

5. **Analyse results:** The model generates reaction rate reports. To generate figures as shown in the case study run:
```bash
cd SOUP/plotting
python plot-fenton.py
```
Note you may need to change the files being plotted.

**Expected behavior:** Fast initial radical burst, followed by substrate consumption. The system naturally handles the extreme concentration gradients typical of radical chemistry.

---

## Inverse Kinetics (Zymonic Acid Chemistry)

**Demonstrates:** Advanced parameter fitting using experimental data, showcased with the zymonic acid system.

### Setup and Execution

1. **Prepare experimental data:** The inverse kinetics module comes pre-configured with zymonic acid experimental data, but can be adapted for other systems.

2. **Run parameter estimation:**
    ```bash
    python inverse_kinetics.py
    ```

3. **Choose your approach:**
    - **Option 1:** Test simulation only (quick verification, ~1 minute)
    - **Option 2:** Standard optimisation (~10-30 minutes, good for initial parameter estimates)
    - **Option 3:** Multi-start optimisation (~8-12 hours, best parameter estimates)
    - **Option 4:** Optimisation with a quick bootstrap, 50 samples (~2 hours, initial uncertainty estimate)

    `python inverse_kinetics.py` currently runs Option 2 by default — edit the `choice` variable in the script to pick a different one. For the full 2245-sample bootstrap, run `seeded_inverse_kinetics.py` separately with your best-fit parameters from Option 2 or 3.

4. **What happens automatically:**
    - **Database creation:** Temporary zymonic reaction database with your rate constants
    - **Simulation engine:** Uses the existing SOUP framework for forward modeling
    - **Optimisation:** Differential evolution with intelligent bounds
    - **Uncertainty analysis:** Bootstrap resampling for confidence intervals
    - **Validation:** Parameter correlation analysis and convergence monitoring

5. **Outputs generated:**
    - **JSON results:** Complete parameter estimates with uncertainties
    - **CSV tables:** Parameter summaries
    - **Plots:** Experimental vs. model comparisons and parameter distributions

### Advanced Usage

For your own chemical systems:

1. **Adapt the experimental data:** Replace the zymonic acid data with your experimental measurements
2. **Modify the network:** Update the reaction network in the database creation function
3. **Adjust bounds:** Set realistic parameter bounds for your chemistry
4. **Run optimisation:** Follow the same workflow for your system

**Expected behavior:** The system finds optimal rate constants that reproduce your experimental data and provides rigorous uncertainty estimates suitable for publication.

---

## Prebiotic Chemistry

**Demonstrates:** Complex multi-pathway reaction networks typical of prebiotic chemistry.

### Setup and Execution

1. **Prepare your database:** Use or create a database with "prebiotic" in the filename. The model applies moderate rate constants suitable for complex networks. Note that to replicate the exact case study you should use `prebiotic.db`.

2. **Configure input files:** Set up your `atmosphere.dat` and `aqueous.dat` files with appropriate initial conditions. Once again note that to replicate the exact case study you should use `prebiotic-aqueous.dat` and `prebiotic-atmosphere.dat`.


3. **Run the simulation:**
    ```bash
    python main.py
    ```
    - Enter your prebiotic database name when prompted  
    - **Recommended duration:** Variable depending on your system (typically 1-100 days)

4. **What happens automatically:**
    - Prebiotic chemistry detection applies balanced rate constants (0.1-10 min⁻¹)
    - Universal robust solver system activates with comprehensive error handling
    - Automatic negative concentration prevention maintains physical realism
    - Henry's Law enforcement for atmospheric species

5. **Advanced features:** For highly complex prebiotic networks, the system may automatically escalate to extreme stiffness handling if needed.

**Expected behavior:** Gradual evolution of complex organic chemistry with multiple competing pathways. The robust solver system ensures stable integration even with intricate reaction networks.

---

## Performance Tips for All Case Studies

### Always Start with Diagnostics
```bash
python diagnose.py
```
This health check identifies potential issues before running long simulations and provides performance estimates.

### Progressive Testing Strategy
1. **Short test:** 1 hour to verify system stability
2. **Medium test:** 1 day to check performance  
3. **Full simulation:** Your target duration

### Monitoring Integration Health
Watch for these positive indicators:

- *"Success with..."*,
- *"All concentrations positive!"*, and
- *"Low/Moderate stiffness. Should be manageable"*

### When to Expect Automatic Escalation
For extremely stiff systems (stiffness ratio > 10¹²), you'll see messages about escalation to conservative methods. This is normal and ensures robust completion.

---

