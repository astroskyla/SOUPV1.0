# The SOUP Model

## An Aqueous Chemical Kinetics Model for Prebiotic Systems

This code simulates the time evolution of chemical species in an aqueous environment, accounting for multiple reaction networks, equilibrium chemistry, and gas exchange with the atmosphere. The model features specialised solver systems and advanced capabilities for handling extreme stiffness and long-term integration.

---

## Overview

**Reaction Kinetics:** Simulates chemical reactions between species using rate laws based on rate constants and concentrations, and solves the resulting differential equations to model how species concentrations change over time.

&nbsp;

**Equilibrium Chemistry:** Represents chemical equilibria in aqueous solutions to capture reversible reactions and species distributions.

&nbsp;

**Atmospheric and Aqueous Interaction:** Initialises species from atmospheric and aqueous input files to simulate gas exchange between the atmosphere and the solution.

---

## Modelling Equations

**Chemical Reactions**

A chemical reaction is represented as:

aA + bB → cC + dD

&nbsp;

**Rate Laws**

The rate law governing a reaction is given by:

&nbsp;

\[
\begin{flalign*}
\text{rate} = k [A]^m [B]^n && \\
\end{flalign*}
\]

- \( a \), \( b \): Stoichiometric coefficients of reactants \( A \) and \( B \).
- \( c \), \( d \): Stoichiometric coefficients of products \( C \) and \( D \).
- \( k \): Rate constant.
- \( m \), \( n \): Reaction orders with respect to reactants \( A \) and \( B \).

&nbsp;

**Differential Equations for Reaction Dynamics**

For a reversible chemical reaction of the form:

aA + bB ↔ cC + dD

&nbsp;

The forward reaction rate is:

\[
\begin{flalign*}
\text{rate}_f = k_f [A]^m [B]^n && \\
\end{flalign*}
\]

where k<sub>f</sub> is the forward rate constant and m and n are reaction orders with respect to A and B.

&nbsp;

The reverse rate is:

\[
\begin{flalign*}
\text{rate}_r = k_r [C]^p [D]^q && \\
\end{flalign*}
\]

where k<sub>r</sub> is the reverse rate constant and p and q are reaction orders with respect to products C and D.

&nbsp;

The differential equations for each species over time (t) are:

\[
\begin{flalign*}
\frac{d[A]}{dt} &= -a \cdot \text{rate}_f + a \cdot \text{rate}_r && \\
\frac{d[B]}{dt} &= -b \cdot \text{rate}_f + b \cdot \text{rate}_r && \\
\frac{d[C]}{dt} &= +c \cdot \text{rate}_f - c \cdot \text{rate}_r && \\
\frac{d[D]}{dt} &= +d \cdot \text{rate}_f - d \cdot \text{rate}_r &&
\end{flalign*}
\]

---

## Inputs

#### The SOUP Database
- Chemical species
- The dependence of the reaction on each reactant (the order)
- The rate constant and/or equilibrium constant for a given reaction

#### The Atmospheric Composition
- Chemical species
- Partial Pressures of each species
- Henry's Law Constant (if applicable)

#### The Aqueous Solution Composition
- Chemical species
- Initial concentrations of each species
- pH of the aqueous solution

---

## Chemistry Type Detection and Processing

The SOUP model detects the chemistry type based on the database name and applies chemistry-specific rate constant calculations and solver optimisations:

### Equilibrium Chemistry
- **Rate Strategy:** Very high rates (up to 5×10⁸ s⁻¹) for fast equilibration
- **Solver:** High-precision Radau method for maximum accuracy

### Fenton Chemistry
- **Rate Strategy:** Base rates around 10⁴ s⁻¹, scaling up to a ceiling of 5×10⁷ s⁻¹ depending on equilibrium constants
- **Solver:** Simple BDF approach with gentle constraint handling that allows natural radical behavior

### Prebiotic Chemistry
- **Rate Strategy:** Moderate rates (0.1-10 min⁻¹) for complex multi-pathway networks
- **Solver:** Universal robust system (see details in the repository)

### Unknown Chemistry Types
- **Fallback:** Any database not matching the above patterns
- **Rate Strategy:** Conservative approach with bulletproof numerical handling
- **Solver:** Universal robust system with intelligent escalation to extreme stiffness handlers

---

## Solver Architecture

The SOUP model features three specialised solver systems that are automatically selected based on your chemistry type:

### Legacy System (Equilibrium Chemistry)
**Designed for maximum precision in equilibrium studies**

- **Solver:** High-precision Radau method
- **Tolerances:** rtol=1e-14, atol=1e-8 
- **Strengths:** Exceptional accuracy for equilibrium calculations
- **Features:** Proven approach for systems requiring precise equilibrium distributions
- **Best for:** Chemical systems where equilibrium accuracy is vital

### Simple Fenton System (Radical Chemistry)
**Optimised for fast radical reactions**

- **Solver:** BDF (Backward Differentiation Formula)
- **Tolerances:** rtol=1e-10, atol=1e-8
- **Features:** 
    - Gentle constraint handling that doesn't interfere with radical kinetics
    - Simple iron mass balance for Fenton systems
- **Best for:** Fast radical processes

### Universal Robust System (Default)
**Bulletproof handling for complex and unknown chemical networks**

- **Solver:** Multi-tiered approach with intelligent escalation
    - Primary: Radau and BDF, tiered by simulation duration
    - Backup: Extreme stiffness handler with conservative implicit methods (Radau, BDF, and LSODA)
- **Features:**
    - Negative starting concentrations are floored to zero before the run begins; mid-run negative excursions are reported as numerical uncertainty rather than corrected
    - Rigorous Henry's Law enforcement for atmospheric species
    - Comprehensive constraint handling and error recovery
    - Smart escalation to extreme stiffness methods when needed
    - Adaptive tolerance selection based on integration duration
- **Best for:** Prebiotic networks, custom databases, highly stiff systems, unknown chemistry types

---

## Advanced Capabilities

### Long-Term Integration System

For simulations longer than one week, the SOUP model automatically switches to specialised long-term integration with adaptive time chunking and performance optimisation:

#### Adaptive Time Chunking Strategy
- **1+ month simulations:** 1-week chunks for optimal memory management
- **1+ week simulations:** 1-day chunks for balanced performance
- **Shorter periods:** 10 adaptive chunks for standard integration

#### Performance Optimization
- **Relaxed tolerances:** rtol=1e-4, atol=1e-7 optimized for speed while maintaining accuracy
- **LSODA solver:** Ideal for mixed stiff/non-stiff behavior in long-term systems
- **Large step sizes:** Allows efficient integration over extended time periods
- **Memory management:** Intelligent sampling to prevent memory overflow during long simulations

#### Real-Time Monitoring
- **Performance tracking:** Live statistics on ODE evaluations and integration speed
- **Progress reporting:** Regular updates on completion percentage and time remaining
- **Numerical health:** Automatic detection and correction of negative concentrations
- **Memory usage:** Tracks data points stored and optimizes sampling frequency

#### Automatic Activation
The long-term solver automatically activates for simulations longer than 604800 seconds (1 week).

### Extreme Stiffness Handling

When chemical systems exhibit extreme mathematical stiffness (rate constant ratios > 10¹²), the SOUP model provides specialised handling through intelligent solver escalation:

#### Stiffness Assessment
The diagnostic tool automatically calculates your system's stiffness ratio and provides guidance:

- **< 10⁶:** Low stiffness - integrates easily with standard methods
- **10⁶ - 10⁸:** Moderate stiffness - manageable with stiff solvers  
- **10⁸ - 10¹²:** High stiffness - requires careful solver selection
- **> 10¹²:** Extreme stiffness - automatically triggers specialised handling

#### Tiered Solver Escalation
When extreme stiffness is detected, the system automatically escalates through increasingly conservative approaches:

1. **Standard Stiff Methods:** High-precision Radau, BDF, and LSODA with tight tolerances
2. **Conservative Implicit:** Relaxed tolerances with limited iterations for stability over precision
3. **Ultra-Conservative:** Very small time steps with maximum stability settings

#### Automatic Stiffness Detection
The system continuously monitors for stiffness indicators and adapts the integration strategy accordingly.

### Inverse Kinetics and Parameter Estimation

The SOUP model includes a parameter estimation system for determining rate constants from experimental data:

#### Advanced Optimisation Methods
- **Differential evolution:** Comprehensive parameter space exploration with literature-constrained bounds
- **Multi-start optimisation:** Multiple independent optimisation runs to avoid local minima

#### Uncertainty Quantification
- **Bootstrap resampling:** Statistical uncertainty estimation with 50-2245 samples
- **95% confidence intervals:** Uncertainty bounds for each parameter
- **Parameter correlation analysis:** Detection of problematic parameter combinations
- **Convergence monitoring:** Rigorous convergence criteria 


#### Outputs
- **JSON results:** Machine-readable parameter estimates and statistics
- **CSV tables:** Spreadsheet-compatible parameter summaries with uncertainties
- **Convergence plots:** Visualisation of parameter distributions and correlations

#### Example Usage
Run the inverse kinetics system — note that `python inverse_kinetics.py` currently runs Option 2 by default (edit the `choice` variable in the script to pick a different one):

- **Option 1:** Test simulation only (~1 minute)
- **Option 2:** Standard optimisation (~10-30 minutes)
- **Option 3:** Multi-start optimisation (~8-12 hours)
- **Option 4:** Optimisation with a quick bootstrap, 50 samples (~2 hours)

For the full 2245-sample bootstrap used in the paper, run `seeded_inverse_kinetics.py` separately once you have best-fit parameters from Option 2 or 3.

---

## Running the Code

For more information see the [Running the Model](RunMe.md#running-the-model) page

---

## Adaptability

This framework is designed to model a wide range of chemical reaction networks. We have made every effort to keep the tool as general and flexible as possible to make customisation and updates easy. We would love to hear your feedback on how we can improve it further!

---

## Future Model Innovations

Navigate to the [Future Innovations](Future.md#future-innovations) page to find out how the model is currently being adapted.

---

## Where to find the Model

Click on the link [here](https://github.com/astroskyla/soup-prototype) to be redirected to the SOUP repository on GitHub

---