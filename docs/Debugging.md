# Debugging the SOUP Model

This guide helps you diagnose and fix common issues that might pop up when running the SOUP Model.

---

## Quick Diagnostic Tool

**Before running long simulations, always run the diagnostic tool first:**

```bash
python diagnose.py
```

This will check your database, initial conditions, and system setup for potential problems.

### What the Diagnostic Shows

The diagnostic tool provides a sort of health check:

- **Database File Check**: Confirms the database file exists and loads correctly
- **Rate Constant Analysis**: Identifies extreme values that could cause stiffness
- **Initial Conditions Check**: Finds negative, zero, or unrealistic concentrations  
- **Henry's Law Validation**: Verifies atmospheric species equilibrium setup
- **Stiffness Assessment**: Calculates the mathematical difficulty of your system
- **Species Balance Analysis**: Checks for elemental mass-balance issues across your reaction network
- **ODE Evaluation Test**: Tests if your equations work at the starting point

---

## Common Issues and Solutions

### **Solver Crashes or "All Solver Attempts Failed"**

**Symptoms:**

- Integration stops early
- "All bulletproof solver attempts failed" message
- Very slow progress then crashes 

**Diagnosis Steps:**

1. Check stiffness ratio in diagnostic output
2. Look for extreme rate constants (>1e8 or <1e-10)
3. Verify initial conditions are reasonable

**Solutions:**
```bash
# 1. Run diagnostic first
python diagnose.py

# 2. If stiffness ratio > 1e10, try shorter time periods
# Instead of 1 year, try 1 week first

# 3. Check for problematic reactions in your database
# Look for rate constants spanning too many orders of magnitude
```

**Prevention:**

- Always test short periods before long simulations
- Use the bulletproof rate constant method (default)
- Ensure your database has physically reasonable rate constants

---

### **Negative Concentrations**

**Symptoms:**

- Warning messages about negative species
- Unphysical results
- Species "going negative" during simulation

**What This Means:**
Negative concentrations are physically impossible therefore the solver is taking steps that overshoot and create unrealistic values.

**Automatic Fixes:**
Negative, NaN, or Inf starting concentrations are floored to zero before the run begins, and any NaN/Inf values left in the final results are replaced with 1e-15 M. Mid-run negative excursions aren't actively corrected — they're reported as numerical uncertainty rather than clamped.

**Manual Check:**
```python
# In your results, check for any remaining negatives:
negative_mask = results < 0
if np.any(negative_mask):
    print("Negative concentrations found!")
    # The system should have prevented this
```

---

### **Very Slow Integration**

**Symptoms:**

- Integration takes hours for short time periods
- Computer becomes unresponsive
- Many tiny time steps reported

**Causes:**

1. **Extreme Stiffness**: Rate constants span too many orders of magnitude
2. **Tight Tolerances**: Solver trying to be too accurate
3. **Poor Initial Conditions**: Starting from unrealistic concentrations

**Solutions:**

**For Short-Term Issues:**
```bash
# Check if your system is extremely stiff
python diagnose.py

# Look for "Stiffness ratio: X.XeYY"
# If ratio > 1e12, expect slow integration
```

**For Long-Term Integration:**
```python
# The system automatically switches to optimized long-term solver
# for simulations > 1 week

# Test progression:
# 1 hour → 1 day → 1 week → 1 month → 1 year
```

---

### **Computer Overheating/Crashes**

**Symptoms:**

- Computer fan running loudly
- System becomes very slow
- Unexpected shutdowns during integration

**Immediate Actions:**

1. **Stop the simulation** (Ctrl+C)
2. **Cool down your computer**
3. **Test with much shorter times first**

**Prevention:**
```bash
# Always test performance first:
python main.py
# Enter end time: 3600  (1 hour test)

# If 1 hour takes more than 1 minute, don't try longer periods
# Your system may be too stiff for long integration
```

---

### **Unrealistic Results**

**Symptoms:**

- Concentrations that don't make chemical sense
- Species not behaving as expected
- Equilibrium species changing when they shouldn't

**Check Henry's Law Species:**
These should NEVER change from their atmospheric equilibrium.

**Validation Steps:**
There's no specific console message to watch for here — instead, check whether the concentration values for your Henry's Law–governed species stay constant across the output timesteps in `species_concentrations.md`.

---

## Understanding System Messages

### **Good Signs:**
```
Success with Radau!  (or BDF, etc. - the solver method name varies)
ODE evaluation successful
```

### **Warning Signs:**
```
WARNING: 3 extremely high rates (>1e8)
Negative concentrations: 5
MODERATELY STIFF - should be manageable
EXTREMELY STIFF - expect solver difficulties
```

###  **Problem Signs:**
```
All universal robust solver attempts failed
Conservative retry also failed
Error in <species>: <exception message>
```

---

## Performance Testing Strategy

### **Step-by-Step Approach:**

**Phase 1: Quick Tests (< 5 minutes each)**
```bash
# Test these durations in order:
python main.py
# 1. End time: 3600     (1 hour)
# 2. End time: 86400    (1 day)  
# 3. End time: 604800   (1 week)
```

**Phase 2: Medium Tests (30-60 minutes each)**
```bash
# Only if Phase 1 succeeds quickly:
python main.py  
# 4. End time: 2592000   (1 month)
# 5. End time: 7776000   (3 months)
```

**Phase 3: Long-Term Test**
```bash
# Only if Phase 2 completes in reasonable time:
python main.py
# 6. End time: 31536000  (1 year)
```

**Performance Benchmarks:**

- 1 hour should complete in < 30 seconds
- 1 day should complete in < 2 minutes
- 1 week should complete in < 30 minutes
- If any step takes 10x longer than expected, stop and investigate

---

## Advanced Debugging

### **Examining Rate Constants**

```python
# After running diagnose.py, look for:
print("Rate constant range: X.XXe-XX to X.XXe+XX")
print("Stiffness ratio: X.XXe+XX")

# Stiffness interpretation:
# < 1e6:  Easy to integrate
# 1e6-1e8: Manageable 
# 1e8-1e12: Difficult but possible
# > 1e12: Very challenging, expect slow integration
```


### **Custom Diagnostic Checks**

Add these checks to your analysis:
```python
# Check concentration ranges
min_conc = np.min(results)
max_conc = np.max(results)
print(f"Concentration range: {min_conc:.2e} to {max_conc:.2e}")

# Check for species that changed significantly
for i, species in enumerate(database_species):
    initial = initial_conditions[i]
    final = results[-1][i]
    change_ratio = abs(final - initial) / max(initial, 1e-12)
    if change_ratio > 10:
        print(f"{species}: {initial:.2e} -> {final:.2e} (changed {change_ratio:.1f}x)")
```

---

## Getting Help

### **Information to Collect:**

Please feel free to get in touch if issues persist. To do this please email skyla.white@unibe.ch. When reporting issues, include:

1. **Diagnostic Output**: Full output from `python diagnose.py`
2. **Error Messages**: Exact text of error messages
3. **System Specs**: Computer type, RAM, operating system
4. **Test Results**: How long different time periods took
5. **Database Info**: Which database you're using
6. **Input Files**: Your atmosphere.dat and aqueous.dat contents

### **Common Solutions Summary:**

| Problem | Quick Fix | Long-term Solution |
|---------|-----------|------------------|
| Solver crashes | Reduce time span | Check database for extreme rates |
| Negative concentrations | Use bulletproof system (default) | Verify initial conditions |
| Very slow | Test shorter periods | Consider system stiffness |
| Computer overheating | Stop simulation, test 1 hour first | Optimize rate constants |
| Unrealistic results | Check Henry's Law species | Validate input data |

---

**Remember: Always run `python diagnose.py` first. It catches most issues before they cause problems!**