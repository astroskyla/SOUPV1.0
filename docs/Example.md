# Example Simulation Run

This page demonstrates an example of how to configure the model, download or create a custom database, run the simulation, and interactively explore the output.

---

## 1. Prepare Required Files

Before running the simulation, make sure you have:

1. **Input files** (`aqueous.dat` and `atmosphere.dat`) in the `SOUP/Inputs` folder.  
2. **A database file** in the `SOUP/Inputs/Databases` folder.

---

### Adjust Input Files

Examples of filled out versions of the input files can be found below:

#### `aqueous.dat`

```dat
# Aqueous Solution Data File
species, concentration (M)
A, 2.0
B, 0.01
C, 0.0
D, 0.0
E, 0.0
F, 0.0
G, 0.0
H, 0.0
I, 0.0
J, 0.0
K, 0.0
L, 0.0
M, 0.0
Q, 0.0
O, 2.52e-4
N, 4.68e-4
S, 5.5e-4


Total Volume (L)
1.0

Total Temperature (K)
298.15

pH
7.0
```

#### `atmosphere.dat`

```dat
# Atmosphere Data File
species, partial_pressure, constant (mol/L atm)
O, 0.21, 1.2e-3
N, 0.78, 6.0e-4
S, 1.0e-3, 5.5e-2

Total pressure (atm)
1.0

Total Temperature (K)
298.15
```


### 1.1 Create a Database

To run the model you will also need a database file. 

- **Option 1**: Use the [master database](http://localhost:5000/) (requires the Flask app to be running locally — see [Database](Database.md))
- **Option 2**: Create your own custom database by selecting specific reactions and exporting the cutom database. More details regarding how to do this can be found by clicking [here](Database.md#database)

Following the download of either a a pre-made or custom database move the downloaded `.db` file into the folder:  
   ```
   SOUP/Inputs/Databases
   ```
This file will now be selectable when running the model.

---

## 2. Running the Model

Once you have placed:

- `aqueous.dat` and `atmosphere.dat` in `SOUP/Inputs`
- your database file in `SOUP/Inputs/Databases`

Run the model:

```bash
python main.py
```

You’ll be prompted to:

1. **Select the database** you just downloaded or created.  
2. **Enter the end time (seconds)** to evolve the system.  
3. **Open system (Henry's law atmospheric replenishment)?** — only asked for prebiotic-type databases.
4. Specify **species to plot** (e.g., `A, B, C, D, E`, or `all`).

`species_concentrations.md` saves automatically once the run finishes. For each species you plot you'll also be asked about a log-scale axis and whether to save the plot as an image.

---

## Outputs

The SOUP framework provides two outputs, a plot of species concentrations over time and a reaction rates report.

### Output Table

| Time (years) | C (M) | Q (M) | A (M) | H (M) | E (M) | N (M) | R (M) | J (M) | 1 (M) | O (M) | F (M) | D (M) | M (M) | L (M) | B (M) | K (M) | G (M) | S (M) | I (M) |
|--------------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| 0            | 0.0000| 0.0000| 2.0000| 0.0000| 0.0000| 0.0005| 0.0000| 0.0000| 0.0000| 0.0003| 0.0000| 0.0000| 0.0000| 0.0000| 0.0100| 0.0000| 0.0000| 0.0001| 0.0000|
| 5.1          | 0.2229| 0.0067| 0.7533| 0.0000| 0.7726| 0.0005| 0.0000| 0.0918| 0.0000| 0.0003| 0.0096| 0.0460| 0.0024| 0.0000| 0.2709| 0.0495| 0.1083| 0.0001| 0.0109|
| 10.2         | 0.2726| 0.0085| 0.1796| 0.0000| 1.0968| 0.0005| 0.0000| 0.3066| 0.0000| 0.0003| 0.0048| 0.0243| 0.0059| 0.0000| 0.0689| 0.3599| 0.1630| 0.0001| 0.0324|
| 15.3         | 0.2390| 0.0091| 0.0665| 0.0000| 1.1796| 0.0005| 0.0000| 0.3494| 0.0000| 0.0003| 0.0012| 0.0078| 0.0064| 0.0000| 0.0107| 0.5429| 0.0538| 0.0001| 0.0377|
| 20.4         | 0.2243| 0.0093| 0.0274| 0.0000| 1.2074| 0.0005| 0.0000| 0.3585| 0.0000| 0.0003| 0.0003| 0.0029| 0.0065| 0.0000| 0.0032| 0.5973| 0.0152| 0.0001| 0.0389|
| 25.5         | 0.2182| 0.0094| 0.0116| 0.0000| 1.2186| 0.0005| 0.0000| 0.3617| 0.0000| 0.0003| 0.0001| 0.0012| 0.0066| 0.0000| 0.0012| 0.6131| 0.0047| 0.0001| 0.0391|
| 30.6         | 0.2157| 0.0094| 0.0049| 0.0000| 1.2234| 0.0005| 0.0000| 0.3630| 0.0000| 0.0003| 0.0000| 0.0005| 0.0066| 0.0000| 0.0005| 0.6183| 0.0017| 0.0001| 0.0392|
| 35.7         | 0.2146| 0.0094| 0.0021| 0.0000| 1.2254| 0.0005| 0.0000| 0.3635| 0.0000| 0.0003| 0.0000| 0.0002| 0.0066| 0.0000| 0.0002| 0.6202| 0.0006| 0.0001| 0.0392|
| 40.8         | 0.2141| 0.0094| 0.0009| 0.0000| 1.2263| 0.0005| 0.0000| 0.3637| 0.0000| 0.0003| 0.0000| 0.0001| 0.0066| 0.0000| 0.0001| 0.6209| 0.0003| 0.0001| 0.0392|
| 45.9         | 0.2139| 0.0094| 0.0004| 0.0000| 1.2267| 0.0005| 0.0000| 0.3638| 0.0000| 0.0003| 0.0000| 0.0000| 0.0066| 0.0000| 0.0000| 0.6212| 0.0001| 0.0001| 0.0392|
| 50.0         | 0.2138| 0.0094| 0.0002| 0.0000| 1.2268| 0.0005| 0.0000| 0.3639| 0.0000| 0.0003| 0.0000| 0.0000| 0.0066| 0.0000| 0.0000| 0.6213| 0.0001| 0.0001| 0.0392|

---

### Output Graph

Use the dropdown or checkboxes below to explore the data.

<div id="speciesPlot" style="height:450px; max-width: 700px; margin: 0 auto;"></div>

<select id="speciesDropdown" onchange="plotSingleSpecies(this.value)" style="margin-top:10px;">
  <option value="C">C</option>
  <option value="Q">Q</option>  
  <option value="A">A</option>
  <option value="H">H</option>
  <option value="E">E</option>
  <option value="N">N</option>
  <option value="R">R</option>
  <option value="J">J</option>
  <option value="1">1</option>
  <option value="O">O</option>
  <option value="F">F</option>
  <option value="D">D</option>
  <option value="M">M</option>
  <option value="L">L</option>
  <option value="B">B</option>
  <option value="K">K</option>
  <option value="G">G</option>
  <option value="S">S</option>
  <option value="I">I</option>
</select>

<div style="margin-top:10px;">
  <label><input type="checkbox" value="C" onchange="toggleSpecies(this)"> C</label>
  <label><input type="checkbox" value="Q" onchange="toggleSpecies(this)"> Q</label>
  <label><input type="checkbox" value="A" onchange="toggleSpecies(this)"> A</label>
  <label><input type="checkbox" value="H" onchange="toggleSpecies(this)"> H</label>
  <label><input type="checkbox" value="E" onchange="toggleSpecies(this)"> E</label>
  <label><input type="checkbox" value="N" onchange="toggleSpecies(this)"> N</label>
  <label><input type="checkbox" value="R" onchange="toggleSpecies(this)"> R</label>
  <label><input type="checkbox" value="J" onchange="toggleSpecies(this)"> J</label>
  <label><input type="checkbox" value="1" onchange="toggleSpecies(this)"> 1</label>
  <label><input type="checkbox" value="O" onchange="toggleSpecies(this)"> O</label>
  <label><input type="checkbox" value="F" onchange="toggleSpecies(this)"> F</label>
  <label><input type="checkbox" value="D" onchange="toggleSpecies(this)"> D</label>
  <label><input type="checkbox" value="M" onchange="toggleSpecies(this)"> M</label>
  <label><input type="checkbox" value="L" onchange="toggleSpecies(this)"> L</label>
  <label><input type="checkbox" value="B" onchange="toggleSpecies(this)"> B</label>
  <label><input type="checkbox" value="K" onchange="toggleSpecies(this)"> K</label>
  <label><input type="checkbox" value="G" onchange="toggleSpecies(this)"> G</label>
  <label><input type="checkbox" value="S" onchange="toggleSpecies(this)"> S</label>
  <label><input type="checkbox" value="I" onchange="toggleSpecies(this)"> I</label>
</div>

<div style="margin-top:5px;">
  <label for="timeSlider">Max Time: <span id="timeValue">50.00</span> seconds</label><br/>
  <input id="timeSlider" type="range" min="1" max="50" step="1" value="50" oninput="updateTimeRange(this.value)" />
</div>

<div style="margin-top:10px;">
  <label>
    <input type="checkbox" id="logToggle" onchange="toggleLogScale(this.checked)"> 
    Logarithmic Y-axis
  </label>
</div>

<style>
  body {
    font-family: Helvetica, Arial, sans-serif;
  }
  #timeSlider {
    -webkit-appearance: none;
    width: 100%;
    height: 8px;
    border-radius: 5px;
    background: #f37b3a; /* orange from palette */
    outline: none;
  }
  #timeSlider::-webkit-slider-thumb {
    -webkit-appearance: none;
    appearance: none;
    width: 20px;
    height: 20px;
    border-radius: 50%;
    background: #000000; /* black thumb */
    cursor: pointer;
    border: 2px solid #f37b3a;
    margin-top: -6px;
  }
  #timeSlider::-moz-range-thumb {
    width: 20px;
    height: 20px;
    border-radius: 50%;
    background: #000000;
    cursor: pointer;
    border: 2px solid #f37b3a;
  }
</style>

<script src="https://cdn.plot.ly/plotly-2.26.0.min.js"></script>
<script>
// Your actual simulation data
var time = [0.0, 1.0204081632653061, 2.0408163265306123, 3.0612244897959187, 4.081632653061225, 5.1020408163265305, 6.122448979591837, 7.142857142857143, 8.16326530612245, 9.183673469387756, 10.204081632653061, 11.224489795918368, 12.244897959183675, 13.26530612244898, 14.285714285714286, 15.306122448979592, 16.3265306122449, 17.346938775510203, 18.367346938775512, 19.387755102040817, 20.408163265306122, 21.42857142857143, 22.448979591836736, 23.46938775510204, 24.48979591836735, 25.510204081632654, 26.53061224489796, 27.551020408163264, 28.571428571428573, 29.591836734693878, 30.612244897959183, 31.63265306122449, 32.6530612244898, 33.673469387755105, 34.69387755102041, 35.714285714285715, 36.734693877551024, 37.755102040816325, 38.775510204081634, 39.79591836734694, 40.816326530612244, 41.83673469387755, 42.85714285714286, 43.87755102040816, 44.89795918367347, 45.91836734693878, 46.93877551020408, 47.95918367346939, 48.9795918367347, 50.0];

var speciesData = {
    'C': [0.000000e+00, 1.371645e-02, 4.203551e-02, 9.276242e-02, 1.599589e-01, 2.229102e-01, 2.643293e-01, 2.826718e-01, 2.855147e-01, 2.805446e-01, 2.726322e-01, 2.642565e-01, 2.564862e-01, 2.496807e-01, 2.438801e-01, 2.389991e-01, 2.349148e-01, 2.315034e-01, 2.286544e-01, 2.262732e-01, 2.242808e-01, 2.226117e-01, 2.212118e-01, 2.200365e-01, 2.190489e-01, 2.182182e-01, 2.175192e-01, 2.169305e-01, 2.164344e-01, 2.160163e-01, 2.156638e-01, 2.153664e-01, 2.151154e-01, 2.149037e-01, 2.147249e-01, 2.145740e-01, 2.144466e-01, 2.143390e-01, 2.142482e-01, 2.141714e-01, 2.141066e-01, 2.140518e-01, 2.140056e-01, 2.139665e-01, 2.139335e-01, 2.139056e-01, 2.138820e-01, 2.138621e-01, 2.138453e-01, 2.138311e-01],
    'Q': [0.000000e+00, 1.807335e-03, 3.404088e-03, 4.768960e-03, 5.874415e-03, 6.716499e-03, 7.331027e-03, 7.775439e-03, 8.102669e-03, 8.350668e-03, 8.543913e-03, 8.697844e-03, 8.822430e-03, 8.924391e-03, 9.008483e-03, 9.078214e-03, 9.136264e-03, 9.184732e-03, 9.225290e-03, 9.259288e-03, 9.287826e-03, 9.311809e-03, 9.331982e-03, 9.348962e-03, 9.363265e-03, 9.375318e-03, 9.385480e-03, 9.394050e-03, 9.401281e-03, 9.407382e-03, 9.412532e-03, 9.416879e-03, 9.420550e-03, 9.423649e-03, 9.426267e-03, 9.428477e-03, 9.430345e-03, 9.431922e-03, 9.433254e-03, 9.434380e-03, 9.435331e-03, 9.436134e-03, 9.436813e-03, 9.437386e-03, 9.437871e-03, 9.438280e-03, 9.438626e-03, 9.438918e-03, 9.439165e-03, 9.439374e-03],
    'A': [2.000000e+00, 1.784196e+00, 1.555526e+00, 1.297111e+00, 1.016677e+00, 7.533338e-01, 5.440452e-01, 3.959560e-01, 2.958399e-01, 2.276597e-01, 1.796306e-01, 1.444068e-01, 1.176295e-01, 9.669654e-02, 7.999734e-02, 6.648403e-02, 5.543839e-02, 4.634470e-02, 3.881845e-02, 3.256470e-02, 2.735234e-02, 2.299742e-02, 1.935184e-02, 1.629524e-02, 1.372917e-02, 1.157262e-02, 9.758634e-03, 8.231679e-03, 6.945553e-03, 5.861720e-03, 4.947974e-03, 4.177341e-03, 3.527216e-03, 2.978612e-03, 2.515580e-03, 2.124700e-03, 1.794680e-03, 1.516008e-03, 1.280671e-03, 1.081912e-03, 9.140307e-04, 7.722228e-04, 6.524320e-04, 5.512353e-04, 4.657432e-04, 3.935161e-04, 3.324942e-04, 2.809378e-04, 2.373779e-04, 2.005736e-04],
    'H': [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
    'E': [0.000000e+00, 1.934079e-01, 3.662118e-01, 5.195885e-01, 6.548428e-01, 7.726323e-01, 8.722550e-01, 9.528340e-01, 1.015174e+00, 1.061981e+00, 1.096759e+00, 1.122752e+00, 1.142492e+00, 1.157784e+00, 1.169860e+00, 1.179553e+00, 1.187435e+00, 1.193904e+00, 1.199252e+00, 1.203696e+00, 1.207402e+00, 1.210501e+00, 1.213098e+00, 1.215277e+00, 1.217109e+00, 1.218649e+00, 1.219946e+00, 1.221038e+00, 1.221959e+00, 1.222735e+00, 1.223390e+00, 1.223942e+00, 1.224408e+00, 1.224802e+00, 1.225134e+00, 1.225414e+00, 1.225651e+00, 1.225851e+00, 1.226020e+00, 1.226163e+00, 1.226283e+00, 1.226385e+00, 1.226471e+00, 1.226544e+00, 1.226605e+00, 1.226657e+00, 1.226701e+00, 1.226738e+00, 1.226769e+00, 1.226796e+00],
    'N': [4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04, 4.680000e-04],
    'R': [0.000000e+00, 4.703302e-07, 2.120538e-06, 5.122141e-06, 9.001587e-06, 1.219312e-05, 1.336933e-05, 1.271270e-05, 1.118510e-05, 9.395579e-06, 7.692237e-06, 6.223037e-06, 5.018038e-06, 4.053203e-06, 3.287412e-06, 2.679693e-06, 2.195329e-06, 1.806866e-06, 1.493220e-06, 1.238352e-06, 1.030043e-06, 8.589228e-07, 7.177387e-07, 6.008210e-07, 5.036944e-07, 4.227945e-07, 3.552594e-07, 2.987745e-07, 2.514564e-07, 2.117642e-07, 1.784312e-07, 1.504114e-07, 1.268391e-07, 1.069947e-07, 9.027906e-08, 7.619195e-08, 6.431513e-08, 5.429838e-08, 4.584789e-08, 3.871698e-08, 3.269830e-08, 2.761749e-08, 2.332776e-08, 1.970549e-08, 1.664649e-08, 1.406293e-08, 1.188077e-08, 1.003751e-08, 8.480439e-09, 7.165063e-09],
    'J': [0.000000e+00, 5.275409e-04, 4.211929e-03, 1.680702e-02, 4.539060e-02, 9.176280e-02, 1.483759e-01, 2.035411e-01, 2.491729e-01, 2.830263e-01, 3.065796e-01, 3.224912e-01, 3.332075e-01, 3.405388e-01, 3.456921e-01, 3.494314e-01, 3.522301e-01, 3.543824e-01, 3.560741e-01, 3.574263e-01, 3.585210e-01, 3.594155e-01, 3.601515e-01, 3.607604e-01, 3.612662e-01, 3.616875e-01, 3.620395e-01, 3.623341e-01, 3.625810e-01, 3.627882e-01, 3.629624e-01, 3.631089e-01, 3.632322e-01, 3.633360e-01, 3.634235e-01, 3.634972e-01, 3.635594e-01, 3.636119e-01, 3.636562e-01, 3.636935e-01, 3.637251e-01, 3.637517e-01, 3.637742e-01, 3.637932e-01, 3.638092e-01, 3.638228e-01, 3.638342e-01, 3.638439e-01, 3.638521e-01, 3.638590e-01],
    '1': [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
    'O': [2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04, 2.520000e-04],
    'F': [0.000000e+00, 2.748329e-04, 1.612493e-03, 4.189001e-03, 7.318560e-03, 9.560656e-03, 1.004687e-02, 9.134819e-03, 7.629848e-03, 6.091785e-03, 4.754675e-03, 3.668436e-03, 2.812008e-03, 2.146376e-03, 1.633298e-03, 1.240179e-03, 9.404682e-04, 7.129270e-04, 5.407515e-04, 4.107774e-04, 3.127913e-04, 2.389422e-04, 1.832443e-04, 1.411661e-04, 1.092965e-04, 8.507856e-05, 6.660131e-05, 5.243902e-05, 4.152874e-05, 3.307775e-05, 2.649418e-05, 2.133492e-05, 1.726780e-05, 1.404250e-05, 1.146997e-05, 9.406542e-06, 7.742604e-06, 6.394073e-06, 5.296035e-06, 4.398092e-06, 3.660878e-06, 3.053464e-06, 2.551391e-06, 2.135201e-06, 1.789327e-06, 1.501243e-06, 1.260823e-06, 1.059835e-06, 8.915598e-07, 7.504899e-07],
    'D': [0.000000e+00, 2.451233e-03, 9.357288e-03, 2.121097e-02, 3.547125e-02, 4.597430e-02, 4.867697e-02, 4.482305e-02, 3.792677e-02, 3.067561e-02, 2.432080e-02, 1.918349e-02, 1.517381e-02, 1.208043e-02, 9.692160e-03, 7.835538e-03, 6.378590e-03, 5.223910e-03, 4.300185e-03, 3.555010e-03, 2.949493e-03, 2.454406e-03, 2.047487e-03, 1.711561e-03, 1.433217e-03, 1.201869e-03, 1.009082e-03, 8.480759e-04, 7.133652e-04, 6.004807e-04, 5.057633e-04, 4.262012e-04, 3.593084e-04, 3.030230e-04, 2.556321e-04, 2.157079e-04, 1.820580e-04, 1.536855e-04, 1.297545e-04, 1.095642e-04, 9.252556e-05, 7.814390e-05, 6.600278e-05, 5.575169e-05, 4.709534e-05, 3.978490e-05, 3.361057e-05, 2.839540e-05, 2.399012e-05, 2.026877e-05],
    'M': [0.000000e+00, 6.907408e-05, 2.581637e-04, 6.784107e-04, 1.406634e-03, 2.384647e-03, 3.423915e-03, 4.337753e-03, 5.039077e-03, 5.532091e-03, 5.862134e-03, 6.078983e-03, 6.222070e-03, 6.318463e-03, 6.385420e-03, 6.433553e-03, 6.469311e-03, 6.496642e-03, 6.518019e-03, 6.535036e-03, 6.548766e-03, 6.559953e-03, 6.569137e-03, 6.576720e-03, 6.583007e-03, 6.588238e-03, 6.592603e-03, 6.596252e-03, 6.599308e-03, 6.601872e-03, 6.604025e-03, 6.605835e-03, 6.607358e-03, 6.608640e-03, 6.609720e-03, 6.610630e-03, 6.611397e-03, 6.612044e-03, 6.612590e-03, 6.613051e-03, 6.613440e-03, 6.613769e-03, 6.614046e-03, 6.614280e-03, 6.614478e-03, 6.614645e-03, 6.614786e-03, 6.614905e-03, 6.615005e-03, 6.615090e-03],
    'L': [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
    'B': [1.000000e-02, 2.914387e-02, 7.400549e-02, 1.477747e-01, 2.273195e-01, 2.709277e-01, 2.594214e-01, 2.107558e-01, 1.533271e-01, 1.044281e-01, 6.891287e-02, 4.523465e-02, 3.011475e-02, 2.059824e-02, 1.456821e-02, 1.066118e-02, 8.046197e-03, 6.230005e-03, 4.921616e-03, 3.947749e-03, 3.202923e-03, 2.620882e-03, 2.158443e-03, 1.786353e-03, 1.484057e-03, 1.236631e-03, 1.032938e-03, 8.644774e-04, 7.246421e-04, 6.082217e-04, 5.110591e-04, 4.298049e-04, 3.617416e-04, 3.046481e-04, 2.567012e-04, 2.163960e-04, 1.824871e-04, 1.539398e-04, 1.298927e-04, 1.096263e-04, 9.253926e-05, 7.812784e-05, 6.596954e-05, 5.570959e-05, 4.704978e-05, 3.973928e-05, 3.356694e-05, 2.835492e-05, 2.395333e-05, 2.023584e-05],
    'K': [0.000000e+00, 8.148268e-05, 1.105387e-03, 5.895663e-03, 1.992859e-02, 4.945962e-02, 9.714163e-02, 1.594748e-01, 2.288402e-01, 2.974517e-01, 3.598681e-01, 4.134352e-01, 4.576119e-01, 4.930855e-01, 5.210824e-01, 5.429457e-01, 5.599203e-01, 5.730678e-01, 5.832519e-01, 5.911550e-01, 5.973069e-01, 6.021147e-01, 6.058893e-01, 6.088674e-01, 6.112292e-01, 6.131122e-01, 6.146210e-01, 6.158362e-01, 6.168197e-01, 6.176193e-01, 6.182721e-01, 6.188071e-01, 6.192473e-01, 6.196106e-01, 6.199113e-01, 6.201609e-01, 6.203686e-01, 6.205416e-01, 6.206862e-01, 6.208071e-01, 6.209084e-01, 6.209933e-01, 6.210647e-01, 6.211246e-01, 6.211750e-01, 6.212174e-01, 6.212531e-01, 6.212832e-01, 6.213086e-01, 6.213300e-01],
    'G': [0.000000e+00, 8.081479e-04, 6.192186e-03, 2.327117e-02, 5.847607e-02, 1.082796e-01, 1.569846e-01, 1.883150e-01, 1.963785e-01, 1.852944e-01, 1.630369e-01, 1.367299e-01, 1.110125e-01, 8.823114e-02, 6.918800e-02, 5.383109e-02, 4.172233e-02, 3.230596e-02, 2.504197e-02, 1.946100e-02, 1.517836e-02, 1.188943e-02, 9.357945e-03, 7.402964e-03, 5.886970e-03, 4.705876e-03, 3.781019e-03, 3.052945e-03, 2.476662e-03, 2.018035e-03, 1.651085e-03, 1.355952e-03, 1.117404e-03, 9.236816e-04, 7.656741e-04, 6.362725e-04, 5.299044e-04, 4.421761e-04, 3.696013e-04, 3.093993e-04, 2.593399e-04, 2.176256e-04, 1.827998e-04, 1.536772e-04, 1.292886e-04, 1.088390e-04, 9.167360e-05, 7.725131e-05, 6.512396e-05, 5.491924e-05],
    'S': [5.500000e-05, 5.500000e-05, 5.500000e-05, 5.500000e-05, 5.500000e-05, 5.500000e-05, 5.499981e-05, 5.565951e-05, 5.718712e-05, 5.897664e-05, 6.067998e-05, 6.214918e-05, 6.335418e-05, 6.431901e-05, 6.508481e-05, 6.569252e-05, 6.617689e-05, 6.656535e-05, 6.687900e-05, 6.713387e-05, 6.734217e-05, 6.751329e-05, 6.765448e-05, 6.777140e-05, 6.786852e-05, 6.794942e-05, 6.801696e-05, 6.807344e-05, 6.812076e-05, 6.816045e-05, 6.819379e-05, 6.822181e-05, 6.824538e-05, 6.826522e-05, 6.828194e-05, 6.829603e-05, 6.830790e-05, 6.831792e-05, 6.832637e-05, 6.833350e-05, 6.833952e-05, 6.834460e-05, 6.834889e-05, 6.835251e-05, 6.835557e-05, 6.835815e-05, 6.836034e-05, 6.836218e-05, 6.836374e-05, 6.836505e-05],
    'I': [0.000000e+00, 4.556512e-05, 5.584870e-04, 2.275595e-03, 5.779692e-03, 1.091873e-02, 1.672473e-02, 2.212413e-02, 2.654209e-02, 2.990402e-02, 3.237773e-02, 3.417918e-02, 3.549246e-02, 3.645462e-02, 3.716306e-02, 3.768662e-02, 3.807440e-02, 3.836198e-02, 3.857539e-02, 3.873388e-02, 3.885169e-02, 3.893938e-02, 3.900478e-02, 3.905367e-02, 3.909035e-02, 3.911796e-02, 3.913884e-02, 3.915470e-02, 3.916682e-02, 3.917612e-02, 3.918331e-02, 3.918891e-02, 3.919328e-02, 3.919673e-02, 3.919946e-02, 3.920163e-02, 3.920337e-02, 3.920478e-02, 3.920592e-02, 3.920684e-02, 3.920760e-02, 3.920822e-02, 3.920874e-02, 3.920916e-02, 3.920951e-02, 3.920980e-02, 3.921005e-02, 3.921025e-02, 3.921042e-02, 3.921056e-02]
};

// Extended color palette for 19 species
const colors = ['#173f3f', '#5c4033', '#f37b3a', '#6c9a8b', '#dde7d8', '#2c5f2d', '#8b4513', '#ff6347', '#4682b4', '#32cd32', '#ff1493', '#00ced1', '#ffd700', '#dc143c', '#9370db', '#ff8c00', '#20b2aa', '#9acd32', '#ff69b4'];
const speciesOrder = ['C', 'Q', 'A', 'H', 'E', 'N', 'R', 'J', '1', 'O', 'F', 'D', 'M', 'L', 'B', 'K', 'G', 'S', 'I'];

var layout = {
    title: '',
    font: {
        family: 'Helvetica, Arial, sans-serif',
        size: 14,
        color: '#000000',
    },
    xaxis: {
        title: {text: "Time (seconds)", standoff: 35},
        range: [0, 50],
        tickfont: {size: 14, color: '#000000'},
        zerolinecolor: '#000000',
        linecolor: '#000000',
        tickcolor: '#000000',
        rangeslider: {visible: false}
    },
    yaxis: {
        title: {text: "Concentration (M)", standoff: 20},
        rangemode: "tozero",
        tickfont: {size: 14, color: '#000000'},
        zerolinecolor: '#000000',
        linecolor: '#000000',
        tickcolor: '#000000',
        automargin: true,
        ticks: 'outside',
        ticklen: 8,
        tickwidth: 2,
        tickson: 'labels',
        ticklabelposition: 'outside left'
    },
    margin: {l: 80, r: 40, t: 20, b: 50},
    hovermode: "closest"
};

var activeTraces = {};
var currentMaxTime = 50;

function getColorForSpecies(name) {
    let idx = speciesOrder.indexOf(name);
    if(idx < 0) return '#000000';
    return colors[idx % colors.length];
}

function filterDataByTime(xArray, yArray) {
    let filteredX = [];
    let filteredY = [];
    for (let i = 0; i < xArray.length; i++) {
        if (xArray[i] <= currentMaxTime) {
            filteredX.push(xArray[i]);
            filteredY.push(yArray[i]);
        }
    }
    return {x: filteredX, y: filteredY};
}

// Make functions globally accessible
window.plotSingleSpecies = function(name) {
    if (!name) return;
    
    activeTraces = {};
    let filtered = filterDataByTime(time, speciesData[name]);
    var trace = {
        x: filtered.x,
        y: filtered.y,
        mode: 'lines', 
        name: name,
        line: {
            width: 3, 
            shape: 'spline', 
            smoothing: 1,
            color: getColorForSpecies(name)
        },
    };
    
    // Use current layout with log scale if enabled
    var currentLayout = {...layout};
    if (isLogScale) {
        currentLayout.yaxis.type = 'log';
        currentLayout.yaxis.title.text = "Concentration (M) - Log Scale";
    }
    
    Plotly.newPlot('speciesPlot', [trace], {...currentLayout, xaxis: {...currentLayout.xaxis, range: [0, currentMaxTime]}});
    document.querySelectorAll('input[type=checkbox]').forEach(cb => cb.checked = false);
    document.getElementById('speciesDropdown').value = name;
};

window.toggleSpecies = function(checkbox) {
    if (checkbox.checked) {
        let filtered = filterDataByTime(time, speciesData[checkbox.value]);
        activeTraces[checkbox.value] = {
            x: filtered.x,
            y: filtered.y,
            mode: 'lines',
            name: checkbox.value,
            line: {
                width: 2,
                shape: 'spline',
                smoothing: 1,
                color: getColorForSpecies(checkbox.value)
            }
        };
    } else {
        delete activeTraces[checkbox.value];
    }
    
    // Use current layout with log scale if enabled
    var currentLayout = {...layout};
    if (isLogScale) {
        currentLayout.yaxis.type = 'log';
        currentLayout.yaxis.title.text = "Concentration (M) - Log Scale";
    }
    
    Plotly.newPlot('speciesPlot', Object.values(activeTraces), {...currentLayout, xaxis: {...currentLayout.xaxis, range: [0, currentMaxTime]}});
    document.getElementById('speciesDropdown').value = '';
};

window.updateTimeRange = function(newMaxTime) {
    currentMaxTime = Number(newMaxTime);
    document.getElementById('timeValue').textContent = currentMaxTime.toFixed(2);

    // Use current layout with log scale if enabled
    var currentLayout = {...layout};
    if (isLogScale) {
        currentLayout.yaxis.type = 'log';
        currentLayout.yaxis.title.text = "Concentration (M) - Log Scale";
    }

    if (Object.keys(activeTraces).length > 0) {
        for (let key in activeTraces) {
            let filtered = filterDataByTime(time, speciesData[key]);
            activeTraces[key].x = filtered.x;
            activeTraces[key].y = filtered.y;
        }
        Plotly.react('speciesPlot', Object.values(activeTraces), {...currentLayout, xaxis: {...currentLayout.xaxis, range: [0, currentMaxTime]}});
    } else {
        let selected = document.getElementById('speciesDropdown').value;
        if (selected) {
            window.plotSingleSpecies(selected);
        }
    }
};

// Wait for Plotly to load and initialize
function initializePlot() {
    if (typeof Plotly === 'undefined') {
        setTimeout(initializePlot, 100);
        return;
    }
    
    // Initialize plot with species A
    window.plotSingleSpecies('A');
}

initializePlot();

// Track log scale state
var isLogScale = false;

window.toggleLogScale = function(useLogScale) {
    isLogScale = useLogScale;
    
    // Update layout
    var newLayout = {...layout};
    if (isLogScale) {
        newLayout.yaxis.type = 'log';
        newLayout.yaxis.title.text = "Concentration (M) - Log Scale";
    } else {
        newLayout.yaxis.type = 'linear';
        newLayout.yaxis.title.text = "Concentration (M)";
    }
    
    // Re-plot current data with new scale
    if (Object.keys(activeTraces).length > 0) {
        // Multi-species plot
        Plotly.react('speciesPlot', Object.values(activeTraces), {...newLayout, xaxis: {...newLayout.xaxis, range: [0, currentMaxTime]}});
    } else {
        // Single species plot
        let selected = document.getElementById('speciesDropdown').value;
        if (selected) {
            let filtered = filterDataByTime(time, speciesData[selected]);
            var trace = {
                x: filtered.x,
                y: filtered.y,
                mode: 'lines', 
                name: selected,
                line: {
                    width: 3, 
                    shape: 'spline', 
                    smoothing: 1,
                    color: getColorForSpecies(selected)
                },
            };
            Plotly.react('speciesPlot', [trace], {...newLayout, xaxis: {...newLayout.xaxis, range: [0, currentMaxTime]}});
        }
    }
};

</script>

### Reaction Rates Report

After running the model, an HTML report is generated summarizing reaction rates at different time steps.

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Chemical Reaction Rates Analysis</title>
    <style>
        :root {
            --primary-light: #DDE7D8;
            --primary-medium: #6C9A8B;
            --primary-dark: #173F3F;
            --accent-brown: #5C4033;
            --accent-orange: #F37B3A;
            --production-green: #27ae60;
            --loss-red: #e74c3c;
        }

        body {
            font-family: 'Segoe UI', -apple-system, BlinkMacSystemFont, sans-serif;
            line-height: 1.6;
            color: var(--primary-dark);
            background: white;
            margin: 0;
            padding: 0;
        }
        
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(92, 64, 51, 0.1);
            padding: 30px;
        }
        
        .header {
            text-align: center;
            margin-bottom: 40px;
            padding-bottom: 20px;
            border-bottom: 3px solid var(--primary-medium);
        }
        
        h1 {
            color: var(--primary-dark);
            font-size: 2.5em;
            margin-bottom: 10px;
        }
        
        .subtitle {
            color: var(--accent-brown);
            font-size: 1.1em;
        }
        
        .time-section {
            margin: 40px 0;
            border: 1px solid var(--primary-light);
            border-radius: 10px;
            overflow: hidden;
            box-shadow: 0 4px 6px rgba(23, 63, 63, 0.1);
        }
        
        .time-header {
            background: linear-gradient(135deg, var(--primary-medium), var(--primary-dark));
            color: white;
            padding: 15px 25px;
            font-size: 1.3em;
            font-weight: 600;
        }
        
        .species-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(350px, 1fr));
            gap: 20px;
            padding: 20px;
            background: #f8faf9;
        }
        
        .species-card {
            background: white;
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 8px rgba(23, 63, 63, 0.08);
            border-left: 4px solid var(--primary-medium);
        }
        
        .species-name {
            font-size: 1.2em;
            font-weight: 600;
            color: var(--primary-dark);
            margin-bottom: 10px;
        }
        
        .rate-info {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 15px;
        }
        
        .concentration {
            color: var(--accent-brown);
            font-size: 0.9em;
        }
        
        .net-rate {
            font-weight: 600;
            font-size: 1.1em;
            padding: 5px 12px;
            border-radius: 20px;
        }
        
        .rate-positive {
            background: #d5f4e6;
            color: var(--production-green);
        }
        
        .rate-negative {
            background: #fadbd8;
            color: var(--loss-red);
        }
        
        .rate-zero {
            background: var(--primary-light);
            color: var(--accent-brown);
        }
        
        .reactions-section {
            margin-top: 15px;
        }
        
        .reaction-item {
            margin: 8px 0;
            padding: 15px;
            border-radius: 8px;
            font-size: 0.9em;
            border-left: 4px solid;
        }
        
        .production {
            background: linear-gradient(90deg, #e8f8f2, #f0fbf5);
            border-left-color: var(--production-green);
        }
        
        .loss {
            background: linear-gradient(90deg, #fdeaea, #fef2f2);
            border-left-color: var(--loss-red);
        }
        
        .reaction-rate {
            font-weight: 600;
            float: right;
            padding: 2px 8px;
            border-radius: 4px;
            font-size: 0.85em;
        }
        
        .production .reaction-rate {
            background: var(--production-green);
            color: white;
        }
        
        .loss .reaction-rate {
            background: var(--loss-red);
            color: white;
        }
        
        .reaction-equation {
            color: var(--accent-brown);
            font-size: 0.85em;
            margin-top: 8px;
            font-family: 'Courier New', monospace;
            background: rgba(255,255,255,0.5);
            padding: 5px 8px;
            border-radius: 4px;
        }
        
        .chemical-equation {
            color: var(--primary-dark);
            font-weight: 500;
            margin-top: 5px;
            font-size: 0.9em;
        }
        
        .summary-stats {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 30px 0;
            padding: 20px;
            background: linear-gradient(135deg, var(--primary-light), #eef4ed);
            border-radius: 10px;
        }
        
        .stat-item {
            text-align: center;
        }
        
        .stat-value {
            font-size: 1.5em;
            font-weight: 600;
            color: var(--primary-dark);
        }
        
        .stat-label {
            color: var(--accent-brown);
            font-size: 0.9em;
        }
        
        .toggle-button {
            background: linear-gradient(135deg, var(--primary-medium), var(--primary-dark));
            color: white;
            border: none;
            padding: 10px 16px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 0.9em;
            margin-top: 10px;
            transition: all 0.3s ease;
        }
        
        .toggle-button:hover {
            background: linear-gradient(135deg, var(--primary-dark), var(--accent-brown));
            transform: translateY(-1px);
        }
        
        .reactions-details {
            display: none;
            margin-top: 15px;
        }
        
        .reactions-details.show {
            display: block;
        }
        
        .arrow {
            color: var(--accent-orange);
            font-weight: bold;
            padding: 0 5px;
        }
    </style>
    <script>
        function toggleReactions(speciesId) {
            const details = document.getElementById(speciesId + '-details');
            const button = document.getElementById(speciesId + '-button');
            
            if (details.classList.contains('show')) {
                details.classList.remove('show');
                button.textContent = button.textContent.replace('Hide', 'Show');
            } else {
                details.classList.add('show');
                button.textContent = button.textContent.replace('Show', 'Hide');
            }
        }
    </script>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Reaction Rates</h1>
            <div class="subtitle">Generated on September 24, 2025 at 10:35 AM</div>
        </div>

        <div class="time-section">
            <div class="time-header">
                Time: 0.00 seconds (0.000 hours)
            </div>
            
            <div class="summary-stats">
                <div class="stat-item">
                    <div class="stat-value">6</div>
                    <div class="stat-label">Active Species</div>
                </div>
                <div class="stat-item">
                    <div class="stat-value">2.00e-01</div>
                    <div class="stat-label">Max Production</div>
                </div>
                <div class="stat-item">
                    <div class="stat-value">2.12e-01</div>
                    <div class="stat-label">Max Loss</div>
                </div>
                <div class="stat-item">
                    <div class="stat-value">4.34e-01</div>
                    <div class="stat-label">Total Activity</div>
                </div>
            </div>
            
            <div class="species-grid">
                <div class="species-card">
                    <div class="species-name">C</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 0.00e+00 M</div>
                        <div class="net-rate rate-positive">+1.00e-02</div>
                    </div>

                    <button class="toggle-button" id="species_0_t_0-button" onclick="toggleReactions('species_0_t_0')">
                        Show Reactions (1 contributing)
                    </button>
                    
                    <div class="reactions-details" id="species_0_t_0-details">
                        <div class="reactions-section">
                            <div class="reaction-item production">
                                <div>
                                    Reaction 1 (forward)
                                    <span class="reaction-rate">+1.00e-02</span>
                                </div>
                                <div class="chemical-equation">A + B <span class="arrow">→</span> 2 B + C</div>
                                <div class="reaction-equation">+ +1.0 * 0.5*A*B</div>
                            </div>
                        </div>
                    </div>
                </div>

                <div class="species-card">
                    <div class="species-name">Q</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 0.00e+00 M</div>
                        <div class="net-rate rate-positive">+1.87e-03</div>
                    </div>

                    <button class="toggle-button" id="species_1_t_0-button" onclick="toggleReactions('species_1_t_0')">
                        Show Reactions (1 contributing)
                    </button>
                    
                    <div class="reactions-details" id="species_1_t_0-details">
                        <div class="reactions-section">
                            <div class="reaction-item production">
                                <div>
                                    Reaction 12 (forward)
                                    <span class="reaction-rate">+1.87e-03</span>
                                </div>
                                <div class="chemical-equation">A + N <span class="arrow">→</span> Q</div>
                                <div class="reaction-equation">+1.0 * 2.0*A*N</div>
                            </div>
                        </div>
                    </div>
                </div>

                <div class="species-card">
                    <div class="species-name">A</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 2.00e+00 M</div>
                        <div class="net-rate rate-negative">-2.12e-01</div>
                    </div>

                    <button class="toggle-button" id="species_2_t_0-button" onclick="toggleReactions('species_2_t_0')">
                        Show Reactions (2 contributing)
                    </button>
                    
                    <div class="reactions-details" id="species_2_t_0-details">
                        <div class="reactions-section">
                            <div class="reaction-item loss">
                                <div>
                                    Reaction 3 (forward)
                                    <span class="reaction-rate">-2.00e-01</span>
                                </div>
                                <div class="chemical-equation">A <span class="arrow">→</span> E</div>
                                <div class="reaction-equation">+ -1.0 * 0.1*A</div>
                            </div>

                            <div class="reaction-item loss">
                                <div>
                                    Reaction 1 (forward)
                                    <span class="reaction-rate">-1.00e-02</span>
                                </div>
                                <div class="chemical-equation">A + B <span class="arrow">→</span> 2 B + C</div>
                                <div class="reaction-equation">-1.0 * 0.5*A*B</div>
                            </div>
                        </div>
                    </div>
                </div>

                <div class="species-card">
                    <div class="species-name">H</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 0.00e+00 M</div>
                        <div class="net-rate rate-zero">+0.00e+00</div>
                    </div>
                </div>

                <div class="species-card">
                    <div class="species-name">E</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 0.00e+00 M</div>
                        <div class="net-rate rate-positive">+2.00e-01</div>
                    </div>

                    <button class="toggle-button" id="species_4_t_0-button" onclick="toggleReactions('species_4_t_0')">
                        Show Reactions (1 contributing)
                    </button>
                    
                    <div class="reactions-details" id="species_4_t_0-details">
                        <div class="reactions-section">
                            <div class="reaction-item production">
                                <div>
                                    Reaction 3 (forward)
                                    <span class="reaction-rate">+2.00e-01</span>
                                </div>
                                <div class="chemical-equation">A <span class="arrow">→</span> E</div>
                                <div class="reaction-equation">+ +1.0 * 0.1*A</div>
                            </div>
                        </div>
                    </div>
                </div>

                <div class="species-card">
                    <div class="species-name">N</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 4.68e-04 M</div>
                        <div class="net-rate rate-zero">+0.00e+00</div>
                    </div>
                </div>

                <div class="species-card">
                    <div class="species-name">R</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 0.00e+00 M</div>
                        <div class="net-rate rate-zero">+0.00e+00</div>
                    </div>
                </div>

                <div class="species-card">
                    <div class="species-name">J</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 0.00e+00 M</div>
                        <div class="net-rate rate-zero">+0.00e+00</div>
                    </div>
                </div>

                <div class="species-card">
                    <div class="species-name">1</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 0.00e+00 M</div>
                        <div class="net-rate rate-zero">-0.00e+00</div>
                    </div>
                </div>

                <div class="species-card">
                    <div class="species-name">O</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 2.52e-04 M</div>
                        <div class="net-rate rate-zero">+0.00e+00</div>
                    </div>
                </div>

                <div class="species-card">
                    <div class="species-name">F</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 0.00e+00 M</div>
                        <div class="net-rate rate-zero">+0.00e+00</div>
                    </div>
                </div>

                <div class="species-card">
                    <div class="species-name">D</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 0.00e+00 M</div>
                        <div class="net-rate rate-zero">+0.00e+00</div>
                    </div>
                </div>

                <div class="species-card">
                    <div class="species-name">M</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 0.00e+00 M</div>
                        <div class="net-rate rate-positive">+3.78e-05</div>
                    </div>

                    <button class="toggle-button" id="species_12_t_0-button" onclick="toggleReactions('species_12_t_0')">
                        Show Reactions (1 contributing)
                    </button>
                    
                    <div class="reactions-details" id="species_12_t_0-details">
                        <div class="reactions-section">
                            <div class="reaction-item production">
                                <div>
                                    Reaction 11 (forward)
                                    <span class="reaction-rate">+3.78e-05</span>
                                </div>
                                <div class="chemical-equation">B + O <span class="arrow">→</span> C + M</div>
                                <div class="reaction-equation">+1.0 * 15.0*B*O</div>
                            </div>
                        </div>
                    </div>
                </div>

                <div class="species-card">
                    <div class="species-name">L</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 0.00e+00 M</div>
                        <div class="net-rate rate-zero">+0.00e+00</div>
                    </div>
                </div>

                <div class="species-card">
                    <div class="species-name">B</div>
                    <div class="rate-info">
                        <div class="concentration">Conc: 1.00e-02 M</div>
                        <div class="net-rate rate-positive">+9.96e-03</div>
                    </div>

                    <button class="toggle-button" id="species_14_t_0-button" onclick="toggleReactions('species_14_t_0')">
                        Show Reactions (1 contributing)
                    </button>
                    
                    <div class="reactions-details" id="species_14_t_0-details">
                        <div class="reactions-section">
                            <div class="reaction-item production">
                                <div>
                                    Reaction 1 (forward)
                                    <span class="reaction-rate">+9.96e-03</span>
                                </div>
                                <div class="chemical-equation">A + B <span class="arrow">→</span> 2 B + C</div>
                                <div class="reaction-equation">+ +1.0 * 0.5*A*B</div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Additional time sections would be added here for multiple time points -->
        
    </div>
</body>
</html>

---
