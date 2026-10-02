# SOUPV1.0.0

## Description
This repository contains the SOUPV1.0.0 model and database. 

##
**Files:**
1. seeded_inverse_kinetics.py | used to produce the bootstrap resampling results from the zymonic acid case study. To run this you need to first run "inverse_kinetics.py" to get best fit parameters. These must then be hardcoded into this script. 
2. inverse_kinetics.py | used to get best fit parameters for the rate constants in the zymonic acid case study.
3. data.py | used to read all input files (apart from the database) required to run the SOUP model.
4. diagnose.py | tries to trouble shoot any potential problems from running the SOUP model. Run this if you are having any issues and it should advise on what the problem is and how to fix the problem.
5. differentials.py | chooses a solving method based on the reaction network. If something is going wrong, this file is probably the culprit. 
6. extreme_stuffness_solver.py | this is the solver used for extrememly stiff systems. 
7. improved_reaction_rates_logger.py | this creates a html report of reaction fluxes over time.
8. long_term_solver.py | this script is for simulations that would otherwise take too long. It uses adaptive time chunking to make the system easier to solve. 
9. main.py | this is where the main model is run from. Again, if you have problems, this file is probably to blame.
10. photoaquation_rate.py | this is the file that generates the prebiotic retrieval case study from the paper.
11. plotting.py | this is the plotting script called in main following the running of the full model. Specific plots generated for the paper are in the plotting folder. 


Note that for some scripts, code was cleaned up and doc strings formatted using Github Copilot. These have been thoroughly checked by SBWhite to ensure information is correct. 

## Setup Instructions

### Prerequisites 

1. **Visual Studio Build Tools**:
    - Download and install the [Visual Studio Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/).
    - During installation, ensure you select the "Desktop development with C++" workload.

### Setup

1. **Clone the repository**:
    ```sh
    git clone https://github.com/astroskyla/soup-prototype.git
    cd soup-prototype
    ```

2. **Install [uv](https://docs.astral.sh/uv/getting-started/installation/)** (if you don't already have it), then create the environment and install all dependencies in one step:
    ```sh
    uv sync
    ```
    This creates a `.venv` and installs everything pinned in `uv.lock` — no separate `requirements.txt` or `setup.bat` needed.

3. **Activate the environment**:

    #### For macOS or Linux:
    ```sh
    source .venv/bin/activate
    ```

    #### For Windows:
    ```sh
    .venv\Scripts\activate
    ```

4. **Launch the website**:

    #### For Windows:
    ```sh
    cd Flask
    start /B python flask-app.py
    cd ..
    start /B mkdocs serve > NUL 2>&1
    ```

    #### For macOS or Linux:
    ```sh
    cd Flask
    nohup python flask-app.py > /dev/null 2>&1 &
    cd ..
    nohup mkdocs serve > /dev/null 2>&1 &
    ```

5. **Navigate to the website** to find out how to run the model and manipulate the database.