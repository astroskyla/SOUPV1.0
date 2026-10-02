# SOUPV1.0.0

## Description
This repository contains the SOUPV1.0.0 model and database. 

Note that for some scripts, code was cleaned up and coments formatted using Github Copilot. These have been thoroughly checked by SBWhite to ensure information is correct. 

## Setup Instructions

### Prerequisites 

1. **Visual Studio Build Tools**:
    - Download and install the [Visual Studio Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/).
    - During installation, ensure you select the "Desktop development with C++" workload.

### Setup

1. **Clone the repository**:
    ```sh
    git clone https://github.com/astroskyla/SOUPV1.0.git
    cd SOUPV1.0
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

5. **Navigate to the website** to find out how to run the model and manipulate the database. The website can be accessed by clicking [here](https://astroskyla.github.io/SOUPV1.0/)**