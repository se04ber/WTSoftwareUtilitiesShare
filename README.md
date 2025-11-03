# WTSoftwareUtilitiesShare

Repository for wind tunnel analysis utilities and deployable notebooks.

## Quick Start Options

### Option 1: Jupyter Server (Recommended for Students)

**Best for:** Using the notebook on a (f.e. the UHH) Jupyter Notebook online server without installing the full package.

**Requirements:** Jupyter Notebook server access

**Steps:**
1. Upload two files to your Jupyter server:
   - `Example_PointConc_Analysis_200325_Analysis_deploy.ipynb`
   - `deploy_config.py`

2. Open the notebook and run the first cell. The notebook will:
   - Automatically install the `windtunnel` package from GitHub
   - Create the necessary folder structure (`Data/InputData/`, `Data/ParameterFiles/`, `Results/`)
   - Download example data

3. **Using local data** (optional):
   - Upload your data files to `Data/InputData/your_folder_name/`
   - Upload your parameter CSV file to `Data/ParameterFiles/`
   - In **Cell 2** (Data Setup and Configuration), set:
     ```python
     USE_GITHUB_EXAMPLE_DATA = False
     DATA_FOLDER_NAME = "your_folder_name"  # Folder name inside InputData
     PARAMETER_FILE_NAME = "your_parameter_file.csv"  # CSV file name in ParameterFiles
     MEASUREMENT_PREFIX = "your_measurement_prefix"  # Prefix of your .ts files
     ```

**Note:** The notebook handles all setup automatically. No manual installation needed for this option.

---

### Option 2: Windows Local Installation

**Best for:** Installing and running locally on Windows in a clean environment with two clicks using helper scripts.

**Requirements:**
- Python >3.6 installed
- (Git installed (or download repository manually as ZIP from GitHub))

**Steps:**

1. **First-time setup** (run once):
   ```
   Run windows_deploy\setup_WTConc2.bat
   ```
   This creates a virtual environment `WTConc2` and installs all dependencies from requirements.txt.

2. **Start Jupyter Notebook** (every time you want to use it):
   ```
   Run windows_deploy\Start_WTConc2_Notebook.bat
   ```

3. Open any notebook f.e. `Example_PointConc_Analysis_200325_Analysis.ipynb` in the opened jupyter folder tree.

**Note:** If Git is not available, download the repository as a ZIP from GitHub, extract it, and the scripts will work the same way. In case the environment is not automatically working you may need to 1) Choose the kernel in the script titled "WTConc2" manually or if something in the installations of the packages went wrong you need to in the worst case change in setup_WTConc2.bat the path PY_CMD manually to the correct python.exe location.

---

### Option 3: Full Development Setup

**Best for:** Editing the project, contributing, or using on Linux/Windows locally with full control.

**Requirements:**
- Python >3.6 installed
- Git installed

**Steps:**

1. **Clone the repository:**
   ```bash
   git clone https://github.com/se04ber/WTSoftwareUtilitiesShare.git
   cd WTSoftwareUtilitiesShare
   ```

2. **Create a virtual environment:**
   ```bash
   python -m venv WTConc2
   ```

3. **Activate the environment:**
   - **Windows:**
     ```bash
     WTConc2\Scripts\activate
     ```
   - **Linux/Mac:**
     ```bash
     source WTConc2/bin/activate
     ```

4. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

5. **Install the package (optional, for development):**
   ```bash
   pip install -e .
   ```

6. **Launch Jupyter:**
   ```bash
   jupyter notebook
   ```

---

## File Structure

Notebooks:
- `Example_PointConc_Analysis_200325_stepByStep.ipynb` - Tutorial notebook, going through each processing,calculation per cell and explaining
- `Example_PointConc_Analysis_200325_Analysis_deploy.ipynb` - Main server deployable analysis notebook
- `Example_PointConc_Analysis_200325_Analysis.ipynb` - Analysis notebook
- `Example_PointConc_Analysis_200325_Map.ipynb` - Maps plotting notebook

Project:
- `windtunnel/` - Main Python package

Deploying:
- `requirements.txt` - Python package dependencies
- `deploy_config.py` - Configuration and setup functions for the jupyter server deploy notebook
- `windows_deploy/` - Windows setup scripts
  - `setup_WTConc2.bat` - Creates virtual environment and installs dependencies
  - `Start_WTConc2_Notebook.bat` - Launches Jupyter Notebook

---

## Troubleshooting and ideas
Contact:
In case there are further problems or ideas feel free to contact me at: Sabrina.ebert@proton.me :)