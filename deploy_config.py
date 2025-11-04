# Windtunnel Deploy Configuration
# Contains all setup functions for GitHub package installation and data downloading

import subprocess
import sys
import os
import urllib.request
from urllib.parse import quote

# Configuration constants
GITHUB_REPO_URL = "https://github.com/se04ber/WTSoftwareUtilitiesShare.git"
GITHUB_DATA_BASE = "https://raw.githubusercontent.com/se04ber/WTSoftwareUtilitiesShare/main/ExampleData/"

def install_windtunnel():
    """Install windtunnel package from GitHub"""
    try:
        print("📦 Installing windtunnel package from GitHub...")
        subprocess.check_call([
            sys.executable, "-m", "pip", "install", 
            f"git+{GITHUB_REPO_URL}"
        ])
        print("✅ windtunnel package installed successfully!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Installation failed: {e}")
        return False

def verify_installation():
    """Verify windtunnel can be imported"""
    try:
        import windtunnel as wt
        print("✅ windtunnel package verified and ready to use!")
        return True
    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False

def setup_folder_structure():
    """Create local folder structure"""
    print("📁 Setting up folder structure...")
    
    base_dir = os.getcwd()
    data_dir = os.path.join(base_dir, "Data")
    input_dir = os.path.join(data_dir, "InputData")
    param_dir = os.path.join(data_dir, "ParameterFiles")
    results_dir = os.path.join(base_dir, "Results")
    
    for directory in [data_dir, input_dir, param_dir, results_dir]:
        os.makedirs(directory, exist_ok=True)
        print(f"✅ Created directory: {directory}")
    
    return base_dir, data_dir, input_dir, param_dir, results_dir

def download_example_data(input_dir, param_dir):
    """Download example data from GitHub repository"""
    try:
        print("📥 Downloading example data from GitHub...")
        
        # Files to download with proper URL encoding
        files_to_download = {
            "InputData/Beispiel%20Umrechnung%20zur%20Kontrolle/": [
                "UBA_GA_02_04_01_000_1_001.txt.ts%230",
                "UBA_GA_02_04_01_000_1_001.txt.ts%231", 
                "UBA_GA_02_04_01_000_1_001.txt.ts%232",
                "UBA_GA_02_04_01_000_1_001.txt.ts%233",
                "UBA_GA_02_04_01_000_1_001.txt.ts%234",
                "UBA_GA_02_04_01_000_1_001.txt.ts%235"
            ],
            "ParameterFiles/": [
                "ambient_conditions_.UBA_GA.csv"
            ]
        }
        
        # Download files
        for folder, files in files_to_download.items():
            # Create the correct folder path
            if "InputData" in folder:
                # Decode the folder name for local path
                folder_name = folder.replace("InputData/", "").replace("%20", " ").replace("%2F", "/")
                folder_path = os.path.join(input_dir, folder_name)
            else:
                folder_path = param_dir  # ParameterFiles goes directly to param_dir
            
            os.makedirs(folder_path, exist_ok=True)
            
            for file in files:
                # Create URL with proper encoding
                url = GITHUB_DATA_BASE + folder + file
                
                # Decode filename for local storage
                local_filename = file.replace("%23", "#").replace("%20", " ")
                file_path = os.path.join(folder_path, local_filename)
                
                if not os.path.exists(file_path):
                    try:
                        urllib.request.urlretrieve(url, file_path)
                        print(f"✅ Downloaded: {local_filename}")
                    except Exception as e:
                        print(f"⚠️  Could not download {local_filename}: {e}")
                        print(f"   URL: {url}")
                else:
                    print(f"✅ Already exists: {local_filename}")
        
        return True
        
    except Exception as e:
        print(f"❌ Error downloading example data: {e}")
        return False

def setup_github_data(input_dir, param_dir, results_dir):
    """Complete setup for GitHub example data"""
    print("🌐 Using GitHub example data")
    
    if download_example_data(input_dir, param_dir):
        print("✅ Example data setup complete!")
        
        # Set paths to downloaded data
        path_dir = os.path.join(os.getcwd(), "Data")
        path = os.path.join(input_dir, "Beispiel Umrechnung zur Kontrolle")
        csv_file = os.path.join(param_dir, "ambient_conditions_.UBA_GA.csv")
        output_path = results_dir + "/"
        namelist = ['UBA_GA_02_04_01_000_1_001']
        
        return path_dir, path, csv_file, output_path, namelist
    else:
        print("⚠️  Using fallback paths")
        path_dir = os.path.join(os.getcwd(), "Data")
        path = os.path.join(input_dir, "Beispiel Umrechnung zur Kontrolle")
        csv_file = os.path.join(param_dir, "ambient_conditions_.UBA_GA.csv")
        output_path = results_dir + "/"
        namelist = ['UBA_GA_02_04_01_000_1_001']
        
        return path_dir, path, csv_file, output_path, namelist

def setup_local_data(input_dir, param_dir, results_dir, data_folder_name, parameter_file_name, measurement_prefix):
    """Setup for local data"""
    print("📁 Using local data")
    print("=" * 60)
    print("📋 INSTRUCTIONS FOR LOCAL DATA:")
    print("1. Put your data files in: Data/InputData/[your_folder]/")
    print("2. Put your parameter file in: Data/ParameterFiles/")
    print("3. Update the configuration variables:")
    print(f"   - DATA_FOLDER_NAME = '{data_folder_name}'")
    print(f"   - PARAMETER_FILE_NAME = '{parameter_file_name}'")
    print(f"   - MEASUREMENT_PREFIX = '{measurement_prefix}'")
    print("=" * 60)
    
    # Set paths based on user configuration
    path_dir = os.path.join(os.getcwd(), "Data")
    path = os.path.join(input_dir, data_folder_name)
    csv_file = os.path.join(param_dir, parameter_file_name)
    output_path = results_dir + "/"
    namelist = [measurement_prefix]
    
    return path_dir, path, csv_file, output_path, namelist

