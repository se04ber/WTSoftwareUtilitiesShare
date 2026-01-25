import numpy as np
import pandas as pd
import csv
import os
import glob
import re
from pathlib import Path

# Note: plotting/STL helpers lazy-import heavy deps (matplotlib, numpy-stl) inside functions.
def load_avg_file(filepath):
    """
    Load and parse an avg file to extract metadata and concentration data.
    
    Args:
        filepath (str): Path to the avg file
        
    Returns:
        dict: Parsed data containing metadata and concentration values
    """
    try:
        with open(filepath, 'r', encoding='utf-8') as file:
            content = file.read()
        
        # Extract filename from path and clean it
        full_filename = Path(filepath).name
        # Extract only the part starting with _avg_ and ending with .txt.ts#0
        filename_match = re.search(r'(_avg_.*?\.txt\.ts#\d+)', full_filename)
        filename = filename_match.group(1) if filename_match else full_filename
        
        # Parse metadata from header
        metadata = {}
        # Extract geometric scale
        scale_match = re.search(r'geometric scale: 1:(\d+\.?\d*)', content)
        if scale_match:
            metadata['geometric_scale'] = float(scale_match.group(1))
        
        # Extract measurement positions (in mm)
        x_match = re.search(r'x \(measurement relativ to source\): ([\d\.-]+) \[mm\]', content)
        y_match = re.search(r'y \(measurement relativ to source\): ([\d\.-]+) \[mm\]', content)
        z_match = re.search(r'z \(measurement relativ to source\): ([\d\.-]+) \[mm\]', content)
        if x_match and y_match and z_match:
            # Convert to full scale (mm to m, then scale up)
            scale = metadata.get('geometric_scale', 200.0)
            metadata['x_fs'] = float(x_match.group(1)) / 1000.0 * scale  # mm to m, then scale
            metadata['y_fs'] = float(y_match.group(1)) / 1000.0 * scale
            metadata['z_fs'] = float(z_match.group(1)) / 1000.0 * scale

        #wtref = re.search(r'\(wtref): ([\d\.-]+) \[m/s\]', content)
        wtref_fs = re.search(r'full scale wtref: ([\d\.-]+)\[m/s\]', content)
        metadata['wtref_fs']  = float(wtref_fs.group(1))
        
        # Find the data line (last non-comment line)
        lines = content.strip().split('\n')
        data_line = None
        for line in reversed(lines):
            if line.strip() and not line.strip().startswith('#'):
                data_line = line.strip()
                break
        
        if data_line:
            # Parse concentration data
            values = data_line.split()
            if len(values) >= 3:
                metadata['c_star'] = float(values[0])
                metadata['net_concentration'] = float(values[1])
                metadata['full_scale_concentration'] = float(values[2])
        
        metadata['filename'] = filename
        return metadata
        
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None

def load_stats_file(filepath):
    """
    Load and parse a stats file to extract metadata and statistical concentration data.
    
    Args:
        filepath (str): Path to the stats file
        
    Returns:
        dict: Parsed data containing metadata and statistical values
    """
    try:
        with open(filepath, 'r', encoding='utf-8') as file:
            content = file.read()
        
        # Extract filename from path and clean it
        full_filename = Path(filepath).name
        # Extract only the part starting with _stats_ and ending with .txt.ts#0
        filename_match = re.search(r'(_stats_.*?\.txt\.ts#\d+)', full_filename)
        filename = filename_match.group(1) if filename_match else full_filename
        
        # Parse metadata from header (same as avg file)
        metadata = {}
        # Extract geometric scale
        scale_match = re.search(r'geometric scale: 1:(\d+\.?\d*)', content)
        if scale_match:
            metadata['geometric_scale'] = float(scale_match.group(1))
        
        # Extract measurement positions (in mm)
        x_match = re.search(r'x \(measurement relativ to source\): ([\d\.-]+) \[mm\]', content)
        y_match = re.search(r'y \(measurement relativ to source\): ([\d\.-]+) \[mm\]', content)
        z_match = re.search(r'z \(measurement relativ to source\): ([\d\.-]+) \[mm\]', content)
        if x_match and y_match and z_match:
            # Convert to full scale (mm to m, then scale up)
            scale = metadata.get('geometric_scale', 200.0)
            metadata['x_fs'] = float(x_match.group(1)) / 1000.0 * scale  # mm to m, then scale
            metadata['y_fs'] = float(y_match.group(1)) / 1000.0 * scale
            metadata['z_fs'] = float(z_match.group(1)) / 1000.0 * scale
        
        #wtref = re.search(r'\(wtref): ([\d\.-]+) \[m/s\]', content)
        wtref_fs = re.search(r'full scale wtref: ([\d\.-]+)\[m/s\]', content)
        metadata['wtref_fs']  = float(wtref_fs.group(1))
        
        # Find the statistical data line (should be the first non-comment line after header)
        lines = content.strip().split('\n')
        data_line = None
        for line in lines:
            if line.strip() and not line.strip().startswith('#'):
                data_line = line.strip()
                break
        
        if data_line:
            # Parse statistical data - should have 12 values
            # Order: means(3), percentiles95(3), percentiles5(3), peak2mean(3)
            values = data_line.split()
            if len(values) >= 12:
                # Means
                metadata['c_star'] = float(values[0])
                metadata['net_concentration'] = float(values[1])
                metadata['full_scale_concentration'] = float(values[2])
                
                # 95th percentiles
                metadata['percentiles95_c_star'] = float(values[3])
                metadata['percentiles95_net_concentration'] = float(values[4])
                metadata['percentiles95_full_scale_concentration'] = float(values[5])
                
                # 5th percentiles
                metadata['percentiles5_c_star'] = float(values[6])
                metadata['percentiles5_net_concentration'] = float(values[7])
                metadata['percentiles5_full_scale_concentration'] = float(values[8])
                
                # Peak2Mean ratios
                metadata['peak2mean_c_star'] = float(values[9])
                metadata['peak2mean_net_concentration'] = float(values[10])
                metadata['peak2mean_full_scale_concentration'] = float(values[11])
        
        metadata['filename'] = filename
        return metadata
        
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None

def combine_to_csv(file_names, base_path, file_type='avg', output_filename='combined_data.csv'):
    """
    Combine specific avg or stats files into a single CSV file.
    
    Args:
        file_names (list): Array of specific filenames to process
        base_path (str): Base directory path containing the files
        file_type (str): Type of files to process ('avg' or 'stats')
        output_filename (str): Name of the output CSV file
        
    Returns:
        pd.DataFrame: Combined data as a pandas DataFrame
    """
    
    # Define subdirectory and load function based on file type
    if file_type.lower() == 'avg':
        #subdirectory = 'Point_Data_avg'
        load_function = load_avg_file
    elif file_type.lower() == 'stats':
        #subdirectory = 'Point_Data_stats'
        load_function = load_stats_file
    else:
        raise ValueError("file_type must be 'avg' or 'stats'")
    
    # Build full file paths from the provided filenames
    file_paths = []
    for filename in file_names:
        # Search for the file in the subdirectory structure
        search_pattern = os.path.join(base_path, '**', filename)
        matches = glob.glob(search_pattern, recursive=True)
        
        if matches:
            file_paths.extend(matches)
        else:
            print(f"Warning: File '{filename}' not found in {base_path}{filename}")
    
    if not file_paths:
        print(f"No valid {file_type} files found from the provided list")
        return None
    
    print(f"Processing {len(file_paths)} {file_type} files")
    
    # Process all files
    all_data = []
    for filepath in file_paths:
        print(f"Processing: {filepath}")
        data = load_function(filepath)
        if data:
            all_data.append(data)
    
    if not all_data:
        print("No valid data found in any files")
        return None
    
    # Create DataFrame
    df = pd.DataFrame(all_data)
    
    # Define column order based on file type
    if file_type.lower() == 'avg':
        columns = [
            'filename', 'x_fs', 'y_fs', 'z_fs', 'wtref_fs',
            'c_star', 'net_concentration', 'full_scale_concentration'
        ]
        # Rename columns to match desired format
        column_mapping = {
            'filename': 'Filename',
            'x_fs': 'X_fs [m]',
            'y_fs': 'Y_fs [m]',
            'z_fs': 'Z_fs [m]',
            'wtref_fs': 'Wtref_fs [m/s]',
            'c_star': 'Avg_c_star [-]',
            'net_concentration': 'Avg_net_concentration [ppmV]',
            'full_scale_concentration': 'Avg_full_scale_concentration [ppmV]'
        }
    else:  # stats
        columns = [
            'filename', 'x_fs', 'y_fs', 'z_fs', 'wtref_fs',
            'c_star', 'net_concentration', 'full_scale_concentration',
            'percentiles95_c_star', 'percentiles95_net_concentration', 'percentiles95_full_scale_concentration',
            'peak2mean_c_star', 'peak2mean_net_concentration', 'peak2mean_full_scale_concentration'
        ]
        # Rename columns to match desired format
        column_mapping = {
            'filename': 'Filename',
            'x_fs': 'X_fs [m]',
            'y_fs': 'Y_fs [m]',
            'z_fs': 'Z_fs [m]',
            'wtref_fs': 'Wtref_fs [m/s]',
            'c_star': 'Avg_c_star [-]',
            'net_concentration': 'Avg_net_concentration [ppmV]',
            'full_scale_concentration': 'Avg_full_scale_concentration [ppmV]',
            'percentiles95_c_star': 'Percentiles 95_cstar',
            'percentiles95_net_concentration': 'Percentiles 95_net_concentration',
            'percentiles95_full_scale_concentration': 'percentiles95_full_scale_concentration',
            'peak2mean_c_star': 'Peak2MeanRatio_cstar',
            'peak2mean_net_concentration': 'Peak2MeanRatio_net_conc',
            'peak2mean_full_scale_concentration': 'Peak2MeanRatio_full_scale_conc'
        }
    
    # Select and rename columns
    df = df[columns].rename(columns=column_mapping)
    
    # Fill NaN values with 0.0
    df = df.fillna(0.0)
    
    # Save to CSV
    df.to_csv(output_filename, index=False)
    print(f"Data saved to {output_filename}")
    
    return df


def calculate_uncertainties(conc_ts_dict, split_factor=1.8, uncertainty_threshold=1e-4, 
                           verbose=True, metrics_to_calculate=None, concentration_types=None,
                           include_abs=True, include_pct=True):
    """
    Calculate measurement uncertainties from repeated measurements.
    
    Args:
        conc_ts_dict (dict): Dictionary mapping filenames to PointConcentration objects
        split_factor (float): Factor for splitting long time series (default: 1.8)
        uncertainty_threshold (float): Threshold for relative deviation calculation (default: 1e-4)
        verbose (bool): Print detailed logging information (default: True)
        metrics_to_calculate (list): List of metrics to calculate uncertainties for. 
                                     Options: ["Mean", "Median", "Peak2Mean", "P95"]
                                     If None, calculates for all (default: None)
        concentration_types (list): List of concentration types to calculate uncertainties for.
                                   Options: ["c_star", "net_concentration", "full_scale_concentration"]
                                   If None, calculates only for c_star (default: None)
        include_abs (bool): Include absolute uncertainties (default: True)
        include_pct (bool): Include percentage uncertainties (default: True)
        
    Returns:
        dict: Dictionary mapping config_name to uncertainty values with statistics column names
    """
    conc_type_map = {
        "c_star": {
            "attr": "c_star",
            "Mean": "Avg_c_star [-]",
            "Median": "Median_c_star",
            "Peak2Mean": "Peak2MeanRatio_cstar",
            "P95": "Percentiles 95_cstar"
        },
        "net_concentration": {
            "attr": "net_concentration",
            "Mean": "Avg_net_concentration [ppmV]",
            "Median": "Median_net_concentration",
            "Peak2Mean": "Peak2MeanRatio_net_conc",
            "P95": "Percentiles 95_net_concentration"
        },
        "full_scale_concentration": {
            "attr": "full_scale_concentration",
            "Mean": "Avg_full_scale_concentration [ppmV]",
            "Median": "Median_full_scale_concentration",
            "Peak2Mean": "Peak2MeanRatio_full_scale_conc",
            "P95": "percentiles95_full_scale_concentration"
        }
    }
    
    if concentration_types is None:
        concentration_types = ["c_star"]
    
    if metrics_to_calculate is None:
        metrics_to_calculate = ["Mean", "Median", "Peak2Mean", "P95"]
    
    def compute_stats(series, conc_type):
        stats = {}
        type_map = conc_type_map[conc_type]
        
        if "Mean" in metrics_to_calculate:
            stats[type_map["Mean"]] = float(np.mean(series))
        if "Median" in metrics_to_calculate:
            stats[type_map["Median"]] = float(np.median(series))
        if "Peak2Mean" in metrics_to_calculate:
            mean = float(np.mean(series))
            stats[type_map["Peak2Mean"]] = float(np.max(series) / mean) if mean != 0 else np.nan
        if "P95" in metrics_to_calculate:
            stats[type_map["P95"]] = float(np.percentile(series, 95))
        return stats
    uncertainty_results = {}
    
    for conc_type in concentration_types:
        config_stats = {}
        split_report = {}
        file_report = {}
        
        if conc_type not in conc_type_map:
            if verbose:
                print(f"⚠️ Unknown concentration type: {conc_type}, skipping")
            continue
        
        attr_name = conc_type_map[conc_type]["attr"]
        
        config_file_lengths = {}
        for file, ts in conc_ts_dict.items():
            if not hasattr(ts, attr_name):
                continue
            series_data = getattr(ts, attr_name)
            if series_data is None:
                continue
            cfg = ts.config_name
            if cfg not in config_file_lengths:
                config_file_lengths[cfg] = []
            config_file_lengths[cfg].append((file, len(series_data)))
        
        for cfg in config_file_lengths:
            split_report[cfg] = []
            file_report[cfg] = {"n_files": 0, "n_samples": 0}
        
        all_metrics = [conc_type_map[conc_type][m] for m in metrics_to_calculate if m in conc_type_map[conc_type]]
        
        for file, ts in conc_ts_dict.items():
            if not hasattr(ts, attr_name):
                continue
            series_data = getattr(ts, attr_name)
            if series_data is None:
                continue
            
            cfg = ts.config_name
            
            if cfg not in config_stats:
                config_stats[cfg] = {m: [] for m in all_metrics}
            
            series = series_data.values if hasattr(series_data, 'values') else series_data
            min_len = min([len(getattr(conc_ts_dict[f], attr_name)) for f in conc_ts_dict 
                          if hasattr(conc_ts_dict[f], attr_name) and getattr(conc_ts_dict[f], attr_name) is not None
                          and conc_ts_dict[f].config_name == cfg])
            
            to_split = len(series) >= split_factor * min_len and len(series) >= 2
            
            if to_split:
                half = len(series) // 2
                for part in [series[:half], series[half:]]:
                    stats = compute_stats(part, conc_type)
                    for m in all_metrics:
                        if m in stats:
                            config_stats[cfg][m].append(stats[m])
                file_report[cfg]["n_samples"] += 2
                split_report[cfg].append((file, len(series), half, len(series) - half))
                if verbose:
                    print(f"🔀 {cfg} ({conc_type}): File {file} (len={len(series)}) split -> halves {half} & {len(series) - half}")
            else:
                stats = compute_stats(series, conc_type)
                for m in all_metrics:
                    if m in stats:
                        config_stats[cfg][m].append(stats[m])
                file_report[cfg]["n_samples"] += 1
                if verbose:
                    print(f"   {cfg} ({conc_type}): File {file} (len={len(series)}) normal processing")
            
            file_report[cfg]["n_files"] += 1
        
        for cfg in config_stats:
            if verbose:
                n_files = file_report[cfg]["n_files"]
                n_samples = file_report[cfg]["n_samples"]
                splits = len(split_report[cfg])
                print(f"✅ Configuration {cfg} ({conc_type}): {n_files} files -> {n_samples} samples (splits: {splits})")
            
            if cfg not in uncertainty_results:
                uncertainty_results[cfg] = {}
            
            for m in all_metrics:
                if m not in config_stats[cfg]:
                    continue
                vals = np.array(config_stats[cfg][m])
                if vals.size >= 2:
                    overall = float(np.mean(vals))
                    deviation_abs = (np.max(vals) - np.min(vals)) / 2.0
                    deviation_pct = (deviation_abs / overall * 100.0) if overall >= uncertainty_threshold else np.nan
                    
                    if include_abs:
                        uncertainty_results[cfg][f"{m}_uncertainty_abs"] = deviation_abs
                    if include_pct:
                        uncertainty_results[cfg][f"{m}_uncertainty_pct"] = deviation_pct
                else:
                    if include_abs:
                        uncertainty_results[cfg][f"{m}_uncertainty_abs"] = np.nan
                    if include_pct:
                        uncertainty_results[cfg][f"{m}_uncertainty_pct"] = np.nan
    
    return uncertainty_results


def combine_to_csv2(file_names, base_path, file_type='avg', output_filename='combined_data.csv',
                    conc_ts_dict=None, config_name_mapping=None, uncertainty_results=None,
                    save_config_names=False, save_uncertainties=False, columns_to_save=None):
    """
    Extended version: Combine files and add config names and uncertainties.
    
    Args:
        file_names (list): Array of specific filenames to process
        base_path (str): Base directory path containing the files
        file_type (str): Type of files to process ('avg' or 'stats')
        output_filename (str): Name of the output CSV file
        conc_ts_dict (dict): Dictionary mapping filenames to PointConcentration objects
        config_name_mapping (dict): Dict mapping filename to config_name (optional if conc_ts_dict provided)
        uncertainty_results (dict): Dict mapping config_name to uncertainty values
        save_config_names (bool): Add ConfigName column
        save_uncertainties (bool): Add uncertainty columns
        columns_to_save (list): List of column names to keep (None = all)
        
    Returns:
        pd.DataFrame: Combined data as a pandas DataFrame
    """
    
    if file_type.lower() == 'avg':
        load_function = load_avg_file
    elif file_type.lower() == 'stats':
        load_function = load_stats_file
    else:
        raise ValueError("file_type must be 'avg' or 'stats'")
    
    file_paths = []
    for filename in file_names:
        search_pattern = os.path.join(base_path, '**', filename)
        matches = glob.glob(search_pattern, recursive=True)
        if matches:
            file_paths.extend(matches)
        else:
            print(f"Warning: File '{filename}' not found in {base_path}{filename}")
    
    if not file_paths:
        print(f"No valid {file_type} files found from the provided list")
        return None
    
    print(f"Processing {len(file_paths)} {file_type} files")
    
    all_data = []
    for filepath in file_paths:
        print(f"Processing: {filepath}")
        data = load_function(filepath)
        if data:
            all_data.append(data)
    
    if not all_data:
        print("No valid data found in any files")
        return None
    
    df = pd.DataFrame(all_data)
    
    if file_type.lower() == 'avg':
        columns = ['filename', 'x_fs', 'y_fs', 'z_fs', 'wtref_fs',
                   'c_star', 'net_concentration', 'full_scale_concentration']
        column_mapping = {
            'filename': 'Filename',
            'x_fs': 'X_fs [m]',
            'y_fs': 'Y_fs [m]',
            'z_fs': 'Z_fs [m]',
            'wtref_fs': 'Wtref_fs [m/s]',
            'c_star': 'Avg_c_star [-]',
            'net_concentration': 'Avg_net_concentration [ppmV]',
            'full_scale_concentration': 'Avg_full_scale_concentration [ppmV]'
        }
    else:
        columns = ['filename', 'x_fs', 'y_fs', 'z_fs', 'wtref_fs',
                   'c_star', 'net_concentration', 'full_scale_concentration',
                   'percentiles95_c_star', 'percentiles95_net_concentration', 'percentiles95_full_scale_concentration',
                   'peak2mean_c_star', 'peak2mean_net_concentration', 'peak2mean_full_scale_concentration']
        column_mapping = {
            'filename': 'Filename',
            'x_fs': 'X_fs [m]',
            'y_fs': 'Y_fs [m]',
            'z_fs': 'Z_fs [m]',
            'wtref_fs': 'Wtref_fs [m/s]',
            'c_star': 'Avg_c_star [-]',
            'net_concentration': 'Avg_net_concentration [ppmV]',
            'full_scale_concentration': 'Avg_full_scale_concentration [ppmV]',
            'percentiles95_c_star': 'Percentiles 95_cstar',
            'percentiles95_net_concentration': 'Percentiles 95_net_concentration',
            'percentiles95_full_scale_concentration': 'percentiles95_full_scale_concentration',
            'peak2mean_c_star': 'Peak2MeanRatio_cstar',
            'peak2mean_net_concentration': 'Peak2MeanRatio_net_conc',
            'peak2mean_full_scale_concentration': 'Peak2MeanRatio_full_scale_conc'
        }
    
    df = df[columns].rename(columns=column_mapping)
    df = df.fillna(0.0)
    
    if save_config_names:
        if config_name_mapping:
            mapping = config_name_mapping
        elif conc_ts_dict:
            mapping = {k: v.config_name for k, v in conc_ts_dict.items()}
        else:
            mapping = {}
        df.insert(1, "ConfigName", df["Filename"].apply(
            lambda x: mapping.get(x.replace("_stats_", "").replace("_avg_", ""), "")))
    
    if save_uncertainties and uncertainty_results:
        for idx, row in df.iterrows():
            fname = row["Filename"].replace("_stats_", "").replace("_avg_", "")
            cfg = None
            if config_name_mapping and fname in config_name_mapping:
                cfg = config_name_mapping[fname]
            elif conc_ts_dict and fname in conc_ts_dict:
                cfg = conc_ts_dict[fname].config_name
            
            if cfg and cfg in uncertainty_results:
                for key, val in uncertainty_results[cfg].items():
                    if key not in df.columns:
                        df[key] = np.nan
                    df.at[idx, key] = val
        
        uncertainty_map = {
            "Avg_c_star [-]": ["Avg_c_star [-]_uncertainty_abs", "Avg_c_star [-]_uncertainty_pct"],
            "Median_c_star": ["Median_c_star_uncertainty_abs", "Median_c_star_uncertainty_pct"],
            "Peak2MeanRatio_cstar": ["Peak2MeanRatio_cstar_uncertainty_abs", "Peak2MeanRatio_cstar_uncertainty_pct"],
            "Percentiles 95_cstar": ["Percentiles 95_cstar_uncertainty_abs", "Percentiles 95_cstar_uncertainty_pct"],
            
            "Avg_net_concentration [ppmV]": ["Avg_net_concentration [ppmV]_uncertainty_abs", "Avg_net_concentration [ppmV]_uncertainty_pct"],
            "Median_net_concentration": ["Median_net_concentration_uncertainty_abs", "Median_net_concentration_uncertainty_pct"],
            "Peak2MeanRatio_net_conc": ["Peak2MeanRatio_net_conc_uncertainty_abs", "Peak2MeanRatio_net_conc_uncertainty_pct"],
            "Percentiles 95_net_concentration": ["Percentiles 95_net_concentration_uncertainty_abs", "Percentiles 95_net_concentration_uncertainty_pct"],
            
            "Avg_full_scale_concentration [ppmV]": ["Avg_full_scale_concentration [ppmV]_uncertainty_abs", "Avg_full_scale_concentration [ppmV]_uncertainty_pct"],
            "Median_full_scale_concentration": ["Median_full_scale_concentration_uncertainty_abs", "Median_full_scale_concentration_uncertainty_pct"],
            "Peak2MeanRatio_full_scale_conc": ["Peak2MeanRatio_full_scale_conc_uncertainty_abs", "Peak2MeanRatio_full_scale_conc_uncertainty_pct"],
            "percentiles95_full_scale_concentration": ["percentiles95_full_scale_concentration_uncertainty_abs", "percentiles95_full_scale_concentration_uncertainty_pct"]
        }
        
        ordered_cols = []
        added_uncertainties = set()
        for col in df.columns:
            if "_uncertainty_" not in col:
                ordered_cols.append(col)
                if col in uncertainty_map:
                    for unc_col in uncertainty_map[col]:
                        if unc_col in df.columns and unc_col not in added_uncertainties:
                            ordered_cols.append(unc_col)
                            added_uncertainties.add(unc_col)
        
        for col in df.columns:
            if "_uncertainty_" in col and col not in added_uncertainties:
                ordered_cols.append(col)
        
        df = df[ordered_cols]
    
    if columns_to_save is not None:
        available_cols = [c for c in columns_to_save if c in df.columns]
        df = df[available_cols]
    
    df.to_csv(output_filename, index=False)
    print(f"Data saved to {output_filename}")
    
    return df


from collections import defaultdict
def load_csv(file_path, columns=["X_fs [m]", "Y_fs [m]", "Z_fs [m]","Avg_net_concentration [ppmV]"], combineCoordinates=True):
    """
    Load data from CSV file and return dict of arrays with column names.
    
    Args:
        file_path (str): Path to the CSV file
        columns (list): Column names to extract (None = all columns)
        combineCoordinates (bool): Combine x,y,z columns into list of tuples
    
    Returns:
        dict: Column names as keys, data arrays as values
    """
    data = defaultdict(list)
    coord_keys = ["X_fs [m]", "Y_fs [m]", "Z_fs [m]"]
    
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file)
        headers = reader.fieldnames
        
        if columns is None:
            columns = headers
        
        for row in reader:
            for col in columns:
                if col in row and row[col]:
                    try:
                        data[col].append(float(row[col]))
                    except ValueError:
                        data[col].append(row[col])
    
    if combineCoordinates and all(key in data for key in coord_keys):
        coords = list(zip(data[coord_keys[0]], data[coord_keys[1]], data[coord_keys[2]]))
        data['coordinates'] = coords
        for key in coord_keys:
            del data[key]
    
    return data


#Functions for plotting locations of measurements in a CAD loaded background  

def stl_to_2d_plot(stl_file, projection='xy', alpha=0.7, color='lightgrey',
                 edgecolor='grey', figsize=(10, 8), toFullScale="False", scaling=1, 
                 toLatLon="False", ref_coords=None, ref_pos=None):
    """
    Convert an STL file to a 2D shaded plot, with option to plot in lat/lon coordinates.
    
    Parameters:
    -----------
    stl_file : str
        Path to the STL file
    projection : str
        Which projection to use: 'xy', 'xz', 'yz', or 'diagonal'
    alpha : float
        Transparency of the shaded regions
    color : str
        Color for the shaded regions
    edgecolor : str
        Color for the edges
    figsize : tuple
        Figure size for the plot
    toFullScale : str
        "True" to scale the mesh by scaling factor / 1000
    scaling : float
        Scaling factor when toFullScale is "True"
    toLatLon : str
        "True" to convert coordinates to lat/lon
    ref_coords : list
        [lat, lon] reference coordinates corresponding to ref_pos position in the mesh
    ref_pos : list
        [x, y] position in the mesh corresponding to ref_coords. Defaults to [0, 0] if not specified.
        
    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from stl import mesh
   
    print(f"Loading STL file: {stl_file}")
    your_mesh = mesh.Mesh.from_file(stl_file)
    
    # Scale mesh to meters if FS
    if toFullScale == "True" and scaling is not None:
        your_mesh.vectors = your_mesh.vectors * scaling / 1000  # [mm] to [m]
    
    # Get min/max values for the mesh (in meters)
    min_coords = your_mesh.min_ * scaling / 1000 if toFullScale != "True" else your_mesh.min_
    max_coords = your_mesh.max_ * scaling / 1000 if toFullScale != "True" else your_mesh.max_
    
    # Convert to lat/lon
    lat_lon_bounds = None
    if toLatLon == "True" and ref_coords is not None:
        # Earth radius in meters
        r_earth = 6371000  # m
        
        #Origin if not provided
        if ref_pos is None:
            ref_pos = [0, 0]     
      
        ref_pos_m = [ref_pos[0] * scaling /1000, ref_pos[1] * scaling /1000] if toFullScale != "True" else ref_pos
            
        x_range = [min_coords[0], max_coords[0]]
        y_range = [min_coords[1], max_coords[1]]
        
        # Calculate lat/lon bounds considering reference position offset
        y_min_offset = y_range[0] - ref_pos_m[1]
        y_max_offset = y_range[1] - ref_pos_m[1]
        
        lat_range = [
            ref_coords[0] + (y_min_offset / r_earth) * (180 / np.pi),  # min lat
            ref_coords[0] + (y_max_offset / r_earth) * (180 / np.pi)   # max lat
        ]
        cos_lat = np.cos(ref_coords[0] * np.pi / 180)
        x_min_offset = x_range[0] - ref_pos_m[0]
        x_max_offset = x_range[1] - ref_pos_m[0]
        
        lon_range = [
            ref_coords[1] + (x_min_offset / (r_earth * cos_lat)) * (180 / np.pi),  # min lon
            ref_coords[1] + (x_max_offset / (r_earth * cos_lat)) * (180 / np.pi)   # max lon
        ]
        
        lat_lon_bounds = {
            'lat_range': lat_range,
            'lon_range': lon_range
        }
        print(f"Lat range: {lat_range}")
        print(f"Lon range: {lon_range}")
        print(f"Reference position [x,y]: {ref_pos}")
        print(f"Reference coordinates [lat,lon]: {ref_coords}")
    
   
    fig, ax = plt.subplots(figsize=figsize)
    
    #Model vertices
    for i, face in enumerate(your_mesh.vectors):
        
        if projection == 'xy':
            # XY projection (top view)
            if toLatLon == "True" and ref_coords is not None:
               
                points = []
                for vertex in face:
                    
                    vert_x = vertex[0]    #/1000?
                    vert_y = vertex[1]   
                   
                    x_offset = vert_x - ref_pos_m[0]
                    y_offset = vert_y - ref_pos_m[1]
                    
                    lat = ref_coords[0] + (y_offset / r_earth) * (180 / np.pi)
                    lon = ref_coords[1] + (x_offset / (r_earth * cos_lat)) * (180 / np.pi)
                    
                    points.append((lon, lat))  # Longitude is x, Latitude is y
            else:
                points = [(vertex[0], vertex[1]) for vertex in face]
        #XZ projection (front view)
        elif projection == 'xz':
            points = [(vertex[0], vertex[2]) for vertex in face]
        # YZ projection (side view)
        elif projection == 'yz':
            points = [(vertex[1], vertex[2]) for vertex in face]
        # Simple diagonal view (custom angle)
        elif projection == 'diagonal':
            angle = np.pi/6
            points = [(vertex[1]*np.cos(angle) + vertex[0]*np.sin(angle),
                     vertex[2] - 0.3*vertex[0]) for vertex in face]
        else:
            raise ValueError("Projections: 'xy', 'xz', 'yz', 'diagonal'")
            
        poly = Polygon(points, closed=True, alpha=alpha,
                      facecolor=color, edgecolor=edgecolor)
        ax.add_patch(poly)
    
    if projection == 'xy':
        if toLatLon == "True" and ref_coords is not None:
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            title = 'Geographic View (Lat/Lon Projection)'
            
            if lat_lon_bounds:
                print("test")
                #ax.set_xlim(lat_lon_bounds['lon_range'])
                #ax.set_ylim(lat_lon_bounds['lat_range'])
        else:
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            title = 'Top View (XY Projection)'
    elif projection == 'xz':
        ax.set_xlabel('X')
        ax.set_ylabel('Z')
        title = 'Front View (XZ Projection)'
    elif projection == 'yz':
        ax.set_xlabel('Y')
        ax.set_ylabel('Z')
        title = 'Side View (YZ Projection)'
    elif projection == 'diagonal':
        ax.set_xlabel('Custom X')
        ax.set_ylabel('Custom Z')
        title = 'Diagonal View'
        
    ax.set_title(title)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.set_aspect('equal')
    ax.autoscale()
    
    return fig, ax

def meters2latlon(points, ref_coords, scaling=1.0, toFullScale="False", ref_pos=None):
    """
    Convert 2D points from meters to latitude/longitude coordinates.
    
    Parameters:
    -----------
    points : array-like
        Array of [x, y] coordinates in meters to convert to lat/lon
        Can be a single point [x, y] or a list of points [[x1, y1], [x2, y2], ...]
    ref_coords : list
        [lat, lon] reference coordinates corresponding to ref_pos position
    scaling : float
        Scaling factor to convert input units to meters (default: 1.0)
    toFullScale : str
        "True" to apply the scaling factor to convert units to meters
    ref_pos : list
        [x, y] position in the input coordinates corresponding to ref_coords
        Defaults to [0, 0] if not specified
        
    Returns:
    --------
    points_latlon : array-like
        Array of [lat, lon] coordinates corresponding to the input points
        Format matches the input format (single point or list of points)
    """
    import numpy as np
    
    # Check if points is a single point or a list of points
    is_single_point = not isinstance(points[0], (list, tuple, np.ndarray))
    # Convert to numpy array for easier processing
    points_array = np.array([points]) if is_single_point else np.array(points)
  
    r_earth = 6371000# m
    
    #Default ref position to origin
    if ref_pos is None:
        ref_pos = [0, 0]
    
    scale_factor = scaling/1000 if toFullScale == "True" else 1.0
    ref_pos_m = np.array(ref_pos) * scale_factor
    
 
    points_latlon = []
    
   
    cos_lat = np.cos(ref_coords[0] * np.pi / 180)
    for i, point in enumerate(points_array):
        # Scale point coordinates to meters if needed
        point_m = point * scale_factor
        
        x_offset = point_m[0] - ref_pos_m[0]
        y_offset = point_m[1] - ref_pos_m[1]
        
        lat = ref_coords[0] + (y_offset / r_earth) * (180 / np.pi)
        lon = ref_coords[1] + (x_offset / (r_earth * cos_lat)) * (180 / np.pi)
        
        points_latlon.append((lat, lon))
    
    if is_single_point:
        return points_latlon[0]
    else:
        return points_latlon
    

def add_crosses(ax, points, values=None, thresholds=None, colors=None, 
                marker='x', size=100, linewidth=2, cmap=None, add_colorbar=True,
                is_latlon=False):
    """
    Add color-coded crosses to the plot based on threshold values,
    with optional handling for geographic lat/lon data.
    
    Parameters:
    -----------
    ax : matplotlib axis
        The axis to add crosses to
    points : list or array
        List of (x, y) coordinates or (lat, lon) coordinates for crosses
    values : list or array, optional
        Values corresponding to each point (for coloring)
    thresholds : list, optional
        Threshold values for color changes
    colors : list, optional
        Colors corresponding to threshold ranges
    marker : str, optional
        Marker style ('x', '+', etc.)
    size : int, optional
        Size of the markers
    linewidth : int, optional
        Width of the marker lines
    cmap : str or colormap, optional
        Custom colormap (if not using thresholds/colors)
    add_colorbar : bool, optional
        Whether to add a colorbar to the plot
    is_latlon : bool, optional
        Whether the points are in (lat, lon) format, which requires special handling
        
    Returns:
    --------
    scatter : matplotlib scatter plot
        The scatter plot object for the crosses
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D
    
    if points is None or len(points) == 0:
        print("No points provided for crosses")
        return None
    #Lat/Lon conversion
    points_array = np.array(points)
    if is_latlon:
        x_coords = points_array[:, 1]  # Longitude as x
        y_coords = points_array[:, 0]  # Latitude as y
    else:
        # Standard x/y coordinates
        x_coords, y_coords = zip(*points)
    
    
    if values is not None:
        if thresholds is not None and colors is not None:
            # Use threshold-based coloring
            point_colors = []
            for value in values:
                for i, threshold in enumerate(thresholds):
                    if i == 0 and value < threshold:
                        point_colors.append(colors[0])
                        break
                    elif i == len(thresholds) - 1 or value < threshold:
                        point_colors.append(colors[i])
                        break
                    
            scatter = ax.scatter(x_coords, y_coords, c=point_colors, marker=marker, 
                                s=size, linewidth=linewidth)
            
            if add_colorbar:
                legend_elements = []
                
                # Add first range
                legend_elements.append(
                    Line2D([0], [0], marker=marker, color=colors[0], markerfacecolor=colors[0],
                          markersize=10, label=f'< {thresholds[0]}')
                )
                
                # Add middle ranges
                for i in range(len(thresholds) - 1):
                    legend_elements.append(
                        Line2D([0], [0], marker=marker, color=colors[i+1], markerfacecolor=colors[i+1],
                              markersize=10, label=f'{thresholds[i]} - {thresholds[i+1]}')
                    )
                
                # Add last range
                legend_elements.append(
                    Line2D([0], [0], marker=marker, color=colors[-1], markerfacecolor=colors[-1],
                          markersize=10, label=f'> {thresholds[-1]}')
                )
                
                ax.legend(handles=legend_elements, title="Value Ranges", 
                         loc='best', framealpha=0.7)
        else:
            # Use continuous coloring with a colormap
            if cmap is None:
                cmap = 'viridis'
                
            norm = Normalize(vmin=min(values), vmax=max(values))
            scatter = ax.scatter(x_coords, y_coords, c=values, marker=marker, 
                                s=size, linewidth=linewidth, cmap=cmap, norm=norm)
            
            if add_colorbar:
                plt.colorbar(scatter, ax=ax, label='Value')
    else:
        # Just add crosses without color coding
        scatter = ax.scatter(x_coords, y_coords, marker=marker, s=size, 
                            linewidth=linewidth, color='red')
    
    # For geographic plots, equal aspect ratio can distort the map at non-equatorial latitudes
    if is_latlon:
        # For lat/lon data, we typically don't want an equal aspect ratio
        # because 1 degree lat doesn't equal 1 degree lon except at the equator
        ax.set_aspect('auto')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
    else:
        # For Cartesian data, equal aspect ratio makes sense
        ax.set_aspect('equal')
        
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.autoscale()
    
    return scatter


def add_crosses_geo(ax, points, values=None, thresholds=None, colors=None, 
                marker='x', size=100, linewidth=2, cmap=None, add_colorbar=True,
                is_latlon=False):
    """
    Add color-coded crosses to the plot based on threshold values,
    with special handling for geographic lat/lon data.
    
    Parameters:
    -----------
    ax : matplotlib axis
        The axis to add crosses to
    points : list or array
        List of (x, y) coordinates or (lat, lon) coordinates for crosses
    values : list or array, optional
        Values corresponding to each point (for coloring)
    thresholds : list, optional
        Threshold values for color changes
    colors : list, optional
        Colors corresponding to threshold ranges
    marker : str, optional
        Marker style ('x', '+', etc.)
    size : int, optional
        Size of the markers
    linewidth : int, optional
        Width of the marker lines
    cmap : str or colormap, optional
        Custom colormap (if not using thresholds/colors)
    add_colorbar : bool, optional
        Whether to add a colorbar to the plot
    is_latlon : bool, optional
        Whether the points are in (lat, lon) format, which requires special handling
        
    Returns:
    --------
    scatter : matplotlib scatter plot
        The scatter plot object for the crosses
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D
    
    if points is None or len(points) == 0:
        print("No points provided for crosses")
        return None
    
    # Handle points format based on whether it's lat/lon
    points_array = np.array(points)
    if is_latlon:
        # For lat/lon data, we need to make sure the order is correct
        # Most GIS uses (lat, lon) format but matplotlib expects (x=lon, y=lat)
        x_coords = points_array[:, 1]  # Longitude as x
        y_coords = points_array[:, 0]  # Latitude as y
    else:
        # Standard x/y coordinates
        x_coords, y_coords = zip(*points)
    
    # Different coloring approaches based on inputs
    if values is not None:
        if thresholds is not None and colors is not None:
            # Use threshold-based coloring
            point_colors = []
            for value in values:
                color_assigned = False
                for i, threshold in enumerate(thresholds):
                    if i == 0 and value < threshold:
                        point_colors.append(colors[0])
                        color_assigned = True
                        break
                    elif i < len(thresholds) - 1 and value >= thresholds[i] and value < thresholds[i+1]:
                        point_colors.append(colors[i+1])
                        color_assigned = True
                        break
                
                # Handle the last range (above the highest threshold)
                if not color_assigned:
                    point_colors.append(colors[-1])
                    
            scatter = ax.scatter(x_coords, y_coords, c=values, marker=marker,  #point_colors
                                s=size, linewidth=linewidth)
            
            if add_colorbar:
                legend_elements = []
                
                # Add first range
                legend_elements.append(
                    Line2D([0], [0], marker=marker, color=colors[0], markerfacecolor=colors[0],
                          markersize=10, label=f'< {thresholds[0]}')
                )
                
                # Add middle ranges
                for i in range(len(thresholds) - 1):
                    legend_elements.append(
                        Line2D([0], [0], marker=marker, color=colors[i+1], markerfacecolor=colors[i+1],
                              markersize=10, label=f'{thresholds[i]} - {thresholds[i+1]}')
                    )
                
                # Add last range
                legend_elements.append(
                    Line2D([0], [0], marker=marker, color=colors[-1], markerfacecolor=colors[-1],
                          markersize=10, label=f'> {thresholds[-1]}')
                )
                
                ax.legend(handles=legend_elements, title="Value Ranges", 
                         loc='best', framealpha=0.7)
        else:
            # Use continuous coloring with a colormap
            if cmap is None:
                cmap = 'viridis'
                
            norm = Normalize(vmin=min(values), vmax=max(values))
            scatter = ax.scatter(x_coords, y_coords, c=values, marker=marker, 
                                s=size, linewidth=linewidth, cmap=cmap, norm=norm)
            
            if add_colorbar:
                plt.colorbar(scatter, ax=ax, label='Value')
    else:
        # Just add crosses without color coding
        scatter = ax.scatter(x_coords, y_coords, marker=marker, s=size, 
                            linewidth=linewidth, color='red')
    
    # For geographic plots, equal aspect ratio can distort the map at non-equatorial latitudes
    if is_latlon:
        # For lat/lon data, we typically don't want an equal aspect ratio
        # because 1 degree lat doesn't equal 1 degree lon except at the equator
        ax.set_aspect('auto')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
    else:
        # For Cartesian data, equal aspect ratio makes sense
        ax.set_aspect('equal')
        
    ax.grid(True, linestyle='--', alpha=0.7)
    
    return scatter

"""
def add_crosses(ax, points, values=None, thresholds=None, colors=None, 
                marker='x', size=100, linewidth=2, cmap=None, add_colorbar=True):
    
    Add color-coded crosses to the plot based on threshold values.
    
    Parameters:
    -----------
    ax : matplotlib axis
        The axis to add crosses to
    points : list or array
        List of (x, y) coordinates for crosses
    values : list or array, optional
        Values corresponding to each point (for coloring)
    thresholds : list, optional
        Threshold values for color changes
    colors : list, optional
        Colors corresponding to threshold ranges
    marker : str, optional
        Marker style ('x', '+', etc.)
    size : int, optional
        Size of the markers
    linewidth : int, optional
        Width of the marker lines
    cmap : str or colormap, optional
        Custom colormap (if not using thresholds/colors)
    add_colorbar : bool, optional
        Whether to add a colorbar to the plot
        
    Returns:
    --------
    scatter : matplotlib scatter plot
        The scatter plot object for the crosses
    
    if points is None or len(points) == 0:
        print("No points provided for crosses")
        return None
    
    x_coords, y_coords = zip(*points)
    
    # Different coloring approaches based on inputs
    if values is not None:
        if thresholds is not None and colors is not None:
            # Use threshold-based coloring
            point_colors = []
            for value in values:
                for i, threshold in enumerate(thresholds):
                    if i == 0 and value < threshold:
                        point_colors.append(colors[0])
                        break
                    elif i == len(thresholds) - 1 or value < threshold:
                        point_colors.append(colors[i])
                        break
                    
            scatter = ax.scatter(x_coords, y_coords, c=point_colors, marker=marker, 
                                s=size, linewidth=linewidth)
            
            # Add a custom legend for threshold colors
            if add_colorbar:
                from matplotlib.lines import Line2D
                legend_elements = []
                
                # Add first range
                legend_elements.append(
                    Line2D([0], [0], marker=marker, color=colors[0], markerfacecolor=colors[0],
                          markersize=10, label=f'< {thresholds[0]}')
                )
                
                # Add middle ranges
                for i in range(len(thresholds) - 1):
                    legend_elements.append(
                        Line2D([0], [0], marker=marker, color=colors[i+1], markerfacecolor=colors[i+1],
                              markersize=10, label=f'{thresholds[i]} - {thresholds[i+1]}')
                    )
                
                # Add last range
                legend_elements.append(
                    Line2D([0], [0], marker=marker, color=colors[-1], markerfacecolor=colors[-1],
                          markersize=10, label=f'> {thresholds[-1]}')
                )
                
                ax.legend(handles=legend_elements, title="Value Ranges", 
                         loc='best', framealpha=0.7)
        else:
            # Use continuous coloring with a colormap
            if cmap is None:
                cmap = 'viridis'
                
            norm = Normalize(vmin=min(values), vmax=max(values))
            scatter = ax.scatter(x_coords, y_coords, c=values, marker=marker, 
                                s=size, linewidth=linewidth, cmap=cmap, norm=norm)
            
            if add_colorbar:
                plt.colorbar(scatter, ax=ax, label='Value')
    else:
        # Just add crosses without color coding
        scatter = ax.scatter(x_coords, y_coords, marker=marker, s=size, 
                            linewidth=linewidth, color='red')
        
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.set_aspect('equal')
    ax.autoscale()
    
    return scatter
"""


def add_crosses_3d(ax, points, values=None, thresholds=None, colors=None,
                  marker='x', size=100, linewidth=2, cmap=None, add_colorbar=True):
    """
    Add color-coded crosses to a 3D plot based on threshold values.
    
    Parameters:
    -----------
    ax : matplotlib 3D axis
        The axis to add crosses to
    points : list or array
        List of (x, y, z) coordinates for crosses
    values : list or array, optional
        Values corresponding to each point (for coloring)
    thresholds : list, optional
        Threshold values for color changes
    colors : list, optional
        Colors corresponding to threshold ranges
    marker : str, optional
        Marker style ('x', '+', etc.)
    size : int, optional
        Size of the markers
    linewidth : int, optional
        Width of the marker lines
    cmap : str or colormap, optional
        Custom colormap (if not using thresholds/colors)
    add_colorbar : bool, optional
        Whether to add a colorbar to the plot
        
    Returns:
    --------
    scatter : matplotlib scatter plot
        The scatter plot object for the crosses
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    if points is None or len(points) == 0:
        print("No points provided for crosses")
        return None
    
    # For 3D points, we need x, y, and z coordinates
    if len(points[0]) == 3:
        x_coords, y_coords, z_coords = zip(*points)
    else:
        # For 2D points, set z to zeros
        x_coords, y_coords = zip(*points)
        z_coords = [0] * len(x_coords)
    
    # Different coloring approaches based on inputs
    if values is not None:
        if thresholds is not None and colors is not None:
            # Use threshold-based coloring
            point_colors = []
            for value in values:
                color_assigned = False
                for i, threshold in enumerate(thresholds):
                    if i == 0 and value < threshold:
                        point_colors.append(colors[0])
                        color_assigned = True
                        break
                    elif i < len(thresholds) - 1 and threshold <= value < thresholds[i+1]:
                        point_colors.append(colors[i+1])
                        color_assigned = True
                        break
                    elif i == len(thresholds) - 1 and value >= threshold:
                        point_colors.append(colors[-1])
                        color_assigned = True
                        break
                # If no threshold matched, use the last color
                if not color_assigned:
                    point_colors.append(colors[-1])
            
            scatter = ax.scatter(x_coords, y_coords, z_coords, c=point_colors, marker=marker,
                               s=size, linewidth=linewidth)
            
            if add_colorbar:
                from matplotlib.lines import Line2D
                legend_elements = []
                
                # Add first range
                legend_elements.append(
                    Line2D([0], [0], marker=marker, color=colors[0], markerfacecolor=colors[0],
                          markersize=10, label=f'< {thresholds[0]}')
                )
                
                # Add middle ranges
                for i in range(len(thresholds) - 1):
                    legend_elements.append(
                        Line2D([0], [0], marker=marker, color=colors[i+1], markerfacecolor=colors[i+1],
                              markersize=10, label=f'{thresholds[i]} - {thresholds[i+1]}')
                    )
                
                # Add last range
                legend_elements.append(
                    Line2D([0], [0], marker=marker, color=colors[-1], markerfacecolor=colors[-1],
                          markersize=10, label=f'> {thresholds[-1]}')
                )
                
                ax.legend(handles=legend_elements, title="Value Ranges",
                        loc='best', framealpha=0.7)
        else:
            # Use continuous coloring with a colormap
            if cmap is None:
                cmap = 'viridis'
            
            norm = Normalize(vmin=min(values), vmax=max(values))
            scatter = ax.scatter(x_coords, y_coords, z_coords, c=values, marker=marker,
                               s=size, linewidth=linewidth, cmap=cmap, norm=norm)
            
            if add_colorbar:
                plt.colorbar(scatter, ax=ax, label='Value')
    else:
        # Just add crosses without color coding
        scatter = ax.scatter(x_coords, y_coords, z_coords, marker=marker, s=size,
                           linewidth=linewidth, color='red')
    
    return scatter


def show_multiple_projections(stl_file, points_dict=None, values_dict=None, 
                             thresholds=None, colors=None):
    """
    Show all three standard projections of an STL file side by side,
    with optional color-coded crosses.
    
    Parameters:
    -----------
    stl_file : str
        Path to the STL file
    points_dict : dict, optional
        Dictionary with projection keys ('xy', 'xz', 'yz') and point coordinates
    values_dict : dict, optional
        Dictionary with projection keys and values for coloring points
    thresholds : list, optional
        Threshold values for color coding
    colors : list, optional
        Colors corresponding to thresholds
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from stl import mesh

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    your_mesh = mesh.Mesh.from_file(stl_file)
    
    projections = ['xy', 'xz', 'yz']
    titles = ['Top View (XY)', 'Front View (XZ)', 'Side View (YZ)']
    
    for i, (proj, title) in enumerate(zip(projections, titles)):
        ax = axes[i]
        
        #Model vertices
        for face in your_mesh.vectors:
            if proj == 'xy':
                points = [(vertex[0], vertex[1]) for vertex in face]
            elif proj == 'xz':
                points = [(vertex[0], vertex[2]) for vertex in face]
            else:  # yz
                points = [(vertex[1], vertex[2]) for vertex in face]
            
            poly = Polygon(points, closed=True, alpha=0.7, 
                          facecolor='lightblue', edgecolor='blue')
            ax.add_patch(poly)
        
        # Add crosses if provided for this projection
        if points_dict is not None and proj in points_dict:
            proj_points = points_dict[proj]
            proj_values = None if values_dict is None else values_dict.get(proj, None)
            
            add_crosses(ax, proj_points, values=proj_values, 
                       thresholds=thresholds, colors=colors)
        
        if proj == 'xy':
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
        elif proj == 'xz':
            ax.set_xlabel('X')
            ax.set_ylabel('Z')
        else:  # yz
            ax.set_xlabel('Y')
            ax.set_ylabel('Z')
            
        ax.set_title(title)
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.set_aspect('equal')
        ax.autoscale()
    
    plt.tight_layout()
    return fig, axes



def plot_stl_3d(stl_file, color_map=None, default_color='gray', 
                x_range=None, y_range=None, z_range=None, 
                clip_polygons=True,elev=30,azim=45,toFullScale="False",scaling=None):
    """
    Create a 3D visualization of an STL file with optional region filtering.
    
    Parameters:
    -----------
    stl_file : str
        Path to the STL file
    color_map : dict, optional
        Dictionary mapping face indices to colors
    default_color : str, optional
        Default color for faces not in the color map
    x_range, y_range, z_range : tuple, optional
        (min, max) ranges to filter polygons (None = no filtering)
    clip_polygons : bool, optional
        If True, clip polygons to the specified ranges
        If False, only show polygons whose centroid is in range
        
    Returns:
    --------
    fig, ax : matplotlib figure and axis
        The figure and 3D axis containing the plot
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from stl import mesh

    your_mesh = mesh.Mesh.from_file(stl_file)
   
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    #Transform to FS
    if toFullScale == "True" and scaling is not None:
        your_mesh.vectors = your_mesh.vectors * scaling / 1000  # [mm] to [m]
    
    # If no color map provided, use default color for all faces
    if color_map is None:
        color_map = {}

    
    mesh_min = your_mesh.points.reshape(-1, 3).min(axis=0)
    mesh_max = your_mesh.points.reshape(-1, 3).max(axis=0)
    #If range overgiven use this instead of mesh range
    x_min = x_range[0] if x_range is not None else mesh_min[0]
    x_max = x_range[1] if x_range is not None else mesh_max[0]
    y_min = y_range[0] if y_range is not None else mesh_min[1]
    y_max = y_range[1] if y_range is not None else mesh_max[1]
    z_min = z_range[0] if z_range is not None else mesh_min[2]
    z_max = z_range[1] if z_range is not None else mesh_max[2]
    
   
    polygons = []
    face_colors = []
    for i, face in enumerate(your_mesh.vectors):
        #Center of of the faces
        centroid = np.mean(face, axis=0)
        # Check if this face should be included based on ranges
        #if (x_range is not None and (centroid[0] < x_min or centroid[0] > x_max)):
        #    continue
        #if (y_range is not None and (centroid[1] < y_min or centroid[1] > y_max)):
        #    continue
        #if (z_range is not None and (centroid[2] < z_min or centroid[2] > z_max)):
        #    continue
        if ((centroid[0] < x_min or centroid[0] > x_max)):
            continue
        if ((centroid[1] < y_min or centroid[1] > y_max)):
            continue
        if ((centroid[2] < z_min or centroid[2] > z_max)):
            continue
        
        polygons.append(face)
        
        if i in color_map:
            face_colors.append(color_map[i])
        else:
            face_colors.append(default_color)
    
    poly = Poly3DCollection(polygons, alpha=0.7)
    poly.set_facecolor(face_colors)
    poly.set_edgecolor('black')
    ax.add_collection3d(poly)
    
    # Set view limits and angle
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)
    ax.view_init(elev=elev, azim=azim)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D STL Viewer')
    
    
    return fig, ax



def add_velocity_field(ax, points, values, components=None, 
                 interpolate=True, grid_density=30, transparencyFactor=0.5, cmap='viridis',
                 add_colorbar=True, show_vectors=False, vector_scale=20, 
                 vector_spacing=5, vector_color='black', vector_width=0.5, show_contour_lines=False,
                 is_latlon=False, zorder=5):
    """
    Add a 2D color-coded velocity field to the plot, with optional vector overlay.
    
    Parameters:
    -----------
    ax : matplotlib axis
        The axis to add the velocity field to
    points : list or array
        List of (x, y) or (lat, lon) coordinates for velocity data points
    values : list or array
        Speed values corresponding to each point (magnitude)
    components : list or array, optional
        List of [u, v] components for each point, required if show_vectors=True
    interpolate : bool, optional
        Whether to interpolate velocities between points (True) or just show at discrete points (False)
    grid_density : int, optional
        Number of grid points in each dimension for interpolation
    transparencyFactor : float, optional
        Transparency of the velocity field, low=high transparency, high=low transparency
    cmap : str or colormap, optional
        Colormap for the velocity field
    add_colorbar : bool, optional
        Whether to add a colorbar to the plot
    show_vectors : bool, optional
        Whether to overlay velocity vectors on the field
    vector_scale : float, optional
        Scaling factor for vector arrows
    vector_spacing : int, optional
        Spacing between vectors when using interpolation (higher = fewer vectors)
    vector_color : str, optional
        Color for the vector arrows
    vector_width : float, optional
        Line width for the vector arrows
    is_latlon : bool, optional
        Whether the points are in (lat, lon) format
    zorder : int, optional
        Z-order for drawing (lower number = drawn earlier/underneath)
        
    Returns:
    --------
    field : matplotlib contourf or scatter
        The field object for the velocities
    vectors : matplotlib quiver or None
        The vector field object if show_vectors=True, otherwise None
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    
    if points is None or len(points) == 0 or values is None or len(values) == 0:
        print("No valid points or values provided")
        return None, None

    points_array = np.array(points)
    values_array = np.array(values)
    
    # Handle coordinate format for plotting
    if is_latlon:
        x_coords = points_array[:, 1]  # Longitude as x
        y_coords = points_array[:, 0]  # Latitude as y
    else:
        # Standard x/y coordinates
        x_coords = points_array[:, 0]
        y_coords = points_array[:, 1]
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    
    field = None
    vectors = None
    if interpolate:
        # Create a regular grid for interpolation
        xi = np.linspace(xlim[0], xlim[1], grid_density)
        yi = np.linspace(ylim[0], ylim[1], grid_density)
        Xi, Yi = np.meshgrid(xi, yi)
        
        # Interpolate the speed values onto the grid
        Zi = griddata((x_coords, y_coords), values_array, (Xi, Yi), method='cubic', fill_value=np.nan)

        # Create a filled contour plot with the velocity field
        if show_contour_lines:
            field = ax.contourf(Xi, Yi, Zi, levels=50, cmap=cmap, alpha=transparencyFactor, zorder=zorder)
        else:
            # Remove visible isolines by setting linewidths=0 and antialiased=True
            field = ax.contourf(Xi, Yi, Zi, levels=50, cmap=cmap, alpha=transparencyFactor, zorder=zorder,
                               linewidths=0, antialiased=True)
        
        #Filled contour plot with the velocity field
       # field = ax.contourf(Xi, Yi, Zi, levels=50, cmap=cmap, alpha=alpha, zorder=zorder)
        
        if add_colorbar:
            cbar = plt.colorbar(field, ax=ax, label='Speed')
        
        #Add vector overlaying
        if show_vectors and components is not None:
            # Prepare vector components
            if is_latlon:
                u_components = np.array([comp[0] for comp in components])
                v_components = np.array([comp[1] for comp in components])
            else:
                u_components = np.array([comp[0] for comp in components])
                v_components = np.array([comp[1] for comp in components])
            
            # Interpolate the vector components onto sparser grid to avoid overcrowding
            step = vector_spacing
            Xi_sparse = Xi[::step, ::step]
            Yi_sparse = Yi[::step, ::step]
            # Interpolate u and v components separately
            Ui = griddata((x_coords, y_coords), u_components, (Xi_sparse, Yi_sparse), method='cubic', fill_value=0)
            Vi = griddata((x_coords, y_coords), v_components, (Xi_sparse, Yi_sparse), method='cubic', fill_value=0)
            
            # Create vector field
            vectors = ax.quiver(Xi_sparse, Yi_sparse, Ui, Vi, 
                              scale=vector_scale, color=vector_color, 
                              width=vector_width, zorder=zorder+1)
    else:
        # For discrete points, use a scatter plot
        field = ax.scatter(x_coords, y_coords, c=values_array, cmap=cmap, 
                         alpha=transparencyFactor, s=50, zorder=zorder)
        
    
        if add_colorbar:
            cbar = plt.colorbar(field, ax=ax, label='Speed')
        
        # Add vector overlay
        if show_vectors and components is not None:
            # Prepare vector components
            if is_latlon:
                u_components = np.array([comp[0] for comp in components])
                v_components = np.array([comp[1] for comp in components])
            else:
                u_components = np.array([comp[0] for comp in components])
                v_components = np.array([comp[1] for comp in components])
            
            # Create vector field at discrete points
            vectors = ax.quiver(x_coords, y_coords, u_components, v_components, 
                              scale=vector_scale, color=vector_color, 
                              width=vector_width, zorder=zorder+1)
    
    return field, vectors


if __name__ == "__main__":
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='Convert STL to 2D shaded plot with optional markers')
    parser.add_argument('stl_file', help='Path to the STL file')
    parser.add_argument('--projection', choices=['xy', 'xz', 'yz', 'all'], default='xy',
                        help='Projection plane (default: xy)')
    parser.add_argument('--output', help='Output file path (default: display only)')
    parser.add_argument('--points', help='CSV file with marker coordinates (x,y,value)')
    parser.add_argument('--thresholds', help='Comma-separated threshold values (e.g., 10,20,30)')
    parser.add_argument('--colors', help='Comma-separated colors (e.g., red,yellow,green,blue)')
    
    args = parser.parse_args()
    
    # Load points and values if provided
    points_dict = {}
    values_dict = {}
    
    if args.points:
        import csv
        with open(args.points, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)  # Skip header row
            
            xy_points = []
            xz_points = []
            yz_points = []
            xy_values = []
            xz_values = []
            yz_values = []
            
            for row in reader:
                if len(row) >= 3:  # At least x, y, z
                    x, y, z = float(row[0]), float(row[1]), float(row[2])
                    
                    # Add points to respective projections
                    xy_points.append((x, y))
                    xz_points.append((x, z))
                    yz_points.append((y, z))
                    
                    # Add values if provided
                    if len(row) >= 4:
                        value = float(row[3])
                        xy_values.append(value)
                        xz_values.append(value)
                        yz_values.append(value)
            
            points_dict = {'xy': xy_points, 'xz': xz_points, 'yz': yz_points}
            if xy_values:
                values_dict = {'xy': xy_values, 'xz': xz_values, 'yz': yz_values}
    
    # Parse thresholds and colors
    thresholds = None
    colors = None
    
    if args.thresholds:
        thresholds = [float(t) for t in args.thresholds.split(',')]
        
    if args.colors:
        colors = args.colors.split(',')
    
    # Generate the plot(s)
    if args.projection == 'all':
        fig, axes = show_multiple_projections(
            args.stl_file, 
            points_dict=points_dict, 
            values_dict=values_dict if values_dict else None,
            thresholds=thresholds,
            colors=colors
        )
    else:
        fig, ax = stl_to_2d_plot(args.stl_file, projection=args.projection)
        
        # Add crosses if points are provided
        if args.projection in points_dict:
            add_crosses(
                ax, 
                points_dict[args.projection], 
                values=values_dict.get(args.projection) if values_dict else None,
                thresholds=thresholds,
                colors=colors
            )
    
    # Save or display the plot
    if args.output:
        plt.savefig(args.output, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {args.output}")
    else:
        plt.show()
    """

    
    

def load_combined_data_from_csv(file_path):
    """Load combined data from CSV file"""
    filenames = []
    points = []
    avg_data = {}
    stats_data = {}
    
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file)
        
        for row in reader:
            try:
                filenames.append(row['Filename'])
                
                x = float(row['X_fs [m]'])
                y = float(row['Y_fs [m]'])
                z = float(row['Z_fs [m]'])
                points.append((x, y, z))
                
                for key, value in row.items():
                    if key.startswith('Avg_'):
                        param_name = key[4:]
                        if param_name not in avg_data:
                            avg_data[param_name] = []
                        avg_data[param_name].append(float(value))
                    elif key.startswith('Stats_'):
                        param_name = key[6:]
                        if param_name not in stats_data:
                            stats_data[param_name] = []
                        stats_data[param_name].append(float(value))
                        
            except (ValueError, KeyError) as e:
                print(f"Skipping invalid row: {e}")
                
    return filenames, points, avg_data, stats_data


def load_data_from_csv(file_path):
    """
    Load data from a CSV file returns points and values lists for printing/or further analysis.
    
    Args:
        filename (str): Path to the CSV file
        
    Returns:
        tuple: (points, values) where:
            - points is a list of (x,y,z) tuples
            - values is a list of concentration values
    """
    points = []
    values = []
    
    with open(file_path, 'r') as file:
        # Skip the header line
        next(file)
        csv_reader = csv.reader(file)
        for row in csv_reader:
            # Check if we have enough elements in the row
            if len(row) >= 4:
                try:
                    # Parse data, assuming the format matches X_fs [m],Y_fs [m],Z_fs [m],C* [-]
                    x = float(row[0])
                    y = float(row[1])
                    z = float(row[2])
                    c = float(row[3])
                    
                    # Add to our lists
                    points.append((x, y, z))
                    values.append(c)
                except ValueError:
                    # Skip rows that can't be converted to float
                    print(f"Skipping invalid row: {row}")
    
    return points, values

    # Example usage:
    # data_dict = load_concentration_data('concentration_data.txt')
    # print(data_dict)
    
    
    
