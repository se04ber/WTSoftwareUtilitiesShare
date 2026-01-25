#!/usr/bin/env python3
# Wind tunnel concentration analysis script
# Configuration: lines 8-95
# Add new plots: see ### ADD NEW PLOTS HERE ### in run_visualizations()

import os
import sys
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import importlib
import windtunnel as wt
from windtunnel.concentration.CompareDatasets import *

# Paths
path_dir = "/home/sabrina/Desktop/Schreibtisch/Arbeit_2025/WTSoftwareUtilitiesShare"
path = f"{path_dir}/ExampleData/InputData/Beispiel Umrechnung zur Kontrolle/"
namelist = ['UBA_GA_02_04_01_000_1_001']
output_path = f"{path_dir}/ExampleData/Results/"
csv_file = f"{path_dir}/ExampleData/ParameterFiles/ambient_conditions_.UBA_GA.csv"
parameters_PerFolder = False  # Read params per folder (True) or per file (False)
# Processing
full_scale = 'ms'  # Scale: 'ms' (model), 'fs' (full), 'nd' (non-dimensional)
applyPostprocessing = True
averageInterval = 60  # seconds
measurementFreq = 0.005  # Hz
averagingColumns = ["net_concentration"]
# Saving
osType = "Linux"  # "Linux" or "Windows"
outputName = None
saveTs = True  # Save time series
saveAvg = True  # Save averages
saveStats = True  # Save statistics
saveCombined = True  # Save combined file
combinedFileName = "combined_file_nora.csv"
base_path = None
saveAll = True  # Override individual save flags
# Uncertainty
calculateUncertainty = True
saveUncertainties = True
saveConfigNames = True
split_factor = 1.8
uncertainty_threshold = 1e-4
verboseUncertainty = True
uncertaintyMetrics = None  # None = all metrics
uncertaintyConcentrationTypes = None  # None = c_star only
includeAbsoluteUncertainty = True
includePercentageUncertainty = True
columnsToSave = None
uncertainty_value = None
uncertainty_representation = "percentage"
# Plotting
plot_validation = True  # Print data validation stats
plot_timeseries = True  # Time series with stats
plot_bandwidth_convergence = True  # Bandwidth vs averaging interval
plot_convergence = True  # Full convergence plots
plot_means = True  # Mean comparison plot
plot_pdf = True  # Probability density function
plot_violinplot = True  # Violin plot
plot_histogram = True  # Histogram
plot_cdf = True  # Cumulative distribution function
plot_power_density = True  # Power density spectrum
plot_comparison = True  # Comprehensive comparison
plot_fluctuations = True  # Concentration fluctuation analysis
# Default ambient conditions (used if CSV not found)
x_source = 0
y_source = 0
z_source = 0
mass_flow_controller = 0.300
calibration_curve = 1.0
calibration_factor = 0
gas_name = 'C12'
gas_factor = 0.5
mol_weight = 29.0
x_measure = 1020
y_measure = 0
z_measure = 5
pressure = 101426.04472
temperature = 23
scale = 400
scaling_factor = 0.5614882
ref_length = 1/400
ref_height = 100/400
full_scale_wtref = 10
full_scale_flow_rate = 0.002
full_scale_temp = 20
full_scale_pressure = 101325

def main():
    if 'windtunnel.concentration.utils' in sys.modules:
        importlib.reload(sys.modules['windtunnel.concentration.utils'])
    print(f"CSV file: {csv_file}")
    if not os.path.exists(csv_file):
        print("CSV file not found")
    uncertainty_results = {}
    conc_ts = {}
    conc_ts.fromkeys(namelist)
    conc_ts_fs = conc_ts
    conc_ts_nd = conc_ts
    dict_conc_ts = conc_ts
    dict_conc_nd = conc_ts
    dict_conc_fs = conc_ts
    data_dict = {}
    data_dict.fromkeys(namelist)
    for name in namelist:
        files = wt.get_files(path, name)
        print(f"files: {files}")
        conc_ts[name] = {}
        conc_ts[name].fromkeys(files)
        # Read ambient conditions per folder if enabled
        if parameters_PerFolder:
            ambient_conditions = wt.PointConcentration.get_ambient_conditions(path=path, name=name, input_file=csv_file)
            if ambient_conditions is not None:
                (x_source, y_source, z_source, x_measure, y_measure, z_measure, pressure, temperature, calibration_curve, mass_flow_controller, calibration_factor, scaling_factor, scale, ref_length, ref_height, gas_name, mol_weight, gas_factor, full_scale_wtref, full_scale_flow_rate, full_scale_temp, full_scale_pressure, config_name) = wt.PointConcentration.read_ambient_conditions(ambient_conditions, name)
        # Process each file
        for file in files:
            # Read ambient conditions per file if enabled
            if not parameters_PerFolder:
                ambient_conditions = wt.PointConcentration.get_ambient_conditions(path=path, name=file, input_file=csv_file)
                if ambient_conditions is not None:
                    (x_source, y_source, z_source, x_measure, y_measure, z_measure, pressure, temperature, calibration_curve, mass_flow_controller, calibration_factor, scaling_factor, scale, ref_length, ref_height, gas_name, mol_weight, gas_factor, full_scale_wtref, full_scale_flow_rate, full_scale_temp, full_scale_pressure, config_name) = wt.PointConcentration.read_ambient_conditions(ambient_conditions, file)
            # Load data and set conditions
            conc_ts[name][file] = wt.PointConcentration.from_file(path + file)
            conc_ts[name][file].ambient_conditions(x_source=x_source, y_source=y_source, z_source=z_source, x_measure=x_measure, y_measure=y_measure, z_measure=z_measure, pressure=pressure, temperature=temperature, calibration_curve=calibration_curve, mass_flow_controller=mass_flow_controller, calibration_factor=calibration_factor, config_name=config_name)
            print("Store information into PointConcentration class objects array")
            conc_ts[name][file].scaling_information(scaling_factor=scaling_factor, scale=scale, ref_length=ref_length, ref_height=ref_height)
            conc_ts[name][file].tracer_information(gas_name=gas_name, mol_weight=mol_weight, gas_factor=gas_factor)
            conc_ts[name][file].full_scale_information(full_scale_wtref=full_scale_wtref, full_scale_flow_rate=full_scale_flow_rate, full_scale_temp=full_scale_temp, full_scale_pressure=full_scale_pressure)
            # Calculate concentrations
            print("Do main calculations")
            conc_ts[name][file].convert_temperature()
            conc_ts[name][file].calc_wtref_mean()
            conc_ts[name][file].calc_model_mass_flow_rate(usingMaxFlowRate="True", applyCalibration="False")
            conc_ts[name][file].calc_net_concentration()
            conc_ts[name][file].calc_c_star()
            conc_ts[name][file].calc_full_scale_concentration()
            # Transform to target scale
            print("Transform scale")
            if full_scale == 'ms':
                dict_conc_ts = conc_ts
            elif full_scale == 'fs':
                dict_conc_ts = conc_ts_fs
                dict_conc_ts[name][file].to_full_scale()
            elif full_scale == 'nd':
                dict_conc_ts = conc_ts_nd
                dict_conc_ts[name][file].to_non_dimensional()
        # Calculate uncertainties
        if calculateUncertainty and saveCombined:
            from windtunnel.concentration.utils import calculate_uncertainties
            if verboseUncertainty:
                metrics_str = ", ".join(uncertaintyMetrics) if uncertaintyMetrics else "all"
                conc_types_str = ", ".join(uncertaintyConcentrationTypes) if uncertaintyConcentrationTypes else "c_star only"
                print(f"Calculating measurement uncertainties for: {metrics_str} ({conc_types_str})")
            try:
                uncertainty_results.update(calculate_uncertainties(conc_ts[name], split_factor=split_factor, uncertainty_threshold=uncertainty_threshold, verbose=verboseUncertainty, metrics_to_calculate=uncertaintyMetrics, concentration_types=uncertaintyConcentrationTypes, include_abs=includeAbsoluteUncertainty, include_pct=includePercentageUncertainty))
            except TypeError as e:
                print(f"Error with new parameters: {e}")
                print("Falling back to basic parameters")
                uncertainty_results.update(calculate_uncertainties(conc_ts[name], split_factor, uncertainty_threshold))
        # Save results
        if saveAll:
            saveTs = True
            saveAvg = True
            saveStats = True
            saveCombined = True
        if saveCombined:
            saveAvg = True
            saveStats = True
        if saveTs or saveAvg or saveStats or saveCombined:
            if osType == "Windows":
                folder = 'Point_Data\\' + name[:name.find('.')] + '\\'
                folder_avg = 'Point_Data_avg\\' + name[:name.find('.')] + '\\'
                folder_stats = 'Point_Data_stats\\' + name[:name.find('.')] + '\\'
            else:
                folder = 'Files/' + 'Point_Data/' + name[:name.find('.')] + '/'
                folder_avg = 'Files/' + 'Point_Data_avg/' + name[:name.find('.')] + '/'
                folder_stats = 'Files/' + 'Point_Data_stats/' + name[:name.find('.')] + '/'
            wt.check_directory(output_path + folder)
            if saveAvg:
                wt.check_directory(output_path + folder_avg)
            if saveStats:
                wt.check_directory(output_path + folder_stats)
            dict_conc_ts[next(iter(conc_ts))][next(iter(conc_ts[next(iter(conc_ts))]))].__check_sum = 8
            # Save all files
            for file in files:
                if saveTs:
                    if full_scale == 'ms':
                        dict_conc_ts[name][file].save2file_ms(file, out_dir=output_path + folder)
                    elif full_scale == 'fs':
                        dict_conc_ts[name][file].save2file_fs(file, out_dir=output_path + folder)
                    elif full_scale == 'nd':
                        dict_conc_ts[name][file].save2file_nd(file, out_dir=output_path + folder)
                    print("Created ts files")
                if saveAvg:
                    dict_conc_ts[name][file].save2file_avg(file, out_dir=output_path + folder_avg)
                    print(f"Created avg files under {output_path + folder_avg}")
                if saveStats:
                    dict_conc_ts[name][file].save2file_fullStats(file, out_dir=output_path + folder_stats)
                    print(f"Created stats files under {output_path + folder_stats}")
            # Save combined file
            if saveCombined:
                from windtunnel.concentration.utils import combine_to_csv
                stats = True
                file_type = "stats" if stats else ""
                file_names = ["_stats_" + file for file in files]
                base_path_local = base_path if base_path else output_path + f"Files/Point_Data_{file_type}/{name[0:-1]}/"
                combine_to_csv(file_names, base_path_local, file_type=file_type, output_filename=output_path + combinedFileName)
                print(output_path + combinedFileName)
                print(f"Created combined file under {output_path + combinedFileName}")
                mapping = {file: conc_ts[name][file].config_name for name in namelist for file in conc_ts[name]}
    return conc_ts, files, namelist

def run_visualizations(conc_ts, files, namelist):
    # Validation stats
    if plot_validation:
        for name in namelist:
            for file in files:
                print(f"\nFile: {file}")
                if hasattr(conc_ts[name][file], "c_star"):
                    print("c_star present")
                print("C_star shape:", conc_ts[name][file].net_concentration.shape)
                print("NaNs present:", np.any(np.isnan(conc_ts[name][file].net_concentration)))
                print("Min/Max:", np.min(conc_ts[name][file].net_concentration), "/", np.max(conc_ts[name][file].net_concentration))
                print(f"Mean: {np.mean(conc_ts[name][file].net_concentration)}")
                print(f"Std: {np.std(conc_ts[name][file].net_concentration)}")
                print(f"Percentiles: {conc_ts[name][file].calc_percentiles(percentiles=[10, 90, 95], var='net_concentration')}")
    # Prepare data
    DataPointsConc = []
    for i in range(len(files)):
        data = conc_ts[namelist[0]][files[i]]
        DataPointsConc.append(data)
    labels = [f"Dataset {i}" for i in range(len(DataPointsConc))]
    # Time series plot
    if plot_timeseries:
        plot_timeseries_with_stats(DataPointsConc, dimensionless=True, labels=labels, color="blue")
    # Bandwidth convergence
    if plot_bandwidth_convergence:
        time_freq = 0.010
        averaging_intervals = [60*0.01, 60*0.1, 60*0.5, 60*1.0, 60*2.0]
        bandwidths = []
        for name in namelist:
            for file in files:
                ts_v_avg = conc_ts[name][file].get_averagedData(name, file, time_freq, averaging_intervals)
                bandwidths.append([np.max(avg) - np.min(avg) for avg in ts_v_avg])
        for i, bandwidth in enumerate(bandwidths):
            print(bandwidth)
            plt.plot(averaging_intervals, bandwidth, "o-", label=f"Dataset {i}")
        plt.xlabel("Averaging Intervals [s]")
        plt.ylabel("Bandwidth [ppmV]")
        plt.legend()
        plt.grid(True)
        plt.show()
    # Full convergence plot
    if plot_convergence:
        plot_type = "both"
        num_points_between = 5
        set_xlog = True
        time_freq = 0.010
        averaging_intervals = [60*0.01, 60*0.1, 60*0.5, 60*1.0, 60*2.0]
        if num_points_between > 0:
            extended_intervals = []
            for i in range(len(averaging_intervals) - 1):
                start = averaging_intervals[i]
                end = averaging_intervals[i + 1]
                points = np.logspace(np.log10(start), np.log10(end), num_points_between + 2)
                extended_intervals.extend(points[:-1])
            extended_intervals.append(averaging_intervals[-1])
            averaging_intervals_use = extended_intervals
        else:
            averaging_intervals_use = averaging_intervals
        if plot_type in ["bandwidth", "both"]:
            bandwidths = []
            for name in namelist:
                for file in files:
                    ts_v_avg = conc_ts[name][file].get_averagedData(name, file, time_freq, averaging_intervals_use, True)
                    bandwidths.append([np.ptp(avg) for avg in ts_v_avg])
            fig1, ax1 = plt.subplots(figsize=(10, 6))
            for i, bandwidth in enumerate(bandwidths):
                ax1.plot(averaging_intervals_use, bandwidth, "o", label=f"Dataset {i}")
            ax1.set_xlabel("Averaging Intervals [s]")
            ax1.set_ylabel("Bandwidth [ppmV]")
            ax1.set_xscale('log')
            ax1.legend()
            ax1.grid(True, which='both', alpha=0.3)
            plt.tight_layout()
            plt.show()
    # Mean comparison
    if plot_means:
        create_means(DataPointsConc, 0.5, dimensionless=False, labels=None, xLabel="Datasets", yLabel="Mean concentration[ppmV]")
    # PDF
    if plot_pdf:
        create_pdf(DataPointsConc, dimensionless="True", labels=None, xLabel="Concentration[-]", yLabel="Density")
    # Violin plot
    if plot_violinplot:
        create_violinplot(DataPointsConc)
    # Histogram
    if plot_histogram:
        create_histogram(DataPointsConc, dimensionless="False", labels=None, xLabel="Datasets", yLabel="Concentration[ppmV]")
    # CDF
    if plot_cdf:
        create_cdf(DataPointsConc, dimensionless=False, labels=None, xLabel="Concentration[-]", yLabel=None)
    # Power density
    if plot_power_density:
        powerDensityPlot(DataPointsConc, dimensionless="False", plot=True, labels=None)
    # Comprehensive comparison
    if plot_comparison:
        functionsForOverview = ["Histogram", "BoxPlot", "Pdf", "Cdf", "Means", "PowerDensity"]
        DataPointsConc_subset = [conc_ts[namelist[0]][files[0]], conc_ts[namelist[0]][files[1]]]
        compare_point_concentrations_3(DataPointsConc_subset, functionsForOverview)
    # Fluctuation analysis
    if plot_fluctuations:
        threshold_type = "ratio"
        threshold_method = "mean"
        intermittency_threshold = 1.5
        for name in namelist:
            for file in files:
                conc_ts[name][file].analyze_concentration_fluctuations(dimensionless="False", intermittency_threshold=intermittency_threshold, threshold_method=threshold_method)
    ### ADD NEW PLOTS HERE ###
    # Example: if plot_my_new_plot:
    #     my_plotting_function(DataPointsConc, param1=value1)
    print("Visualization complete")

if __name__ == "__main__":
    conc_ts, files, namelist = main()
    run_visualizations(conc_ts, files, namelist)
    print("Analysis complete")
