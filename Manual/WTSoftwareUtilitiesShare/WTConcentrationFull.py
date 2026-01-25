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
from windtunnel.pipelines.concentration_point import run_point_concentration

#!/usr/bin/env python3
# ============================================================
# Wind tunnel point concentration analysis
# Configuration section (edit here)
# ============================================================

# Paths and filenames
path_dir="/home/sabrina/Desktop/Schreibtisch/Arbeit_2025/WTSoftwareUtilitiesShare"   # base path
path=f"{path_dir}/Data/InputData/Beispiel Umrechnung zur Kontrolle/"                # raw input data
namelist=['UBA_GA_02_04_01_000_1_001']                                              # prefix of raw input files (can be multiple)
output_path=f"{path_dir}/Data/Results/"                                            # output folder
# Ambient / parameter file (if not found, default values below are used)
csv_file=f"{path_dir}/Data/ParameterFiles/ambient_conditions_.UBA_GA.csv"
parameters_PerFolder=False 
                                                        # True: one parameter set per folder, False: per file
# Scaling and processing
# nd = non-dimensional, ms = model scale, fs = full scale
full_scale='ms'
# Legacy: apply down-averaging after calculation
applyPostprocessing=True
averageInterval=60          # s
measurementFreq=0.005       # Hz
averagingColumns=["net_concentration"]

# Output and saving
# Output: time series, averages, statistics, combined file
osType="Linux"
outputName=None
saveTs=True
saveAvg=True
saveStats=True
saveCombined=True
combinedFileName="combined_file_test.csv"
base_path=None
saveAll=True
"""
# Uncertainty analysis (Legacy)
calculateUncertainty=True
saveUncertainties=True
saveConfigNames=True
split_factor=1.8
uncertainty_threshold=1e-4
verboseUncertainty=True
uncertaintyMetrics=None
uncertaintyConcentrationTypes=None
includeAbsoluteUncertainty=True
includePercentageUncertainty=True
# Combined file options
columnsToSave=None
# Legacy uncertainty (for plotting)
uncertainty_value=None
uncertainty_representation="percentage"  # or "absoluteValue"
"""
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

# Configuration end
# ============================================================

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
    if "windtunnel.concentration.utils" in sys.modules:
        importlib.reload(sys.modules["windtunnel.concentration.utils"])

    result = run_point_concentration(
        path=path,
        namelist=namelist,
        csv_file=csv_file,
        parameters_per_folder=parameters_PerFolder,
        full_scale=full_scale,
        save_all=saveAll,
        save_ts=saveTs,
        save_avg=saveAvg,
        save_stats=saveStats,
        save_combined=saveCombined,
        os_type=osType,
        output_path=output_path,
        combined_filename=combinedFileName,
    )

    # Keep return signature stable for existing downstream notebook cells.
    # `files` was previously the last computed list; we expose the first entry.
    first_name = namelist[0] if namelist else None
    files = result.files_by_name.get(first_name, []) if first_name else []
    return result.conc_ts, files, result.namelist

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
