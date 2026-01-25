"""
Central place to choose which input files are used for regression tests.

If you want to switch test datasets, change only the names here (and copy the
matching files into `tests/fixtures/`).
"""

from __future__ import annotations


# ---- Concentration (UBA_GA) -------------------------------------------------

UBA_GA_PREFIX = "UBA_GA_02_04_01_000_1_001"

# Explicit list so tests are stable and easy to tweak.
UBA_GA_FILES = [
    f"{UBA_GA_PREFIX}.txt.ts#0",
    f"{UBA_GA_PREFIX}.txt.ts#1",
    f"{UBA_GA_PREFIX}.txt.ts#2",
    f"{UBA_GA_PREFIX}.txt.ts#3",
    f"{UBA_GA_PREFIX}.txt.ts#4",
    f"{UBA_GA_PREFIX}.txt.ts#5",
]

UBA_GA_PARAMETER_CSV = "ambient_conditions_.UBA_GA.csv"
UBA_GA_EXPECTED_COMBINED_CSV = "combined_data.csv"
UBA_GA_CONTROL_EXCEL = "Beispiel Umrechnung zur Kontrolle.xlsx"


# ---- Flow (vertical profile) ------------------------------------------------

FLOW_NAME = "UBA_BL_Rep2025_UW_01"

# Keep it fast in CI: pick a representative subset (bottom/mid/top).
# Index is derived from the numeric suffix (000001 -> 0, etc.).
FLOW_TS_FILES = [
    f"{FLOW_NAME}.000001.txt",
    f"{FLOW_NAME}.000012.txt",
    f"{FLOW_NAME}.000024.txt",
]

FLOW_WTREF_FILE = f"{FLOW_NAME}_wtref.txt"
FLOW_TURB_GOLDEN = f"{FLOW_NAME}_turb.txt"
FLOW_SPECTRA_GOLDEN = f"spectra_{FLOW_NAME}.000024.txt"

# These match `windtunnel/example_data_analysis.py` defaults that produced the golden file.
FLOW_SCALE_FACTOR = 1.0
FLOW_SCALE = 400.0
FLOW_DATA_ND = 0  # 0 => nondimensionalise() is called in legacy script


