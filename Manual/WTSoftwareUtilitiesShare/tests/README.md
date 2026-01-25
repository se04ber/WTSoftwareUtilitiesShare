# Tests

This folder contains **unit tests** (fast guardrails) and **regression tests**
(compare against known-good example outputs).

## Why this exists

This project is used by students and non-IT users. To keep the code safe to
refactor, we added automated tests that:

- Re-run the **UBA_GA concentration** pipeline and compare against a golden CSV
  and a hand-made Excel control table.
- Re-run the **flow (vertical profile)** pipeline and compare against saved
  reference outputs used by plots (turb table + spectra).

## Folder structure

- **`tests/unit/`**: fast checks for formulas and config branches
- **`tests/regression/`**: golden-file regression tests (bigger but still CI-friendly)
- **`tests/support/`**: shared helpers and configuration
- **`tests/fixtures/`**: copied input data used by tests (easy to find)

## Changing which files are used for testing

Edit **one file**:

- `tests/support/config.py`

That’s the only place where the tested filenames are defined. After changing
the names, copy the matching files into `tests/fixtures/` using the same folder
structure.

## What is tested right now

### Concentration (UBA_GA)

- **Regression (golden files + control table)**:
  - `tests/regression/concentration/test_concentration_uba_ga_regression.py`
  - Compares computed results vs `tests/fixtures/concentration/expected/combined_data.csv`
  - Cross-checks key intermediate values vs
    `tests/fixtures/concentration/excel/Beispiel Umrechnung zur Kontrolle.xlsx`
  - Also checks saving `_avg_`/`_stats_` files + combining to CSV works end-to-end
- **Formulas + output files**:
  - `tests/unit/concentration/test_pointconcentration_formulas.py`
  - Guards `net_concentration`, `c_star`, full-scale concentration, time scaling, and save-file creation
- **Config branch coverage**:
  - `tests/unit/concentration/test_pointconcentration_config_paths.py`
  - Covers `calc_model_mass_flow_rate()` branches (max-flow, calibration, legacy path)
- **Multiple prefixes (no overwrites)**:
  - `tests/unit/concentration/test_multi_prefix_processing.py`
- **`parameters_PerFolder=True` path**:
  - `tests/unit/concentration/test_parameters_per_folder.py`

### Flow (vertical profile)

- **Regression against saved “turb table” (plot input)**:
  - `tests/regression/flow/test_flow_vertical_profile_regression.py`
  - Uses `tests/fixtures/flow/*` and matches rows in `*_turb.txt` for a small
    representative subset of time series files.
- **Regression against saved spectra file (plot input)**:
  - `tests/regression/flow/test_flow_spectra_regression.py`
- **Regression guardrails for fitted parameters**:
  - `tests/regression/flow/test_alpha_and_z0_regression.py`
- **Basic data contracts (shape/ranges)**:
  - `tests/unit/flow/test_flow_profile_contracts.py`

### Utilities / performance

- `tests/unit/common/test_utils_and_uncertainty.py`
  - Checks `get_files()` prefix matching + sorting
  - Checks `calculate_uncertainties()` output keys are stable (shape/contract)
- `tests/unit/common/test_performance_smoke.py`
  - Soft runtime budget for a small representative workload
  - Can be overridden with `WT_TEST_MAX_SECONDS`

## Running tests locally

From `Manual/WTSoftwareUtilitiesShare/`:

```bash
python -m pip install -r requirements.txt -r requirements-dev.txt
python -m pytest
```

## CI

GitHub Actions runs pytest on every push/PR:

- `.github/workflows/tests.yml`


