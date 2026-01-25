# Tests

This folder contains **unit tests** (fast, formula-level guardrails) and
**regression tests** (end-to-end checks against known-good example outputs).

## Why this exists

This project is used by students and non-IT users. To keep the code safe to
refactor, we added automated tests that:

- Re-run the **UBA_GA concentration** pipeline and compare results to a golden CSV + the original Excel control table.
- Re-run a **flow (vertical profile)** pipeline and compare intermediate results to the saved `*_turb.txt` reference (before plotting).

## Folder structure (SOTA)

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

- **End-to-end regression**:
  - `tests/regression/test_concentration_uba_ga_regression.py`
  - Compares computed results vs `tests/fixtures/concentration/expected/combined_data.csv`
  - Cross-checks intermediate values vs `tests/fixtures/concentration/excel/Beispiel Umrechnung zur Kontrolle.xlsx`
- **Intermediate formulas + output files**:
  - `tests/unit/test_pointconcentration_formulas.py`
  - Locks down `net_concentration`, `c_star`, full-scale concentration, time scaling, and save-file creation
- **Config branch coverage**:
  - `tests/unit/test_pointconcentration_config_paths.py`
  - Covers the different `calc_model_mass_flow_rate()` branches (max-flow, calibration, legacy path)

### Flow (vertical profile)

- **Regression against saved intermediate results (before plots)**:
  - `tests/regression/test_flow_vertical_profile_regression.py`
  - Uses `tests/fixtures/flow/*` and matches rows in `*_turb.txt` for a small representative subset of time series files.

### Utilities

- `tests/unit/test_utils_and_uncertainty.py`
  - Checks `get_files()` sorting/selection
  - Checks `calculate_uncertainties()` output keys are stable (shape contract)

---

## TODO / missing tests (planned next steps)

This is the current shortlist of missing coverage to add before larger refactors.

### Flow: more outputs used by plots

Right now we regression-test the **pre-plot “turb table”** (`*_turb.txt`) for a small subset
of files. We still need tests for other plot-relevant computed outputs, e.g.:

- **Spectra** (power spectral density): intermediate arrays used by spectra plots
- **Wavelet transform** outputs (if used by `plot_wavelet`)
- **Convergence** metrics used by convergence plots
- **Reference spectra alignment** (matching the correct reference dataset by height)

Goal: tests should validate the **computed data that feeds plots**, not the plot images.

### Concentration: more settings / configuration paths

Right now we cover UBA_GA with a single parameter CSV shape and `parameters_PerFolder=False`.
We still need tests for:

- **Ambient conditions per folder vs per file** (`parameters_PerFolder=True/False`)
- **Parameter CSV variants**:
  - “one column per file” (current)
  - “one column for folder / config” (if supported)
  - multiple configs in one file (multiple columns / names)
- **Multiple `namelist` entries**:
  - read and process multiple prefixes in one run
  - ensure outputs are combined correctly and don’t overwrite each other

### Performance guardrail

Add a simple performance test that measures **processing time** for a known workload
(e.g. a small fixed subset of UBA_GA or flow files) so we can detect large slowdowns
when refactoring. This should be a **soft threshold** (CI-safe) and mainly used as a
trend/alert, not a strict benchmark.

Status: implemented as `tests/unit/common/test_performance_smoke.py` (budget can be
overridden with `WT_TEST_MAX_SECONDS`).

### Test organization (further cleanup)

We should further separate tests by domain for fast navigation:

- `tests/unit/concentration/` and `tests/unit/flow/`
- `tests/regression/concentration/` and `tests/regression/flow/`

## Running tests locally

From the repo root:

```bash
pip install -r requirements.txt -r requirements-dev.txt
pytest
```

## CI

GitHub Actions runs pytest on every push/PR:

- `.github/workflows/tests.yml`


