from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from support.concentration_uba_ga import get_paths as _paths
from support.concentration_uba_ga import list_uba_ga_files
from support.concentration_uba_ga import process_one_file as _process_one_file


def _build_expected_df(expected_combined_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(expected_combined_csv)
    df = df.copy()
    # Normalize windows-style paths to bare filename; tests run on Linux in CI.
    df["Filename"] = (
        df["Filename"].astype(str).str.replace("\\\\", "/", regex=True).str.split("/").str[-1]
    )
    # The stored golden file references saved avg files (prefix `_avg_`).
    # For comparisons we normalize to the raw input filename.
    df["Filename"] = df["Filename"].astype(str).str.replace(r"^_avg_", "", regex=True)
    return df.sort_values("Filename").reset_index(drop=True)


def _build_actual_df(files: list[str], pcs: dict[str, object]) -> pd.DataFrame:
    rows = []
    for file in files:
        pc = pcs[file]
        rows.append(
            {
                "Filename": file,
                "X_fs [m]": float(pc.x) / 1000.0 * float(pc.scale),
                "Y_fs [m]": float(pc.y) / 1000.0 * float(pc.scale),
                "Z_fs [m]": float(pc.z) / 1000.0 * float(pc.scale),
                "Avg_c_star [-]": float(np.nanmean(pc.c_star)),
                "Avg_net_concentration [ppmV]": float(np.nanmean(pc.net_concentration)),
                "Avg_full_scale_concentration [ppmV]": float(
                    np.nanmean(pc.full_scale_concentration)
                ),
            }
        )

    return pd.DataFrame(rows).sort_values("Filename").reset_index(drop=True)


def _load_excel_control_table(excel_path: Path) -> pd.DataFrame:
    """
    The sheet is a hand-made calculation table. We locate the header row that
    starts with "File name" and then parse the rows below it as the actual table.
    """
    raw = pd.read_excel(excel_path, sheet_name=0, header=None)

    header_row_idx = None
    for i in range(len(raw)):
        first = str(raw.iloc[i, 0]).strip()
        second = str(raw.iloc[i, 1]).strip()
        if first == "File name" and second == "Run number":
            header_row_idx = i
            break
    if header_row_idx is None:
        raise RuntimeError(
            "Could not find the concentration table header row ('File name' / 'Run number') in the Excel sheet."
        )

    header = raw.iloc[header_row_idx].tolist()
    table = raw.iloc[header_row_idx + 1 :].copy()
    table.columns = header
    table = table.dropna(how="all")
    table = table[table["File name"].notna()].copy()
    table["File name"] = table["File name"].astype(str).str.strip()
    return table


def test_uba_ga_pipeline_matches_expected_combined_csv():
    p = _paths()

    assert p.input_dir.exists()
    assert p.parameter_csv.exists()
    assert p.expected_combined_csv.exists()

    files = list_uba_ga_files(p)
    pcs = {f: _process_one_file(p, f) for f in files}

    expected = _build_expected_df(p.expected_combined_csv)
    actual = _build_actual_df(files, pcs)

    assert list(actual["Filename"]) == list(expected["Filename"])

    numeric_cols = [
        "X_fs [m]",
        "Y_fs [m]",
        "Z_fs [m]",
        "Avg_c_star [-]",
        "Avg_net_concentration [ppmV]",
        "Avg_full_scale_concentration [ppmV]",
    ]
    # The golden CSV stores values rounded to ~4 decimals, so we compare with a
    # small absolute tolerance.
    for col in numeric_cols:
        if col == "Avg_c_star [-]":
            np.testing.assert_allclose(
                actual[col].to_numpy(), expected[col].to_numpy(), rtol=0.0, atol=5e-5
            )
        else:
            np.testing.assert_allclose(
                actual[col].to_numpy(), expected[col].to_numpy(), rtol=0.0, atol=5e-4
            )


def test_excel_control_table_matches_key_intermediate_values():
    p = _paths()
    assert p.excel_control.exists()

    files = list_uba_ga_files(p)
    pcs = {f: _process_one_file(p, f) for f in files}

    table = _load_excel_control_table(p.excel_control)
    by_file = table.set_index("File name")

    for file in files:
        pc = pcs[file]

        # Excel: C_mean = measured concentration reduced by background concentration.
        excel_c_mean = float(by_file.loc[file, "C_mean [ppmV]"])
        np.testing.assert_allclose(
            float(np.nanmean(pc.net_concentration)), excel_c_mean, rtol=1e-6, atol=1e-6
        )

        # Excel: Q_standard is the flow rate under standard conditions.
        #
        # Note: the current implementation multiplies `open_rate` by 10 internally.
        # The control table stores the flow rate already in [l/h] (Q_ambient).
        excel_q_ambient = float(by_file.loc[file, "Q_ambient [l/h]"])
        np.testing.assert_allclose(
            float(pc.mass_flow_rate) / 10.0, excel_q_ambient, rtol=0.0, atol=1e-6
        )


def test_saving_avg_and_stats_and_combining_roundtrip(tmp_path: Path):
    p = _paths()
    files = list_uba_ga_files(p)
    pcs = {f: _process_one_file(p, f) for f in files}

    avg_dir = tmp_path / "Point_Data_avg" / "UBA_GA_02_04_01_000_1_00"
    stats_dir = tmp_path / "Point_Data_stats" / "UBA_GA_02_04_01_000_1_00"
    avg_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)

    for file in files:
        pcs[file].save2file_avg(file, out_dir=str(avg_dir) + "/")
        pcs[file].save2file_fullStats(file, out_dir=str(stats_dir) + "/")

    from windtunnel.concentration.utils import combine_to_csv

    avg_file_names = ["_avg_" + f for f in files]
    out_csv = tmp_path / "combined_data.csv"
    df = combine_to_csv(
        avg_file_names, base_path=str(tmp_path), file_type="avg", output_filename=str(out_csv)
    )
    assert df is not None
    assert out_csv.exists()
    assert len(df) == len(files)


