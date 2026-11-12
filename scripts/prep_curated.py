"""Format conversion only: extract per-visit LEDD and NSD-ISS stage from the PPMI Curated Data Cut
(.xlsx) into a CSV that scripts/ppmi_load.py (frozen in prereg-v1.0.0, CSV-only) can read.

Columns copied unchanged as text: PATNO, EVENT_ID, visit_date, LEDD, NSD_STAGE. No rows are
filtered or transformed here. The analysis logic stays in ppmi_load.py: LEDD attached by EVENT_ID;
NSD-ISS = latest dated stageable record, integer stage.
Run: python -I scripts/prep_curated.py <curated.xlsx> <out.csv>
"""
import sys
from pathlib import Path

import pandas as pd

COLS = ["PATNO", "EVENT_ID", "visit_date", "LEDD", "NSD_STAGE"]


def main(src, out):
    xl = pd.ExcelFile(src)
    sheet = xl.sheet_names[0]  # the data sheet is named by the cut date; others are dictionary/info
    df = pd.read_excel(xl, sheet, usecols=COLS, dtype=str)
    out = Path(out)
    root_data = Path(__file__).resolve().parents[1] / "data"
    if root_data not in out.resolve().parents:
        sys.exit("output must be inside data/")
    df.to_csv(out, index=False)
    print(f"sheet '{sheet}': {len(df)} rows x {len(COLS)} columns -> {out.name}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
