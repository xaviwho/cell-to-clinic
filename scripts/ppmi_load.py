"""Build tidy clinical tables from a PPMI data cut (run only after the OSF registration is approved).

Analysis-relevant logic, fixed here before any clinical data is seen:
- OFF state (O1): untreated at the exam (treated == 0) OR state recorded as "OFF". ON state
  (sensitivity 6): treated == 1 AND state "ON". Repeated exams of the same state on the same date
  are averaged.
- Time: years since PD diagnosis = (exam date - diagnosis date) / 365.25. Sensitivity 9 uses years
  since the PPMI baseline (BL) exam date.
- Age at PD onset = (symptom-onset date - birth date) / 365.25. Disease duration at PPMI baseline =
  (BL date - diagnosis date) / 365.25. A donor missing either is excluded from models that need it.
- DaTscan putamen SBR (O2) = mean of left and right; scans missing either side are dropped.
- LEDD (sensitivity 6) is attached to the UPDRS III exam with the same EVENT_ID.
- SAA (O3) = result of the earliest dated assay with a positive or negative result.
- NSD-ISS (O3) = stage at the latest dated record with a stage, coded numerically as the integer
  stage (subtype letters A/B dropped, e.g. "2B" -> 2).
- Month/year-only dates are set to day 15.
Only PATNOs in --patnos (the FOUNDIN-PD PD donors) are kept; all other participants are dropped on read.
File patterns and column names: config/ppmi_columns.json.
"""
from __future__ import annotations

import glob, json, re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def load_colmap(path=ROOT / "config" / "ppmi_columns.json"):
    with open(path) as f:
        return json.load(f)


def parse_date(s, fmts, day15):
    if pd.isna(s):
        return pd.NaT
    s = str(s).strip()
    for f in fmts:
        try:
            d = datetime.strptime(s, f)
            if f == "%m/%Y":
                d = d.replace(day=day15)
            return pd.Timestamp(d)
        except ValueError:
            continue
    return pd.NaT


def _read(ppmi_dir, spec, patnos, cm):
    files = sorted(glob.glob(str(Path(ppmi_dir) / spec["file"])))
    if not files:
        raise FileNotFoundError(f"no file matches {spec['file']}")
    df = pd.concat([pd.read_csv(f, low_memory=False) for f in files], ignore_index=True)
    cols = [v for k, v in spec.items() if k not in ("file",) and isinstance(v, str) and k not in
            ("baseline_event",)]
    miss = [c for c in cols if c not in df.columns]
    if miss:
        raise KeyError(f"{spec['file']}: columns not found {miss} (fix names in config/ppmi_columns.json)")
    df = df[df[spec["PATNO"]].astype(str).isin(patnos)].copy()
    df["PATNO"] = df[spec["PATNO"]].astype(str)
    for k in ("date", "birth", "dx_date", "onset_date"):
        if k in spec:
            df[k] = df[spec[k]].map(lambda s: parse_date(s, cm["date_formats"], cm["date_day_if_month_only"]))
    return df


def yrs(a, b):
    return (a - b).dt.days / 365.25


def nsd_numeric(v):
    m = re.match(r"^\s*(\d+)", str(v))
    return float(m.group(1)) if m else np.nan


def build(ppmi_dir, patnos, cm=None):
    cm = cm or load_colmap()
    T = cm["tables"]
    patnos = set(map(str, patnos))
    dem = _read(ppmi_dir, T["demographics"], patnos, cm)
    dem = dem.groupby("PATNO").first()
    dx = _read(ppmi_dir, T["diagnosis"], patnos, cm).groupby("PATNO")[["dx_date", "onset_date"]].first()
    u3 = _read(ppmi_dir, T["updrs3"], patnos, cm)
    bl = u3[u3[T["baseline"]["event"]].astype(str) == T["baseline"]["baseline_event"]].groupby("PATNO")["date"].min()

    donors = pd.DataFrame(index=sorted(patnos))
    donors["sex_female"] = (dem[T["demographics"]["sex"]] == T["demographics"]["sex_female_value"]).astype(float).where(
        dem[T["demographics"]["sex"]].notna())
    donors["birth"] = dem["birth"]
    donors = donors.join(dx).join(bl.rename("bl_date"))
    donors["age_onset"] = yrs(donors["onset_date"], donors["birth"])
    donors["dur_baseline"] = yrs(donors["bl_date"], donors["dx_date"])

    # UPDRS III OFF / ON
    s = T["updrs3"]
    u3["total"] = pd.to_numeric(u3[s["total"]], errors="coerce")
    tr = pd.to_numeric(u3[s["treated"]], errors="coerce")
    st = u3[s["state"]].astype(str).str.upper().str.strip()
    u3["state"] = np.where((tr == 0) | (st == "OFF"), "OFF", np.where((tr == 1) & (st == "ON"), "ON", "other"))
    u3["event"] = u3[s["event"]].astype(str)
    u3 = u3.dropna(subset=["total", "date"])
    upd3 = u3[u3.state != "other"].groupby(["PATNO", "date", "event", "state"], as_index=False)["total"].mean()
    led = _read(ppmi_dir, T["ledd"], patnos, cm)
    led["ledd"] = pd.to_numeric(led[T["ledd"]["ledd"]], errors="coerce")
    led["event"] = led[T["ledd"]["event"]].astype(str)
    upd3 = upd3.merge(led.groupby(["PATNO", "event"], as_index=False)["ledd"].mean(), on=["PATNO", "event"], how="left")

    u2 = _read(ppmi_dir, T["updrs2"], patnos, cm)
    u2["total"] = pd.to_numeric(u2[T["updrs2"]["total"]], errors="coerce")
    upd2 = u2.dropna(subset=["total", "date"]).groupby(["PATNO", "date"], as_index=False)["total"].mean()

    ds = _read(ppmi_dir, T["datscan"], patnos, cm)
    ds["total"] = (pd.to_numeric(ds[T["datscan"]["put_l"]], errors="coerce") +
                   pd.to_numeric(ds[T["datscan"]["put_r"]], errors="coerce")) / 2
    dat = ds.dropna(subset=["total", "date"]).groupby(["PATNO", "date"], as_index=False)["total"].mean()

    for v in (upd3, upd2, dat):
        v["t_dx"] = yrs(v["date"], v["PATNO"].map(donors["dx_date"]))
        v["t_bl"] = yrs(v["date"], v["PATNO"].map(donors["bl_date"]))

    sa = _read(ppmi_dir, T["saa"], patnos, cm)
    r = sa[T["saa"]["result"]].astype(str).str.strip()
    sa["saa_pos"] = np.where(r.isin(T["saa"]["positive_values"]), 1.0,
                             np.where(r.isin(T["saa"]["negative_values"]), 0.0, np.nan))
    sa = sa.dropna(subset=["saa_pos", "date"]).sort_values("date")
    donors["saa_pos"] = sa.groupby("PATNO")["saa_pos"].first()

    nd = _read(ppmi_dir, T["nsd_iss"], patnos, cm)
    nd["stage"] = nd[T["nsd_iss"]["stage"]].map(nsd_numeric)
    nd = nd.dropna(subset=["stage", "date"]).sort_values("date")
    donors["nsd_iss"] = nd.groupby("PATNO")["stage"].last()
    donors.index.name = "PATNO"
    return {"donors": donors.reset_index(), "updrs3": upd3, "updrs2": upd2, "datscan": dat}
