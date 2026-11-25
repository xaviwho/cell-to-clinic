"""Fig. 5: Cell-to-Clinic II forest plot (osf.io/4pzxa). Aggregates only (results/c2c2_models_da.json).
Writes paper/figures/fig5_c2c2_forest.svg (hand-built SVG: matplotlib is blocked by Windows App Control).

Same conventions as Fig. 3: partial r with 90% CI; shaded equivalence region +/-0.40; class encoded twice
(colour + marker shape) and written on every row. Supported = ink square (as in Fig. 3); equivalent =
slot 1 circle; inconclusive = slot 2 diamond; not tested = muted open circle.
Sensitivity rows are registered sensitivity analyses of the primary, not separate tests.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "results" / "c2c2_models_da.json"
OUT = ROOT / "paper" / "figures" / "fig5_c2c2_forest.svg"
INK, INK2, MUTED, GRID, SURFACE, BAND = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb", "#f0efec"
STYLE = {"supported": (INK, "s", "Supported"), "equivalent": ("#2a78d6", "o", "Equivalent"),
         "inconclusive": ("#eb6834", "D", "Inconclusive"), "not_tested": (MUTED, "open", "Not tested")}
BOUND = 0.40


def rows(r):
    s = r["sensitivity"]
    R = [("h", "Primary  H1-DA: motor progression (OFF-state UPDRS III slope)"),
         ("r", "SNCA in DA neurons (D2)", r["H1_DA_primary"], "p", True),
         ("h", "Sensitivity analyses of the primary"),
         ("r", "Old DA label (TH ≥ 1 UMI)", s["S1_old_label"]["H1_DA_primary"], "p", False),
         ("r", "Composition-adjusted", s["S2_composition_adjusted"]["H1_DA_primary"], "p", False),
         ("r", "Idiopathic PD only", s["S3_idiopathic"]["H1_DA_primary"], "p", False),
         ("r", "Five-level batch covariate", s["S4_batch"]["H1_DA_primary"], "p", False),
         ("r", "First 5 years after diagnosis", s["S6_first_5y"]["H1_DA_primary"], "p", False),
         ("h", "Secondary family (BH over 6 tests)")]
    sec = r["secondary"]
    for key, label in [("H1_D4", "Lysosome (D4): UPDRS III slope"), ("H2_D2", "SNCA (D2): putamen SBR slope"),
                       ("H2_D4", "Lysosome (D4): putamen SBR slope"), ("H3_SAA_D2", "SNCA (D2): SAA status"),
                       ("H3_SAA_D4", "Lysosome (D4): SAA status"), ("H3_NSDISS_D2", "SNCA (D2): NSD-ISS stage"),
                       ("H3_NSDISS_D4", "Lysosome (D4): NSD-ISS stage")]:
        R.append(("r", label, sec[key], "q", False))
    R += [("h", "H4  Bulk RNA, day 65 (single test)"),
          ("r", "OXPHOS, composition-adjusted: UPDRS III slope", r["H4_bulk_oxphos"], "p", False),
          ("h", "Not tested: failed the reliability gate (descriptive)")]
    for k, label in [("D1", "DA fraction (D1): UPDRS III slope"), ("D3", "OXPHOS in DA neurons (D3): UPDRS III slope")]:
        R.append(("r", label, dict(r["not_tested_descriptive_H1"][k], **{"class": "not_tested", "why": "G-icc"}), "p", False))
    return R


def cls(v):
    c = v.get("class", "not_tested")
    return "not_tested" if c == "dropped_by_A0-5" else c


def note(v, stat):
    c = v.get("class")
    if c == "dropped_by_A0-5":
        return "not tested (A0-5)"
    if c == "not_tested":
        return f"not tested ({v['why']})"
    val = v.get(stat, v.get("p"))
    return f"{STYLE[c][2].lower()}, {stat} = {val:.3f}" if val < 0.1 else f"{STYLE[c][2].lower()}, {stat} = {val:.2f}"


def marker(kind, x, y, col, big=False):
    rr = 6.5 if big else 5
    if kind == "o":
        return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{rr}" fill="{col}" stroke="{SURFACE}" stroke-width="2"/>'
    if kind == "open":
        return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{rr - 0.5}" fill="{SURFACE}" stroke="{col}" stroke-width="2"/>'
    if kind == "s":
        return f'<rect x="{x - rr:.1f}" y="{y - rr:.1f}" width="{2 * rr}" height="{2 * rr}" rx="1.5" fill="{col}" stroke="{SURFACE}" stroke-width="2"/>'
    d = rr + 1.5
    return f'<polygon points="{x:.1f},{y - d:.1f} {x + d:.1f},{y:.1f} {x:.1f},{y + d:.1f} {x - d:.1f},{y:.1f}" fill="{col}" stroke="{SURFACE}" stroke-width="2"/>'


def main():
    r = json.load(open(SRC))
    R = rows(r)
    W, top, rh, hh = 820, 30, 22, 30
    x0, x1 = 300, 590
    H = top + sum(hh if t[0] == "h" else rh for t in R) + 95
    X = lambda v: x0 + (v + 1) / 2 * (x1 - x0)
    y_end = top + sum(hh if t[0] == "h" else rh for t in R) + 6
    e = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
         f'font-family="DejaVu Sans, Helvetica, Arial, sans-serif">', f'<rect width="{W}" height="{H}" fill="{SURFACE}"/>',
         f'<rect x="{X(-BOUND):.1f}" y="{top}" width="{X(BOUND) - X(-BOUND):.1f}" height="{y_end - top}" fill="{BAND}"/>']
    for b in (-BOUND, BOUND):
        e.append(f'<line x1="{X(b):.1f}" x2="{X(b):.1f}" y1="{top}" y2="{y_end}" stroke="{MUTED}" stroke-width="1" stroke-dasharray="3 3"/>')
    e.append(f'<line x1="{X(0):.1f}" x2="{X(0):.1f}" y1="{top}" y2="{y_end}" stroke="{INK2}" stroke-width="1"/>')
    y = top
    for t in R:
        if t[0] == "h":
            y += hh
            e.append(f'<text x="16" y="{y - 9}" font-size="11.5" font-weight="bold" fill="{INK}">{t[1]}</text>')
            continue
        _, label, v, stat, big = t
        yc = y + rh / 2
        c = cls(v); col, mk, _ = STYLE[c]
        lo, hi = v["ci90"]
        bold = ' font-weight="bold"' if big else ""
        e.append(f'<text x="{x0 - 10}" y="{yc + 4:.1f}" font-size="10.5" fill="{INK}" text-anchor="end"{bold}>{label}</text>')
        e.append(f'<line x1="{X(lo):.1f}" x2="{X(hi):.1f}" y1="{yc:.1f}" y2="{yc:.1f}" stroke="{col}" stroke-width="2" stroke-linecap="round"/>')
        e.append(marker(mk, X(v["r"]), yc, col, big))
        e.append(f'<text x="{x1 + 14}" y="{yc + 4:.1f}" font-size="10" fill="{INK2}">n = {v["n_donors"]} · {note(v, stat)}</text>')
        y += rh
    e.append(f'<line x1="{x0}" x2="{x1}" y1="{y_end}" y2="{y_end}" stroke="{MUTED}" stroke-width="1"/>')
    for tv in (-0.8, -0.4, 0.0, 0.4, 0.8):
        e.append(f'<text x="{X(tv):.1f}" y="{y_end + 16}" font-size="10" fill="{INK2}" text-anchor="middle">' + f"{tv:.1f}".replace("-", "−") + "</text>")
    e.append(f'<text x="{(x0 + x1) / 2}" y="{y_end + 34}" font-size="10.5" fill="{INK2}" text-anchor="middle">'
             'Partial correlation r with 90% CI (phenotype per SD)</text>')
    lx, ly = 16, y_end + 62
    for k in ("supported", "equivalent", "inconclusive", "not_tested"):
        col, mk, name = STYLE[k]
        e.append(marker(mk, lx + 6, ly - 4, col))
        e.append(f'<text x="{lx + 18}" y="{ly}" font-size="10" fill="{INK}">{name}</text>')
        lx += 22 + 7 * len(name) + 18
    e.append(f'<rect x="{lx}" y="{ly - 9}" width="16" height="10" fill="{BAND}" stroke="{MUTED}" stroke-width="0.6" stroke-dasharray="2 2"/>')
    e.append(f'<text x="{lx + 22}" y="{ly}" font-size="10" fill="{INK}">Equivalence region (±0.40)</text>')
    e.append(f'<text x="16" y="{H - 10}" font-size="9" fill="{MUTED}">Cell-to-Clinic II (osf.io/4pzxa), 47 PD donors. '
             'Primary and H4 are single tests (p); secondary q = BH. Sensitivity rows re-estimate the primary.</text>')
    e.append("</svg>")
    OUT.write_text(chr(10).join(e), encoding="utf-8")
    print("written", OUT)


if __name__ == "__main__":
    main()
