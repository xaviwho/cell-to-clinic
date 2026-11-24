"""A4: design guide v2. EXPLORATORY. Aggregates only.

Combines the across-batch variance components measured in Jerber (A3: line, batch and per-cell noise)
with the clinical-slope reliabilities measured in PPMI for the FOUNDIN donors (E4: lambda = 0.64 for
OFF-state UPDRS III, 0.54 for putamen SBR), using cell2clinic.design.

"Typical gene": the median gene of A3 is summarised by components fitted to its ICC curve
(ICC_inf = 0.64, ICC_500 = 0.54): line 1, batch 1/0.64 - 1, cell noise from the 500-cell point.

Outputs: results/design_guide_v2.json; paper/figures/fig4_design.svg
Figure: donors needed (80% power, alpha 0.05, 6 covariates) for true r = 0.3 vs runs per line,
one line per phenotype, 500 DA cells per run (P1: 3,000 cells), lambda = 0.64. Log y axis.
Palette: reference categorical slots 1-5, validated (CVD adjacent dE >= 9.1; normal >= 19.6);
contrast WARN for slots 3-5 -> every line is direct-labelled and has its own marker shape.
"""
import json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, "D:/cell2clinic/src")
from cell2clinic.design import reliability_from_components, n_donors  # noqa: E402

LAMBDA = {"UPDRS III (OFF)": 0.64, "Putamen SBR": 0.54}
N_FOUNDIN = 47
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
SERIES = [("P2da", "SNCA in DA neurons", "#2a78d6", "o"),
          ("P4da", "Lysosome score in DA neurons", "#eb6834", "D"),
          ("P1", "DA-neuron fraction", "#1baf7a", "s"),
          ("gene", "Typical single gene, DA neurons", "#eda100", "^"),
          ("P3da", "OXPHOS score in DA neurons", "#e87ba4", "v")]


def components():
    a3 = json.load(open(ROOT / "data/derived/jerber/a3_noise.json"))
    comp = {k: (a3[k]["sigma_line2"], a3[k]["sigma_batch2"], a3[k]["c_per_cell"]) for k in ("P1", "P2da", "P3da", "P4da")}
    g = a3["genes_random1000"]["icc_by_cells_median"]
    i_inf, i_500 = g["inf"], g["500"]
    b = 1 / i_inf - 1
    comp["gene"] = (1.0, b, 500 * (1 / i_500 - 1 - b))
    return comp, a3


def cells_for(k, m):
    return 3000.0 if k == "P1" and np.isfinite(m) else m


def table(comp):
    out, seen = [], set()
    for k, (L, B, C) in comp.items():
        for oc, lam in LAMBDA.items():
            for k_runs in (1, 2, 4):
                for m in (100, 500, np.inf):
                    if (k, oc, k_runs, cells_for(k, m)) in seen:
                        continue
                    seen.add((k, oc, k_runs, cells_for(k, m)))
                    rel = reliability_from_components(L, B, C, k_runs, cells_for(k, m))
                    row = {"phenotype": k, "outcome": oc, "runs": k_runs, "cells": "inf" if not np.isfinite(m) else int(cells_for(k, m)),
                           "reliability": round(float(rel), 3)}
                    for r in (0.2, 0.3, 0.4):
                        n = float(n_donors(r * np.sqrt(rel * lam)))
                        row[f"N_r{r}"] = None if not np.isfinite(n) else int(np.ceil(n))
                    out.append(row)
    return out


def figure(comp):
    """Hand-built SVG (matplotlib is blocked by Windows App Control on this machine)."""
    W, H = 760, 470
    x0, x1, y0, y1 = 78, 500, 62, 372
    runs = np.arange(1, 7)
    curves = []
    for k, label, col, mk in SERIES:
        L, B, C = comp[k]
        N = [float(n_donors(0.3 * np.sqrt(reliability_from_components(L, B, C, r, cells_for(k, 500)) * 0.64))) for r in runs]
        curves.append((label, col, mk, N))
    lo, hi = 30, 10 ** np.ceil(np.log10(max(max(c[3]) for c in curves)))
    X = lambda r: x0 + (r - 1) / 5 * (x1 - x0)
    Y = lambda n: y1 - (np.log10(n) - np.log10(lo)) / (np.log10(hi) - np.log10(lo)) * (y1 - y0)
    e = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
         f'font-family="DejaVu Sans, Helvetica, Arial, sans-serif">', f'<rect width="{W}" height="{H}" fill="{SURFACE}"/>',
         f'<text x="{x0 - 60}" y="26" font-size="14" fill="{INK}">Repeating differentiations cuts the donors needed, '
         f'but every phenotype still needs 150+ donors</text>']
    for t in [30, 50, 100, 200, 500, 1000, 2000, 5000, 10000]:
        if lo <= t <= hi:
            e.append(f'<line x1="{x0}" x2="{x1}" y1="{Y(t):.1f}" y2="{Y(t):.1f}" stroke="{GRID}" stroke-width="1"/>')
            e.append(f'<text x="{x0 - 8}" y="{Y(t) + 4:.1f}" font-size="11" fill="{INK2}" text-anchor="end">{t:,}</text>')
    for r in runs:
        e.append(f'<text x="{X(r):.1f}" y="{y1 + 18}" font-size="11" fill="{INK2}" text-anchor="middle">{r}</text>')
    e.append(f'<line x1="{x0}" x2="{x1}" y1="{y1}" y2="{y1}" stroke="{MUTED}" stroke-width="1"/>')
    e.append(f'<text x="{(x0 + x1) / 2}" y="{y1 + 40}" font-size="11.5" fill="{INK2}" text-anchor="middle">'
             'Independent differentiation runs per line (each in its own batch)</text>')
    e.append(f'<text transform="translate({x0 - 52},{(y0 + y1) / 2}) rotate(-90)" font-size="11.5" fill="{INK2}" '
             'text-anchor="middle">Donors needed (80% power, true r = 0.3)</text>')
    yf = Y(N_FOUNDIN)
    e.append(f'<line x1="{x0}" x2="{x1}" y1="{yf:.1f}" y2="{yf:.1f}" stroke="{INK2}" stroke-width="1.2" stroke-dasharray="5 4"/>')
    e.append(f'<text x="{x0 + 4}" y="{yf + 15:.1f}" font-size="10.5" fill="{INK2}">FOUNDIN-PD with PPMI follow-up: {N_FOUNDIN} PD donors</text>')

    def marker(mk, x, y, col):
        r = 5
        pts = {"D": [(x, y - r - 1), (x + r + 1, y), (x, y + r + 1), (x - r - 1, y)],
               "^": [(x, y - r - 1), (x + r + 1, y + r), (x - r - 1, y + r)],
               "v": [(x, y + r + 1), (x + r + 1, y - r), (x - r - 1, y - r)]}
        st = f'fill="{col}" stroke="{SURFACE}" stroke-width="2"'
        if mk == "o":
            return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" {st}/>'
        if mk == "s":
            return f'<rect x="{x - r:.1f}" y="{y - r:.1f}" width="{2 * r}" height="{2 * r}" rx="1.5" {st}/>'
        return '<polygon points="' + " ".join(f"{a:.1f},{b:.1f}" for a, b in pts[mk]) + f'" {st}/>'

    labels = []
    for label, col, mk, N in curves:
        d = " ".join(f"{'M' if i == 0 else 'L'}{X(r):.1f},{Y(n):.1f}" for i, (r, n) in enumerate(zip(runs, N)))
        e.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="2"/>')
        e += [marker(mk, X(r), Y(n), col) for r, n in zip(runs, N)]
        labels.append([Y(N[-1]), label, N[0], N[-1], col, mk])
    labels.sort()
    ys = [l[0] for l in labels]
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + 30)
    for (yy, label, n1, n6, col, mk), ly in zip(labels, ys):
        e.append(f'<line x1="{x1 + 6}" x2="{x1 + 22}" y1="{yy:.1f}" y2="{ly:.1f}" stroke="{MUTED}" stroke-width="0.8"/>')
        e.append(marker(mk, x1 + 32, ly, col))
        e.append(f'<text x="{x1 + 44}" y="{ly - 2:.1f}" font-size="11" fill="{INK}">{label}</text>')
        e.append(f'<text x="{x1 + 44}" y="{ly + 12:.1f}" font-size="10" fill="{INK2}">{n1:,.0f} donors (1 run) → {n6:,.0f} (6 runs)</text>')
    e.append(f'<text x="{x0 - 60}" y="{H - 30}" font-size="9.5" fill="{MUTED}">Across-batch variance components from Jerber et al. 2021 '
             '(199 lines, 16 pools; A3). Slope reliability λ = 0.64 (OFF-state UPDRS III, PPMI; E4).</text>')
    e.append(f'<text x="{x0 - 60}" y="{H - 15}" font-size="9.5" fill="{MUTED}">500 DA cells per run (DA fraction: 3,000 cells). '
             'α = 0.05 two-sided, 6 covariates. Log scale.</text>')
    e.append("</svg>")
    out = ROOT / "paper/figures/fig4_design.svg"
    out.write_text(chr(10).join(e), encoding="utf-8")


def main():
    comp, a3 = components()
    res = {"components_line_batch_cell": {k: [float(x) for x in v] for k, v in comp.items()},
           "lambda": LAMBDA, "table": table(comp), "source": "jerber-a3-v0.1.0 + E4 lambda"}
    json.dump(res, open(ROOT / "results/design_guide_v2.json", "w"), indent=1)
    figure(comp)
    for row in res["table"]:
        if row["outcome"].startswith("UPDRS") and row["cells"] in (500, 3000):
            print(row)


if __name__ == "__main__":
    main()
