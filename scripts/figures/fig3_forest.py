"""Fig. 3: equivalence forest plot of the 10 confirmatory tests (+ 2 reported-only H1 rows).
Reads aggregate results only (data/derived/models.json). Writes paper/figures/fig3_forest.{png,pdf}.
Palette: reference categorical slots 1-2 (validated: CVD dE 24.7, normal dE 33.6); grey = not tested.
Class is encoded twice (colour + marker shape) and labelled in text on each row.
"""
import json, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "data" / "derived" / "models.json"
OUT = ROOT / "paper" / "figures"
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
STYLE = {"equivalent": ("#2a78d6", "o", "Equivalent"),
         "inconclusive": ("#eb6834", "D", "Inconclusive"),
         "supported": ("#0b0b0b", "s", "Supported"),
         "dropped_by_A0-5": (MUTED, "o", "Not tested (dropped by A0-5)")}
NAMES = {"P2": "SNCA (P2)", "P3": "OXPHOS (P3)", "P4": "Lysosome (P4)"}
FAMILIES = [("H1", "H1  Motor progression\n(OFF-state UPDRS III slope)"),
            ("H2", "H2  DaTscan putamen decline"),
            ("H3_SAA", "H3  CSF α-syn SAA status"),
            ("H3_NSDISS", "H3  NSD-ISS stage")]
BOUND = 0.40


def main():
    r = json.load(open(SRC))
    rows = []
    for fam, title in FAMILIES:
        rows.append(("header", title, None))
        for ph in ("P2", "P3", "P4"):
            v = r[fam].get(ph)
            if v and "r" in v:
                rows.append(("row", ph, v))
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.edgecolor": MUTED,
                         "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK})
    fig, ax = plt.subplots(figsize=(7.6, 5.2), dpi=300)
    fig.patch.set_facecolor(SURFACE); ax.set_facecolor(SURFACE)
    ax.axvspan(-BOUND, BOUND, color="#f0efec", zorder=0, lw=0, label="Equivalence region (±0.40)")
    for x in (-BOUND, BOUND):
        ax.axvline(x, color=MUTED, lw=0.8, ls=(0, (3, 3)), zorder=1)
    ax.axvline(0, color=INK2, lw=0.8, zorder=1)
    y = 0; yt, yl = [], []
    seen = set()
    for kind, name, v in rows:
        y -= 1
        if kind == "header":
            ax.text(-0.98, y, name, ha="left", va="center", fontsize=8.5, fontweight="bold", color=INK)
            continue
        cls = v["class"]
        col, mk, lab = STYLE[cls]
        lo, hi = v["ci90"]
        filled = cls != "dropped_by_A0-5"
        ax.plot([lo, hi], [y, y], color=col, lw=2, solid_capstyle="round", zorder=2)
        ax.plot(v["r"], y, mk, ms=6.5, mfc=col if filled else SURFACE, mec=col, mew=1.6, zorder=3,
                label=lab if lab not in seen else None)
        seen.add(lab)
        q = v.get("q")
        right = (f"{lab.split(' (')[0].lower()}, q = {q:.2f}" if q is not None else "not tested (A0-5)")
        ax.text(1.03, y, f"n = {v['n_donors']}", ha="left", va="center", fontsize=7.5, color=INK2,
                transform=ax.get_yaxis_transform(), clip_on=False)
        ax.text(1.13, y, right, ha="left", va="center", fontsize=7.5, color=INK2,
                transform=ax.get_yaxis_transform(), clip_on=False)
        yt.append(y); yl.append("    " + NAMES[name])
    ax.set_yticks(yt); ax.set_yticklabels(yl)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(y - 0.7, 0.3)
    ax.set_xlim(-1.0, 1.0)
    ax.set_xticks([-0.8, -0.4, 0, 0.4, 0.8])
    ax.set_xlabel("Partial correlation r with 90% CI (phenotype per SD)")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.grid(False)
    h, l = ax.get_legend_handles_labels()
    order = ["Equivalent", "Inconclusive", "Supported", "Not tested (dropped by A0-5)", "Equivalence region (±0.40)"]
    hl = sorted([x for x in zip(h, l)], key=lambda x: order.index(x[1]))
    leg = ax.legend([x[0] for x in hl], [x[1] for x in hl], loc="lower center", bbox_to_anchor=(0.62, -0.2),
                    ncol=4, frameon=False, fontsize=7.5, handletextpad=0.4, columnspacing=1.4)
    for t in leg.get_texts():
        t.set_color(INK)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig3_forest.{ext}", facecolor=SURFACE, bbox_inches="tight")
    print("written fig3_forest.png/.pdf")


if __name__ == "__main__":
    sys.exit(main())
