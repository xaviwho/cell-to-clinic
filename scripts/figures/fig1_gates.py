"""Fig. 1: the gated framework. (A) flow of phenotypes through the gates to the hypotheses;
(B) each gate's measured value against its pass threshold.
Reads aggregate outputs only. Writes paper/figures/fig1_gates.{png,pdf}.
Gate outcome is a state, so it uses the reference status palette (good / warning / critical),
always paired with a marker shape and a text label (never colour alone).
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "data" / "derived"
OUT = ROOT / "paper" / "figures"
INK, INK2, MUTED, SURFACE, BOX = "#0b0b0b", "#52514e", "#8a8984", "#fcfcfb", "#f0efec"
STATE = {"pass": ("#0ca30c", "o", "pass"), "partial": ("#fab219", "D", "partial"), "fail": ("#d03b3b", "X", "fail")}


def load():
    scrn = json.load(open(D / "scrn" / "gate.json"))
    g1 = json.load(open(D / "day0_gate.json")); g2 = json.load(open(D / "day0_gate_v2.json"))
    pre = json.load(open(D / "phenotypes" / "gates_prelink.json")); mod = json.load(open(D / "models.json"))
    a01 = mod["A0_1"]
    rows = [
        ("A0-6  Reproduction: r with published %TH+", scrn["gate_primary"]["pearson_r"], 0.80, scrn["gate_verdict"]),
        ("A0-4  Positive control v0.1.0: CV AUC", g1["primary"]["auc_mean"], 0.75, g1["verdict"]),
        ("A0-4  Positive control v0.2.0: CV AUC", g2["primary_published_outcome"]["cv_auc_mean"], 0.75, g2["verdict"]),
    ]
    names = {"P1": "DA-neuron fraction", "P2": "SNCA", "P3": "OXPHOS", "P4": "Lysosome"}
    for p in ("P1", "P2", "P3", "P4"):
        v = pre["A0_2"][p]
        rows.append((f"A0-2  Stability ICC: {names[p]} ({p})", v["icc"], pre["A0_2_threshold"], "pass" if v["pass"] else "fail"))
    for p in ("P3", "P4"):
        v = mod["A0_5"][p]
        rows.append((f"A0-5  Coherence, genes same sign: {names[p]}", v["share_same_sign"], 0.80, "pass" if v["pass"] else "fail"))
    rows.append((f"A0-1  Linkage: PD donors meeting rule ({a01['n_meeting_rule']}/{a01['n_pd_donors']})",
                 a01["n_meeting_rule"] / a01["n_pd_donors"], 30 / a01["n_pd_donors"], "pass" if a01["pass"] else "fail"))
    return rows, pre, mod


def box(ax, x, y, w, h, text, weight="normal", edge=MUTED, fs=7.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.015", fc=BOX, ec=edge, lw=0.8))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=INK, fontweight=weight, linespacing=1.3)


def arrow(ax, x0, y0, x1, y1, text=None):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.9, mutation_scale=8))
    if text:
        ax.text((x0 + x1) / 2 + 0.015, (y0 + y1) / 2, text, ha="left", va="center", fontsize=6.8, color="#d03b3b")


def main():
    rows, pre, mod = load()
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
    fig = plt.figure(figsize=(7.4, 7.6), dpi=300, facecolor=SURFACE)
    axA = fig.add_axes([0.02, 0.50, 0.96, 0.47]); axA.set_axis_off(); axA.set_xlim(0, 1); axA.set_ylim(0, 1)
    axA.text(0.0, 1.0, "A", fontsize=11, fontweight="bold", va="top")
    # pipeline-validity gates (left column)
    axA.text(0.04, 0.95, "Pipeline validity (before linkage)", fontsize=7.8, color=INK2, va="top")
    box(axA, 0.04, 0.74, 0.28, 0.16, "A0-6 Reproduction\nr = 0.999 vs published %TH+\nPASS")
    box(axA, 0.04, 0.52, 0.28, 0.17, "A0-4 Positive control\nv0.1.0 AUC 0.57: FAIL\nv0.2.0 AUC 0.72: PARTIAL")
    box(axA, 0.04, 0.30, 0.28, 0.17, "A0-3 Cell sufficiency\n78 / 79 lines kept\n(no PD donor lost)")
    box(axA, 0.04, 0.08, 0.28, 0.17, "A0-1 Clinical linkage\n44 / 47 PD donors\nPASS (≥ 30)")
    # phenotype flow (right)
    axA.text(0.40, 0.95, "Phenotype flow (47 PD donors)", fontsize=7.8, color=INK2, va="top")
    box(axA, 0.40, 0.78, 0.56, 0.12, "4 registered phenotypes\nP1 DA fraction · P2 SNCA · P3 OXPHOS · P4 Lysosome", weight="bold")
    arrow(axA, 0.68, 0.78, 0.68, 0.67, "P1 dropped (ICC 0.33 < 0.5)")
    box(axA, 0.40, 0.55, 0.56, 0.12, "A0-2 Technical stability (ICC ≥ 0.5)\n3 phenotypes kept: P2, P3, P4")
    arrow(axA, 0.52, 0.55, 0.52, 0.42)
    arrow(axA, 0.84, 0.55, 0.84, 0.42, None)
    box(axA, 0.40, 0.27, 0.25, 0.15, "H2 DaTscan, H3 SAA,\nH3 NSD-ISS\nP2, P3, P4 tested")
    box(axA, 0.71, 0.27, 0.25, 0.15, "A0-5 Composite coherence\n(≥ 80% genes same sign)\nP3 75%, P4 62%")
    arrow(axA, 0.835, 0.27, 0.835, 0.17, "P3, P4 dropped")
    box(axA, 0.71, 0.04, 0.25, 0.13, "H1 Motor progression\nP2 tested")
    box(axA, 0.40, 0.04, 0.25, 0.13, "Result: 10 tests\n0 supported · 6 equivalent\n4 inconclusive", weight="bold")
    arrow(axA, 0.525, 0.27, 0.525, 0.17)
    arrow(axA, 0.71, 0.105, 0.65, 0.105)

    axB = fig.add_axes([0.47, 0.06, 0.38, 0.38], facecolor=SURFACE)
    fig.text(0.02, 0.45, "B", fontsize=11, fontweight="bold", va="top")
    n = len(rows)
    seen = set()
    for i, (lab, val, thr, st) in enumerate(rows):
        y = n - i
        col, mk, txt = STATE[st]
        axB.plot([thr, thr], [y - 0.32, y + 0.32], color=INK, lw=1.4, solid_capstyle="butt", zorder=2)
        if lab.startswith("A0-4"):
            axB.plot([0.65, 0.65], [y - 0.22, y + 0.22], color=MUTED, lw=1.2, zorder=2)
        axB.plot([min(val, thr), max(val, thr)], [y, y], color=MUTED, lw=0.8, zorder=1)
        axB.plot(val, y, mk, ms=6.5, mfc=col, mec=INK if st == "partial" else col, mew=0.6, zorder=3,
                 label=txt if txt not in seen else None)
        seen.add(txt)
        axB.text(-0.03, y, lab, ha="right", va="center", fontsize=7, color=INK, transform=axB.get_yaxis_transform())
        axB.text(1.03, y, f"{val:.3f}  {txt}" if val > 0.99 else f"{val:.2f}  {txt}", ha="left", va="center", fontsize=7, color=INK2, transform=axB.get_yaxis_transform())
    axB.plot([], [], "|", color=INK, ms=8, mew=1.4, label="pass threshold")
    axB.plot([], [], "|", color=MUTED, ms=6, mew=1.2, label="partial threshold")
    axB.set_xlim(0, 1.02); axB.set_ylim(0.3, n + 0.7)
    axB.set_yticks([]); axB.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    axB.set_xlabel("Measured value (each row against its own threshold)", color=INK2)
    for s in ("top", "right", "left"):
        axB.spines[s].set_visible(False)
    axB.spines["bottom"].set_color(MUTED); axB.tick_params(colors=INK2)
    h, l = axB.get_legend_handles_labels()
    order = ["pass", "partial", "fail", "pass threshold", "partial threshold"]
    hl = sorted(zip(h, l), key=lambda x: order.index(x[1]))
    axB.legend([x[0] for x in hl], [x[1] for x in hl], loc="upper center", bbox_to_anchor=(0.3, -0.14), ncol=5,
               frameon=False, fontsize=7, handletextpad=0.3, columnspacing=1.0)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig1_gates.{ext}", facecolor=SURFACE, bbox_inches="tight")
    print("written fig1_gates.png/.pdf")


if __name__ == "__main__":
    main()
