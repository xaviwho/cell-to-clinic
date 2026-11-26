"""Fig. 1: system model of the Cell-to-Clinic framework (architecture only, no results).
Hand-built SVG (matplotlib is blocked by Windows App Control); convert with svg2pdf.py.

Layout (U-shaped pipeline): top tiers = cellular and clinical branches, which never meet
phenotypes and outcomes; a horizontal "registration lock" line; bottom tier = linkage ->
decision -> reported result, flowing right to left. Every block cites the Section II equations
it implements. Equation numbers are read from the LaTeX .aux (--aux), so the figure cannot drift
from the text; without --aux the script refuses to run.

Sized for a two-column IEEE figure* (7.16 in): canvas 1080 px wide; at that scale 17 px text
prints at about 8 pt. Grayscale-safe: state is carried by line style (solid flow, dashed
gate-failure / exploratory paths), not colour alone.
Usage: python -I scripts/figures/fig1_system.py --aux <build>/main.aux
"""
import argparse, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "paper" / "figures" / "fig1_system.svg"
W, H = 1080, 600
FS, FT = 17, 18                     # body / title font size (px)
INK, INK2, LINE = "#111111", "#444444", "#222222"
FILL = {"input": "#ffffff", "proc": "#f1f1f1", "gate": "#e6eef8", "decide": "#f1f1f1", "out": "#ffffff", "sink": "#ffffff"}
EDGE = {"gate": "#2b5797"}
COLS = [20, 265, 510, 755]          # x of the four columns
CW = 215                            # column width
TA, TB, TC = (40, 240), (270, 375), (425, 575)   # tier y-ranges
LOCK_Y = 400


def eqnums(aux):
    lab = dict(re.findall(r"\\newlabel\{eq:([^}]+)\}\{\{(\d+)\}", Path(aux).read_text(encoding="utf-8", errors="ignore")))
    need = ["linkset", "qc", "score", "dalabel", "calib", "pseudobulk", "geneset", "resid", "enet", "icc", "icc_ci",
            "a01", "eqtl", "pi1", "lmm", "reml", "satt", "partialr", "tost", "bh", "class", "a05", "within"]
    miss = [k for k in need if k not in lab]
    if miss:
        raise SystemExit(f"labels missing from aux: {miss}")
    return lab


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def block(x, y0, y1, title, lines, kind, dashed=False, w=CW):
    e = EDGE.get(kind, LINE)
    dash = ' stroke-dasharray="6 4"' if dashed else ""
    out = [f'<rect x="{x}" y="{y0}" width="{w}" height="{y1 - y0}" rx="6" fill="{FILL[kind]}" stroke="{e}" stroke-width="1.6"{dash}/>']
    out.append(f'<text x="{x + w / 2}" y="{y0 + 26}" font-size="{FT}" font-weight="bold" fill="{INK}" text-anchor="middle">{esc(title)}</text>')
    for i, s in enumerate(lines):
        out.append(f'<text x="{x + w / 2}" y="{y0 + 52 + 22 * i}" font-size="{FS}" fill="{INK2}" text-anchor="middle">{esc(s)}</text>')
    return out


def path(pts, dashed=False, head=True):
    d = "M " + " L ".join(f"{x} {y}" for x, y in pts)
    dash = ' stroke-dasharray="6 4"' if dashed else ""
    mk = ' marker-end="url(#ah)"' if head else ""
    return f'<path d="{d}" fill="none" stroke="{LINE}" stroke-width="1.6"{dash}{mk}/>'


def label(x, y, s, anchor="start", size=16, style=""):
    return f'<text x="{x}" y="{y}" font-size="{size}" fill="{INK}" text-anchor="{anchor}"{style}>{esc(s)}</text>'


def build(q):
    r = lambda *ks: "(" + ")–(".join(q[k] for k in ks) + ")" if len(ks) > 1 else f"({q[ks[0]]})"
    e = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
         f'font-family="Arial, Helvetica, sans-serif">',
         '<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
         f'<path d="M 0 0 L 10 5 L 0 10 z" fill="{LINE}"/></marker></defs>',
         f'<rect width="{W}" height="{H}" fill="#ffffff"/>']
    c0, c1, c2, c3 = COLS
    # tier A: cellular branch
    e += block(c0, *TA, "FOUNDIN-PD", ["iPSC-derived neurons", "day-65 scRNA-seq", "bulk RNA, day 0/65", "SNP genotypes"], "input")
    e += block(c1, *TA, "Single-cell processing", [f"QC filter {r('qc')}", f"module scores {r('score')}",
                                                  f"DA-neuron label {r('dalabel')}", f"threshold t* {r('calib')}"], "proc")
    e += block(c2, *TA, "Phenotype construction", [f"pseudobulk {r('pseudobulk')}", f"gene-set score {r('geneset')}",
                                                  f"composition adj. {r('resid')}", "z-score per model"], "proc")
    e += block(c3, *TA, "Pre-link gates", ["A0-6 reproduction", f"A0-4 pos. control {r('enet')}",
                                          f"A0-2, G-icc {r('icc', 'icc_ci')}", "A0-3, G-cells", f"G-eQTL {r('eqtl', 'pi1')}"], "gate")
    # tier B: clinical branch + sink
    e += block(c0, *TB, "PPMI", ["MDS-UPDRS III", "DaTscan, SAA, NSD-ISS"], "input")
    e += block(c1, *TB, "Outcome derivation", ["OFF state, time since dx", "SBR, SAA, NSD-ISS stage"], "proc")
    e += block(c2, *TB, f"A0-1 clinical gate {r('a01')}", ["≥ 30 donors with ≥ 3", "OFF visits over ≥ 2 y"], "gate")
    e += block(c3, *TB, "Not tested", ["phenotype removed;", "estimate descriptive"], "sink", dashed=True)
    # tier C: linkage (right to left)
    e += block(c3, *TC, "Linkage models", [f"LMM, REML {r('lmm', 'reml')}", f"Satterthwaite df {r('satt')}",
                                          "OLS for H3", f"A0-5 coherence {r('a05')}"], "proc")
    e += block(c2, *TC, "Decision rule", [f"partial r {r('partialr')}", f"TOST, Δ = 0.40 {r('tost')}",
                                         f"BH q-values {r('bh')}", f"class {r('class')}"], "decide")
    e += block(c1, *TC, "Reported result", ["supported", "equivalent", "inconclusive"], "out")
    e += block(c0, *TC, "Exploratory", ["post hoc, labelled:", f"batch decomposition {r('within')}"], "sink", dashed=True)
    # flow arrows, tier A
    ya = (TA[0] + TA[1]) / 2
    for a, b in ((c0, c1), (c1, c2), (c2, c3)):
        e.append(path([(a + CW, ya), (b - 2, ya)]))
    yb = (TB[0] + TB[1]) / 2
    for a, b in ((c0, c1), (c1, c2)):
        e.append(path([(a + CW, yb), (b - 2, yb)]))
    # linkable set: gate block -> right gutter -> down -> linkage block
    gx = c3 + CW + 40
    yc = (TC[0] + TC[1]) / 2
    e.append(path([(c3 + CW, ya), (gx, ya), (gx, yc), (c3 + CW + 2, yc)]))
    e.append(label(gx + 6, ya + 60, "linkable", size=16))
    e.append(label(gx + 6, ya + 80, f"set {r('linkset')}", size=16))
    # clinical gate -> linkage block (through the column gap)
    gap = c3 - 15
    e.append(path([(c2 + CW, yb), (gap, yb), (gap, TC[0] + 30), (c3 - 2, TC[0] + 30)]))
    # gate failures -> not tested (dashed)
    e.append(path([(c3 + 60, TA[1]), (c3 + 60, TB[0] - 2)], dashed=True))
    e.append(label(c3 + 68, TA[1] + 21, "fail", size=16))
    e.append(path([(c3 + 160, TC[0]), (c3 + 160, TB[1] + 2)], dashed=True))
    e.append(label(c3 + 168, TC[0] - 6, "A0-5 fail", size=16))
    # tier C flow, right to left
    for a, b in ((c3, c2), (c2, c1)):
        e.append(path([(a, yc), (b + CW + 2, yc)]))
    e.append(path([(c1, yc), (c0 + CW + 2, yc)], dashed=True))
    # registration lock
    e.append(f'<line x1="10" x2="{W - 10}" y1="{LOCK_Y}" y2="{LOCK_Y}" stroke="{INK}" stroke-width="2" stroke-dasharray="2 4"/>')
    txt = "Registration lock: rules, models and code frozen (OSF, git tag) before linkage"
    e.append(f'<rect x="24" y="{LOCK_Y - 11}" width="{7.7 * len(txt)}" height="22" fill="#ffffff"/>')
    e.append(label(30, LOCK_Y + 6, txt, size=16, style=' font-style="italic"'))
    e.append("</svg>")
    return "\n".join(e)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aux", required=True)
    a = ap.parse_args()
    OUT.write_text(build(eqnums(a.aux)), encoding="utf-8")
    print("written", OUT)


if __name__ == "__main__":
    main()
