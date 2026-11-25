"""Convert hand-built SVG figures to PDF for LaTeX with headless Microsoft Edge (vector output).
matplotlib/cairo are blocked by Windows App Control here; Edge's print-to-PDF is allowed.
Each SVG is wrapped in an HTML page whose @page size equals the SVG's width/height (px), margin 0.
Usage: python -I scripts/figures/svg2pdf.py paper/figures/fig4_design.svg [...]
"""
import re, subprocess, sys, tempfile
from pathlib import Path

EDGE = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"


def convert(svg):
    svg = Path(svg).resolve()
    txt = svg.read_text(encoding="utf-8")
    w, h = (float(re.search(fr'{k}="([\d.]+)"', txt).group(1)) for k in ("width", "height"))
    html = (f"<!doctype html><html><head><meta charset='utf-8'><style>@page {{ size: {w}px {h}px; margin: 0 }}"
            f" html,body {{ margin:0; padding:0 }} svg {{ display:block }}</style></head><body>{txt}</body></html>")
    with tempfile.TemporaryDirectory() as tmp:
        page = Path(tmp) / "page.html"
        page.write_text(html, encoding="utf-8")
        pdf = svg.with_suffix(".pdf")
        subprocess.run([EDGE, "--headless", "--disable-gpu", "--no-pdf-header-footer",
                        f"--print-to-pdf={pdf}", page.as_uri()], check=True, capture_output=True)
    head = pdf.read_bytes()[:5]
    if head != b"%PDF-":
        raise SystemExit(f"not a PDF: {pdf}")
    print("written", pdf, pdf.stat().st_size, "bytes")


if __name__ == "__main__":
    for f in sys.argv[1:]:
        convert(f)
