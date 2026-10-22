"""Build the typeset PDF report from the single Markdown source.

The same source feeds the live document (display formulas as ```latex blocks,
figures as FIG:<key> placeholders). Here the formulas become display math, the
placeholders become the generated figures, and pandoc + XeLaTeX typeset it.

    python docs/report/build_pdf.py            # -> docs/report/report.pdf
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

FIGURES = {
    "attitude": "outputs/sensitivity/figures/attitude_characterisation.png",
    "gantt": "outputs/part1/figures/gantt_docking_windows.pdf",
    "anatomy": "outputs/part1/figures/near_miss_pass_anatomy.png",
    "geometry": "outputs/sensitivity/figures/geometry_why_only_wellington.png",
    "relaxation": "outputs/sensitivity/figures/relaxation_curves.png",
    "mapzoom": "outputs/sensitivity/figures/feasibility_map_ltan_zoom.png",
    "trade": "outputs/part2/figures/strategy_trade_sun.png",
    "trajectory": "outputs/part2/figures/approach_trajectory.png",
    "timeline": "outputs/part2/figures/part2_timeline.png",
    "tltrade": "outputs/part2/figures/timeline_trade.png",
    "mc": "outputs/part2/figures/monte_carlo.png",
}

FRONT = """---
title: "Docking-Window Analysis and RPOD Final-Approach Design"
date: "October 2026"
documentclass: article
classoption: [a4paper]
fontsize: 10pt
geometry: [top=2.1cm, bottom=2.3cm, left=2.2cm, right=2.2cm]
fontfamily: fontspec
mainfont: TeX Gyre Pagella
sansfont: Inter
monofont: DejaVu Sans Mono
mathfont: Latin Modern Math
linestretch: 1.08
colorlinks: true
linkcolor: "accent"
urlcolor: "accent"
header-includes: |
  ```{=latex}
  \\usepackage{xcolor}
  \\definecolor{accent}{HTML}{1F5FAD}
  \\definecolor{ink}{HTML}{111111}
  \\definecolor{muted}{HTML}{5A5A55}
  \\usepackage{etoolbox}
  \\usepackage{newunicodechar}
  \\newunicodechar{¹}{\\textsuperscript{1}}
  \\newunicodechar{²}{\\textsuperscript{2}}
  \\newunicodechar{³}{\\textsuperscript{3}}
  \\newunicodechar{⁴}{\\textsuperscript{4}}
  \\newunicodechar{⁵}{\\textsuperscript{5}}
  \\newunicodechar{⁷}{\\textsuperscript{7}}
  \\newunicodechar{⁺}{\\textsuperscript{+}}
  \\newunicodechar{⁻}{\\textsuperscript{−}}
  \\newunicodechar{₀}{\\textsubscript{0}}
  \\newunicodechar{₁}{\\textsubscript{1}}
  \\newunicodechar{₂}{\\textsubscript{2}}
  \\newunicodechar{₃}{\\textsubscript{3}}
  \\newunicodechar{₄}{\\textsubscript{4}}
  \\newfontfamily\\headingfont{Inter}[UprightFont=* SemiBold, BoldFont=* Bold]
  \\usepackage[nobottomtitles*]{titlesec}
  \\renewcommand{\\bottomtitlespace}{0.12\\textheight}
  \\titleformat{\\section}{\\headingfont\\large\\color{ink}}{}{0pt}{}
  \\titlespacing*{\\section}{0pt}{2.0ex plus 0.6ex minus 0.3ex}{0.9ex}
  \\titleformat{\\subsection}{\\headingfont\\normalsize\\color{ink}}{}{0pt}{}
  \\usepackage{caption}
  \\captionsetup{font={small,color=muted}, labelfont={bf,color=ink}, labelsep=period, skip=4pt}
  \\usepackage{float}
  \\floatplacement{figure}{H}
  \\AtBeginEnvironment{longtable}{\\small\\renewcommand{\\arraystretch}{1.18}}
  \\setlength{\\LTpre}{6pt}
  \\setlength{\\LTpost}{8pt}
  \\usepackage{enumitem}
  \\setlist[itemize]{leftmargin=1.3em, itemsep=2pt, topsep=3pt}
  \\usepackage{fancyhdr}
  \\pagestyle{fancy}
  \\fancyhf{}
  \\renewcommand{\\headrulewidth}{0pt}
  \\fancyfoot[C]{\\color{muted}\\footnotesize\\thepage}
  \\fancypagestyle{plain}{\\fancyhf{}\\renewcommand{\\headrulewidth}{0pt}\\fancyfoot[C]{\\color{muted}\\footnotesize\\thepage}}
  \\usepackage{titling}
  \\setlength{\\droptitle}{-3.5em}
  \\pretitle{\\begin{flushleft}\\headingfont\\LARGE\\color{ink}}
  \\posttitle{\\par\\end{flushleft}\\vspace{-0.6em}}
  \\preauthor{}
  \\postauthor{}
  \\predate{\\begin{flushleft}\\vspace{0.4em}\\color{muted}\\small}
  \\postdate{\\par\\end{flushleft}\\vspace{-0.4em}\\noindent\\color{muted}\\rule{\\linewidth}{0.4pt}\\vspace{0.2em}}
  \\renewcommand{\\maketitlehookb}{\\begin{flushleft}\\color{muted}\\large\\headingfont\\mdseries Rendezvous with a 550 km sun-synchronous target, 10–20 July 2026\\end{flushleft}\\vspace{-0.8em}}
  \\setlength{\\parskip}{0.45em}
  \\setlength{\\parindent}{0pt}
  ```
---

"""


def to_pandoc_markdown(src: str) -> str:
    def latex_block(m):
        body = m.group(1).strip()
        return "$$\n" + body + "\n$$"
    out = re.sub(r"```latex\n(.*?)```", latex_block, src, flags=re.S)

    def fig(m):
        key = m.group(1)
        if key not in FIGURES:
            raise KeyError(f"unknown figure placeholder FIG:{key}")
        p = ROOT / FIGURES[key]
        if not p.exists():
            raise FileNotFoundError(p)
        return f"({p.as_posix()})"
    out = re.sub(r"\(FIG:([a-z0-9_]+)\)", fig, out)
    return FRONT + out


def main():
    src = (HERE / "report_source.md").read_text(encoding="utf-8")
    md = to_pandoc_markdown(src)
    build = HERE / "build"
    build.mkdir(exist_ok=True)
    (build / "report.md").write_text(md, encoding="utf-8")
    pdf = HERE / "report.pdf"
    cmd = ["pandoc", str(build / "report.md"), "-o", str(pdf), "--pdf-engine=xelatex",
           "--shift-heading-level-by=-1", "-f", "markdown",
           f"--resource-path={ROOT}"]
    tex = subprocess.run(cmd[:3] + [str(build / "report.tex")] + cmd[4:] + ["-s"], capture_output=True, text=True)
    if tex.returncode:
        sys.stderr.write(tex.stderr)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        sys.stderr.write(r.stderr)
        return r.returncode
    if r.stderr.strip():
        sys.stderr.write(r.stderr)
    print(pdf)
    return 0


if __name__ == "__main__":
    if not shutil.which("pandoc"):
        sys.exit("pandoc not found")
    sys.exit(main())
