"""Matplotlib style for the NHFOAM figures: all text typeset by LaTeX
(Computer Modern), so the figures match a LaTeX paper.

Needs latex, dvipng and the cm-super fonts (Ubuntu: texlive-latex-extra,
dvipng, cm-super). Without them the figures fall back to matplotlib's own
Computer Modern fonts, so the scripts still run.
"""
import shutil
import subprocess

import matplotlib as mpl


def _have_latex():
    if not (shutil.which("latex") and shutil.which("dvipng") and shutil.which("kpsewhich")):
        return False
    try:
        out = subprocess.run(["kpsewhich", "type1ec.sty"], capture_output=True, text=True)
        return bool(out.stdout.strip())
    except OSError:
        return False


def use_latex(fontsize=10):
    """Set the rcParams; returns True if real LaTeX is used."""
    usetex = _have_latex()
    mpl.rcParams.update({
        "text.usetex": usetex,
        "text.latex.preamble": r"\usepackage{amsmath}",
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"] if usetex else ["cmr10"],
        "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True,
        "font.size": fontsize,
        "axes.titlesize": fontsize,
        "axes.labelsize": fontsize,
        "legend.fontsize": fontsize - 1,
        "xtick.labelsize": fontsize - 1,
        "ytick.labelsize": fontsize - 1,
        "figure.titlesize": fontsize + 1,
    })
    return usetex


def save(fig, png, dpi=150):
    """Save a PNG (for the README) and a vector PDF (for papers)."""
    fig.savefig(png, dpi=dpi)
    fig.savefig(png.rsplit(".", 1)[0] + ".pdf")
