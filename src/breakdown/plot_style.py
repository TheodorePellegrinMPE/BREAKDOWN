"""Shared look of all figures: sizes for journal columns, fonts, a colour blind safe palette and saving.

Use `with plotStyle():` around the code that draws a figure, so that nothing is changed for other matplotlib users.
Every figure is saved twice, as a vector pdf for the manuscript and as a 300 dpi png for slides.
"""
import os
import re
import zlib
from contextlib import contextmanager
import matplotlib.pyplot as plt

# Okabe and Ito colour blind safe palette, without yellow and black which do not work for thin lines
COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#8C6D31', '#666666']
NEUTRAL = '#B0B0B0'
# Line styles that stay distinguishable in grey scale
LINE_STYLES = ['-', '--', '-.', ':', (0, (5, 1)), (0, (3, 1, 1, 1))]
# Journal column widths in inches: one column, one and a half columns, two columns
WIDTHS = {'single': 3.4, 'medium': 5.2, 'double': 7.0}

STYLE = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
    'mathtext.fontset': 'dejavusans',
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7.5,
    'axes.linewidth': 0.8,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'lines.linewidth': 1.4,
    'legend.frameon': False,
    'axes.axisbelow': True,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.03,
}


# Only fonts and file settings, for the charts that keep their own layout and colours
TEXT_STYLE = {key: STYLE[key] for key in ('font.family', 'font.sans-serif', 'mathtext.fontset', 'pdf.fonttype',
                                          'ps.fonttype', 'savefig.dpi', 'savefig.bbox', 'savefig.pad_inches')}


@contextmanager
def textStyle():
    """Apply only the fonts and the file settings inside a with block."""
    with plt.rc_context(TEXT_STYLE):
        yield


def savePdfAndPng(folder, name, fig=None):
    """Save the current (or given) figure as pdf and png, close it, and return the path of the pdf."""
    return saveFigure(fig or plt.gcf(), folder, name)


@contextmanager
def plotStyle():
    """Apply the figure style inside a with block."""
    with plt.rc_context(STYLE):
        yield


def newFigure(width='single', height=None, **kwargs):
    """Make a figure with the width of a journal column.

    Args:
        width (str or float): 'single', 'medium', 'double' or a width in inches.
        height (float or None): Height in inches, 0.75 times the width if None.
        **kwargs: Passed to plt.subplots (nrows, ncols, sharex, ...).

    Returns:
        Tuple[Figure, Axes]: The figure and its axes.
    """
    inches = WIDTHS.get(width, width)
    return plt.subplots(figsize=(inches, height or 0.75 * inches), **kwargs)


def saveFigure(fig, folder, name):
    """Save a figure as pdf and png and close it.

    Args:
        fig (Figure): The figure.
        folder (str): Folder for the files, made if needed.
        name (str): File name without extension.

    Returns:
        str: Path of the pdf.
    """
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, f'{name}.pdf')
    fig.savefig(path)
    fig.savefig(os.path.join(folder, f'{name}.png'))
    plt.close(fig)
    return path


def rangeFormulaLabel(carbons, hydrogens, others=None):
    """Formula label for legends where the number of H atoms (and of other atoms) varies, like C_6H_{6-7}O_{0-1}.
    Every element that is present gets a subscript, a range when more than one number was seen.

    Args:
        carbons (int): Number of carbon atoms.
        hydrogens (list): Observed numbers of hydrogen atoms.
        others (dict or None): Element symbol to the observed numbers of atoms of that element, for example
            {'O': [0, 1]}. Elements that were never seen (all zero) are left out.

    Returns:
        str: Mathtext label.
    """
    def part(symbol, counts):
        counts = list(counts)
        if not counts or max(counts) == 0:
            return ''
        low, high = min(counts), max(counts)
        return f'{symbol}_{{{low}}}' if low == high else f'{symbol}_{{{low}\u2013{high}}}'
    body = part('C', [carbons]) + part('H', hydrogens)
    for symbol in sorted((others or {}), key=lambda e: e):
        body += part(symbol, others[symbol])
    return f'$\\mathrm{{{body}}}$'


def fractionYLimits(ax, floorStep=10, margin=4):
    """Set the y axis of a percentage plot to 0-102, or when all data is high, to start lower than the lowest
    lower band edge. A truncated axis is marked by the tick labels, so mention it in the caption.

    Args:
        ax (Axes): Axes with fill_between bands drawn on it.
        floorStep (int): Round the lower limit down to a multiple of this.
        margin (float): Space in percent under the lowest band edge.
    """
    import numpy as np
    lowest = min((np.min(c.get_paths()[0].vertices[:, 1]) for c in ax.collections if c.get_paths()), default=0)
    lines = [np.min(l.get_ydata()) for l in ax.lines if len(l.get_ydata()) > 2]
    lowest = min([lowest] + lines) if ax.collections else (min(lines) if lines else 0)
    bottom = 0 if lowest < 50 else max(0, floorStep * int((lowest - margin) // floorStep))
    ax.set_ylim(bottom, 102)
    if bottom > 0:
        ax.text(0.99, 0.02, f'y axis starts at {bottom} %', transform=ax.transAxes, ha='right', va='bottom',
                fontsize=6.5, color='#555555')
    return bottom


def styleFor(labels):
    """Colour and line style for every label, the same for the same label in every figure.

    The choice comes from a hash of the label, so a molecule keeps its look whichever molecules it is plotted
    with. Two labels in one figure never get the same look: the second one moves on to the next free one.
    The colour and the line style are tied together, so the curves also differ in grey scale.

    Args:
        labels (list): Names of the series in a figure.

    Returns:
        dict: Label to (colour, line style).
    """
    size = 6
    taken = set()
    styles = {}
    for label in labels:
        index = zlib.crc32(str(label).encode()) % size
        for _ in range(size):
            if index not in taken:
                break
            index = (index + 1) % size
        taken.add(index)
        styles[label] = (COLORS[index], LINE_STYLES[index])
    return styles


def ciProxy(text='95 % CI'):
    """Legend handle for a shaded confidence band, so the figure says what the shading is."""
    from matplotlib.patches import Patch
    return Patch(facecolor='#777777', alpha=0.3, linewidth=0, label=text)


def mixedLegend(ax, curveHandles, showBands, **kwargs):
    """Legend of the curves, with an entry for the shaded band when bands are drawn."""
    handles, labels = ax.get_legend_handles_labels()
    if showBands:
        handles.append(ciProxy())
        labels.append('95 % CI')
    return ax.legend(handles, labels, **kwargs)

# With more curves than this, bands of different curves overlap into an unreadable mix and only lines are drawn
MAX_CURVES_WITH_BANDS = 2
