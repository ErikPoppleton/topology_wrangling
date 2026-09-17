"""Atom colouring: the Jmol element scheme, and partial-charge colour maps.

The element table is Jmol's own, transcribed from the constants its colour
documentation is built from (https://jmol.sourceforge.net/jscolors/), so a
drawing made here matches what a 3D viewer would show.

Label colours are chosen by contrast, not by a fixed guess: white text on
Jmol's nitrogen blue is readable, black text on its sulphur yellow is, and a
threshold picked by eye gets one of the two wrong.
"""

from __future__ import annotations

from .errors import TopologyError

# Jmol's default element colours, by symbol.  Generated from the table at
# https://jmol.sourceforge.net/jscolors/ (elements 1-109).
JMOL_COLORS = {
    "H": "#FFFFFF", "He": "#D9FFFF", "Li": "#CC80FF", "Be": "#C2FF00",
    "B": "#FFB5B5", "C": "#909090", "N": "#3050F8", "O": "#FF0D0D",
    "F": "#90E050", "Ne": "#B3E3F5", "Na": "#AB5CF2", "Mg": "#8AFF00",
    "Al": "#BFA6A6", "Si": "#F0C8A0", "P": "#FF8000", "S": "#FFFF30",
    "Cl": "#1FF01F", "Ar": "#80D1E3", "K": "#8F40D4", "Ca": "#3DFF00",
    "Sc": "#E6E6E6", "Ti": "#BFC2C7", "V": "#A6A6AB", "Cr": "#8A99C7",
    "Mn": "#9C7AC7", "Fe": "#E06633", "Co": "#F090A0", "Ni": "#50D050",
    "Cu": "#C88033", "Zn": "#7D80B0", "Ga": "#C28F8F", "Ge": "#668F8F",
    "As": "#BD80E3", "Se": "#FFA100", "Br": "#A62929", "Kr": "#5CB8D1",
    "Rb": "#702EB0", "Sr": "#00FF00", "Y": "#94FFFF", "Zr": "#94E0E0",
    "Nb": "#73C2C9", "Mo": "#54B5B5", "Tc": "#3B9E9E", "Ru": "#248F8F",
    "Rh": "#0A7D8C", "Pd": "#006985", "Ag": "#C0C0C0", "Cd": "#FFD98F",
    "In": "#A67573", "Sn": "#668080", "Sb": "#9E63B5", "Te": "#D47A00",
    "I": "#940094", "Xe": "#429EB0", "Cs": "#57178F", "Ba": "#00C900",
    "La": "#70D4FF", "Ce": "#FFFFC7", "Pr": "#D9FFC7", "Nd": "#C7FFC7",
    "Pm": "#A3FFC7", "Sm": "#8FFFC7", "Eu": "#61FFC7", "Gd": "#45FFC7",
    "Tb": "#30FFC7", "Dy": "#1FFFC7", "Ho": "#00FF9C", "Er": "#00E675",
    "Tm": "#00D452", "Yb": "#00BF38", "Lu": "#00AB24", "Hf": "#4DC2FF",
    "Ta": "#4DA6FF", "W": "#2194D6", "Re": "#267DAB", "Os": "#266696",
    "Ir": "#175487", "Pt": "#D0D0E0", "Au": "#FFD123", "Hg": "#B8B8D0",
    "Tl": "#A6544D", "Pb": "#575961", "Bi": "#9E4FB5", "Po": "#AB5C00",
    "At": "#754F45", "Rn": "#428296", "Fr": "#420066", "Ra": "#007D00",
    "Ac": "#70ABFA", "Th": "#00BAFF", "Pa": "#00A1FF", "U": "#008FFF",
    "Np": "#0080FF", "Pu": "#006BFF", "Am": "#545CF2", "Cm": "#785CE3",
    "Bk": "#8A4FE3", "Cf": "#A136D4", "Es": "#B31FD4", "Fm": "#B31FBA",
    "Md": "#B30DA6", "No": "#BD0D87", "Lr": "#C70066", "Rf": "#CC0059",
    "Db": "#D1004F", "Sg": "#D90045", "Bh": "#E00038", "Hs": "#E6002E",
    "Mt": "#EB0026",
}

# What to use for an atom whose element we could not work out.
UNKNOWN_COLOR = "#FF1493"


def element_color(symbol):
    """The Jmol colour for an element symbol, as '#rrggbb'."""
    return JMOL_COLORS.get(symbol, UNKNOWN_COLOR)


def to_rgb(color):
    """'#rrggbb' -> (r, g, b) as 0-255 ints."""
    text = color.lstrip("#")
    if len(text) != 6:
        raise TopologyError("not a #rrggbb colour: %r" % color)
    return tuple(int(text[i:i + 2], 16) for i in (0, 2, 4))


def relative_luminance(color):
    """WCAG relative luminance of a colour, 0 (black) to 1 (white)."""
    channels = []
    for value in to_rgb(color):
        srgb = value / 255.0
        channels.append(srgb / 12.92 if srgb <= 0.04045
                        else ((srgb + 0.055) / 1.055) ** 2.4)
    r, g, b = channels
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def text_color(background):
    """Black or white, whichever contrasts better with the given fill.

    Uses the WCAG contrast ratio rather than a brightness threshold, so the
    choice is the one that is actually easier to read.
    """
    luminance = relative_luminance(background)
    against_black = (luminance + 0.05) / 0.05
    against_white = 1.05 / (luminance + 0.05)
    return "#000000" if against_black >= against_white else "#ffffff"


def value_colors(values, cmap="viridis", vmin=None, vmax=None,
                 symmetric=False):
    """Map per-atom numbers to colours with a matplotlib colour map.

    Returns (colours, vmin, vmax).  ``symmetric`` centres the range on zero,
    which is what a signed quantity like a partial charge wants: the middle of
    a diverging map then lands on a neutral atom and the sign is read straight
    off the colour.  A quantity with no meaningful zero crossing -- an LJ well
    depth, say -- is better served by a sequential map over its own range.
    """
    try:
        from matplotlib import colormaps
        from matplotlib.colors import Normalize, to_hex
    except ImportError:
        raise TopologyError(
            "colouring by a numeric property needs matplotlib; install it, or "
            "use --color-by element"
        )
    values = [float(v) for v in values]
    if symmetric:
        limit = max((abs(v) for v in values), default=0.0) or 1.0
        vmin = -limit if vmin is None else vmin
        vmax = limit if vmax is None else vmax
    else:
        vmin = min(values, default=0.0) if vmin is None else vmin
        vmax = max(values, default=1.0) if vmax is None else vmax
        if vmax - vmin < 1e-12:      # every atom the same: keep the bar sane
            vmin, vmax = vmin - 0.5, vmax + 0.5
    try:
        mapping = colormaps[cmap]
    except KeyError:
        raise TopologyError(
            "unknown matplotlib colour map %r; try one of viridis, plasma, "
            "magma, coolwarm, bwr, RdBu" % cmap
        )
    normalize = Normalize(vmin=vmin, vmax=vmax)
    return [to_hex(mapping(normalize(v))) for v in values], vmin, vmax


def charge_colors(charges, cmap="coolwarm", vmin=None, vmax=None):
    """Partial charges, on a diverging map centred on neutral."""
    return value_colors(charges, cmap, vmin, vmax, symmetric=True)
