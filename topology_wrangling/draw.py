"""Drawing a topology as a 2D graph.

There are plenty of 3D molecular viewers; this makes the flat, labelled
picture a topology file actually calls for -- every atom a circle carrying its
unique name, every bond an edge, laid out so the connectivity is readable.

A drawing is built as a Scene of plain primitives and then written out by one
of two backends: SVG, which needs nothing but the standard library, and PNG
through matplotlib.  Keeping the scene separate from the backends is what
stops the two outputs from drifting apart.
"""

from __future__ import annotations

import math

from . import colors as color_module
from .chem import element_of, is_hydrogen
from .errors import TopologyError
from .layout import force_directed, relax_overlaps
from .skeleton import skeletal

# How much smaller a hydrogen is drawn than a heavy atom.  Hydrogens outnumber
# everything else in an all-atom topology, and shrinking them lets the heavy
# atom skeleton carry the picture.
HYDROGEN_SCALE = 0.68

# Matplotlib wants a list of candidates; the SVG gets the same stack as a
# CSS font-family string.
FONT_STACK = ["Helvetica", "Arial", "DejaVu Sans", "sans-serif"]

class Property:
    """A per-atom number that can be painted onto the drawing."""

    __slots__ = ("name", "label", "unit", "cmap", "symmetric", "needs_ff")

    def __init__(self, name, label, unit, cmap, symmetric, needs_ff=False):
        self.name = name
        self.label = label
        self.unit = unit
        self.cmap = cmap            # the map that suits this quantity
        self.symmetric = symmetric  # centre the scale on zero?
        self.needs_ff = needs_ff    # does it come from a force field file?

    def caption(self):
        return "%s (%s)" % (self.label, self.unit) if self.unit else self.label


# A signed quantity gets a diverging map centred on zero, so the sign reads
# straight off the colour; one with no zero crossing gets a sequential map.
PROPERTIES = {
    "charge": Property("charge", "partial charge", "e", "coolwarm", True),
    "epsilon": Property("epsilon", "LJ epsilon", "kJ/mol", "viridis", False,
                        needs_ff=True),
}

BACKGROUND = "#ffffff"
BOND_COLOR = "#555555"
OUTLINE = "#333333"


class Scene:
    """A drawing, as primitives in pixel coordinates with y running down."""

    def __init__(self, width, height, background=BACKGROUND):
        self.width = width
        self.height = height
        self.background = background
        self.lines = []     # (x1, y1, x2, y2, color, width)
        self.circles = []   # (x, y, r, fill, stroke, stroke_width)
        self.texts = []     # (x, y, text, size, color, anchor)
        self.rects = []     # (x, y, w, h, fill)

    def line(self, x1, y1, x2, y2, color=BOND_COLOR, width=2.0):
        self.lines.append((x1, y1, x2, y2, color, width))

    def circle(self, x, y, r, fill, stroke=OUTLINE, stroke_width=1.0):
        self.circles.append((x, y, r, fill, stroke, stroke_width))

    def text(self, x, y, text, size, color="#000000", anchor="middle"):
        self.texts.append((x, y, text, size, color, anchor))

    def rect(self, x, y, w, h, fill):
        self.rects.append((x, y, w, h, fill))

    # -- output ----------------------------------------------------------

    def to_svg(self):
        out = [
            '<?xml version="1.0" encoding="UTF-8"?>',
            '<svg xmlns="http://www.w3.org/2000/svg" width="%g" height="%g" '
            'viewBox="0 0 %g %g">' % (self.width, self.height,
                                      self.width, self.height),
            '<rect width="100%%" height="100%%" fill="%s"/>' % self.background,
            '<g font-family="Helvetica,Arial,sans-serif">',
        ]
        for x, y, w, h, fill in self.rects:
            out.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" '
                       'fill="%s"/>' % (x, y, w, h, fill))
        for x1, y1, x2, y2, color, width in self.lines:
            out.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" '
                       'stroke="%s" stroke-width="%.2f" stroke-linecap="round"/>'
                       % (x1, y1, x2, y2, color, width))
        for x, y, r, fill, stroke, stroke_width in self.circles:
            out.append('<circle cx="%.2f" cy="%.2f" r="%.2f" fill="%s" '
                       'stroke="%s" stroke-width="%.2f"/>'
                       % (x, y, r, fill, stroke, stroke_width))
        for x, y, text, size, color, anchor in self.texts:
            out.append('<text x="%.2f" y="%.2f" font-size="%.2f" fill="%s" '
                       'text-anchor="%s" dominant-baseline="central">%s</text>'
                       % (x, y, size, color, anchor, _escape(text)))
        out.append("</g></svg>")
        return "\n".join(out) + "\n"

    def to_png(self, path, dpi=200):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.patches import Circle, Rectangle
        except ImportError:
            raise TopologyError(
                "writing .png needs matplotlib; install it, or write .svg "
                "instead (which needs nothing)"
            )
        figure = plt.figure(figsize=(self.width / 100.0, self.height / 100.0),
                            dpi=dpi)
        axes = figure.add_axes([0, 0, 1, 1])
        axes.set_xlim(0, self.width)
        axes.set_ylim(self.height, 0)          # y runs down, as in the SVG
        axes.set_axis_off()
        axes.set_facecolor(self.background)
        figure.patch.set_facecolor(self.background)

        # The figure is width/100 inches across and spans `width` scene units,
        # so one scene unit is 1/100 inch.  Matplotlib measures text and line
        # widths in points (1/72 inch), and that conversion does not depend on
        # dpi -- scaling it by dpi would make labels grow with resolution and
        # silently disagree with the SVG.
        scale = 72.0 / 100.0
        for x, y, w, h, fill in self.rects:
            axes.add_patch(Rectangle((x, y), w, h, facecolor=fill, lw=0))
        for x1, y1, x2, y2, color, width in self.lines:
            axes.plot([x1, x2], [y1, y2], color=color, lw=width * scale,
                      solid_capstyle="round", zorder=1)
        for x, y, r, fill, stroke, stroke_width in self.circles:
            axes.add_patch(Circle((x, y), r, facecolor=fill, edgecolor=stroke,
                                  lw=stroke_width * scale, zorder=2))
        for x, y, text, size, color, anchor in self.texts:
            axes.text(x, y, text, fontsize=size * scale, color=color,
                      ha={"middle": "center", "start": "left",
                          "end": "right"}[anchor],
                      va="center_baseline", zorder=3, family=FONT_STACK)
        figure.savefig(path, dpi=dpi, facecolor=self.background)
        plt.close(figure)

    def save(self, path, dpi=200):
        """Write the scene, choosing the backend from the file extension."""
        lower = path.lower()
        if lower.endswith(".svg"):
            with open(path, "w") as handle:
                handle.write(self.to_svg())
        elif lower.endswith(".png"):
            self.to_png(path, dpi)
        else:
            raise TopologyError(
                "do not know how to write %r; use a .svg or .png name" % path
            )


def _escape(text):
    return (text.replace("&", "&amp;").replace("<", "&lt;")
                .replace(">", "&gt;"))


def draw_topology(itp, label="name", color_by="element", cmap=None,
                  value_range=(None, None), atomtypes=None,
                  hide_hydrogens=False,
                  layout="skeletal", branch="compact", iterations=400, seed=0,
                  node_radius=24.0, max_size=5000.0, title=True,
                  colorbar=True, extend="auto"):
    """Build a Scene from a topology.

    label        what to write in each atom: name, type, nr or none
    color_by     element (the Jmol scheme), or a name from PROPERTIES painted
                 with a matplotlib colour map
    cmap         overrides the map that property would use by default
    atomtypes    {type name: AtomType}, for a property read from a force field
    title        draw the molecule name and atom count above the structure
    colorbar     draw the colour bar beneath it; turn both off for a panel
                 that will share one bar with its neighbours
    extend       mark the ends of the bar the scale does not reach; "auto"
                 works it out from the atoms actually drawn
    layout       skeletal (structural-formula angles) or force (spring model)
    branch       how a skeletal layout aims branches: compact, spread or auto
    node_radius  drawn radius of a heavy atom, in pixels

    The canvas is sized from the atoms rather than the other way round: a long
    lipid needs more room than a water, and fitting either into one fixed
    canvas would shrink its labels past reading.  max_size caps the result.
    """
    atoms = itp.atoms
    keep = _visible_atoms(atoms, hide_hydrogens)
    if not keep:
        raise TopologyError("nothing left to draw once hydrogens are hidden")

    new_index = {old: i for i, old in enumerate(keep)}
    edges = [(new_index[i] + 1, new_index[j] + 1)
             for i, j in itp.connectivity()
             if i in new_index and j in new_index]

    hydrogens = {i + 1 for i, old in enumerate(keep)
                 if is_hydrogen(atoms[old - 1])}
    if layout == "skeletal":
        positions = skeletal(len(keep), edges, hydrogens, branch)
    elif layout == "force":
        positions = force_directed(len(keep), edges, iterations, seed)
    else:
        raise TopologyError(
            "unknown layout %r; use skeletal or force" % layout)
    fills, legend = _fills(itp, keep, color_by, cmap, value_range, atomtypes,
                           extend)
    labels = [_label_for(atoms[old - 1], label) for old in keep]

    if not colorbar:
        legend = None
    return _compose(itp, keep, positions, edges, fills, labels, legend,
                    hydrogens, node_radius, max_size, title)


def _visible_atoms(atoms, hide_hydrogens):
    """The 1-based indices of the atoms to draw."""
    if not hide_hydrogens:
        return list(range(1, len(atoms) + 1))
    return [i for i, atom in enumerate(atoms, start=1) if not is_hydrogen(atom)]


def _label_for(atom, label):
    if label == "name":
        return atom.name
    if label == "type":
        return atom.type
    if label == "nr":
        return str(atom.nr)
    if label == "none":
        return ""
    raise TopologyError("unknown label %r; use name, type, nr or none" % label)


def _fills(itp, keep, color_by, cmap, value_range, atomtypes, extend="auto"):
    """The fill colour of every drawn atom, and any legend it implies."""
    atoms = itp.atoms
    if color_by == "element":
        return ([color_module.element_color(element_of(atoms[old - 1]))
                 for old in keep], None)

    prop = PROPERTIES.get(color_by)
    if prop is None:
        raise TopologyError(
            "unknown colouring %r; use element, %s"
            % (color_by, " or ".join(sorted(PROPERTIES)))
        )
    values = _values(itp, keep, prop, atomtypes)
    cmap = cmap or prop.cmap
    fills, vmin, vmax = color_module.value_colors(
        values, cmap, value_range[0], value_range[1], prop.symmetric)
    if extend == "auto":
        extend = _extend_for(values, vmin, vmax)
    elif extend not in EXTEND_CHOICES:
        raise TopologyError("extend must be one of %s"
                            % ", ".join(EXTEND_CHOICES))
    return fills, (prop.caption(), cmap, vmin, vmax, extend)


def _extend_for(values, vmin, vmax, tolerance=1e-9):
    """Which ends of the scale the drawn atoms actually run past."""
    below = any(v < vmin - tolerance for v in values)
    above = any(v > vmax + tolerance for v in values)
    if below and above:
        return "both"
    return "min" if below else "max" if above else "neither"


def _values(itp, keep, prop, atomtypes):
    """One number per drawn atom, from the molecule or from a force field."""
    atoms = itp.atoms
    if prop.name == "charge":
        # Deliberately the molecule's own per-atom charge, not the force
        # field's default for the type: the two routinely differ.
        values = []
        for old in keep:
            try:
                values.append(float(atoms[old - 1].charge))
            except (TypeError, ValueError):
                raise TopologyError(
                    "atom %d (%s) has no numeric charge, so it cannot be "
                    "coloured by charge" % (old, atoms[old - 1].name))
        return values

    if not atomtypes:
        raise TopologyError(
            "colouring by %s needs the force field that defines the atom "
            "types; give it with --ff" % prop.name)
    missing = sorted({atoms[old - 1].type for old in keep} - set(atomtypes))
    if missing:
        raise TopologyError(
            "these atom types are not in the force field given: %s"
            % ", ".join(missing))
    values = []
    for old in keep:
        entry = atomtypes[atoms[old - 1].type]
        value = getattr(entry, prop.name, None)
        if value is None:
            raise TopologyError(
                "atom type %s has no %s in the force field"
                % (entry.name, prop.name))
        values.append(value)
    return values


def _compose(itp, keep, positions, edges, fills, labels, legend, hydrogens,
             node_radius, max_size, title):
    """Turn laid-out atoms into a Scene, sized so nothing collides."""
    xs = [x for x, _ in positions]
    ys = [y for _, y in positions]
    span_x = max(max(xs) - min(xs), 1e-6)
    span_y = max(max(ys) - min(ys), 1e-6)

    # A radius derived from the typical bond length keeps circles apart at any
    # molecule size, instead of a fixed guess that only suits one of them.
    radius_units = _radius_units(positions, edges, len(keep), hydrogens)

    radii_units = [radius_units * (HYDROGEN_SCALE if i + 1 in hydrogens else 1.0)
                   for i in range(len(keep))]
    # Separate any atoms whose circles (and so whose labels) would collide,
    # then re-measure, since the relaxation moves things.  A skeletal layout
    # usually has nothing to fix, so this is a no-op on the common path.
    if _has_overlap(positions, radii_units):
        positions = relax_overlaps(positions, radii_units)
    xs = [x for x, _ in positions]
    ys = [y for _, y in positions]
    span_x = max(max(xs) - min(xs), 1e-6)
    span_y = max(max(ys) - min(ys), 1e-6)
    margin_units = radius_units * 2.2

    # Scale so that a heavy atom comes out at the requested radius, then pull
    # the whole drawing back if that would make an unwieldy canvas.
    scale = node_radius / radius_units
    header = 46.0 if title else 0.0
    footer = 58.0 if legend else 0.0
    if max_size:
        longest = max((span_x + 2 * margin_units) * scale,
                      (span_y + 2 * margin_units) * scale + header + footer)
        if longest > max_size:
            scale *= max_size / longest
    width = (span_x + 2 * margin_units) * scale
    height = (span_y + 2 * margin_units) * scale + header + footer

    def place(point):
        x, y = point
        return ((x - min(xs) + margin_units) * scale,
                (y - min(ys) + margin_units) * scale + header)

    scene = Scene(round(width, 2), round(height, 2))
    radius = radius_units * scale
    bond_width = max(1.2, radius * 0.22)

    for i, j in edges:
        x1, y1 = place(positions[i - 1])
        x2, y2 = place(positions[j - 1])
        scene.line(x1, y1, x2, y2, BOND_COLOR, bond_width)

    for position, r_units, fill, text in zip(positions, radii_units, fills,
                                             labels):
        x, y = place(position)
        r = r_units * scale
        scene.circle(x, y, r, fill, OUTLINE, max(0.6, radius * 0.06))
        if text:
            scene.text(x, y, text, _font_size(text, r),
                       color_module.text_color(fill))

    if title:
        name = itp.molecule_name or itp.residue_names[0]
        total = _total_charge(itp)
        scene.text(14, header / 2,
                   "%s  -  %d atoms%s" % (name, itp.natoms,
                                          _charge_caption(total)),
                   20, "#222222", anchor="start")
    if legend:
        _draw_colorbar(scene, legend, width, height - footer + 14)
    return scene


def _radius_units(positions, edges, natoms, hydrogens=()):
    """A third of a typical skeleton bond, so bonded circles nearly touch.

    Hydrogen spokes are deliberately shorter than skeleton bonds, so they are
    left out of the measurement: sizing from them would shrink every heavy
    atom in an all-atom topology.
    """
    heavy = [(i, j) for i, j in edges
             if i not in hydrogens and j not in hydrogens]
    lengths = sorted(math.dist(positions[i - 1], positions[j - 1])
                     for i, j in (heavy or edges))
    if lengths:
        typical = lengths[len(lengths) // 2]
    else:
        typical = 2.0 / max(math.sqrt(natoms), 1.0)
    return max(typical * 0.34, 1e-3)


def _has_overlap(positions, radii, padding=1.02):
    """Whether any two circles collide; cheap enough to check before relaxing."""
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            if math.dist(positions[i], positions[j]) < \
                    (radii[i] + radii[j]) * padding:
                return True
    return False


def _font_size(text, radius):
    """The largest size that still fits the label inside its circle."""
    # Helvetica capitals and digits run about 0.68 em wide.
    fits_width = 1.7 * radius / (0.68 * max(len(text), 1))
    return max(min(fits_width, radius * 1.05), radius * 0.34)


def _total_charge(itp):
    """The net charge, rounded past the noise of summing decimal charges."""
    total = 0.0
    for atom in itp.atoms:
        try:
            total += float(atom.charge)
        except (TypeError, ValueError):
            return None
    return round(total, 6) + 0.0    # + 0.0 turns a rounded -0.0 into 0.0


def _charge_caption(total):
    if total is None:
        return ""
    if total == 0:
        return ", neutral"
    return ", net charge %+g" % total


BAR_STEPS = 64          # fewest blocks a gradient is ever drawn in
BAR_THICKNESS = 14.0


def _bar_steps(length):
    """How many blocks to draw a gradient of this length in.

    The gradient is a row of solid rectangles, so the block count sets how
    smooth it is.  Scaling it with the bar's length keeps every block well
    under a scene unit, which stays smooth however high the output resolution
    is turned up -- a fixed count bands visibly once a bar is rendered large.
    """
    return max(BAR_STEPS, int(round(length * 2)))

# How a colour bar says its range does not cover every value.  Atoms outside
# it are drawn in the end colour, so the end label has to admit as much: a
# scale topped at the nitrogens keeps the carbons apart, but the phosphorus
# sitting on it must not read as though it were 0.92.
EXTEND_CHOICES = ("auto", "neither", "min", "max", "both")


def _bar_labels(vmin, vmax, extend):
    """The two end labels, marked where the scale does not reach."""
    low, high = "%.3g" % vmin, "%.3g" % vmax
    if extend in ("min", "both"):
        low = "<" + low
    if extend in ("max", "both"):
        high = ">" + high
    return low, high


def _bar_fills(cmap, vmin, vmax, steps=BAR_STEPS):
    """The colour map sampled evenly across its range."""
    fills, _, _ = color_module.value_colors(
        [vmin + (vmax - vmin) * i / (steps - 1.0) for i in range(steps)],
        cmap, vmin, vmax,
    )
    return fills


def _draw_colorbar(scene, legend, width, top):
    """The colour bar that sits under a molecule."""
    caption, cmap, vmin, vmax, extend = legend
    bar_width = min(340.0, width * 0.5)
    left = (width - bar_width) / 2
    fills = _bar_fills(cmap, vmin, vmax, _bar_steps(bar_width))
    step = bar_width / len(fills)
    for i, fill in enumerate(fills):
        scene.rect(left + i * step, top, step + 0.5, BAR_THICKNESS, fill)
    low, high = _bar_labels(vmin, vmax, extend)
    scene.text(left - 8, top + BAR_THICKNESS / 2, low, 13, "#222222",
               anchor="end")
    scene.text(left + bar_width + 8, top + BAR_THICKNESS / 2, high, 13,
               "#222222", anchor="start")
    # The colour map's name is an implementation detail of the figure, not
    # something a reader of it needs; the caption says what is being shown.
    scene.text(width / 2, top + 30, caption, 13, "#222222")


def colorbar_scene(color_by, cmap=None, value_range=(None, None), length=360.0,
                   thickness=18.0, orientation="horizontal", font=14.0,
                   extend="neither"):
    """A Scene holding nothing but a colour bar.

    Wanted when several molecules go into one figure: the panels are drawn
    bare and share a single bar, instead of each carrying its own.  The range
    must be given, since there is no molecule here to take it from.
    """
    prop = PROPERTIES.get(color_by)
    if prop is None:
        raise TopologyError(
            "unknown property %r; use %s"
            % (color_by, " or ".join(sorted(PROPERTIES))))
    vmin, vmax = value_range
    if vmin is None or vmax is None:
        raise TopologyError(
            "a standalone colour bar needs its range: pass both ends")
    if vmax <= vmin:
        raise TopologyError("the colour bar range must increase: got %g to %g"
                            % (vmin, vmax))
    if orientation not in ("horizontal", "vertical"):
        raise TopologyError("orientation must be horizontal or vertical")
    if extend not in EXTEND_CHOICES or extend == "auto":
        raise TopologyError(
            "a standalone colour bar needs an explicit extend: %s"
            % ", ".join(c for c in EXTEND_CHOICES if c != "auto"))

    cmap = cmap or prop.cmap
    fills = _bar_fills(cmap, vmin, vmax, _bar_steps(length))
    caption = prop.caption()
    low, high = _bar_labels(vmin, vmax, extend)
    pad = font * 0.75

    if orientation == "horizontal":
        # end labels either side of the bar, caption centred beneath it
        side = max(len(low), len(high)) * font * 0.62 + pad
        width = length + 2 * side
        height = thickness + font * 2.4
        scene = Scene(round(width, 2), round(height, 2))
        step = length / len(fills)
        for i, fill in enumerate(fills):
            scene.rect(side + i * step, 0.0, step + 0.5, thickness, fill)
        scene.text(side - pad, thickness / 2, low, font, "#222222", anchor="end")
        scene.text(side + length + pad, thickness / 2, high, font, "#222222",
                   anchor="start")
        scene.text(width / 2, thickness + font * 1.2, caption, font, "#222222")
        return scene

    # vertical: caption above, low at the bottom, high at the top, both to the
    # right of the bar.  No rotated text, which keeps every backend honest.
    label = max(len(low), len(high)) * font * 0.62 + pad
    width = max(thickness + label, len(caption) * font * 0.62)
    height = length + font * 1.8
    scene = Scene(round(width, 2), round(height, 2))
    top = font * 1.8
    step = length / len(fills)
    for i, fill in enumerate(fills):
        # the first block is the low end, which belongs at the bottom
        y = top + length - (i + 1) * step
        scene.rect(0.0, y, thickness, step + 0.5, fill)
    scene.text(thickness + pad, top + length, low, font, "#222222",
               anchor="start")
    scene.text(thickness + pad, top, high, font, "#222222", anchor="start")
    scene.text(0.0, font * 0.7, caption, font, "#222222", anchor="start")
    return scene
