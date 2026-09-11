"""The permutation that the tools pass around, and its on-disk form.

Every way of deciding a new atom order -- from a bond graph, from atom-name
matching against another file, from a map written by an earlier run -- produces
an AtomOrder, and everything that applies an order consumes one.  That is what
lets `update_group_ordering --write-map` and a later `--read-map` share a code path.
"""

from __future__ import annotations

from .errors import TopologyError


class AtomOrder:
    """A permutation of the atoms of one molecule.

    ``old_for_new[i]`` is the 1-based original index of the atom that ends up
    at 1-based position ``i + 1``.  ``names`` holds the atom names in the
    ORIGINAL order, so it can be checked against a coordinate file that has not
    been permuted yet.
    """

    __slots__ = ("old_for_new", "names", "resname")

    def __init__(self, old_for_new, names=None, resname=None):
        old_for_new = list(old_for_new)
        if sorted(old_for_new) != list(range(1, len(old_for_new) + 1)):
            raise TopologyError(
                "the atom order is not a permutation of 1..%d" % len(old_for_new)
            )
        self.old_for_new = old_for_new
        self.names = list(names) if names is not None else [""] * len(old_for_new)
        if len(self.names) != len(old_for_new):
            raise TopologyError("the atom order and its name list differ in length")
        self.resname = resname

    @classmethod
    def identity(cls, natoms, names=None, resname=None):
        return cls(range(1, natoms + 1), names, resname)

    @property
    def natoms(self):
        return len(self.old_for_new)

    @property
    def new_for_old(self):
        """{original 1-based index: new 1-based index}."""
        return {old: pos for pos, old in enumerate(self.old_for_new, start=1)}

    @property
    def new_names(self):
        """The atom names in the NEW order."""
        return [self.names[old - 1] for old in self.old_for_new]

    def inverse(self):
        """The order that undoes this one."""
        inv = [0] * self.natoms
        for pos, old in enumerate(self.old_for_new, start=1):
            inv[old - 1] = pos
        return AtomOrder(inv, self.new_names, self.resname)

    def moved(self):
        """How many atoms change index."""
        return sum(1 for pos, old in enumerate(self.old_for_new, start=1)
                   if pos != old)

    def apply(self, items):
        """Permute any sequence of one atom per position."""
        if len(items) != self.natoms:
            raise TopologyError(
                "cannot permute %d items with an order for %d atoms"
                % (len(items), self.natoms)
            )
        return [items[old - 1] for old in self.old_for_new]

    # -- the map file ----------------------------------------------------

    def write_map(self, path, source=None):
        with open(path, "w") as handle:
            handle.write("; atom index map written by topology_wrangling\n")
            if source:
                handle.write("; source %s\n" % source)
            handle.write("; resname %s\n" % (self.resname or "-"))
            handle.write("; natoms %d\n" % self.natoms)
            handle.write(";  new    old  name\n")
            for pos, old in enumerate(self.old_for_new, start=1):
                handle.write("%6d %6d  %s\n" % (pos, old, self.names[old - 1]))

    @classmethod
    def read_map(cls, path):
        """Read a map file written by write_map."""
        rows, resname = [], None
        with open(path) as handle:
            for raw in handle:
                line = raw.strip()
                if line.startswith(";"):
                    tok = line.lstrip(";").split()
                    if len(tok) >= 2 and tok[0] == "resname" and tok[1] != "-":
                        resname = tok[1]
                    continue
                if not line:
                    continue
                tok = line.split()
                if len(tok) < 2:
                    raise TopologyError(
                        "%s: malformed map line: %r" % (path, raw.rstrip())
                    )
                rows.append((int(tok[0]), int(tok[1]),
                             tok[2] if len(tok) > 2 else ""))
        if not rows:
            raise TopologyError("%s: the map file holds no entries" % path)
        rows.sort()
        if [new for new, _, _ in rows] != list(range(1, len(rows) + 1)):
            raise TopologyError(
                "%s: the map does not list new indices 1..N exactly once" % path
            )
        old_for_new = [old for _, old, _ in rows]
        names = [""] * len(rows)
        for _, old, name in rows:
            if not 1 <= old <= len(rows):
                raise TopologyError(
                    "%s: old index %d is outside 1..%d" % (path, old, len(rows))
                )
            names[old - 1] = name
        return cls(old_for_new, names, resname)
