"""Matching the residues of a coordinate file to topologies.

A .gro written by one tool and an .itp from a force field often disagree about
names as well as order: CHARMM-GUI writes TIP3/OH2 where the topology says
TP3/O, and SOD where it says Na+.  An AliasTable records those disagreements,
and apply_topologies walks a frame, re-orders every molecule it has a topology
for and rewrites its names, leaving everything else exactly where it is.

Alias files use one directive per line:

    residue TIP3 TP3        ; the .gro calls TP3 "TIP3"
    atom    TP3  O   OH2    ; within TP3, the .gro calls the atom O "OH2"
"""

from __future__ import annotations

from .errors import TopologyError
from .ordering import order_from_names


class AliasTable:
    """Name disagreements between a coordinate file and a set of topologies."""

    def __init__(self):
        self.residues = {}   # gro residue name -> itp residue name
        self.atoms = {}      # itp residue name -> {itp atom name: gro atom name}

    def add_residue(self, gro_name, itp_name):
        self.residues[gro_name] = itp_name

    def add_atom(self, resname, itp_name, gro_name):
        self.atoms.setdefault(resname, {})[itp_name] = gro_name

    def atoms_for(self, itp_resname):
        return self.atoms.get(itp_resname, {})

    def gro_name_for(self, itp_resname):
        """The residue name the coordinate file uses for this topology."""
        for gro_name, itp_name in self.residues.items():
            if itp_name == itp_resname:
                return gro_name
        return itp_resname

    # -- parsing ---------------------------------------------------------

    def add_residue_spec(self, spec):
        """`GRO=ITP`, as given on the command line."""
        if "=" not in spec:
            raise TopologyError("--residue wants GRO_NAME=ITP_NAME, got %r" % spec)
        gro_name, itp_name = spec.split("=", 1)
        self.add_residue(gro_name.strip(), itp_name.strip())

    def add_atom_spec(self, spec, default_resname=None):
        """`[RES:]ITP_NAME=GRO_NAME`, as given on the command line."""
        resname = default_resname
        body = spec
        if ":" in spec.split("=", 1)[0]:
            resname, body = spec.split(":", 1)
            resname = resname.strip()
        if "=" not in body:
            raise TopologyError(
                "--alias wants [RESIDUE:]ITP_NAME=GRO_NAME, got %r" % spec
            )
        itp_name, gro_name = body.split("=", 1)
        if resname is None:
            raise TopologyError(
                "--alias %r does not say which residue it applies to; write "
                "RESIDUE:ITP_NAME=GRO_NAME" % spec
            )
        self.add_atom(resname, itp_name.strip(), gro_name.strip())

    def load(self, path):
        """Read an alias file."""
        with open(path) as handle:
            for lineno, raw in enumerate(handle, start=1):
                line = raw.split(";", 1)[0].split("#", 1)[0].strip()
                if not line:
                    continue
                tok = line.split()
                where = "%s line %d" % (path, lineno)
                if tok[0] == "residue" and len(tok) == 3:
                    self.add_residue(tok[1], tok[2])
                elif tok[0] == "atom" and len(tok) == 4:
                    self.add_atom(tok[1], tok[2], tok[3])
                else:
                    raise TopologyError(
                        "%s: expected 'residue GRO ITP' or "
                        "'atom RESIDUE ITP_NAME GRO_NAME', got %r"
                        % (where, raw.strip())
                    )
        return self


class ResidueTopology:
    """One topology, paired with the residue name the coordinates use for it."""

    __slots__ = ("itp", "itp_resname", "gro_resname", "names", "aliases")

    def __init__(self, itp, itp_resname, gro_resname, aliases=None):
        self.itp = itp
        self.itp_resname = itp_resname
        self.gro_resname = gro_resname
        self.names = itp.atom_names
        self.aliases = aliases or {}

    @property
    def size(self):
        return len(self.names)

    @classmethod
    def from_itp(cls, itp, table=None, resname=None, gro_resname=None):
        table = table or AliasTable()
        itp_resname = itp.residue_name(resname)
        if gro_resname is None:
            gro_resname = table.gro_name_for(itp_resname)
        return cls(itp, itp_resname, gro_resname, table.atoms_for(itp_resname))


def apply_topologies(frame, topologies, rename=True, strict=False,
                     out_resnames=None, report=None):
    """Re-order and rename every molecule of the frame that has a topology.

    Returns {residue name as the coordinates spell it: molecules touched}.  Residues with no
    topology are left untouched; with ``strict`` they are an error instead.
    """
    out_resnames = out_resnames or {}
    by_resname = {}
    for top in topologies:
        if top.gro_resname in by_resname:
            raise TopologyError(
                "two topologies both claim the residue %s in the coordinates; "
                "separate them with --residue GRO_NAME=ITP_NAME" % top.gro_resname
            )
        by_resname[top.gro_resname] = top

    missing = [name for name in frame.resnames if name not in by_resname]
    if missing and strict:
        raise TopologyError(
            "no topology for residue(s) %s in %s"
            % (", ".join(missing), frame._where())
        )
    if missing and report:
        report("no topology for %s; left untouched" % ", ".join(missing))

    counts = {}
    out = []
    i, total = 0, len(frame.atoms)
    orders = {}   # (resname, atom-name tuple) -> AtomOrder
    while i < total:
        atom = frame.atoms[i]
        gro_resname = atom.resname          # read before any renaming
        top = by_resname.get(gro_resname)
        if top is None:
            out.append(atom)
            i += 1
            continue

        block = frame.atoms[i:i + top.size]
        resid = block[0].resid
        if len(block) != top.size or any(
            a.resname != gro_resname or a.resid != resid for a in block
        ):
            raise TopologyError(
                "%s: the %s residue starting at atom %d is not %d contiguous "
                "atoms of a single residue"
                % (frame._where(), gro_resname, i + 1, top.size)
            )

        # The permutation depends only on the sequence of atom names, which
        # is the same for every molecule of a residue in practice, so solve it
        # once per distinct sequence instead of once per molecule.
        key = (gro_resname, tuple(a.name for a in block))
        order = orders.get(key)
        if order is None:
            order = order_from_names(
                top.names, [a.name for a in block], top.aliases,
                where="%s: %s at atom %d: "
                      % (frame._where(), gro_resname, i + 1),
            )
            orders[key] = order

        block = order.apply(block)
        if rename:
            out_resname = out_resnames.get(top.itp_resname, top.itp_resname)
            for position, a in enumerate(block):
                a.name = top.names[position]
                a.resname = out_resname
        out.extend(block)
        counts[gro_resname] = counts.get(gro_resname, 0) + 1
        i += top.size

    frame.atoms = out
    frame.renumber()
    if not counts:
        raise TopologyError(
            "%s: none of the topologies matched a residue in the coordinates "
            "(the file has %s)" % (frame._where(), ", ".join(frame.resnames))
        )
    return counts
