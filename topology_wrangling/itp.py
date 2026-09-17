"""Reading, editing and writing GROMACS .itp topologies.

An ItpFile keeps every physical line of the input, tagged with the directive
it belongs to, so directives this library knows nothing about survive a
round-trip unchanged.  On top of that it offers typed views -- the [ atoms ]
records, the rows of the index-carrying directives, the bond graph -- that the
tools work with.
"""

from __future__ import annotations

from .chem import bond_graph
from .errors import TopologyError


class DirectiveSpec:
    """How many leading columns of a directive's rows are atom indices.

    ``kind`` is one of:
      "fixed"   -- the first ``ncols`` tokens are atom indices
      "all"     -- every token is an atom index ([ exclusions ])
      "vsiten"  -- site, funct, then a variable-length list ([ virtual_sitesn ])
    """

    __slots__ = ("name", "ncols", "kind")

    def __init__(self, name, ncols, kind="fixed"):
        self.name = name
        self.ncols = ncols
        self.kind = kind

    def index_positions(self, tokens):
        """The token positions that hold atom indices, for one row."""
        if self.kind == "all":
            return range(len(tokens))
        if self.kind == "vsiten":
            # site funct a1 a2 ... : columns 0 and 2..end
            return [0] + list(range(2, len(tokens)))
        return range(min(self.ncols, len(tokens)))

    def sort_key(self, tokens):
        """Key that orders the rows of a block by their atom indices."""
        positions = list(self.index_positions(tokens))
        return tuple(int(tokens[i]) for i in positions)

    def format_row(self, tokens):
        """Right-align the atom columns, then the function type, then the rest."""
        if self.kind != "fixed":
            return "".join("%6s" % t for t in tokens) + " "
        ncols = self.ncols
        parts = ["%6s" % t for t in tokens[:ncols]]
        parts.extend("%6s" % t for t in tokens[ncols:ncols + 1])      # funct
        parts.extend("%14s" % t for t in tokens[ncols + 1:])          # parameters
        return "".join(parts) + " "


def _spec(name, ncols, kind="fixed"):
    return name, DirectiveSpec(name, ncols, kind)


# Directives whose leading columns are atom indices.  Everything after those
# columns (function type, parameters) is copied verbatim.  Adding support for
# a new directive is one entry here.
DIRECTIVES = dict(
    [
        _spec("bonds", 2),
        _spec("pairs", 2),
        _spec("pairs_nb", 2),
        _spec("angles", 3),
        _spec("dihedrals", 4),
        _spec("constraints", 2),
        _spec("settles", 1),
        _spec("position_restraints", 1),
        _spec("distance_restraints", 2),
        _spec("dihedral_restraints", 4),
        _spec("orientation_restraints", 2),
        _spec("angle_restraints", 4),
        _spec("angle_restraints_z", 2),
        _spec("virtual_sites1", 1),
        _spec("virtual_sites2", 3),
        _spec("virtual_sites3", 4),
        _spec("virtual_sites4", 5),
        _spec("virtual_sitesn", None, "vsiten"),
        _spec("exclusions", None, "all"),
    ]
)


class Atom:
    """One record of the [ atoms ] directive."""

    __slots__ = ("nr", "type", "resnr", "residue", "name", "cgnr", "charge",
                 "mass", "extra", "comment")

    def __init__(self, nr, type, resnr, residue, name, cgnr, charge=None,
                 mass=None, extra=(), comment=""):
        self.nr = nr
        self.type = type
        self.resnr = resnr
        self.residue = residue
        self.name = name
        self.cgnr = cgnr
        self.charge = charge
        self.mass = mass
        self.extra = list(extra)
        self.comment = comment

    def __repr__(self):
        return "Atom(%d, %r, %r)" % (self.nr, self.name, self.type)

    def format(self, nr=None, cgnr=None):
        """The .itp text for this record, optionally with new nr / cgnr."""
        fields = [
            "%6d" % (self.nr if nr is None else nr),
            "%11s" % self.type,
            "%7s" % self.resnr,
            "%7s" % self.residue,
            "%7s" % self.name,
            "%7s" % (self.cgnr if cgnr is None else cgnr),
            "%11s" % (self.charge if self.charge is not None else ""),
        ]
        if self.mass is not None:
            fields.append("%11s" % ("%.4f" % self.mass).rstrip("0").rstrip("."))
        fields.extend("%9s" % t for t in self.extra)
        line = "".join(fields)
        if self.comment:
            line += "   " + self.comment
        return line


# The particle-type column of [ atomtypes ]: A(tom), S(hell), V/D (virtual).
PTYPE_CODES = ("A", "S", "V", "D")


class AtomType:
    """One record of the [ atomtypes ] directive of a force field.

    Note that its ``charge`` is the force field's default for the type, which
    a molecule's [ atoms ] records routinely override -- Slipids gives the
    type NTL a charge of -0.60, while the nitrogen of DOPC carries +0.20.  For
    a per-atom charge, read the molecule, not the force field.
    """

    __slots__ = ("name", "atomic_number", "mass", "charge", "ptype", "sigma",
                 "epsilon", "comment")

    def __init__(self, name, atomic_number=None, mass=None, charge=None,
                 ptype="A", sigma=None, epsilon=None, comment=""):
        self.name = name
        self.atomic_number = atomic_number
        self.mass = mass
        self.charge = charge
        self.ptype = ptype
        self.sigma = sigma          # nm
        self.epsilon = epsilon      # kJ/mol
        self.comment = comment

    def __repr__(self):
        return "AtomType(%r, sigma=%r, epsilon=%r)" % (self.name, self.sigma,
                                                       self.epsilon)


class Line:
    """One physical line of the file, tagged with the directive it is inside."""

    __slots__ = ("raw", "directive", "data", "comment")

    def __init__(self, raw, directive, data=None, comment=""):
        self.raw = raw
        self.directive = directive
        self.data = data          # list of tokens, for data lines only
        self.comment = comment    # trailing ';...' comment, if any

    def is_data(self):
        return self.data is not None


def split_comment(text):
    """Return (code, comment); the comment keeps its leading ';'."""
    idx = text.find(";")
    if idx == -1:
        return text, ""
    return text[:idx], text[idx:]


def split_text(text):
    """Split .itp text into lines without inventing or losing a final one."""
    raws = text.split("\n")
    if raws and raws[-1] == "" and text.endswith("\n"):
        raws.pop()          # the newline that ends the last line, not a line
    return raws


def parse_lines(raws):
    """Turn raw .itp lines (without their newlines) into Line objects."""
    lines = []
    directive = None
    for raw in raws:
        code, comment = split_comment(raw)
        stripped = code.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            directive = stripped[1:-1].strip().lower()
            lines.append(Line(raw, directive))
            continue
        entry = Line(raw, directive)
        if stripped and not stripped.startswith("#"):
            entry.data = stripped.split()
            entry.comment = comment
        lines.append(entry)
    return lines


class ItpFile:
    """A parsed .itp topology."""

    def __init__(self, lines, path=None):
        self.lines = lines
        self.path = path
        self._atoms = None
        self._atomtypes = None

    # -- construction ----------------------------------------------------

    @classmethod
    def read(cls, path):
        with open(path) as handle:
            return cls.from_text(handle.read(), path)

    @classmethod
    def from_text(cls, text, path=None):
        return cls(parse_lines(split_text(text)), path)

    def write(self, path):
        with open(path, "w") as handle:
            handle.write("\n".join(line.raw for line in self.lines) + "\n")

    def text(self):
        return "\n".join(line.raw for line in self.lines) + "\n"

    # -- typed views -----------------------------------------------------

    @property
    def atoms(self):
        """The [ atoms ] records, in file order."""
        if self._atoms is None:
            self._atoms = self._read_atoms()
        return self._atoms

    def _read_atoms(self):
        atoms = []
        for entry in self.lines:
            if entry.directive != "atoms" or not entry.is_data():
                continue
            tok = entry.data
            if len(tok) < 6:
                raise TopologyError(
                    "%s: malformed [ atoms ] line: %r" % (self._where(), entry.raw)
                )
            try:
                mass = float(tok[7]) if len(tok) > 7 else None
            except ValueError:
                mass = None
            atoms.append(
                Atom(
                    nr=int(tok[0]),
                    type=tok[1],
                    resnr=tok[2],
                    residue=tok[3],
                    name=tok[4],
                    cgnr=tok[5],
                    charge=tok[6] if len(tok) > 6 else "",
                    mass=mass,
                    extra=tok[8:],
                    comment=entry.comment,
                )
            )
        if not atoms:
            raise TopologyError("%s: no [ atoms ] directive found" % self._where())
        if [a.nr for a in atoms] != list(range(1, len(atoms) + 1)):
            raise TopologyError(
                "%s: [ atoms ] is not numbered 1..N consecutively" % self._where()
            )
        return atoms

    def _where(self):
        return self.path or "<topology>"

    @property
    def atomtypes(self):
        """The [ atomtypes ] records of a force field, by type name."""
        if self._atomtypes is None:
            self._atomtypes = self._read_atomtypes()
        return self._atomtypes

    def _read_atomtypes(self):
        """Parse [ atomtypes ], whose columns vary between force fields.

        GROMACS accepts the directive with or without an atomic number and
        with or without a bonded type, so counting columns from the left is
        unreliable.  The particle-type column is a lone A, S, V or D, and
        everything else is placed relative to it: sigma and epsilon follow,
        mass and charge precede.
        """
        types = {}
        for entry in self.lines:
            if entry.directive != "atomtypes" or not entry.is_data():
                continue
            tok = entry.data
            where = [i for i, t in enumerate(tok)
                     if t.upper() in PTYPE_CODES and len(t) == 1]
            if not where or len(tok) < 4:
                raise TopologyError(
                    "%s: cannot read [ atomtypes ] line (no particle-type "
                    "column): %r" % (self._where(), entry.raw))
            at = where[0]
            if at < 2 or at + 2 >= len(tok):
                raise TopologyError(
                    "%s: malformed [ atomtypes ] line: %r"
                    % (self._where(), entry.raw))
            name = tok[0]
            # Between the name and the mass sit an optional bonded type and an
            # optional atomic number, in that order; the number is whichever
            # of them parses as one.
            atomic_number = None
            for token in tok[1:max(at - 2, 1)]:
                atomic_number = _maybe_int(token)
                if atomic_number is not None:
                    break
            atom_type = AtomType(
                name=name,
                atomic_number=atomic_number,
                mass=_maybe_float(tok[at - 2]),
                charge=_maybe_float(tok[at - 1]),
                ptype=tok[at].upper(),
                sigma=_maybe_float(tok[at + 1]),
                epsilon=_maybe_float(tok[at + 2]),
                comment=entry.comment,
            )
            types[name] = atom_type
        if not types:
            raise TopologyError("%s: no [ atomtypes ] directive found"
                                % self._where())
        return types

    @property
    def natoms(self):
        return len(self.atoms)

    @property
    def atom_names(self):
        return [a.name for a in self.atoms]

    @property
    def molecule_name(self):
        """The name from [ moleculetype ], or None."""
        for entry in self.lines:
            if entry.directive == "moleculetype" and entry.is_data():
                return entry.data[0]
        return None

    @property
    def residue_names(self):
        """The distinct residue names in [ atoms ], in order of first use."""
        seen = []
        for atom in self.atoms:
            if atom.residue not in seen:
                seen.append(atom.residue)
        return seen

    def residue_name(self, requested=None):
        """The single residue this topology describes.

        Raises when the topology spans several residues and the caller did not
        say which one it means.
        """
        names = self.residue_names
        if requested is not None:
            if requested not in names:
                raise TopologyError(
                    "%s: no residue named %s (has %s)"
                    % (self._where(), requested, ", ".join(names))
                )
            return requested
        if len(names) > 1:
            raise TopologyError(
                "%s: the topology covers several residue names (%s); say which "
                "one to use" % (self._where(), ", ".join(names))
            )
        return names[0]

    def rows(self, directive):
        """The parsed rows of one directive, as (tokens, comment) pairs."""
        return [
            (entry.data, entry.comment)
            for entry in self.lines
            if entry.directive == directive and entry.is_data()
        ]

    def bonds(self):
        """The (i, j) atom index pairs of [ bonds ]."""
        return [(int(t[0]), int(t[1])) for t, _ in self.rows("bonds")]

    def bond_graph(self):
        """Adjacency lists (1-based) built from [ bonds ]."""
        return bond_graph(self.bonds(), self.natoms)

    def connectivity(self):
        """Every atom pair held at a fixed distance, not just [ bonds ].

        Rigid water has no [ bonds ] at all -- its geometry lives in
        [ settles ], and other molecules put bonds in [ constraints ] -- so a
        drawing or a graph traversal that used [ bonds ] alone would show such
        a molecule as unconnected dots.  A [ settles ] row names the first atom
        of a three-site water, whose other two sites follow it.
        """
        pairs = self.bonds()
        pairs.extend((int(t[0]), int(t[1])) for t, _ in self.rows("constraints"))
        for tokens, _ in self.rows("settles"):
            first = int(tokens[0])
            if first + 2 <= self.natoms:
                pairs.extend([(first, first + 1), (first, first + 2)])
        return pairs

    # -- editing ---------------------------------------------------------

    def permute(self, order, cgnr_mode="renumber", sort_rows=True):
        """Return a new ItpFile with the atoms in the order given by AtomOrder.

        Only numbers change.  Functional forms, parameters and the order of the
        atoms *inside* a bonded interaction are left untouched, so impropers
        and other order-sensitive terms keep their meaning.
        """
        if order.natoms != self.natoms:
            raise TopologyError(
                "%s: the atom order covers %d atoms but the topology has %d"
                % (self._where(), order.natoms, self.natoms)
            )
        new_index = order.new_for_old
        atoms = self.atoms
        cgnr_for = _cgnr_mapper(atoms, order, cgnr_mode)

        out = []
        i, n = 0, len(self.lines)
        while i < n:
            entry = self.lines[i]
            directive = entry.directive

            # [ atoms ]: emit the whole block in the new order.
            if directive == "atoms" and entry.is_data():
                while i < n and self.lines[i].directive == "atoms" \
                        and self.lines[i].is_data():
                    i += 1
                for pos, old in enumerate(order.old_for_new, start=1):
                    out.append(atoms[old - 1].format(pos, cgnr_for(old, pos)))
                continue

            # Index-carrying directives: re-index, optionally re-sort.
            spec = DIRECTIVES.get(directive)
            if spec is not None and entry.is_data():
                block = []
                while i < n and self.lines[i].directive == directive \
                        and self.lines[i].is_data():
                    block.append(self.lines[i])
                    i += 1
                rows = [
                    (_remap(item.data, spec, new_index), item.comment)
                    for item in block
                ]
                if sort_rows:
                    rows.sort(key=lambda row: spec.sort_key(row[0]))
                for tokens, comment in rows:
                    text = spec.format_row(tokens)
                    if comment:
                        text += "  " + comment
                    out.append(text)
                continue

            out.append(entry.raw)
            i += 1

        return ItpFile(parse_lines(out), self.path)


def _maybe_float(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _maybe_int(text):
    try:
        return int(text)
    except (TypeError, ValueError):
        return None


def _remap(tokens, spec, new_index):
    out = list(tokens)
    for i in spec.index_positions(tokens):
        out[i] = str(new_index[int(tokens[i])])
    return out


def _cgnr_mapper(atoms, order, mode):
    """Return a function (old index, new position) -> charge group number."""
    if mode == "keep":
        return lambda old, pos: atoms[old - 1].cgnr
    if mode == "per-atom":
        return lambda old, pos: str(pos)
    if mode != "renumber":
        raise TopologyError("unknown charge-group mode %r" % mode)
    # Renumber the existing groups so they still ascend in the new order.
    cgnr_map = {}
    for old in order.old_for_new:
        cg = atoms[old - 1].cgnr
        if cg not in cgnr_map:
            cgnr_map[cg] = str(len(cgnr_map) + 1)
    return lambda old, pos: cgnr_map[atoms[old - 1].cgnr]
