"""Reading, editing and writing GROMACS .gro coordinate files.

The fixed-column layout of an atom line is

    cols  0:5   residue number
    cols  5:10  residue name
    cols 10:15  atom name
    cols 15:20  atom number
    cols 20:    coordinates, and optionally velocities

A GroAtom keeps the identity columns as parsed fields and everything from
column 20 on as the verbatim text of the input, so re-ordering and renaming
never perturb a coordinate.  Scripts that need numbers -- a visualiser, say --
get them from the ``position`` and ``velocity`` properties, which parse that
text on demand.
"""

from __future__ import annotations

from .errors import TopologyError

RESID_COLS = (0, 5)
RESNAME_COLS = (5, 10)
ATOMNAME_COLS = (10, 15)
ATOMNR_COLS = (15, 20)
COORD_COL = 20

# GROMACS wraps the residue and atom number columns at this value.
GRO_NUMBER_WRAP = 100000


class GroAtom:
    """One atom line of a .gro file."""

    __slots__ = ("resid", "resname", "name", "nr", "tail")

    def __init__(self, resid, resname, name, nr, tail):
        self.resid = resid
        self.resname = resname
        self.name = name
        self.nr = nr
        self.tail = tail          # cols 20: verbatim, coordinates onwards

    def __repr__(self):
        return "GroAtom(%s%s, %s, %s)" % (self.resid, self.resname, self.name,
                                          self.nr)

    @classmethod
    def parse(cls, line, where=""):
        if len(line) < COORD_COL:
            raise TopologyError("%sshort .gro atom line: %r" % (where, line))
        return cls(
            resid=line[RESID_COLS[0]:RESID_COLS[1]].strip(),
            resname=line[RESNAME_COLS[0]:RESNAME_COLS[1]].strip(),
            name=line[ATOMNAME_COLS[0]:ATOMNAME_COLS[1]].strip(),
            nr=line[ATOMNR_COLS[0]:ATOMNR_COLS[1]].strip(),
            tail=line[COORD_COL:],
        )

    def format(self, name_left=False, resname_left=True):
        """The .gro text for this atom.

        The justification flags reproduce the convention of the file the atom
        came from; see GroFrame.justifies_left.
        """
        return (
            "%5s" % _wrap(self.resid)
            + _name_field(self.resname, resname_left, "residue")
            + _name_field(self.name, name_left, "atom")
            + "%5s" % _wrap(self.nr)
            + self.tail
        )

    # -- numbers, parsed on demand ---------------------------------------

    @property
    def position(self):
        """(x, y, z) in nm."""
        return self._numbers()[0]

    @property
    def velocity(self):
        """(vx, vy, vz) in nm/ps, or None when the file carries no velocities."""
        return self._numbers()[1]

    def _numbers(self):
        values = _split_fixed(self.tail)
        if len(values) == 3:
            return tuple(values), None
        if len(values) == 6:
            return tuple(values[:3]), tuple(values[3:])
        raise TopologyError(
            "cannot read coordinates from %r: expected 3 or 6 numbers, got %d"
            % (self.tail, len(values))
        )


def _split_fixed(tail):
    """Read the coordinate fields of a .gro line.

    .gro is a fixed-width format, but the width varies with the precision the
    writer chose, so we split on whitespace first and only fall back to equal
    width chunks when the fields have run together.
    """
    text = tail.rstrip()
    tokens = text.split()
    if len(tokens) in (3, 6):
        try:
            return [float(t) for t in tokens]
        except ValueError:
            pass
    for count in (3, 6):
        if len(text) % count == 0:
            width = len(text) // count
            chunks = [text[i:i + width] for i in range(0, len(text), width)]
            try:
                return [float(c) for c in chunks]
            except ValueError:
                continue
    raise TopologyError("cannot read coordinates from %r" % tail)


def _wrap(value):
    """Numbers in the .gro index columns wrap at 100000, as GROMACS does."""
    try:
        return "%d" % (int(value) % GRO_NUMBER_WRAP)
    except (TypeError, ValueError):
        return value


def _name_field(name, justify_left, what="atom"):
    if len(name) > 5:
        raise TopologyError(
            "%s name %r does not fit the .gro's 5-column field" % (what, name)
        )
    return ("%-5s" if justify_left else "%5s") % name


def _vote_justification(raw, votes):
    """Tally how a line justifies its residue and atom name fields."""
    for cols in (RESNAME_COLS, ATOMNAME_COLS):
        field = raw[cols[0]:cols[1]]
        if len(field.strip()) >= 5 or not field.strip():
            continue
        if field[0] != " ":
            votes[cols][0] += 1
        elif field[-1] != " ":
            votes[cols][1] += 1


class GroFrame:
    """One frame of a .gro file: title, atoms and box."""

    def __init__(self, title, atoms, box, path=None, justify=None):
        self.title = title
        self.atoms = atoms
        self.box = box
        self.path = path
        # {(lo, hi): left-justified?} for the residue and atom name fields.
        self.justify = justify or {}

    @classmethod
    def read(cls, path):
        with open(path) as handle:
            title = handle.readline().rstrip("\n")
            count_line = handle.readline()
            if not count_line:
                raise TopologyError("%s: truncated, no atom count" % path)
            try:
                natoms = int(count_line.split()[0])
            except (ValueError, IndexError):
                raise TopologyError(
                    "%s: the second line is not an atom count: %r"
                    % (path, count_line.rstrip())
                )
            atoms = []
            votes = {RESNAME_COLS: [0, 0], ATOMNAME_COLS: [0, 0]}
            for i in range(natoms):
                raw = handle.readline()
                if not raw:
                    raise TopologyError(
                        "%s: the file ends after %d of %d atoms"
                        % (path, len(atoms), natoms)
                    )
                raw = raw.rstrip("\n")
                atoms.append(GroAtom.parse(raw, "%s: atom %d: " % (path, i + 1)))
                _vote_justification(raw, votes)
            box = handle.readline().rstrip("\n")
            trailing = handle.read()
        if not box.strip():
            raise TopologyError("%s: missing box line" % path)
        if trailing.strip():
            raise TopologyError(
                "%s: extra content after the box line - multi-frame .gro files "
                "are not supported" % path
            )
        frame = cls(title, atoms, box, path)
        frame.justify = {cols: left > right for cols, (left, right) in votes.items()}
        return frame

    def write(self, path):
        name_left = self.justifies_left(ATOMNAME_COLS)
        res_left = self.justifies_left(RESNAME_COLS)
        with open(path, "w") as handle:
            handle.write(self.title + "\n")
            handle.write("%5d\n" % len(self.atoms))
            for atom in self.atoms:
                handle.write(atom.format(name_left, res_left) + "\n")
            handle.write(self.box + "\n")

    def __len__(self):
        return len(self.atoms)

    def __iter__(self):
        return iter(self.atoms)

    def __getitem__(self, index):
        return self.atoms[index]

    @property
    def resnames(self):
        """The distinct residue names, in order of first appearance."""
        seen = []
        known = set()
        for atom in self.atoms:
            if atom.resname not in known:
                known.add(atom.resname)
                seen.append(atom.resname)
        return seen

    def justifies_left(self, cols=ATOMNAME_COLS):
        """Whether this file left-justifies the name in the given columns.

        Recorded while parsing and reused on output, so a file keeps the
        convention its writer used.  Only names shorter than the 5-column field
        carry any information, so those are the only ones that vote.
        """
        return self.justify.get(tuple(cols), cols == ATOMNAME_COLS)

    def renumber(self, start=1):
        """Give the atoms consecutive numbers, wrapping as GROMACS does."""
        for i, atom in enumerate(self.atoms, start=start):
            atom.nr = str(i % GRO_NUMBER_WRAP)

    def iter_residue_blocks(self, resname, size):
        """Yield (start index, block) for every ``size``-atom ``resname`` molecule.

        A block must be ``size`` contiguous atoms that all carry the residue
        name and the same residue number; anything else is an error rather than
        a silent mis-ordering.  Atoms of other residues are skipped, so the
        caller can walk the whole frame.
        """
        i, total = 0, len(self.atoms)
        while i < total:
            if self.atoms[i].resname != resname:
                i += 1
                continue
            block = self.atoms[i:i + size]
            resid = block[0].resid
            if len(block) != size or any(
                a.resname != resname or a.resid != resid for a in block
            ):
                raise TopologyError(
                    "%s: the %s residue starting at atom %d is not %d "
                    "contiguous atoms of a single residue"
                    % (self._where(), resname, i + 1, size)
                )
            yield i, block
            i += size

    def _where(self):
        return self.path or "<coordinates>"
