"""Writing GROMACS .ndx index files."""

from __future__ import annotations

# GROMACS writes 15 indices to a line; any layout reads back, but matching it
# keeps the files diffable against ones gmx make_ndx wrote.
PER_LINE = 15


def format_ndx(groups):
    """The .ndx text for {group name: 1-based atom indices}."""
    lines = []
    for name, indices in groups.items():
        lines.append("[ %s ]" % name)
        indices = list(indices)
        for k in range(0, len(indices), PER_LINE):
            lines.append(" ".join("%d" % i for i in indices[k:k + PER_LINE]))
        lines.append("")
    return "".join(line + "\n" for line in lines)


def write_ndx(path, groups, append=False):
    """Write (or, with ``append``, add) groups of 1-based atom indices."""
    with open(path, "a" if append else "w") as handle:
        handle.write(format_ndx(groups))
