"""Periodic boxes and the minimum-image convention, rectangular or triclinic.

GROMACS stores a box as three vectors in lower-triangular form,

    a = (ax,  0,  0)
    b = (bx, by,  0)
    c = (cx, cy, cz)

and a .gro box line lists them as ``ax by cz`` (rectangular) or
``ax by cz ay az bx bz cx cy`` (triclinic), where ``ay``, ``az`` and ``bz``
must be zero.  A zero-length diagonal means that direction is not periodic.

Everything works on numpy arrays of displacements, shape ``(..., 3)``, so a
whole distance matrix is one call.
"""

from __future__ import annotations

import itertools

import numpy as np

from .errors import TopologyError


class Box:
    """A periodic box, as three lower-triangular vectors in nm (rows)."""

    __slots__ = ("vectors", "triclinic", "_periodic", "_safe2", "_shifts")

    def __init__(self, vectors):
        vectors = np.array(vectors, dtype=float).reshape(3, 3)
        a, b, c = vectors
        if a[1] or a[2] or b[2]:
            raise TopologyError(
                "box vectors must be lower-triangular, as GROMACS writes them "
                "(a = (ax, 0, 0), b = (bx, by, 0)); got a = %s, b = %s"
                % (a.tolist(), b.tolist())
            )
        self.vectors = vectors
        self.triclinic = bool(b[0] or c[0] or c[1])
        # Reduce along c, then b, then a: each vector is the only one with a
        # component along its own axis, so this order never undoes a step.
        self._periodic = [(axis, v) for axis, v in ((2, c), (1, b), (0, a))
                          if v[axis] > 0]
        # Every non-zero lattice vector is at least as long as the shortest
        # box height, so a vector no longer than half of it is already the
        # minimum image.  Only longer ones need the neighbour search.
        heights = [v[axis] for axis, v in self._periodic]
        self._safe2 = (min(heights) / 2) ** 2 if heights else np.inf
        # The 26 neighbouring lattice translations, and zero.
        periodic = np.array([v for _, v in self._periodic]).reshape(-1, 3)
        if len(periodic):
            steps = np.array(list(itertools.product((-1, 0, 1),
                                                    repeat=len(periodic))))
            self._shifts = steps @ periodic
        else:
            self._shifts = np.zeros((1, 3))

    @classmethod
    def parse(cls, line, where=""):
        """A Box from the last line of a .gro file."""
        try:
            values = [float(v) for v in line.split()]
        except ValueError:
            raise TopologyError("%scannot read the box line %r" % (where, line))
        if len(values) == 3:
            values += [0.0] * 6
        if len(values) != 9:
            raise TopologyError(
                "%sthe box line has %d numbers; expected 3 or 9: %r"
                % (where, len(values), line)
            )
        ax, by, cz, ay, az, bx, bz, cx, cy = values
        try:
            return cls([(ax, ay, az), (bx, by, bz), (cx, cy, cz)])
        except TopologyError as err:
            raise TopologyError("%s%s" % (where, err))

    def minimum_image(self, d):
        """The shortest periodic images of displacements ``d``, shape (..., 3).

        Exact for any box GROMACS accepts: after the triangular reduction the
        shortest image is the reduced vector or one of its 26 neighbours, which
        are searched only for the vectors long enough for that to matter.
        """
        d = np.array(d, dtype=float)
        for axis, v in self._periodic:
            d -= np.round(d[..., axis] / v[axis])[..., None] * v
        if not self.triclinic:
            return d
        far = np.einsum("...i,...i->...", d, d) > self._safe2
        if far.any():
            candidates = d[far][:, None, :] + self._shifts      # (m, 27, 3)
            best = np.einsum("mki,mki->mk", candidates, candidates).argmin(1)
            d[far] = candidates[np.arange(len(best)), best]
        return d

    def distance2(self, p, q):
        """Squared minimum-image distances between positions ``p`` and ``q``.

        Broadcasts: positions of shape (n, 1, 3) and (1, m, 3) give the
        (n, m) matrix.
        """
        d = self.minimum_image(np.asarray(p, dtype=float)
                               - np.asarray(q, dtype=float))
        return np.einsum("...i,...i->...", d, d)
