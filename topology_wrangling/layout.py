"""Laying a molecular graph out in two dimensions.

A Fruchterman-Reingold spring embedding: bonded atoms attract, every pair
repels, and the whole thing cools over a fixed number of iterations.  For the
sparse, mostly tree-shaped graphs a topology gives, that untangles into
something that reads like a structural formula.

The layout is deterministic for a given seed, so re-running a drawing gives the
same picture.
"""

from __future__ import annotations

import math
from collections import defaultdict

from .errors import TopologyError


def _numpy():
    try:
        import numpy
    except ImportError:
        raise TopologyError(
            "the force-directed layout needs numpy; install it, or use the "
            "default --layout skeletal, which needs nothing"
        )
    return numpy


def force_directed(natoms, edges, iterations=400, seed=0, weights=None):
    """Return [(x, y), ...] for atoms 1..natoms, in a unit-ish box centred on 0.

    ``edges`` are 1-based (i, j) pairs.  ``weights`` optionally gives a
    per-edge multiplier on the attraction, so a shorter bond can be asked for.
    """
    np = _numpy()
    if natoms < 1:
        raise TopologyError("cannot lay out a molecule with no atoms")
    if natoms == 1:
        return [(0.0, 0.0)]

    rng = np.random.default_rng(seed)
    # Start on a circle with a little jitter: a symmetric start makes the
    # embedding of a symmetric molecule come out symmetric, and the jitter
    # stops exactly superimposed atoms from giving a zero-length force.
    angles = np.linspace(0, 2 * math.pi, natoms, endpoint=False)
    pos = np.column_stack([np.cos(angles), np.sin(angles)])
    pos += rng.normal(scale=0.01, size=pos.shape)

    if edges:
        index = np.array([(i - 1, j - 1) for i, j in edges], dtype=int)
        for a, b in index.reshape(-1, 2):
            if not (0 <= a < natoms and 0 <= b < natoms):
                raise TopologyError("a bond refers to an atom outside 1..%d"
                                    % natoms)
    else:
        index = np.empty((0, 2), dtype=int)
    strength = np.ones(len(index)) if weights is None \
        else np.asarray(weights, dtype=float)

    area = 1.0
    k = math.sqrt(area / natoms)          # the natural spring length
    temperature = 0.1 * math.sqrt(area)
    cooling = temperature / (iterations + 1)

    for _ in range(iterations):
        delta = pos[:, None, :] - pos[None, :, :]           # (n, n, 2)
        distance = np.linalg.norm(delta, axis=-1)
        np.fill_diagonal(distance, np.inf)                  # no self-force
        distance = np.maximum(distance, 1e-9)

        # Repulsion between every pair: k^2 / d, along the separation.
        force = (delta / distance[..., None]) * (k * k / distance)[..., None]
        displacement = force.sum(axis=1)

        # Attraction along bonds: d^2 / k, pulling the two ends together.
        if len(index):
            a, b = index[:, 0], index[:, 1]
            edge = pos[a] - pos[b]
            length = np.maximum(np.linalg.norm(edge, axis=-1), 1e-9)
            pull = (edge / length[:, None]) * (length * length / k
                                               * strength)[:, None]
            np.add.at(displacement, a, -pull)
            np.add.at(displacement, b, pull)

        # Move, but never further than the temperature allows.
        step = np.maximum(np.linalg.norm(displacement, axis=-1), 1e-9)
        pos += displacement / step[:, None] * np.minimum(step, temperature)[:, None]
        temperature -= cooling

    pos -= pos.mean(axis=0)
    pos = _align(np, pos)
    span = np.abs(pos).max()
    if span > 0:
        pos /= span
    return [tuple(float(v) for v in row) for row in pos]


def _align(np, pos):
    """Rotate so the molecule's long axis is horizontal.

    A spring embedding lands in whatever orientation it started from, which
    wastes a landscape canvas as often as not.  Turning the principal axis of
    the atoms onto x costs nothing and makes the drawing fit its page.
    """
    if len(pos) < 3:
        return pos
    # The principal axis is the eigenvector of the covariance with the larger
    # eigenvalue; eigh returns them in ascending order.
    covariance = np.cov(pos.T)
    if not np.all(np.isfinite(covariance)):
        return pos
    _, vectors = np.linalg.eigh(covariance)
    axis = vectors[:, -1]
    angle = math.atan2(axis[1], axis[0])
    cos, sin = math.cos(-angle), math.sin(-angle)
    return pos @ np.array([[cos, -sin], [sin, cos]]).T


def relax_overlaps(positions, radii, iterations=60, padding=1.06):
    """Push overlapping atoms apart without undoing the layout.

    A layout only knows about points, so it will happily put two atoms close
    enough that their labels collide.  This nudges any pair closer than the sum
    of their radii apart again, in small steps so the shape the layout found
    survives.

    Pure standard library, and bucketed by position so it stays cheap on a big
    molecule: only atoms in neighbouring cells can possibly overlap.
    """
    pos = [[float(x), float(y)] for x, y in positions]
    count = len(pos)
    if count < 2:
        return [(x, y) for x, y in pos]
    reach = [float(r) * padding for r in radii]
    if len(reach) != count:
        raise TopologyError("got %d radii for %d atoms" % (len(reach), count))
    cell = max(max(reach) * 2.0, 1e-9)

    for _ in range(iterations):
        buckets = defaultdict(list)
        for i, (x, y) in enumerate(pos):
            buckets[(int(x // cell), int(y // cell))].append(i)

        shift = [[0.0, 0.0] for _ in range(count)]
        touched = False
        for (cx, cy), members in buckets.items():
            nearby = []
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    nearby.extend(buckets.get((cx + dx, cy + dy), ()))
            for i in members:
                xi, yi = pos[i]
                for j in nearby:
                    if j <= i:
                        continue
                    dx = xi - pos[j][0]
                    dy = yi - pos[j][1]
                    want = reach[i] + reach[j]
                    distance = math.hypot(dx, dy)
                    if distance >= want:
                        continue
                    touched = True
                    if distance < 1e-9:
                        # exactly superimposed: pick a direction, any direction
                        dx, dy, distance = 1.0, 0.0, 1.0
                    push = (want - distance) * 0.5 / distance
                    shift[i][0] += dx * push
                    shift[i][1] += dy * push
                    shift[j][0] -= dx * push
                    shift[j][1] -= dy * push
        if not touched:
            break
        # damped, so a crowded atom does not shoot out of the drawing
        for i in range(count):
            pos[i][0] += shift[i][0] * 0.5
            pos[i][1] += shift[i][1] * 0.5

    mean_x = sum(p[0] for p in pos) / count
    mean_y = sum(p[1] for p in pos) / count
    return [(p[0] - mean_x, p[1] - mean_y) for p in pos]
