"""ChemDraw-style skeletal layout.

A spring embedding treats a molecule as a blob of charges and gives back
something that sprawls; a chemist draws the same molecule with fixed angles,
and that is what makes a structure readable.  This lays a topology out the way
a structural formula is drawn:

  * a chain zigzags -- consecutive bonds sit at a 120 degree interior angle,
    so bond directions alternate 30 degrees either side of the chain axis,
  * a branching atom fans its bonds evenly (three bonds 120 degrees apart,
    four bonds 90 degrees apart),
  * a ring is a regular polygon, and a ring fused to one already placed is
    reflected across the bond they share,
  * hydrogens hang off their heavy atom on short spokes, in whatever
    directions the heavy-atom skeleton leaves free, instead of taking part in
    the chain.

Every bond is the same length, so the drawing is compact and the picture of a
long lipid reads as a long lipid rather than a sprawl.  Pure standard library:
no numpy needed for this path.
"""

from __future__ import annotations

import math
from collections import defaultdict, deque

from .errors import TopologyError

BOND = 1.0                          # every skeleton bond, in layout units
PENDANT_BOND = 0.62                 # a hydrogen sits closer in
CHAIN_TURN = math.radians(60.0)     # 180 - 120: the zigzag of a chain
START_ANGLE = math.radians(30.0)    # so a chain runs left to right
COMPONENT_GAP = 1.2                 # between disconnected molecules


BRANCH_STRATEGIES = ("auto", "spread", "compact")


def skeletal(natoms, edges, pendant=(), branch="compact"):
    """Return [(x, y), ...] for atoms 1..natoms, laid out as a structure.

    ``edges`` are 1-based (i, j) pairs; ``pendant`` names atoms (hydrogens,
    normally) that should hang off their neighbour rather than join the chain.

    Where a branch should go is a judgement call.  "compact" (the default)
    lets a branch carry on forwards, which packs a lipid's two tails roughly
    parallel and is the tightest arrangement for most molecules; "spread"
    sends the longest branch away from everything drawn so far, putting the
    tails in a line instead.  "auto" builds both and keeps whichever scores
    better by _score.
    """
    if natoms < 1:
        raise TopologyError("cannot lay out a molecule with no atoms")
    if branch not in BRANCH_STRATEGIES:
        raise TopologyError("unknown branch strategy %r; use %s"
                            % (branch, ", ".join(BRANCH_STRATEGIES)))
    if branch == "spread":
        return _layout(natoms, edges, pendant, True)
    if branch == "compact":
        return _layout(natoms, edges, pendant, False)
    best, best_score = None, None
    for outward in (True, False):
        points = _layout(natoms, edges, pendant, outward)
        score = _score(points, edges)
        if best_score is None or score < best_score:
            best, best_score = points, score
    return best


def _score(points, edges):
    """The drawing's cost, in the order of priority the tools care about.

    Clashes first -- two atoms on top of each other is a misreading, not an
    aesthetic complaint -- then crossing bonds, then the area of page the
    drawing takes for the atoms it holds.  Area is what turns two splayed
    lipid tails into two parallel ones: parallel chains pack into a fraction
    of the box that a V does, and nothing else in the score prefers the V.
    Area is per atom, so the number means the same thing for a water and for
    a lipid, and a drawing that wastes its page scores badly whatever its size.
    """
    bonded = {(min(i, j), max(i, j)) for i, j in edges}
    position = {i + 1: point for i, point in enumerate(points)}
    return _quality(position, bonded)


def _count_crossings(points, edges):
    """How many pairs of bonds cross, which is what makes a drawing misread."""
    segments = [(points[i - 1], points[j - 1]) for i, j in edges]
    total = 0
    for a in range(len(segments)):
        p1, p2 = segments[a]
        for b in range(a + 1, len(segments)):
            p3, p4 = segments[b]
            if {p1, p2} & {p3, p4}:      # bonds sharing an atom never count
                continue
            if _crosses(p1, p2, p3, p4):
                total += 1
    return total


def _crosses(p1, p2, p3, p4):
    d1 = _side(p3, p4, p1)
    d2 = _side(p3, p4, p2)
    d3 = _side(p1, p2, p3)
    d4 = _side(p1, p2, p4)
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0))


def _side(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _layout(natoms, edges, pendant, outward):

    adjacency = _adjacency(natoms, edges)
    pendants = {a for a in pendant if len(adjacency[a]) == 1}
    skeleton = [a for a in range(1, natoms + 1) if a not in pendants]
    backbone = {a: [n for n in adjacency[a] if n not in pendants]
                for a in skeleton}

    rings = find_rings(backbone)
    ring_of = defaultdict(list)
    for ring in rings:
        for atom in ring:
            ring_of[atom].append(ring)

    position = {}
    offset = 0.0
    for component in _components(backbone):
        local = _place_component(component, backbone, ring_of, rings,
                                 outward)
        _unfold(local, backbone, ring_of)
        _compact(local, backbone, ring_of)
        width = _shift_into_place(local, position, offset)
        offset += width + COMPONENT_GAP

    _place_pendants(adjacency, pendants, position)
    missing = [a for a in range(1, natoms + 1) if a not in position]
    if missing:
        raise TopologyError("could not place atom(s) %s"
                            % ", ".join(str(a) for a in missing[:5]))
    return _normalised(_orient([position[a] for a in range(1, natoms + 1)]))


# --------------------------------------------------------------------------
# the graph
# --------------------------------------------------------------------------

def _adjacency(natoms, edges):
    adjacency = {a: [] for a in range(1, natoms + 1)}
    seen = set()
    for i, j in edges:
        if i == j:
            continue
        for a in (i, j):
            if not 1 <= a <= natoms:
                raise TopologyError(
                    "a bond refers to atom %d, outside 1..%d" % (a, natoms))
        key = (min(i, j), max(i, j))
        if key in seen:
            continue
        seen.add(key)
        adjacency[i].append(j)
        adjacency[j].append(i)
    return adjacency


def _components(adjacency):
    """Connected groups of atoms, each in a stable order."""
    seen = set()
    out = []
    for start in sorted(adjacency):
        if start in seen:
            continue
        group, queue = [], deque([start])
        seen.add(start)
        while queue:
            atom = queue.popleft()
            group.append(atom)
            for neighbour in adjacency[atom]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    queue.append(neighbour)
        out.append(group)
    return out


def find_rings(adjacency):
    """A smallest-set-of-smallest-rings style cycle basis.

    For every bond, the shortest path between its ends that avoids the bond
    itself closes a ring.  Collecting those, shortest first, and keeping only
    the ones that cover a bond not yet in a ring gives the small rings a
    chemist would draw, without a full ring-perception algorithm.
    """
    bonds = {(min(a, b), max(a, b))
             for a in adjacency for b in adjacency[a]}
    candidates = []
    for a, b in sorted(bonds):
        path = _shortest_path(adjacency, a, b, banned=(a, b))
        if path:
            candidates.append(path)
    candidates.sort(key=len)

    basis_size = len(bonds) - len(adjacency) + len(_components(adjacency))
    rings, covered = [], set()
    for ring in candidates:
        if len(rings) >= basis_size:
            break
        ring_bonds = {(min(x, y), max(x, y)) for x, y in _ring_bonds(ring)}
        if ring_bonds - covered:
            rings.append(ring)
            covered |= ring_bonds
    return rings


def _ring_bonds(ring):
    return [(ring[i], ring[(i + 1) % len(ring)]) for i in range(len(ring))]


def _shortest_path(adjacency, start, goal, banned=None):
    """Shortest path start->goal, optionally forbidding one bond."""
    banned = {min(banned), max(banned)} if banned else None
    previous = {start: None}
    queue = deque([start])
    while queue:
        atom = queue.popleft()
        for neighbour in adjacency[atom]:
            if banned and {atom, neighbour} == banned:
                continue
            if neighbour in previous:
                continue
            previous[neighbour] = atom
            if neighbour == goal:
                path, node = [], goal
                while node is not None:
                    path.append(node)
                    node = previous[node]
                return path
            queue.append(neighbour)
    return None


def _farthest(adjacency, start, within):
    """The atom of `within` farthest from start, by bond count."""
    distance = {start: 0}
    queue = deque([start])
    far, best = start, 0
    while queue:
        atom = queue.popleft()
        for neighbour in adjacency[atom]:
            if neighbour in distance or neighbour not in within:
                continue
            distance[neighbour] = distance[atom] + 1
            if distance[neighbour] > best:
                far, best = neighbour, distance[neighbour]
            queue.append(neighbour)
    return far


def _choose_root(component, adjacency):
    """One end of the molecule's longest chain, so the backbone runs across."""
    within = set(component)
    first = _farthest(adjacency, component[0], within)
    return _farthest(adjacency, first, within)


# --------------------------------------------------------------------------
# placing the skeleton
# --------------------------------------------------------------------------

def _place_component(component, adjacency, ring_of, rings, outward=True):
    position = {}
    placed_rings = set()
    root = _choose_root(component, adjacency)
    position[root] = (0.0, 0.0)

    # (atom, direction we arrived from, zigzag sign)
    queue = deque([(root, None, 1)])
    while queue:
        atom, arrival, sign = queue.popleft()

        for ring in ring_of.get(atom, ()):
            key = id(ring)
            if key in placed_rings:
                continue
            placed_rings.add(key)
            grown = _place_ring(ring, adjacency, position, atom, arrival)
            for member in grown:
                if member != atom:
                    queue.append((member, None, 1))

        pending = [n for n in adjacency[atom] if n not in position]
        if not pending:
            continue
        directions, signs = _directions(atom, adjacency, position, arrival,
                                        sign, len(pending))
        pending = _order_by_subtree(atom, pending, adjacency, position)
        if outward:
            directions, signs = _assign_outward(atom, position, directions,
                                                signs)
        x, y = position[atom]
        for neighbour, angle, next_sign in zip(pending, directions, signs):
            position[neighbour] = (x + BOND * math.cos(angle),
                                   y + BOND * math.sin(angle))
            queue.append((neighbour, angle, next_sign))
    return position


def _directions(atom, adjacency, position, arrival, sign, count):
    """Where this atom's unplaced bonds should point, and the next zigzag sign.

    A two-bonded atom continues the zigzag; anything more fans its bonds evenly
    around the atom, which is what makes a branch point look like one.
    """
    degree = len(adjacency[atom])
    if arrival is None:
        used = [_angle_to(position[atom], position[n])
                for n in adjacency[atom] if n in position]
        if not used:
            start = START_ANGLE if degree <= 2 else math.pi / 2
            step = 2 * math.pi / max(degree, 1)
            angles = [start + i * step for i in range(count)]
            return angles, [-1 if i % 2 == 0 else 1 for i in range(count)]
        return _fill_gaps(used, count), [1] * count

    if degree == 2:
        # a chain: turn 60 degrees, alternating, for a 120 degree interior
        angle = arrival + sign * CHAIN_TURN
        return [angle], [-sign]

    back = arrival + math.pi
    step = 2 * math.pi / degree
    free = [back + i * step for i in range(1, degree)]
    free = [a for a in free][:count] if len(free) >= count else \
        _fill_gaps([back], count)
    # keep zigzagging along whichever bond carries on the chain
    signs = [-sign] * len(free)
    return free, signs


def _assign_outward(atom, position, directions, signs):
    """Send the biggest branch into the emptiest direction.

    The callers hand over branches largest first.  Pointing the largest away
    from everything drawn so far keeps two long tails from curling back into
    each other and into the head group they both hang off.
    """
    if len(directions) < 2 or len(position) < 3:
        return directions, signs
    here = position[atom]
    others = [p for a, p in position.items() if a != atom]
    crowd = (sum(p[0] for p in others) / len(others),
             sum(p[1] for p in others) / len(others))
    away = math.atan2(here[1] - crowd[1], here[0] - crowd[0])
    paired = sorted(zip(directions, signs),
                    key=lambda ds: abs(_wrap(ds[0] - away)))
    return [d for d, _ in paired], [s for _, s in paired]


def _wrap(angle):
    """An angle folded into -pi..pi."""
    return (angle + math.pi) % (2 * math.pi) - math.pi


def _order_by_subtree(atom, pending, adjacency, position):
    """Biggest branch first, so the main chain gets the straightest direction."""
    sizes = {}
    for neighbour in pending:
        sizes[neighbour] = _subtree_size(neighbour, atom, adjacency, position)
    return sorted(pending, key=lambda n: (-sizes[n], n))


def _subtree_size(start, blocked, adjacency, position, limit=400):
    seen = {blocked, start}
    queue = deque([start])
    count = 0
    while queue and count < limit:
        atom = queue.popleft()
        count += 1
        for neighbour in adjacency[atom]:
            if neighbour not in seen and neighbour not in position:
                seen.add(neighbour)
                queue.append(neighbour)
    return count


def _fill_gaps(used, count):
    """`count` directions in the widest gaps left by the used ones."""
    if not used:
        step = 2 * math.pi / max(count, 1)
        return [START_ANGLE + i * step for i in range(count)]
    angles = sorted(a % (2 * math.pi) for a in used)
    gaps = []
    for i, angle in enumerate(angles):
        nxt = angles[(i + 1) % len(angles)]
        size = (nxt - angle) % (2 * math.pi)
        gaps.append([size or 2 * math.pi, angle])
    out = []
    for _ in range(count):
        gaps.sort(reverse=True)
        size, start = gaps[0]
        middle = start + size / 2
        out.append(middle)
        gaps[0] = [size / 2, start]
        gaps.append([size / 2, middle])
    return out


def _angle_to(origin, point):
    return math.atan2(point[1] - origin[1], point[0] - origin[0])


# --------------------------------------------------------------------------
# rings
# --------------------------------------------------------------------------

def _place_ring(ring, adjacency, position, anchor, arrival):
    """Place a ring as a regular polygon; return the atoms it placed."""
    order = _ring_order(ring, adjacency)
    if order is None:
        return []
    size = len(order)
    radius = BOND / (2 * math.sin(math.pi / size))

    shared = [a for a in order if a in position]
    centre = _ring_centre(order, position, shared, anchor, arrival, radius)

    # Walk the ring from the anchor, both windings; keep whichever agrees with
    # the atoms already on the page (a fused ring must reuse its shared bond).
    start = order.index(anchor) if anchor in order else 0
    best, best_error = None, None
    for winding in (1, -1):
        rotated = [order[(start + winding * i) % size] for i in range(size)]
        first = position.get(rotated[0], _on_circle(
            centre, radius, _angle_to(centre, position[anchor])
            if anchor in position else math.pi / 2))
        base = _angle_to(centre, first)
        trial = {atom: _on_circle(centre, radius, base + i * 2 * math.pi / size)
                 for i, atom in enumerate(rotated)}
        error = sum(math.dist(trial[a], position[a]) for a in shared
                    if a in trial)
        if best_error is None or error < best_error:
            best, best_error = trial, error

    placed = []
    for atom, point in best.items():
        if atom not in position:
            position[atom] = point
            placed.append(atom)
    return placed


def _ring_centre(order, position, shared, anchor, arrival, radius):
    """Where the polygon's centre goes."""
    fused = _fused_bond(order, shared)
    if fused:
        a, b = (position[fused[0]], position[fused[1]])
        middle = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
        half = math.dist(a, b) / 2
        height = math.sqrt(max(radius * radius - half * half, 0.0))
        along = _angle_to(a, b)
        normal = along + math.pi / 2
        options = [_on_circle(middle, height, normal),
                   _on_circle(middle, height, normal + math.pi)]
        # grow away from whatever is already drawn
        others = [p for atom, p in position.items() if atom not in order]
        if others:
            crowd = (sum(p[0] for p in others) / len(others),
                     sum(p[1] for p in others) / len(others))
            return max(options, key=lambda c: math.dist(c, crowd))
        return options[0]

    if anchor in position:
        outward = arrival if arrival is not None else math.pi / 2
        return _on_circle(position[anchor], radius, outward)
    return (0.0, 0.0)


def _fused_bond(order, shared):
    """Two ring atoms that are already placed and adjacent in this ring."""
    placed = set(shared)
    for a, b in _ring_bonds(order):
        if a in placed and b in placed:
            return (a, b)
    return None


def _ring_order(ring, adjacency):
    """The ring's atoms in the order they are bonded."""
    members = set(ring)
    start = ring[0]
    order = [start]
    previous = None
    current = start
    while len(order) < len(ring):
        nxt = [n for n in adjacency[current]
               if n in members and n != previous and n not in order]
        if not nxt:
            return None
        previous, current = current, nxt[0]
        order.append(current)
    return order


def _on_circle(centre, radius, angle):
    return (centre[0] + radius * math.cos(angle),
            centre[1] + radius * math.sin(angle))


# --------------------------------------------------------------------------
# hydrogens, packing, scaling
# --------------------------------------------------------------------------

def _place_pendants(adjacency, pendants, position):
    """Hang each hydrogen off its heavy atom, in a direction nothing else uses."""
    hosts = defaultdict(list)
    for atom in sorted(pendants):
        hosts[adjacency[atom][0]].append(atom)
    for host, attached in hosts.items():
        if host not in position:
            continue
        used = [_angle_to(position[host], position[n])
                for n in adjacency[host]
                if n in position and n not in pendants]
        x, y = position[host]
        for atom, angle in zip(attached, _fill_gaps(used, len(attached))):
            position[atom] = (x + PENDANT_BOND * math.cos(angle),
                              y + PENDANT_BOND * math.sin(angle))


def _shift_into_place(local, position, offset):
    """Move one component clear of the ones already laid out."""
    if not local:
        return 0.0
    xs = [p[0] for p in local.values()]
    ys = [p[1] for p in local.values()]
    dx = offset - min(xs)
    middle = (min(ys) + max(ys)) / 2
    for atom, (x, y) in local.items():
        position[atom] = (x + dx, y - middle)
    return max(xs) - min(xs)


def _orient(points, steps=180):
    """Turn the whole drawing to the angle whose bounding box is smallest.

    A rotation is a rigid motion, so this costs nothing in geometry: every
    bond length and angle is untouched.  It is worth doing because the layout
    has no reason to end up axis-aligned -- two parallel acyl chains running
    diagonally need a box twice the area of the same two lying flat -- and the
    box is exactly the page the figure is given.  Landscape is preferred at
    equal area, since figures are wider than they are tall.
    """
    if len(points) < 2:
        return points
    best, best_area = points, None
    for step in range(steps):
        angle = math.pi * step / steps
        cos, sin = math.cos(angle), math.sin(angle)
        turned = [(x * cos - y * sin, x * sin + y * cos) for x, y in points]
        xs = [p[0] for p in turned]
        ys = [p[1] for p in turned]
        width = max(xs) - min(xs)
        height = max(ys) - min(ys)
        area = max(width, 1e-9) * max(height, 1e-9)
        if best_area is None or area < best_area - 1e-9:
            best, best_area = turned, area
    xs = [p[0] for p in best]
    ys = [p[1] for p in best]
    if (max(ys) - min(ys)) > (max(xs) - min(xs)):
        best = [(-y, x) for x, y in best]      # stand it up the other way
    return best


def _normalised(points):
    """Centre the drawing and scale it into a unit-ish box, shape unchanged."""
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    cx = (min(xs) + max(xs)) / 2
    cy = (min(ys) + max(ys)) / 2
    points = [(x - cx, y - cy) for x, y in points]
    span = max(max(abs(x) for x, _ in points),
               max(abs(y) for _, y in points), 1e-9)
    return [(x / span, y / span) for x, y in points]


# --------------------------------------------------------------------------
# unfolding a branch that landed on top of something
# --------------------------------------------------------------------------

CLASH = 1.15 * BOND        # non-bonded atoms closer than this are too close
UNFOLD_PASSES = 12


def _compact(position, adjacency, ring_of, budget=2500):
    """Fold branches inwards until the drawing stops getting smaller.

    Unfolding only looks at bonds between atoms that clash, so once nothing
    clashes it stops -- and a clash-free drawing can still waste most of its
    page, which is what two lipid tails splayed into a V do.  This tries the
    same flip on every bond that has a branch beyond it, biggest branch first
    since those move the most, and keeps whatever shrinks the page.
    """
    bonded = {(min(a, b), max(a, b))
              for a in adjacency for b in adjacency[a] if a in position}
    best = _quality(position, bonded)
    candidates = _flippable(position, adjacency, ring_of, branch_only=True)
    spent = 0
    improved = True
    while improved and spent < budget:
        improved = False
        for bond, subtree in candidates:
            if spent >= budget:
                break
            spent += 1
            _reflect(position, bond, subtree)
            score = _quality(position, bonded)
            if score < best - 1e-9:
                best = score
                improved = True
            else:
                _reflect(position, bond, subtree)   # undo; it did not help
    return best


def _flippable(position, adjacency, ring_of, branch_only=False):
    """Every bond with a branch beyond it, biggest branch first.

    ``branch_only`` keeps to bonds that touch a branch point.  Reflecting one
    of those swings a whole substituent round, which is a move a chemist makes.
    Reflecting a bond in the middle of a chain instead reverses the chain and
    folds it into a hairpin -- that can shrink the page, since the fold tucks
    into space the rest of the molecule already occupies, but an acyl chain
    drawn doubled back reads as wrong.  So compaction is allowed the first kind
    and not the second; clash removal is allowed both, because a clash is worse.
    """
    out = []
    seen = set()
    for a in position:
        for b in adjacency[a]:
            if b not in position:
                continue
            key = (min(a, b), max(a, b))
            if key in seen or _in_a_ring((a, b), ring_of):
                continue
            seen.add(key)
            if branch_only and len(adjacency[a]) < 3 and len(adjacency[b]) < 3:
                continue
            subtree = _beyond((a, b), adjacency, position)
            if subtree and len(subtree) < len(position):
                out.append(((a, b), subtree))
    out.sort(key=lambda item: -len(item[1]))
    return out


def _unfold(position, adjacency, ring_of):
    """Flip branches that clash until nothing improves.

    Which way a bond turns is decided locally as the walk passes through it,
    but whether that choice clashes is only known once the rest of the
    molecule is down -- a lipid's head group folds back onto the glycerol it
    hangs off, and nothing at the moment of placing it could have said so.

    So: find the atoms that are too close, and for each such pair look at the
    bonds on the path between them.  Reflecting everything beyond one of those
    bonds across its own axis swings that branch to the other side while
    leaving every bond length and every angle exactly as it was, because a
    reflection is a rigid motion.  Keep the flip that helps most, and repeat.
    """
    bonded = {(min(a, b), max(a, b))
              for a in adjacency for b in adjacency[a] if a in position}
    best = _quality(position, bonded)
    for _ in range(UNFOLD_PASSES):
        candidate = _best_flip(position, adjacency, ring_of, bonded, best)
        if candidate is None:
            return
        score, bond, subtree = candidate
        _reflect(position, bond, subtree)
        best = score


def _best_flip(position, adjacency, ring_of, bonded, current):
    """The single bond flip that most improves the drawing, if any does."""
    best = None
    for bond in _relieving_bonds(position, adjacency, ring_of, bonded):
        subtree = _beyond(bond, adjacency, position)
        if subtree is None:
            continue
        _reflect(position, bond, subtree)
        score = _quality(position, bonded)
        _reflect(position, bond, subtree)        # a reflection undoes itself
        if score < current - 1e-9 and (best is None or score < best[0]):
            best = (score, bond, subtree)
    return best


def _relieving_bonds(position, adjacency, ring_of, bonded):
    """Bonds that could separate a clashing pair: those on the path between.

    A bond inside a ring is not a candidate -- reflecting across it would tear
    the ring open.
    """
    out = []
    seen = set()
    for a, b in _clashes(position, bonded):
        path = _shortest_path(adjacency, a, b)
        if not path:
            continue
        for i in range(len(path) - 1):
            bond = (path[i], path[i + 1])
            key = (min(bond), max(bond))
            if key in seen:
                continue
            seen.add(key)
            if _in_a_ring(bond, ring_of):
                continue
            out.append(bond)
    return out


def _in_a_ring(bond, ring_of):
    a, b = bond
    for ring in ring_of.get(a, ()):
        if b in ring:
            return True
    return False


def _beyond(bond, adjacency, position):
    """Everything reachable past a bond; None when it is not a bridge."""
    a, b = bond
    seen = {a, b}
    queue = deque([b])
    group = [b]
    while queue:
        atom = queue.popleft()
        for neighbour in adjacency[atom]:
            if neighbour in seen or neighbour not in position:
                continue
            if neighbour == a:                   # a cycle: not a bridge
                return None
            seen.add(neighbour)
            group.append(neighbour)
            queue.append(neighbour)
    return group


def _reflect(position, bond, subtree):
    """Mirror a subtree in the line through its bond, a rigid motion."""
    a, b = bond
    ax, ay = position[a]
    bx, by = position[b]
    dx, dy = bx - ax, by - ay
    length = math.hypot(dx, dy)
    if length < 1e-12:
        return
    dx, dy = dx / length, dy / length
    for atom in subtree:
        px, py = position[atom]
        vx, vy = px - ax, py - ay
        along = vx * dx + vy * dy
        position[atom] = (ax + 2 * along * dx - vx, ay + 2 * along * dy - vy)


def _clashes(position, bonded):
    """Non-bonded atoms sitting closer than a bond length."""
    cell = CLASH
    buckets = defaultdict(list)
    for atom, (x, y) in position.items():
        buckets[(int(x // cell), int(y // cell))].append(atom)
    out = []
    for (cx, cy), members in buckets.items():
        nearby = []
        for ox in (-1, 0, 1):
            for oy in (-1, 0, 1):
                nearby.extend(buckets.get((cx + ox, cy + oy), ()))
        for a in members:
            for b in nearby:
                if b <= a or (a, b) in bonded:
                    continue
                if math.dist(position[a], position[b]) < CLASH:
                    out.append((a, b))
    return out


# Clashes outrank crossings, crossings outrank a roomy page, and the page is
# measured per atom so the weights mean the same thing for any molecule.
#
# There is deliberately no term for the shape of the page.  A straight alkane
# is thin, and it is also exactly how one is drawn; penalising that folds the
# chain in half to square the page up, which costs area and looks wrong.  Wide
# and thin is only a problem when it is also empty, and area already says so.
CLASH_WEIGHT = 400.0
CROSSING_WEIGHT = 40.0


def _quality(position, bonded):
    """Lower is better.  See _score for what the terms mean."""
    penalty = 0.0
    for a, b in _clashes(position, bonded):
        gap = CLASH - math.dist(position[a], position[b])
        penalty += gap * gap

    order = sorted(position)
    points = [position[a] for a in order]
    index = {a: i for i, a in enumerate(order)}
    edges = [(index[a] + 1, index[b] + 1) for a, b in bonded]

    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    width = max(max(xs) - min(xs), BOND)
    height = max(max(ys) - min(ys), BOND)
    area_per_atom = width * height / (len(points) * BOND * BOND)

    return (CLASH_WEIGHT * penalty
            + CROSSING_WEIGHT * _count_crossings(points, edges)
            + area_per_atom)
