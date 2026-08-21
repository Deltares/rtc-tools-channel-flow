from math import atan2, pi, sqrt

TWO_PI = 2.0 * pi


class DeadEndError(Exception):
    pass


def enclosing_segments(point, lines, return_lines=False):
    lines = _split_lines(lines)

    # Get all segments in all lines
    line_segments = [
        (s, i) for i, points in enumerate(lines) for s in zip(points[:-1], points[1:])
    ]

    # TODO: We are assuming unique segments here. It should be possible for a
    # segment to be part of more than one line.
    segment_to_line = {s: i for s, i in line_segments}

    # All segments twice, with start and end points reversed.
    segment_to_line.update({tuple(reversed(s)): i for s, i in line_segments})

    # Process on segments only. Mapping to original line number will be done
    # later.
    segments = [s for s, i in line_segments]

    # Some segments might be zero length with the start and ending point
    # equal. Get rid of them, as we cannot derive any direction from them (and
    # do not need them anyway).
    segments = [(s, e) for s, e in segments if s != e]

    segment_dict = {}

    for s in segments:
        # Original segment
        segment_dict.setdefault(s[0], []).append(s)

        # Flipped segment
        flip_s = tuple(reversed(s))
        segment_dict.setdefault(flip_s[0], []).append(flip_s)

    # Find closest segment/point
    sorted_segments = sorted(
        segments, key=lambda x: _distance_point_to_segment(point, x)
    )

    # Determine in which way we should traverse the closest segment. If we
    # cannot figure out the direction, we pick the next closest segment
    # instead.
    for closest_segment in sorted_segments:
        (x1, y1), (x2, y2) = closest_segment
        (x0, y0) = point
        angle_start = atan2((y1 - y0).item(), (x1 - x0).item())
        angle_end = atan2((y2 - y0).item(), (x2 - x0).item())
        angle_diff = (angle_end - angle_start) % TWO_PI

        if angle_diff > pi:
            closest_segment = tuple(reversed(closest_segment))

        if angle_diff != 0:
            break

    enclosing_segments = []
    enclosing_segments.append(closest_segment)

    # Start adding segments until we get back to our start point again.
    start_point = enclosing_segments[0][0]
    prev_segment = closest_segment
    prev_end_point = closest_segment[1]

    while prev_end_point != start_point:
        results = segment_dict[prev_end_point]

        # Remove the reverse of the last added segment in the possibilities
        results = [x for x in results if x != tuple(reversed(prev_segment))]

        if not results:
            raise DeadEndError(
                "Could not find another segment starting from {}".format(prev_end_point)
            )

        # Pick the segment which goes most counter-clockwise
        cur_segment_angle = _segment_angle(prev_segment)
        results = sorted(
            results,
            key=lambda x: (_segment_angle(x) - cur_segment_angle + pi) % TWO_PI,
            reverse=True,
        )  # noqa B023
        next_segment = results[0]

        enclosing_segments.append(next_segment)

        prev_end_point = next_segment[1]
        prev_segment = next_segment

    # Select lines from the enclosing segments, based on their source line
    enclosing_lines = []
    source_lines = []

    i_start, i_stop = 0, None  # [i_start, i_stop)

    prev_line_no = segment_to_line[enclosing_segments[0]]

    for i, s in enumerate(enclosing_segments):
        line_no = segment_to_line[s]

        if line_no != prev_line_no:
            i_stop = i
            enclosing_lines.append(enclosing_segments[i_start:i_stop])
            source_lines.append(prev_line_no)
            i_start = i_stop
            prev_line_no = line_no

    # Last slice
    enclosing_lines.append(enclosing_segments[i_start:])
    source_lines.append(line_no)

    # Check if last line and first line belong to the same source line. If so, merge them
    if source_lines[0] == source_lines[-1]:
        enclosing_lines[-1].extend(enclosing_lines.pop(0))

    # Convert segments to points
    for i, line in enumerate(enclosing_lines):
        points = [s[0] for s in line]
        points.append(line[-1][1])
        enclosing_lines[i] = points

    if return_lines:
        return enclosing_segments, enclosing_lines
    else:
        return enclosing_segments


def _segment_angle(segment):
    """
    Calculates the angle of the segment, with the start point as the origin.
    The returned angle is in the range [0, 2*pi)
    """
    (x1, y1), (x2, y2) = segment
    return atan2(y2 - y1, x2 - x1) % TWO_PI


def _split_lines(lines):
    """
    Input is a list of lines, each consisting of a list of points: [[(x1, y1),
    (x2, y2), ...], ...]

    The output is a similar structure, but with intersection points between
    all lines added.
    """

    new_lines = []

    for points_i in lines:
        new_line = [points_i[0]]

        for segment in zip(points_i[:-1], points_i[1:]):
            # Note that we also check other segments of the current line, as we
            # have to account for the possibility that the line intersects itself.

            other_segments = (
                s for points_j in lines for s in zip(points_j[:-1], points_j[1:])
            )

            new_points = _split_segment(segment, other_segments)

            new_line.extend(new_points[1:])  # Start point already was in the list

        new_lines.append(new_line)

    return new_lines


def _split_segment(segment, other_segments):
    """
    If segment goes from (x_start, y_start) to (x_end, y_end), this function
    will return a list [(x_start, y_start), (x2, y2), ..., (x_end, y_end)].
    Each new point is a point where the segment intersects with any of the
    other segments.
    """

    def _sort_distance_points(a, b):
        # Just for sorting, so no need to take the sqrt
        xa, ya = a
        xb, yb = b

        return (xb - xa) ** 2 + (yb - ya) ** 2

    start_point, end_point = segment

    intersection_points = set()

    for s in other_segments:
        p = _segment_intersection(segment, s)

        if p is not None:
            intersection_points.add(p)

    intersection_points.add(start_point)
    intersection_points.add(end_point)

    ret = sorted(
        intersection_points, key=lambda x: _sort_distance_points(start_point, x)
    )

    return ret


def _point_on_segment(point, segment):
    """
    Checks whether a point lies on a particular segment.
    """
    x0, y0 = point

    (x1, y1), (x2, y2) = segment

    x_min, x_max = min(x1, x2), max(x1, x2)
    y_min, y_max = min(y1, y2), max(y1, y2)

    return (x_min <= x0 <= x_max) and (y_min <= y1 <= y_max)


def _general_equation_form(segment):
    """
    Calculates the general equation form (ax + by + c) of the line of which
    this segment is part.
    """
    (x1, y1), (x2, y2) = segment

    a = y1 - y2
    b = x2 - x1
    c = x1 * y2 - x2 * y1

    return a, b, c


def _distance_point_to_segment(point, segment):
    """
    Calculates the distance from a point to a line segment.
    """
    x0, y0 = point
    (x1, y1), (x2, y2) = segment

    # Zero-length segments sometimes occur in the segment list of e.g. contour
    # plots. Although "segment" is then somewhat of a misnomer, we do want to
    # handle such cases transparently. We therefore terminate early, to avoid
    # divisions by zero.
    if (x1, y1) == (x2, y2):
        return sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)

    # Formulate segment as ax + by + c = 0.
    a, b, c = _general_equation_form(segment)

    # Calculate position of closest point from point to line through segment
    # See https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line for proof/derivation
    xt = (b * (b * x0 - a * y0) - a * c) / (a**2 + b**2)
    yt = (a * (-b * x0 + a * y0) - b * c) / (a**2 + b**2)

    # We check whether this point is on the segment. If not, we take the
    # minimum distance to either of the end points as the distance from the
    # point to this segment.
    if _point_on_segment((xt, yt), segment):
        distance = sqrt(((xt - x0) ** 2).item() + ((yt - y0) ** 2).item())
    else:
        d_1 = sqrt((x1 - x0) ** 2 + ((y1 - y0) ** 2).item())
        d_2 = sqrt((x2 - x0) ** 2 + ((y2 - y0) ** 2).item())
        distance = min(d_1, d_2)

    return distance


def _segment_intersection(segment_1, segment_2):
    """
    Calculates the point of intersection between two segments. If lines are
    parallel, or if the intersection is not on both segments, None is
    returned.
    """

    # Intersection of two segments
    a1, b1, c1 = _general_equation_form(segment_1)
    a2, b2, c2 = _general_equation_form(segment_2)

    div = a1 * b2 - a2 * b1
    if div != 0:
        yt = (c1 * a2 - c2 * a1) / div
        xt = (b1 * c2 - c1 * b2) / div
    else:
        # lines are parallel, and do not intersect
        return None

    if _point_on_segment((xt, yt), segment_1) and _point_on_segment(
        (xt, yt), segment_2
    ):
        return (xt, yt)
    else:
        return None
