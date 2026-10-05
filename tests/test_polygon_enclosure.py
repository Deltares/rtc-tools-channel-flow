from rtctools_channel_flow.polygon_enclosure import _split_lines

def test_split_lines_adds_intersection_point():
    lines = [
        [(0.0, 0.0), (2.0, 2.0)],  # diagonal
        [(0.0, 2.0), (2.0, 0.0)],  # crossing diagonal
    ]

    result = _split_lines(lines)

    assert result == [
        [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)],
        [(0.0, 2.0), (1.0, 1.0), (2.0, 0.0)],
    ]
    
def test_split_lines_no_intersection():
    lines = [
        [(0.0, 0.0), (1.0, 0.0)],
        [(0.0, 1.0), (1.0, 1.0)],
    ]

    result = _split_lines(lines)

    assert result == lines
    
def test_split_lines_multiple_intersections():
    lines = [
        [(0.0, 0.0), (4.0, 0.0)],  # horizontal
        [(1.0, -1.0), (1.0, 1.0)],  # vertical
        [(3.0, -1.0), (3.0, 1.0)],  # vertical
    ]

    result = _split_lines(lines)

    assert result[0] == [
        (0.0, 0.0),
        (1.0, 0.0),
        (3.0, 0.0),
        (4.0, 0.0),
    ]