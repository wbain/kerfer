import svgelements
import math
from typing import Optional
import logging

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(levelname)s %(message)s")
for h in logging.getLogger().handlers:
    if h.formatter is not None:
        h.formatter.default_msec_format = '%s.%03d'
logger = logging.getLogger(__name__)


##########################################################################################
# Analyze winding

def remove_dupes(points: list[svgelements.Point]) -> list[svgelements.Point]:
    """
    Remove duplicate consecutive points
    Args:
        points (list): A list of svgelements.Point objects.
    Returns:
        list: A list of svgelements.Point objects with consecutive duplicates removed.
    """
    unique_points = []
    if points:
        unique_points.append(points[0])
        for i in range(1, len(points)):
            if points[i] != points[i-1]:
                unique_points.append(points[i])
    return unique_points


def get_signed_area(points: list[svgelements.Point]) -> float:
    """
    Calculates the signed area of a polygon defined by a list of points.
    A positive area generally indicates counter-clockwise winding,
    while a negative area indicates clockwise winding.
    Args:
        points (list): A list of svgelements.Point objects defining the polygon vertices.
    Returns:
        float: The signed area of the polygon.
    """
    area = 0.0
    for i in range(len(points)):
        p1 = points[i]
        p2 = points[(i + 1) % len(points)]  # Wrap around to the first point
        if p1.x is not None and p1.y is not None and p2.x is not None and p2.y is not None:
            area += (p1.x * p2.y) - (p2.x * p1.y)
    return area / 2.0


def get_subpath_area(subpath: svgelements.Subpath):
    """
    Calculates the signed area of a subpath.
    Args:
        subpath (svgelements.Subpath): The subpath to analyze.
    Returns:
        float or None: The signed area of the subpath, or None if undetermined.
    """
    # Extract points from the subpath segments
    # logger.debug(f"  Subpath")  # : {subpath.d()}
    points = []
    for segment in subpath:
        if isinstance(segment, (svgelements.Line, svgelements.Close)):
            points.append(segment.start)
            points.append(segment.end)
        elif isinstance(segment, (svgelements.CubicBezier, svgelements.QuadraticBezier, svgelements.Arc)):
            # For curved segments, a more complex integration or
            # approximation (e.g., by sampling points) would be needed
            # for accurate signed area calculation.
            # For simplicity, we'll just consider the start and end points here.
            points.append(segment.start)
            points.append(segment.end)
    
    unique_points = remove_dupes(points)

    if len(unique_points) >= 3:  # A polygon needs at least 3 unique points
        signed_area = get_signed_area(unique_points)
        return signed_area
    else:
        return None
    

def calculate_is_subpath_clockwise(subpath: svgelements.Subpath) -> Optional[bool]:
    """
    Calculates if the subpath is clockwise.
    Args:
        subpath (svgelements.Subpath): The subpath to analyze.
    Returns:
        True if clockwise, False if counter-clockwise, None if undetermined.
    """
    signed_area = get_subpath_area(subpath)
    if signed_area is None:
        logger.warning("    Subpath has too few points to determine winding order.")
        return None
    else:
        return signed_area < 0


def calculate_is_path_clockwise(path: svgelements.Path) -> Optional[bool]:
    """
    Calculates if the first subpath in the path is clockwise.
    Args:
        path (svgelements.Path): The path to analyze.
    Returns:
        True if clockwise, False if counter-clockwise, None if undetermined.
    """
    # logger.debug(f"Path")  # : {path.d()}
    for subpath in path.as_subpaths():  # really only looks at the *first* subpath
        return calculate_is_subpath_clockwise(subpath)


##########################################################################################
# Opening/Closing

def show_path_open_closed(path: svgelements.Path):
    """
    Shows whether subpaths in the path are open or closed.
    Args:
        path (svgelements.Path): The path object to process.
    """
    logger.debug(f"  Showing open/closed for path \"{path.id}\"")
    subpath_idx = -1
    for subpath in path.as_subpaths():
        subpath_idx += 1

        # iterate to find start and final segment of subpath
        subpath_start_seg = None
        for seg in subpath:
            if subpath_start_seg is None:
                subpath_start_seg = seg
            subpath_final_seg = seg

        if isinstance(subpath_final_seg, svgelements.Close):
            logger.debug(f"    Subpath {subpath_idx} closed")
            # logger.debug(f"        {subpath.d()}")
        elif isinstance(subpath_final_seg, svgelements.Line):
            logger.debug(f"    Subpath {subpath_idx} open")
            # logger.debug(f"      {subpath.d()}")
        else:
            logger.warning(f"    Subpath {subpath_idx} final segment is neither a line or a close!")


def show_svg_open_closed(svg: svgelements.SVG):
    """
    Shows whether subpaths in the SVG are open or closed.
    Args:
        svg (svgelements.SVG): The SVG object to process.
    """
    logger.debug("Showing open/closed for SVG")
    for element in svg.elements():
        if isinstance(element, svgelements.Path):
            show_path_open_closed(element)


def open_path(path: svgelements.Path):
    """
    Opens all subpaths in the SVG, in-place, by replacing close commands with line segments.
    Args:
        path (svgelements.Path): The path object to process.
    """
    logger.debug(f"  Opening path \"{path.id}\"")
    subpath_idx = -1
    for subpath in path.as_subpaths():
        subpath_idx += 1

        found_close = False

        subpath_seg_idx = -1
        for segment in subpath:
            subpath_seg_idx += 1

            if isinstance(segment, svgelements.Close):
                found_close = True
                logger.debug(f"    Opening subpath {subpath_idx}")
                # replace close with line segment
                start = segment.start
                end = segment.end
                path_seg_idx = subpath.index_to_path_index(subpath_seg_idx)
                del path[path_seg_idx]
                path.insert(path_seg_idx, svgelements.Line(start, end))
            
        if not found_close:
            logger.debug(f"    Subpath {subpath_idx} already open")


def open_svg(svg: svgelements.SVG):
    """
    Opens all subpaths in the SVG, in-place, by replacing close commands with line segments.
    Args:
        svg (svgelements.SVG): The SVG object to process.
    """
    logger.debug("Opening SVG")
    for element in svg.elements():
        if isinstance(element, svgelements.Path):
            open_path(element)


def close_path(path: svgelements.Path):
    """
    Closes all subpaths in the SVG, in-place, by replacing final line segments with close commands.
    Args:
        path (svgelements.Path): The path object to process.
    """
    logger.debug(f"  Closing path \"{path.id}\"")
    # logger.debug(f"  {path.d()}")
    subpath_idx = -1
    for subpath in path.as_subpaths():
        subpath_idx += 1

        logger.debug(f"    Closing subpath {subpath_idx}")
        # logger.debug(f"      {subpath.d()}")

        num_segs = len(subpath) # includes initial move
        seg_end_path_idx = subpath.index_to_path_index(num_segs - 1)

        # iterate to find start and final segments of subpath
        subpath_start_seg = None
        for seg in subpath:
            # logger.debug(f'  seg: {seg}, is line? {isinstance(seg, svgelements.Line)}')
            if subpath_start_seg is None:
                subpath_start_seg = seg
            subpath_final_seg = seg

        assert subpath_final_seg == path[seg_end_path_idx], "Final segment mismatch!"

        # logger.debug(f'subpath start seg: {subpath_start_seg}')
        # logger.debug(f'subpath final seg: {subpath_final_seg}')

        if not isinstance(subpath_start_seg, svgelements.Move):
            logger.warning(f"      Subpath {subpath_idx} doesn't start with a move!")
            continue

        if isinstance(subpath_final_seg, svgelements.Close):
            logger.debug(f"      Subpath {subpath_idx} already closed")
            continue
        
        if not isinstance(subpath_final_seg, svgelements.Line):
            logger.warning(f"      Subpath {subpath_idx} final segment is neither a line nor a close!")
            continue

        final_seg_start_point = subpath_start_seg.end     # end coordinate of the initial move item
        final_seg_end_point = subpath_final_seg.end

        if final_seg_end_point != final_seg_start_point:
            logger.warning(f"      Subpath {subpath_idx} end point does not match start point, cannot close!")
            continue

        try:
            del path[seg_end_path_idx]
            path.insert(seg_end_path_idx, svgelements.Close())
        except Exception as e:
            logger.error(f"      Failed to close subpath {subpath_idx}! {e}")
            # logger.error(f"        {subpath.d()}")



def close_svg(svg: svgelements.SVG):
    """
    Closes all subpaths in the SVG, in-place, by replacing final line segments with close commands.
    Args:
        svg (svgelements.SVG): The SVG object to process.
    """
    logger.debug(f"Closing SVG")
    for element in svg.elements():
        if isinstance(element, svgelements.Path):
            close_path(element)


##########################################################################################
# Offsetting

def offset_endpoints(p1: svgelements.Point, p2: svgelements.Point, offset_distance: float):
    """
    Offsets a line segment by a given perpendicular distance.

    Args:
        p1: The first point of the segment.
        p2: The second point of the segment.
        offset_distance (float): The perpendicular distance to offset the line.

    Returns:
        tuple: A tuple containing the (x, y) coordinates of the two new points
               forming the offset line segment.
    """
    if (not p1.y) or (not p2.y):
        logger.warning("offset_endpoints: at least one of the points has no y value!")
        return p1, p2
    
    # x1, y1 = p1
    # x2, y2 = p2
    x1 = p1.x
    y1 = p1.y
    x2 = p2.x
    y2 = p2.y

    # Calculate the vector representing the line segment
    dx = x2 - x1
    dy = y2 - y1

    # Calculate the length of the line segment
    length = math.sqrt(dx**2 + dy**2)

    if length == 0:  # Handle the case of a zero-length segment (a point)
        return p1, p2

    # Calculate the normalized perpendicular vector
    # (rotated 90 degrees clockwise for one side, counter-clockwise for the other)
    # For clockwise: (dy, -dx) / length
    # For counter-clockwise: (-dy, dx) / length
    # We'll use the counter-clockwise direction for positive offset_distance
    # and clockwise for negative offset_distance
    perp_dx = -dy / length
    perp_dy = dx / length

    # Calculate the offset vector
    offset_vec_x = perp_dx * offset_distance
    offset_vec_y = perp_dy * offset_distance

    # Calculate the new points
    new_p1_x = x1 + offset_vec_x
    new_p1_y = y1 + offset_vec_y
    new_p2_x = x2 + offset_vec_x
    new_p2_y = y2 + offset_vec_y

    return svgelements.Point(new_p1_x, new_p1_y), svgelements.Point(new_p2_x, new_p2_y)


def line_segment_intersection(p1: svgelements.Point, p2: svgelements.Point, p3: svgelements.Point, p4: svgelements.Point):
    """
    Calculates the intersection point of two line segments. First segment is p1 to p2, second is p3 to p4.

    Args:
        p1: Start point of the first line segment (x1, y1).
        p2: End point of the first line segment (x2, y2).
        p3: Start point of the second line segment (x3, y3).
        p4: End point of the second line segment (x4, y4).

    Returns:
        tuple or None: The intersection point (x, y) if it exists, otherwise None.
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4

    if (x1 is None) or (x2 is None) or (x3 is None) or (x4 is None) or \
            (y1 is None) or (y2 is None) or (y3 is None) or (y4 is None):
        logger.warning("      line_segment_intersection: at least one of the coordinates is None!")
        return None

    # Calculate the denominator of the intersection formula
    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

    # If the denominator is zero, lines are parallel or collinear
    if denom == 0:
        logger.warning("      lines are parallel or collinear")
        return None

    # Calculate t and u values for the parametric equations
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
    u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom

    intersection_x = x1 + t * (x2 - x1)
    intersection_y = y1 + t * (y2 - y1)
    return (intersection_x, intersection_y)


def offset_line_segment(seg: svgelements.Linear, offset_distance: float) -> svgelements.Line:
    """
    Offsets a line segment (or a close) by a given perpendicular distance. Returns a new line segment.
    Args:
        seg (svgelements.Linear): The line segment to offset.
        offset_distance (float): The perpendicular distance to offset the line.
    Returns:
        svgelements.Line: A new line segment representing the offset line.
    """
    endpoints = offset_endpoints(seg.start, seg.end, offset_distance)
    return svgelements.Line(endpoints[0], endpoints[1])


def offset_subpath(subpath: svgelements.Subpath, offset_dist: float):
    """
    Offsets a subpath by a given distance. Modifies the subpath in-place.
    Args:
        subpath (svgelements.Subpath): The subpath to offset.
        offset_dist (float): The distance to offset the subpath.
    """
    # First, go around the loop and offset each line segment
    off_segs = []
    for segment in subpath:
        if isinstance(segment, svgelements.Line):
            # logger.debug(f"      Offsetting line segment")
            # logger.debug(f"        {segment.start.x} {segment.start.y} -> {segment.end.x} {segment.end.y}")
            off_seg = offset_line_segment(segment, offset_dist)
            # logger.debug(f"        {off_seg.start.x} {off_seg.start.y} -> {off_seg.end.x} {off_seg.end.y}")
            off_segs.append(off_seg)

    # Next, find the intersections of each offset segment with the next
    intersections = []
    seg_idx = -1
    for segment in off_segs:
        seg_idx += 1

        next_seg = off_segs[(seg_idx + 1) % len(off_segs)]
        intersection_point = line_segment_intersection(segment.start, segment.end, next_seg.start, next_seg.end)

        if intersection_point is not None:
            intersections.append(intersection_point)
        else:
            intersections.append(segment.end)

    # Finally, update the original subpath segments to use the intersection points
    seg_end_intersection_idx = 0
    for segment in subpath:
        if isinstance(segment, svgelements.Move):
            if segment.end is not None:
                # initial move command goes to final intersection point
                segment.end.x = intersections[-1][0]
                segment.end.y = intersections[-1][1]
        elif isinstance(segment, svgelements.Line):
             # index of previous intersection
            seg_start_intersection_idx = (seg_end_intersection_idx + len(intersections) - 1) % len(intersections)
            
            if segment.start is not None:
                segment.start.x = intersections[seg_start_intersection_idx][0]
                segment.start.y = intersections[seg_start_intersection_idx][1]

            if segment.end is not None:
                segment.end.x = intersections[seg_end_intersection_idx][0]
                segment.end.y = intersections[seg_end_intersection_idx][1]

            seg_end_intersection_idx += 1


def offset_path(path: svgelements.Path, offset_dist: float):
    """
    Offsets a path by a given distance. Modifies the subpaths in-place.
    Args:
        path (svgelements.Path): The path to offset.
        offset_dist (float): The distance to offset the path.
    """
    is_clockwise = calculate_is_path_clockwise(path)
    dist = offset_dist if is_clockwise else -offset_dist
    logger.debug(f"  Offsetting path \"{path.id}\": {'clockwise' if is_clockwise else 'counter-clockwise'} => {dist}") # : {path.d()}
    for subpath in path.as_subpaths():
        offset_subpath(subpath, dist)


def offset_svg(svg: svgelements.SVG, offset_dist: float):
    """
    Offsets all the paths in an SVG by a given distance. Modifies the paths in-place.
    Args:
        svg (svgelements.SVG): The SVG to offset.
        offset_dist (float): The distance to offset the SVG.
    """
    logger.debug(f"Offsetting SVG")
    for element in svg.elements():
        if isinstance(element, svgelements.Path):
            offset_path(element, offset_dist)

               
##########################################################################################

import argparse
import sys
from pathlib import Path


def _print_summary(svg: svgelements.SVG):
    """Print a short summary of paths/subpaths."""
    for element in svg.elements():
        if isinstance(element, svgelements.Path):
            print(f"Path ID \"{element.id}\"")
            subpath_idx = -1
            for subpath in element.as_subpaths():
                is_clockwise = calculate_is_subpath_clockwise(subpath)
                subpath_idx += 1
                is_closed = isinstance(subpath[-1], svgelements.Close)
                num_segs = len(subpath)
                if isinstance(subpath[0], svgelements.Move):
                    num_segs -= 1  # don't count initial move
                print(f"  Subpath {subpath_idx}: {num_segs} segments {"clockwise" if is_clockwise else "counter-clockwise" if is_clockwise is not None else "undetermined"}, {'closed' if is_closed else 'open'}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="kerfer — SVG path offsetting to account for kerf")
    parser.add_argument("-i", "--input", help="Input SVG file path")
    parser.add_argument("-o", "--output", help="Output SVG file path (writes modified SVG)")
    parser.add_argument("--open", dest="do_open", action="store_true", help="Open closed subpaths (replace Close with Line)")
    parser.add_argument("--close", dest="do_close", action="store_true", help="Close open subpaths (add Close where endpoints match)")
    parser.add_argument("--offset", type=float, help="Perpendicular offset distance (in same units as SVG)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging")

    args = parser.parse_args(argv)

    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    if args.input:
        input_path = Path(args.input)
        if not input_path.exists():
            logger.error(f"Input file does not exist: {input_path}")
            return 2

        logger.info(f"Loading SVG from {input_path}")
        try:
            svg = svgelements.SVG.parse(str(input_path))
        except Exception as e:
            logger.error(f"Failed to parse SVG '{input_path}': {e}")
            return 3
    else:
        # Create a dummy SVG file for testing
        with open("test_output/test_a.svg", "w") as f:
            f.write("""
            <svg width="15" height="7" xmlns="http://www.w3.org/2000/svg">
            <path d="M 1,1 l 5,5 l 0,-5 l -5,0 m 4,1 l 0,2 l -2,-2 l 2,0" id="triangle" fill="blue" />
            <path d="M 7,1 l 5,0 l 0,5 l -5,0 z m 1,1 l 0,1 l 1,0 l 0,-1 z m 2,2 l 0,1 l 1,0 l 0,-1 z" id="square" fill="red" />
            </svg>
            """)
        svg = svgelements.SVG.parse("test_output/test_a.svg")

        _print_summary(svg)

        print()
        print("Testing open...")
        open_svg(svg)
        svg.write_xml("test_output/test_b_opened.svg")
        _print_summary(svg)
        
        offset_dist = 0.25
        print()
        print(f"Testing offset {offset_dist}...")
        offset_svg(svg, offset_dist)
        svg.write_xml(f"test_output/test_c_offset_{offset_dist}.svg")
        _print_summary(svg)

        print()
        print("Testing close...")
        close_svg(svg)
        svg.write_xml("test_output/test_d_closed.svg")
        _print_summary(svg)
        return 0
    
    # Perform operations in a sensible order: open -> offset -> close
    if args.do_open:
        logger.info("Opening subpaths (in-place)")
        open_svg(svg)

    if args.offset is not None:
        try:
            offset_value = float(args.offset)
            logger.info(f"Offsetting SVG by {offset_value}")
            offset_svg(svg, offset_value)
        except Exception as e:
            logger.error(f"Offset failed: {e}")
            return 4

    if args.do_close:
        logger.info("Closing subpaths (in-place)")
        close_svg(svg)

    # If output path specified, write modified SVG; otherwise print a summary
    if args.output:
        out_path = Path(args.output)
        try:
            svg.write_xml(str(out_path))
            logger.info(f"Wrote output SVG to {out_path}")
        except Exception as e:
            logger.error(f"Failed to write output SVG '{out_path}': {e}")
            return 5

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
