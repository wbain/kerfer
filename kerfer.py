import re
import svgelements
import math
from typing import Optional
import logging
import argparse
from pathlib import Path

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
# Subpath operations

def break_apart(path: svgelements.Path) -> list[svgelements.Path]:
    """
    Breaks apart a path into its constituent subpaths.
    Args:
        path (svgelements.Path): The path to break apart.
    Returns:
        list[svgelements.Path]: A list of subpaths.
    """
    subpaths = []
    subpath_idx = 0
    for subpath in path.as_subpaths():
        new_path = svgelements.Path(subpath.d())
        if path.id:
            new_path.id = f"{path.id}_{subpath_idx}"
        subpaths.append(new_path)
        subpath_idx += 1
    return subpaths


def break_apart_svg(svg: svgelements.SVG) -> tuple[int, int]:
    """
    Breaks apart all paths in the SVG into their constituent subpaths.
    Args:
        svg (svgelements.SVG): The SVG to break apart.
    Returns:
        int: The number of new subpaths created.
    """
    new_paths = []
    paths_to_delete = []
    for element in svg.elements():
        if isinstance(element, svgelements.Path):
            paths = break_apart(element)
            new_paths.extend(paths)
            paths_to_delete.append(element)

    for path in paths_to_delete:
        svg.remove(path)

    for path in new_paths:
        svg.append(path)

    return len(paths_to_delete), len(new_paths)


def bounding_box_contains(container: svgelements.Path, contained: svgelements.Path) -> bool:
    """
    Checks if the bounding box of the contained path is completely within the bounding box of the container path.
    Args:
        container (svgelements.Path): The container path.
        contained (svgelements.Path): The path to check.
    Returns:
        bool: True if the contained path's bounding box is within the container's, False otherwise.
    """
    container_bbox = container.bbox()
    contained_bbox = contained.bbox()

    if (container_bbox is None) or (contained_bbox is None):
        logger.warning("bounding_box_contains: at least one of the bounding boxes is None!")
        return False

    return (contained_bbox[0] >= container_bbox[0] and
            contained_bbox[2] <= container_bbox[2] and
            contained_bbox[1] >= container_bbox[1] and
            contained_bbox[3] <= container_bbox[3])


def nest_paths(container: svgelements.Path, contained: svgelements.Path):
    """
    Combines two paths if the contained path's bounding box is completely within the container path's bounding box.
    Args:
        container (svgelements.Path): The container path.
        contained (svgelements.Path): The path to combine.
    """
    # Nest contained into container if bounding box is within container's
    for segment in contained:
        container.append(segment)


def nest_svg(svg: svgelements.SVG) -> int:
    """
    Combines subpaths within the SVG if their bounding boxes are completely within the bounding box of another subpath.
    Args:
        svg (svgelements.SVG): The SVG to nest.
    Returns:
        int: The number of paths nested.
    """
    paths_to_delete = []
    for element in svg.elements():
        if isinstance(element, svgelements.Path):
            for other_element in svg.elements():
                if isinstance(other_element, svgelements.Path) and (element != other_element):
                    if bounding_box_contains(other_element, element):
                        logger.debug(f"  Nesting path \"{element.id}\" into path \"{other_element.id}\"")
                        nest_paths(other_element, element)
                        paths_to_delete.append(element)
                        break

    for path in paths_to_delete:
        svg.remove(path)

    return len(paths_to_delete)


##########################################################################################
# Uniquifying/Simplifying

def zero_cull_path(path: svgelements.Path) -> int:
    """
    Removes zero-length segments from a subpaths in the path, in-place.
    Args:
        path (svgelements.Path): The path object to process.
    Returns:
        int: The number of zero-length segments removed.
    """
    logger.debug(f"  Deduplicating path \"{path.id}\"")
    path_segs_to_delete = []
    subpath_idx = -1
    for subpath in path.as_subpaths():
        subpath_idx += 1

        subpath_seg_idx = -1
        for segment in subpath:
            subpath_seg_idx += 1

            if isinstance(segment, svgelements.Move):
                continue

            if segment.start == segment.end:
                if isinstance(segment, svgelements.Line):
                    path_segs_to_delete.append(subpath.index_to_path_index(subpath_seg_idx))
                elif isinstance(segment, svgelements.Close):
                    path_segs_to_delete.append(subpath.index_to_path_index(subpath_seg_idx - 1))
    
    for seg_idx in reversed(path_segs_to_delete):
        del path[seg_idx]

    return len(path_segs_to_delete)


def zero_cull_svg(svg: svgelements.SVG) -> int:
    """
    Removes zero-length segments from all subpaths in the SVG, in-place.
    Args:
        svg (svgelements.SVG): The SVG object to process.
    Returns:
        int: The number of zero-length segments removed.
    """
    logger.debug("Uniquifying SVG")
    num_culled = 0
    for element in svg.elements():
        if isinstance(element, svgelements.Path):
            num_culled += zero_cull_path(element)

    return num_culled


def are_collinear(seg1: svgelements.Linear, seg2: svgelements.Linear) -> bool:
    """
    Checks if two linear segments are collinear.
    Args:
        seg1 (svgelements.Linear): The first segment.
        seg2 (svgelements.Linear): The second segment.
    Returns:
        bool: True if the segments are collinear, False otherwise.
    """
    if (not seg1.start) or (not seg1.end) or (not seg2.start) or (not seg2.end):
        logger.warning("are_collinear: at least one of the segments has no start or end point!")
        return False
    
    if (not seg1.start.x) or (not seg1.start.y) or (not seg1.end.x) or (not seg1.end.y) or \
            (not seg2.start.x) or (not seg2.start.y) or (not seg2.end.x) or (not seg2.end.y):
        logger.warning("are_collinear: at least one of the segment endpoints has no x or y value!")
        return False
    
    # Calculate vectors from start to end for both segments
    v1_x = seg1.end.x - seg1.start.x
    v1_y = seg1.end.y - seg1.start.y
    v2_x = seg2.end.x - seg2.start.x
    v2_y = seg2.end.y - seg2.start.y

    # Calculate the cross product of the two vectors
    cross_product = v1_x * v2_y - v1_y * v2_x

    # If the cross product is close to zero, the segments are collinear
    return math.isclose(cross_product, 0, rel_tol=1e-6)


def simplify_path(path: svgelements.Path) -> int:
    """
    Simplifies a path by removing collinear segments, in-place.
    Args:
        path (svgelements.Path): The path object to process.
    Returns:
        int: The number of collinear segments removed.
    """
    logger.debug(f"  Simplifying path \"{path.id}\"")
    path_segs_to_delete = []
    subpath_idx = -1
    for subpath in path.as_subpaths():
        subpath_idx += 1

        subpath_seg_idx = -1
        for segment in subpath:
            subpath_seg_idx += 1

            if isinstance(segment, svgelements.Move):
                # logger.debug(f"      Skipping initial move segment")
                continue
            elif subpath_seg_idx == 1:
                # logger.debug(f"      Skipping first segment after initial move")
                continue

            prev_seg = subpath[subpath_seg_idx - 1]

            if are_collinear(segment, prev_seg):
                segment.start = prev_seg.start    # extend current segment backward
                # Make sure the previous gets removed later
                path_segs_to_delete.append(subpath.index_to_path_index(subpath_seg_idx - 1))
    
    for seg_idx in reversed(path_segs_to_delete):
        del path[seg_idx]

    return len(path_segs_to_delete)


def simplify_svg(svg: svgelements.SVG) -> int:
    """
    Simplifies all paths in the SVG by removing collinear segments, in-place.
    Args:
        svg (svgelements.SVG): The SVG object to process.
    Returns:
        int: The number of collinear segments removed.
    """
    logger.debug("Simplifying SVG")
    num_simplified = 0
    for element in svg.elements():
        if isinstance(element, svgelements.Path):
            num_simplified += simplify_path(element)

    return num_simplified


##########################################################################################
# Linifying/Closing

def linify_path(path: svgelements.Path) -> tuple[int, int]:
    """
    Linifies all subpaths in the SVG, in-place, by replacing close commands with line segments.
    Args:
        path (svgelements.Path): The path object to process.
    Returns:
        tuple[int, int]: The number of close commands that were replaced with lines, and the total number of subpaths.
    """
    logger.debug(f"  Linifying path \"{path.id}\"")
    num_total = 0
    num_replaced = 0
    subpath_idx = -1
    for subpath in path.as_subpaths():
        num_total += 1
        subpath_idx += 1

        found_close = False

        subpath_seg_idx = -1
        for segment in subpath:
            subpath_seg_idx += 1

            if isinstance(segment, svgelements.Close):
                found_close = True
                logger.debug(f"    Linifying subpath {subpath_idx}")
                # replace close with line segment
                start = segment.start
                end = segment.end
                path_seg_idx = subpath.index_to_path_index(subpath_seg_idx)
                del path[path_seg_idx]
                path.insert(path_seg_idx, svgelements.Line(start, end))
                num_replaced += 1
            
        if not found_close:
            logger.debug(f"    Subpath {subpath_idx} already ends with line segment")

    return num_replaced, num_total
            

def linify_svg(svg: svgelements.SVG) -> tuple[int, int]:
    """
    Linifies all subpaths in the SVG, in-place, by replacing close commands with line segments.
    Args:
        svg (svgelements.SVG): The SVG object to process.
    Returns:
        tuple[int, int]: The number of close commands that were replaced with lines, and the total number of subpaths.
    """
    logger.debug("Linifying SVG")
    num_total = 0
    num_replaced = 0
    for element in svg.elements():
        if isinstance(element, svgelements.Path):
            r, t = linify_path(element)
            num_replaced += r
            num_total += t

    return num_replaced, num_total


def close_path(path: svgelements.Path) -> tuple[int, int]:
    """
    Closes all subpaths in the SVG, in-place, by replacing final line segments with close commands.
    Args:
        path (svgelements.Path): The path object to process.
    Returns:
        tuple[int, int]: The number of line segments that were replaced with close commands, and the total number of subpaths.
    """
    logger.debug(f"  Closing path \"{path.id}\"")
    num_total = 0
    num_replaced = 0
    subpath_idx = -1
    for subpath in path.as_subpaths():
        subpath_idx += 1
        num_total += 1

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
            num_replaced += 1
        except Exception as e:
            logger.error(f"      Failed to close subpath {subpath_idx}! {e}")
            # logger.error(f"        {subpath.d()}")

    return num_replaced, num_total


def close_svg(svg: svgelements.SVG) -> tuple[int, int]:
    """
    Closes all subpaths in the SVG, in-place, by replacing final line segments with close commands.
    Args:
        svg (svgelements.SVG): The SVG object to process.
    Returns:
        tuple[int, int]: The number of line segments that were replaced with close commands, and the total number of subpaths.
    """
    logger.debug(f"Closing SVG")
    num_total = 0
    num_replaced = 0
    for element in svg.elements():
        if isinstance(element, svgelements.Path):
            r, t = close_path(element)
            num_replaced += r
            num_total += t

    return num_replaced, num_total


##########################################################################################
# Dilating

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


def dilate_subpath(subpath: svgelements.Subpath, offset_dist: float):
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


def dilate_path(path: svgelements.Path, offset_dist: float):
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
        dilate_subpath(subpath, dist)


def dilate_svg(svg: svgelements.SVG, offset_dist: float):
    """
    Offsets all the paths in an SVG by a given distance. Modifies the paths in-place.
    Args:
        svg (svgelements.SVG): The SVG to offset.
        offset_dist (float): The distance to offset the SVG.
    """
    logger.debug(f"Offsetting SVG")
    for element in svg.elements():
        if isinstance(element, svgelements.Path):
            dilate_path(element, offset_dist)


##########################################################################################
# Miscellaneous

def copy_and_group_all_paths(svg: svgelements.SVG, group_name: str):
    """
    Copies and groups all paths in an SVG into a single group with the given name.
    Args:
        svg (svgelements.SVG): The SVG to process.
        group_name (str): The name of the group to create.
    Returns:
        svgelements.Group: The new group containing all paths.
    """
    logger.debug(f"Copying and grouping all paths in SVG into group \"{group_name}\"")
    new_group = svgelements.Group()
    new_group.id = group_name
    for element in svg.elements():
        if isinstance(element, svgelements.Path):
            # new_group.append(element.copy())
            new_group.append(svgelements.Path(element))
    return new_group


def generate_unique_output_path(input_path: Path, dilate: float, units: str) -> str:
    """
    Generates a unique output file path based on the input file path.
    Args:
        input_path (Path): The input file path.
        dilate (float): The dilation distance used in the output filename.
        units (str): The units used in the output filename.
    Returns:
        str: The unique output file path.
    """
    input_path_str = str(input_path)

    dilated = f"{dilate}{units}" if dilate is not None else "0"

    output_path_str = f"{input_path.stem}_d{dilated}.svg"
    output_path = input_path.with_name(output_path_str)

    counter = 1
    while output_path.exists():
        copy_suffix = "copy" if counter == 1 else f"copy{counter}"
        output_path_str = f"{input_path.stem}_d{dilated}_{copy_suffix}.svg"
        output_path = input_path.with_name(output_path_str)
        counter += 1

    return str(output_path)


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


##########################################################################################
# Main

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="kerfer — SVG path offsetting to account for kerf")
    parser.add_argument("-i", "--input", help="Input SVG file path")
    parser.add_argument("-o", "--output", help="Optional: Output SVG file path (writes modified SVG)")
    parser.add_argument("-b", "--break", dest="do_break", action="store_true", help="Break apart subpaths into separate paths")
    parser.add_argument("-n", "--nest", dest="do_nest", action="store_true", help="Nest subpaths into parent paths")
    parser.add_argument("-l", "--line", dest="do_line", action="store_true", help="Open closed subpaths (replace each Close with a *Line*)")
    parser.add_argument("-z", "--zero_cull", dest="do_zero_cull", action="store_true", help="Removes *zero-length* segments from paths")
    parser.add_argument("-s", "--simplify", dest="do_simplify", action="store_true", help="Remove unnecessary points from paths to *simplify* them")
    parser.add_argument("-d", "--dilate", type=float, help="Perpendicular offset *dilation* distance (in same units as SVG)")
    parser.add_argument("-c", "--close", dest="do_close", action="store_true", help="Close open subpaths (add *Close* where endpoints match)")
    parser.add_argument("-r", "--rebreak", dest="do_rebreak", action="store_true", help="Rebreak subpaths into separate paths")
    parser.add_argument("-a", "--all", dest="do_all", action="store_true", help="Default if no other processing specified: Do *all* the steps - line, zero-cull, simplify, dilate, close")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging")

    args = parser.parse_args(argv)

    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    if not args.do_break and not args.do_nest and not args.do_line and not args.do_zero_cull and not args.do_simplify and not args.dilate and not args.do_close:
        args.do_all = True

    if args.do_all:
        args.do_break = True
        args.do_nest = True
        args.do_line = True
        args.do_zero_cull = True
        args.do_simplify = True
        if not args.dilate:
            args.dilate = 0.15
        args.do_close = True
        args.do_rebreak = True

    if args.input:
        input_path = Path(args.input)
        if not input_path.exists():
            logger.error(f"Input file does not exist: {input_path}")
            return 2
        
        # Open the input SVG file as text and search it for the width attribute to parse the units at the end of the value
        units = find_units(input_path)  

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
            <path d="M 1,1 l 2,2 l 2,2 l 1,1 l 0,-5 l -3,0 l -2,0 m 4,1 l 0,2 l -1,-1 l -1,-1 l 2,0" id="triangle" fill="blue" />
            <path d="M 7,1 l 5,0 l 0,2 l 0,3 l -5,0 z m 1,1 l 0,1 l 0,0 l 1,0 l 0,-1 z m 2,2 l 0,1 l 1,0 l 0,-1 z" id="square" fill="red" />
            </svg>
            """)
        svg = svgelements.SVG.parse("test_output/test_a.svg")

        # _print_summary(svg)

        print()
        print("Testing break apart...")
        num_old, num_new = break_apart_svg(svg)
        logger.info(f"  Broke {num_old} paths into {num_new} subpaths")
        svg.write_xml("test_output/test_b_broken.svg")
        # _print_summary(svg)

        print()
        print("Testing nest...")
        num_nested = nest_svg(svg)
        logger.info(f"  Nested {num_nested} paths into others")
        svg.write_xml("test_output/test_c_nested.svg")
        # _print_summary(svg)

        print()
        print("Testing linify...")
        num_replaced, num_total = linify_svg(svg)
        logger.info(f"  Replaced {num_replaced} close commands with lines in {num_total} subpaths")
        svg.write_xml("test_output/test_d_linified.svg")
        # _print_summary(svg)
        
        print()
        print("Testing zero-cull...")
        num_culled = zero_cull_svg(svg)
        logger.info(f"  Removed {num_culled} zero-length segments")
        svg.write_xml("test_output/test_e_zero_culled.svg")
        # _print_summary(svg)

        print()
        print("Testing simplify...")
        num_simplified = simplify_svg(svg)
        logger.info(f"  Removed {num_simplified} collinear segments")
        svg.write_xml("test_output/test_f_simplified.svg")
        # _print_summary(svg)

        print()
        print(f"Testing dilate {args.dilate}...")
        dilate_svg(svg, args.dilate)
        svg.write_xml(f"test_output/test_g_dilated_{args.dilate}.svg")
        # _print_summary(svg)

        print()
        print("Testing close...")
        num_replaced, num_total = close_svg(svg)
        logger.info(f"  Replaced {num_replaced} line segments with close commands in {num_total} subpaths")
        svg.write_xml("test_output/test_h_closed.svg")
        # _print_summary(svg)
        return 0
    
    if not args.output:
        args.output = generate_unique_output_path(input_path, args.dilate, units)

    group = copy_and_group_all_paths(svg, "original_paths")

    # Perform operations in a sensible order: break -> line -> zero_cull -> simplify -> offset -> close -> rebreak

    if args.do_break:
        logger.info("Breaking apart subpaths")
        num_old, num_new = break_apart_svg(svg)
        logger.info(f"  Broke {num_old} paths into {num_new} subpaths")

        # Assign unique IDs to each path
        path_idx = 0
        for element in svg.elements():
            if isinstance(element, svgelements.Path):
                if not element.id:
                    element.id = f"path_{path_idx}"
                path_idx += 1

        for element in svg.elements():
            if isinstance(element, svgelements.Path):
                logger.debug(f"    Path ID \"{element.id}\": {len(list(element.as_subpaths()))} subpaths")

    if args.do_nest:
        logger.info("Nesting subpaths")
        num_nested = nest_svg(svg)
        logger.info(f"  Nested {num_nested} paths into others")

        for element in svg.elements():
            if isinstance(element, svgelements.Path):
                logger.debug(f"    Path ID \"{element.id}\": {len(list(element.as_subpaths()))} subpaths")

    if args.do_line:
        logger.info("Linifying subpaths (in-place)")
        num_replaced, num_total = linify_svg(svg)
        logger.info(f"  Replaced {num_replaced} close commands with lines in {num_total} subpaths")

    if args.do_zero_cull:
        logger.info("Zero-culling subpaths (in-place)")
        num_culled = zero_cull_svg(svg)
        logger.info(f"  Removed {num_culled} zero-length segments")

    if args.do_simplify:
        logger.info("Simplifying subpaths (in-place)")
        num_simplified = simplify_svg(svg)
        logger.info(f"  Removed {num_simplified} collinear segments")

    if args.dilate is not None:
        try:
            offset_value = float(args.dilate)
            logger.info(f"Offsetting SVG by {offset_value}")
            dilate_svg(svg, offset_value)
        except Exception as e:
            logger.error(f"  Offset failed: {e}")
            return 4

    if args.do_close:
        logger.info("Closing subpaths (in-place)")
        num_replaced, num_total = close_svg(svg)
        logger.info(f"  Replaced {num_replaced} line segments with close commands in {num_total} subpaths")

    stroke_width = str(args.dilate) if args.dilate else "1"
    color_fill = svgelements.Color("none")
    color_outer = svgelements.Color("#0000ff")  # blue
    color_inner = svgelements.Color("#ff0000")  # red
    color_orig = svgelements.Color("#00ff00")   # green

    if args.do_rebreak:
        logger.info("Rebreaking subpaths into separate paths")
        num_old, num_new = break_apart_svg(svg)
        logger.info(f"  Rebroke {num_old} paths into {num_new} subpaths")

        for element in svg.elements():
            if isinstance(element, svgelements.Path):
                logger.info(f"    Path ID \"{element.id}\": {len(list(element.as_subpaths()))} subpaths")

    # Assign styles to outer and inner paths
    for element in svg.elements():
        if isinstance(element, svgelements.Path):
            if isinstance(element, svgelements.GraphicObject):
                element.fill = color_fill
                element.stroke_width = stroke_width
                if isinstance(element.id, str) and element.id.endswith("_0"):
                    element.stroke = color_outer
                else:
                    element.stroke = color_inner
            else:
                logger.warning(f"Element ID \"{element.id}\" is not a GraphicObject, cannot assign style!")

    # Assign original style to paths in original group
    for element in group:
        if isinstance(element, svgelements.Path):
            element.fill = color_fill
            element.stroke_width = stroke_width
            element.stroke = color_orig

    svg.append(group)

    # For each path in SVG, find its transform matrix and apply it to the path data, then remove the transform attribute
    # Doesn't work: At this point, paths only have identity transforms.
    # It appears that transforms get added automatically when file is written out.
    # for element in svg.elements():
    #     if isinstance(element, svgelements.Path):
    #         if element.transform is not None:
    #             logger.info(f"Applying transform to path \"{element.id}\"")
    #             try:
    #                 logger.info(f"  Transform before: {element.transform}")
    #                 element.reify()
    #                 logger.info(f"  Transform after: {element.transform}")
    #             except Exception as e:
    #                 logger.error(f"  Failed to apply transform to path \"{element.id}\": {e}")

    # If output path specified, write modified SVG
    if args.output:
        temp_out_file = args.output + ".tmp"
        try:
            temp_out_file_path = Path(temp_out_file)
            svg.write_xml(temp_out_file_path)
            logger.info(f"Wrote temporary output SVG to {temp_out_file_path}")
        except Exception as e:
            logger.error(f"Failed to write output SVG '{temp_out_file_path}': {e}")
            return 5

        # Open the resulting SVG file as text and insert the correct units into the width and height attributes at the end of each
        out_path = Path(args.output)
        inject_units(units, temp_out_file_path, out_path)

        # Remove the temporary output file
        try:
            temp_out_file_path.unlink()
            logger.info(f"Removed temporary file {temp_out_file_path}")
        except Exception as e:
            logger.warning(f"Failed to remove temporary file {temp_out_file_path}: {e}")

    return 0


def find_units(input_path: Path) -> str:
    """
    Find the units of the width attribute in an SVG file by searching line-by-line.
    Args:
        input_path (Path): The path to the SVG file.
    Returns:
        str: The units of the width attribute, or empty string if not found or unitless.
    """
    width_units = ''
    
    try:
        with open(input_path, 'r', encoding='utf-8') as f:
            for line in f:
                # Search for width attribute in this line
                width_attrib_match = re.search(r'width="([^"]+)"', line)
                if width_attrib_match:
                    width_val_str = width_attrib_match.group(1)
                    logger.info(f"  SVG width value string: {width_val_str}")
                    
                    # Extract the units from the width value
                    width_units_match = re.search(r'(\d+(\.\d+)?)\s*([a-zA-Z]+)', width_val_str)
                    if width_units_match:
                        width_units = width_units_match.group(3)
                        logger.info(f"  SVG width units: {width_units}")
                        return width_units
                    else:
                        logger.info("  SVG is unitless.")
                    break
        
        if width_attrib_match is None:
            logger.warning("  Unable to find width attribute in SVG file.")
        elif width_units == '':
            logger.info("  SVG is unitless.")
    except Exception as e:
        logger.error(f"  Error reading SVG file: {e}")
    
    return width_units


def inject_units(units: str, in_path: Path, out_path: Path):
    """
    Inject the correct units into the width and height attributes of an SVG file by processing line-by-line.
    Args:
        units (str): The units to inject (e.g., 'mm', 'in', 'px').
        in_path (Path): The path to the input SVG file.
        out_path (Path): The path to the output SVG file.
    """
    width_updated = False
    height_updated = False
    
    try:
        with open(in_path, 'r', encoding='utf-8') as infile, \
             open(out_path, 'w', encoding='utf-8') as outfile:
            
            for line in infile:
                modified_line = line
                
                # Process width attribute if present in this line
                if not width_updated and re.search(r'width="[^"]+"', line):
                    width_match = re.search(r'width="([^"]+)"', line)
                    if width_match:
                        width_val_str = width_match.group(1)
                        logger.info(f"  SVG width value string: {width_val_str}")
                        
                        # Inject units into width value
                        new_width_val_str = f"{width_val_str}{units}"
                        modified_line = re.sub(
                            r'width="[^"]+"',
                            f'width="{new_width_val_str}"',
                            line
                        )
                        logger.info(f"  Updated SVG width units to {units}")
                        width_updated = True
                    else:
                        logger.warning("  Unable to find width value in width attribute of SVG")
                
                # Process height attribute if present in this line
                if not height_updated and re.search(r'height="[^"]+"', modified_line):
                    height_match = re.search(r'height="([^"]+)"', modified_line)
                    if height_match:
                        height_val_str = height_match.group(1)
                        logger.info(f"  SVG height value string: {height_val_str}")
                        
                        # Inject units into height value
                        new_height_val_str = f"{height_val_str}{units}"
                        modified_line = re.sub(
                            r'height="[^"]+"',
                            f'height="{new_height_val_str}"',
                            modified_line
                        )
                        logger.info(f"  Updated SVG height units to {units}")
                        height_updated = True
                    else:
                        logger.warning("  Unable to find height value in height attribute of SVG")
                
                outfile.write(modified_line)
        
        if not width_updated:
            logger.warning("  Unable to find width attribute of SVG")
        if not height_updated:
            logger.warning("  Unable to find height attribute of SVG")
    
    except Exception as e:
        logger.error(f"  Error processing SVG file: {e}")


if __name__ == "__main__":
    raise SystemExit(main())
