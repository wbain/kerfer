from dataclasses import dataclass, field
import re
import math
from typing import List, Optional
import logging
import argparse
import pathlib

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(levelname)s %(message)s")
for h in logging.getLogger().handlers:
    if h.formatter is not None:
        h.formatter.default_msec_format = '%s.%03d'
logger = logging.getLogger(__name__)


# Data class that represents an x, y coordinate pair
@dataclass
class Point:
    x: float
    y: float


# Data class to represent a simple shape as a list of points
@dataclass
class SimpleShape:
    points: List[Point] = field(default_factory=list)
    attributes: dict[str, str] = field(default_factory=dict)


# Data class to represent a compound shape as a list of simple shapes,
# where the first one is the main shape and the subsequent ones are holes
@dataclass
class CompoundShape:
    component_shapes: List[SimpleShape]
    attributes: dict[str, str] = field(default_factory=dict)


##########################################################################################
# Subpath operations

def compute_bounding_box(shape: SimpleShape) -> Optional[tuple[float, float, float, float]]:
    """
    Computes the bounding box of a shape.
    Args:
        shape (SimpleShape): The shape to compute the bounding box for.
    Returns:
        tuple: A tuple containing the (min_x, min_y, max_x, max_y) coordinates of the bounding box.
    """
    if not shape.points:
        return None

    min_x = min(p.x for p in shape.points)
    max_x = max(p.x for p in shape.points)
    min_y = min(p.y for p in shape.points)
    max_y = max(p.y for p in shape.points)

    return (min_x, min_y, max_x, max_y)


def bounding_box_contains(container: SimpleShape, contained: SimpleShape) -> bool:
    """
    Checks if the bounding box of the contained path is completely within the bounding box of the container path.
    Args:
        container (SimpleShape): The container shape.
        contained (SimpleShape): The shape to check.
    Returns:
        bool: True if the contained shape's bounding box is within the container's, False otherwise.
    """
    container_bbox = compute_bounding_box(container)
    contained_bbox = compute_bounding_box(contained)

    if (container_bbox is None) or (contained_bbox is None):
        logger.warning("bounding_box_contains: at least one of the bounding boxes is None!")
        return False

    return (contained_bbox[0] >= container_bbox[0] and
            contained_bbox[2] <= container_bbox[2] and
            contained_bbox[1] >= container_bbox[1] and
            contained_bbox[3] <= container_bbox[3])


##########################################################################################
# Uniquifying/Simplifying

def zero_cull(shape: SimpleShape) -> int:
    """
    Removes redundant points, in-place.
    Args:
        shape (SimpleShape): The shape to process.
    Returns:
        int: The number of redundant points removed.
    """
    logger.debug(f"  Deduplicating points")
    points_to_delete_indicies = []
    point_idx = -1
    for point in shape.points:
        point_idx += 1

        next_idx = (point_idx + 1) % len(shape.points)
        next_point = shape.points[next_idx]

        if point == next_point:
            points_to_delete_indicies.append(point_idx)

    for pt_idx in reversed(points_to_delete_indicies):
        del shape.points[pt_idx]

    return len(points_to_delete_indicies)


def are_collinear(a: Point, b: Point, c: Point) -> bool:
    """
    Checks if three points are collinear.
    Args:
        a (Point): The first point.
        b (Point): The second point.
        c (Point): The third point.
    Returns:
        bool: True if the points are collinear, False otherwise.
    """
    # Calculate the cross product of vectors AB and BC
    # 0 the cross product is close to zero, the points are collinear    
    # Calculate vectors from start to end for both segments
    v1_x = b.x - a.x
    v1_y = b.y - a.y
    v2_x = c.x - b.x
    v2_y = c.y - b.y

    # Calculate the cross product of the two vectors
    cross_product = v1_x * v2_y - v1_y * v2_x

    # If the cross product is close to zero, the segments are collinear
    return math.isclose(cross_product, 0, rel_tol=1e-6)


def simplify_shape(shape: SimpleShape) -> int:
    """
    Simplifies a shape by removing collinear segments, in-place.
    Args:
        shape (SimpleShape): The shape to process.
    Returns:
        int: The number of collinear segments removed.
    """
    logger.debug(f"  Simplifying shape")
    points_to_delete_indices = []
    point_idx = -1
    for point in shape.points:
        point_idx += 1

        prev_idx = (point_idx - 1) % len(shape.points)
        next_idx = (point_idx + 1) % len(shape.points)

        prev_point = shape.points[prev_idx]
        next_point = shape.points[next_idx]

        if are_collinear(prev_point, point, next_point):
            points_to_delete_indices.append(point_idx)

    for pt_idx in reversed(points_to_delete_indices):
        del shape.points[pt_idx]

    return len(points_to_delete_indices)


##########################################################################################
# Dilating

def line_segment_intersection(p1: Point, p2: Point, p3: Point, p4: Point):
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
    x1 = p1.x
    y1 = p1.y
    x2 = p2.x
    y2 = p2.y
    x3 = p3.x
    y3 = p3.y
    x4 = p4.x
    y4 = p4.y

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


def is_contained(shape: SimpleShape, shapes: List[SimpleShape]) -> bool:
    """
    Determines if a simple shape is contained within another shape (i.e., is a hole).
    Args:
        shape (SimpleShape): The shape to check.
        shapes (List[SimpleShape]): The list of shapes to check against.
    Returns:
        bool: True if the shape is contained (a hole), False otherwise.
    """
    for other_shape in shapes:
        if shape != other_shape:
            if bounding_box_contains(other_shape, shape):
                return True
            
    return False


def offset_endpoints(p1: Point, p2: Point, offset_distance: float):
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

    return Point(new_p1_x, new_p1_y), Point(new_p2_x, new_p2_y)


def dilate_shape(shape: SimpleShape, offset_dist: float):
    """
    Offsets a simple shape by a given distance. Modifies the shape in-place.
    Args:
        shape (SimpleShape): The shape to offset.
        offset_dist (float): The distance to offset the shape.
    """
    if len(shape.points) < 2:
        return

    # First, go around the loop and offset each line segment
    off_segs = []
    for i in range(len(shape.points)):
        p1 = shape.points[i]
        p2 = shape.points[(i + 1) % len(shape.points)]
        offset_p1, offset_p2 = offset_endpoints(
            Point(p1.x, p1.y),
            Point(p2.x, p2.y),
            offset_dist
        )
        off_segs.append((offset_p1, offset_p2))

    # Next, find the intersections of each offset segment with the next
    intersections = []
    for i in range(len(off_segs)):
        current_seg = off_segs[i]
        next_seg = off_segs[(i + 1) % len(off_segs)]
        intersection_point = line_segment_intersection(
            current_seg[1], current_seg[0],
            next_seg[0], next_seg[1]
        )

        if intersection_point is not None:
            intersections.append(intersection_point)
        else:
            intersections.append((current_seg[1].x, current_seg[1].y))

    # Finally, update the original shape points to use the intersection points
    for i in range(len(shape.points)):
        shape.points[i].x = intersections[i][0]
        shape.points[i].y = intersections[i][1]


##########################################################################################
# Miscellaneous

def generate_unique_output_path(input_path: pathlib.Path, dilate: float, units: str) -> str:
    """
    Generates a unique output file path based on the input file path.
    Args:
        input_path (pathlib.Path): The input file path.
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


##########################################################################################
# Parsing

def parse_svg_file(input_path: pathlib.Path) -> tuple[str, List[str], str]:
    """
    Parses an SVG file and returns its preamble, list of paths, and postamble.
    Args:
        input_path (pathlib.Path): The path to the SVG file.
    Returns:
        tuple[str, List[str], str]: A tuple containing the preamble, list of paths, and postamble of the SVG.
    """
    with open(input_path, "r") as f:
        content = f.read()

        # Extract preamble
        path_match = re.search(r"<path ", content)
        if path_match:
            preamble = content[:path_match.start()]
        else:
            preamble = ""

        # Extract paths
        path_node_strs = re.findall(r"<path[^>]*>", content)

        # Extract postamble (or just fake it)
        postamble = "</svg>"  # This is a placeholder; in practice, extract it more carefully

    return preamble, path_node_strs, postamble


def parse_path_node_string(path_str: str) -> tuple[List[str], dict[str, str]]:
    """
    Parses a path string and returns a list of path segments and a dictionary of attributes.
    Args:
        path_str (str): The path string to parse.

    Returns:
        tuple[List[str], dict[str, str]]: A tuple containing the list of path segments and the attributes dictionary.
    """
    # Extract attributes
    attributes = {}
    attr_match = re.search(r'<path\s+([^>]+)>', path_str)
    if attr_match:
        attr_str = attr_match.group(1)
        attr_list = re.findall(r'(\w+)="([^"]+)"', attr_str)
        for attr_name, attr_value in attr_list:
            if not attr_name == "d":
                attributes[attr_name] = attr_value

    # Extract path data
    path_data = re.search(r'd="([^"]+)"', path_str)
    if path_data:
        path_data_str = path_data.group(1)
        path_segments = re.findall(r'([MLCZmlcz][^MLCZmlcz]*)', path_data_str)
    else:
        path_segments = []

    return path_segments, attributes
    

def extract_simple_shapes(path_data_strs: List[str]) -> List[SimpleShape]:
    """
    Extracts simple shapes from a list of path data strings. Resulting simple shapes do not have redundant closing points.
    Args:
        path_data_strs (List[str]): The list of path data strings to parse.
    Returns:
        List[SimpleShape]: A list of simple shapes extracted from the path data strings.
    """
    shapes = []
    current_shape = None

    for path_data_str in path_data_strs:
        # Extract the command and parameters
        command = path_data_str[0]
        parameters = path_data_str[1:]

        # Process the command and parameters
        if command in "Mm":
            # Move to (M) or move to relative (m)
            if current_shape:
                if current_shape.points[-1].x == current_shape.points[0].x and current_shape.points[-1].y == current_shape.points[0].y:
                    current_shape.points.pop()    # Remove the last point if it's a duplicate of the first point
                shapes.append(current_shape)    # Store the previous path, if any
            current_shape = SimpleShape()
            x, y = map(float, parameters.split())
            current_shape.points.append(Point(x, y))
        elif command in "Ll":
            # Line to (L) or line to relative (l)
            x, y = map(float, parameters.split())
            if current_shape:
                current_shape.points.append(Point(x, y))
        elif command in "Cc":
            # Curve to (C) or curve to relative (c)
            # For simplicity, we'll just ignore the control points and treat it as a line to the end point
            x, y = map(float, parameters.split()[-2:])
            if current_shape:
                current_shape.points.append(Point(x, y))
        elif command in "Zz":
            # Close path (Z) or close path relative (z)
            # Ignore closing segments
            pass
        else:
            # Unknown command
            logger.warning(f"Unknown command: {command}")

    return shapes


def to_svg_path(simple_shape: SimpleShape) -> str:
    """
    Converts a simple shape to an SVG path string. Assumes simple shape does not have a redundant closing point.
    Args:
        simple_shape (SimpleShape): The simple shape to convert.
    Returns:
        str: The SVG path string representing the simple shape.
    """
    path_data = "\"M"
    for i, point in enumerate(simple_shape.points):
        if i > 0:
            path_data += " L"
        path_data += f" {point.x},{point.y}"

    path_data += " Z"   # auto close the shape
    path_data += "\""

    # Build the attributes string
    attributes = ""
    for key, value in simple_shape.attributes.items():
        attributes += f' {key}="{value}"'

    return f'<path d={path_data}{attributes} />\n'


##########################################################################################
# Main

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="kerfer — SVG path offsetting to account for kerf")
    parser.add_argument("-i", "--input", help="Input SVG file path")
    parser.add_argument("-o", "--output", help="Optional: Output SVG file path (writes modified SVG)")
    parser.add_argument("-z", "--zero_cull", dest="do_zero_cull", action="store_true", help="Remove *zero-length* segments from paths")
    parser.add_argument("-s", "--simplify", dest="do_simplify", action="store_true", help="Remove unnecessary points from paths to *simplify* them")
    parser.add_argument("-d", "--dilate", type=float, help="Perpendicular offset *dilation* distance (in same units as SVG)")
    parser.add_argument("-a", "--all", dest="do_all", action="store_true", help="Default if no other processing specified except dilation: Do *all* the steps - break, nest, line, zero-cull, simplify, dilate, close, rebreak")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug loggin output")

    args = parser.parse_args(argv)

    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    if not args.do_zero_cull and not args.do_simplify:
        args.do_all = True

    if not args.dilate:
        args.dilate = 0.15

    if args.do_all:
        args.do_zero_cull = True
        args.do_simplify = True

    if args.input:
        input_path = pathlib.Path(args.input)
        if not input_path.exists():
            logger.error(f"Input file does not exist: {input_path}")
            return 2

        # Open the input SVG file as text and search it for the width attribute to parse the units at the end of the value
        units = find_units(input_path)  

        logger.info(f"Parsing SVG from {input_path}")
        preamble, path_node_strs, postamble = parse_svg_file(input_path)
        simple_shapes = []
        for path_node_str in path_node_strs:
            path_data_strs, attrs = parse_path_node_string(path_node_str)
            path_simple_shapes = extract_simple_shapes(path_data_strs)
            simple_shapes.extend(path_simple_shapes)
    
    if not args.output:
        args.output = generate_unique_output_path(input_path, args.dilate, units)

    if args.do_zero_cull:
        logger.info("Zero-culling segments")
        num_culled = 0
        for simple_shape in simple_shapes:
            num_culled += zero_cull(simple_shape)
        logger.info(f"  Removed {num_culled} redundant points")

    if args.do_simplify:
        logger.info("Simplifying subpaths (in-place)")
        num_simplified = 0
        for simple_shape in simple_shapes:
            num_simplified += simplify_shape(simple_shape)
        logger.info(f"  Removed {num_simplified} collinear segments")

    if args.dilate is not None:
        try:
            offset_value = float(args.dilate)
            logger.info(f"Dilating SVG by {offset_value}")
            for simple_shape in simple_shapes:
                dilate_shape(simple_shape, offset_value)
        except Exception as e:
            logger.error(f"  Dilation failed: {e}")
            return 4

    stroke_width = f"{args.dilate}{units}"
    stroke_opacity = "0.5"
    color_fill  = "none"
    color_inner = "#0000ff"  # blue
    color_outer = "#ff0000"  # red
    color_orig  = "#00ff00"   # green

    # Assign styles to outer and inner paths
    for simple_shape in simple_shapes:
        simple_shape.attributes["stroke-width"] = stroke_width
        simple_shape.attributes["stroke-opacity"] = stroke_opacity
        simple_shape.attributes["fill"] = color_fill
        if is_contained(simple_shape, simple_shapes):
            simple_shape.attributes["stroke"] = color_inner
        else:
            simple_shape.attributes["stroke"] = color_outer

    # If output path specified, write modified SVG
    if args.output:
        # Open a new output file for writing
        with open(args.output, "w") as f:
            f.write(preamble)

            for simple_shape in simple_shapes:
                path_str = to_svg_path(simple_shape)
                f.write(path_str)

            f.write("group id='original_paths'>\n")
            for path_node_str in path_node_strs:
                f.write(path_node_str)
            f.write("</group>\n")

            f.write(postamble)

    return 0


def find_units(input_path: pathlib.Path) -> str:
    """
    Find the units of the width attribute in an SVG file by searching line-by-line.
    Args:
        input_path (pathlib.Path): The path to the SVG file.
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


if __name__ == "__main__":
    raise SystemExit(main())
