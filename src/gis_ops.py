"""
This module provides functions to perform various GIS operations such as extracting edges from polygons, extending lines, moving lines to points, and splitting polygons with multiple lines. It uses the Shapely library for geometric operations and supports WKB (Well-Known Binary) format for input geometries.
"""
from shapely.ops import split
from shapely.affinity import scale, translate
from shapely.geometry import Polygon, LineString, Point
from shapely.wkb import loads as wkbloads


def extract_intersecting_edges(polygon: bytes, line: bytes) -> list[LineString]:
    """ Extracts the edges of a polygon that intersect with a given line.
    
    Arguments
    ---------
    polygon : bytes
        The Well-Known-Binary (WKB) representation of the polygon.
    line : bytes
        The WKB representation of the line.
    
    Returns
    -------
    list[LineString] : A list of LineString objects representing the intersecting edges.
    """
    polygon, line = wkbloads(polygon), wkbloads(line)
    edges = []
    boundary_coords = list(polygon.exterior.coords)

    for i in range(len(boundary_coords) - 1):
        edge = LineString([boundary_coords[i], boundary_coords[i + 1]])
        if edge.intersects(line):
            edges.append(edge)
    
    return sort_by_centroid(edges)

def sort_by_centroid(edges: list[LineString]) -> list[LineString]:
    """ Sorts a list of LineString objects by their centroid coordinates.
    
    Arguments
    ---------
    edges : [LineString]
        A list of LineString objects to be sorted.
    
    Returns
    -------
    list[LineString] : A sorted list of LineString objects based on their centroid coordinates.
    """
    centroids = [g.centroid for g in edges]
    sorted_cents = sorted(centroids, key=lambda point: (-point.y, point.x))
    sorted_geom = []
    for i in sorted_cents:
        sorted_geom.append(edges[centroids.index(i)])
    return sorted_geom

def extend_line(line:LineString, extend_length:float) -> LineString:
    """ Extends a line by a specified length on both ends.

    Arguments
    ---------
    line : LineString
        The LineString object to be extended.
    extend_length : float
        The length by which to extend the line on both ends.
    
    Returns
    -------
    LineString : A new LineString object that has been extended by the specified length.
    """

    centroid = line.centroid
    factor = (line.length + 2 * extend_length) / line.length
    extended_line = scale(line, xfact=factor, yfact=factor, origin=centroid)
    return extended_line

def move_line2point(line:LineString, point:Point) -> LineString:
    """ Moves a line to a point.
    
    Arguments
    ---------
    line : LineString
        The LineString object to be moved.
    point : Point
        The Point object to which the line will be moved.
    
    Returns
    -------
    LineString : A new LineString object that has been moved to the specified point.
    """
    line_cent = line.centroid
    dx = point.x - line_cent.x
    dy = point.y - line_cent.y
    moved_line = translate(line, xoff=dx, yoff=dy)
    return moved_line

def move_lines_to_points(lines: list[LineString], points: list[Point]) -> list[LineString]:
    """ Moves lines to the given points. 

    Arguments
    ---------
    lines : list[LineString]
        A list of LineString objects to be moved.
    points : list[Point]
        A list of Point objects to which the lines will be moved.

    Returns
    -------
    list[LineString] : A list of moved LineString objects.
    """

    return [move_line2point(line, point) for line, point in zip(lines, points)]


def multisplit(polygon : Polygon, lines: list[LineString]) -> list[Polygon]:
    """ Split a polygon by multiple lines.

    Arguments
    ---------
    polygon : Polygon
        The polygon to be split.
    lines : list[LineString]
        A list of LineString objects to split the polygon.

    Returns
    -------
    list[Polygon]: A list of resulting polygons after the split.
    """
    
    collection = []
    temp = split(polygon, lines[0])
    for geom in temp.geoms:
        if geom.intersection(lines[1]).is_empty:
            collection.append(geom)
        else:
            geom2 = split(geom, lines[1])
            collection.append(geom2.geoms[0])
            collection.append(geom2.geoms[1])
    return sort_by_centroid(collection)
