import unittest
from safebridge import gis_ops

from shapely.geometry import Polygon, LineString, Point
from shapely.wkb import dumps as wkb_dumps
from random import randint
import logging
logging.disable(logging.INFO)
logging.disable(logging.ERROR)

class TestGisOps(unittest.TestCase):
    def setUp(self):
        """Set up a simple polygon and line for testing."""
        self.polygon = wkb_dumps(Polygon([(0, 0), (2, 0), (2, 2), (0, 2), (0, 0)]))
        self.buffer = Polygon([(0, 0), (2, 0), (2, 2), (0, 2), (0, 0)]).buffer(6)
        self.line = wkb_dumps(LineString([(1, -1), (1, 3)]))
        self.points = [Point(1,0.5), Point(1, 1.5)]

    def test_extract_intersecting_edges(self):
        """Test extracting edges of a polygon that intersect with a line."""
        edges = gis_ops.extract_intersecting_edges(self.polygon, self.line)
        self.assertEqual(len(edges), 2)
        self.assertTrue(all(isinstance(edge, LineString) for edge in edges))

    def test_extract_intersecting_edges_no_intersection(self):
        """Test extracting edges when there is no intersection."""
        line_no_intersection = wkb_dumps(LineString([(3, 3), (4, 4)]))
        edges = gis_ops.extract_intersecting_edges(self.polygon, line_no_intersection)
        self.assertEqual(len(edges), 0)

    def test_extend_line(self):
        """Test extending a line to intersect with a polygon."""
        edges = gis_ops.extract_intersecting_edges(self.polygon, self.line)

        extended_lines = [gis_ops.extend_line(edge, 1e3) for edge in edges]
        self.assertTrue(extended_lines == [LineString([(1002, 2), (-1000, 2)]), LineString([(-1000, 0), (1002, 0)])])

    def test_move_lines_to_points(self):
        """Test moving lines to points."""
        edges = gis_ops.extract_intersecting_edges(self.polygon, self.line)
        moved_lines = gis_ops.move_lines_to_points(edges, self.points)
        
        self.assertEqual(moved_lines, [LineString([(2,0.5), (0,0.5)]), LineString([(-0, 1.5), (2, 1.5)])])
        
    
    def test_multisplit(self):
        split_lines = [gis_ops.extend_line(edge, 1e3) for edge in gis_ops.extract_intersecting_edges(self.polygon, self.line)]
        multi_split = gis_ops.multisplit(self.buffer, split_lines)
        self.assertTrue(isinstance(multi_split, list))
        self.assertTrue(len(multi_split)  == 3)


if __name__ == '__main__':
    unittest.main()