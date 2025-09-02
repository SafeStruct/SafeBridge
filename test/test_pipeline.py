import unittest
from safebridge.pipeline import DBPipeline, DBQueries
from safebridge.data import BridgeDamage, Deck, Axis, Support, Ascending, Descending
from safebridge.damage_assessment import DamageAssessment
import logging
logging.disable(logging.INFO)
logging.disable(logging.ERROR)
       

class TestPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        deck = Deck(source_file='./examples/toy_data/deck.shp')
        axis = Axis(source_file='./examples/toy_data/axis.shp')
        support = Support(source_file='./examples/toy_data/support.shp', source_projection="EPSG:3857")
        ascending = Ascending(
            source_file='./examples/toy_data/ascending_data.csv',
            lat_field="lat",
            lon_field="lon",
            orbit_azimuth=45.0,
            incidence_angle=30.0
        )
        descending = Descending(
            source_file='./examples/toy_data/descending_data.csv',
            lat_field="lat",
            lon_field="lon",
            orbit_azimuth=135.0,
            incidence_angle=60.0
        )
        bridgedamage = BridgeDamage(deck=deck, axis=axis, support=support, ascending=ascending, descending=descending)
        cls.damage = DamageAssessment(deck=deck, axis=axis, support=support, ascending=ascending, descending=descending)
        cls.damage.db.setup()
        cls.pipeline = DBPipeline(bridgedamage, cls.damage.db.con)
        cls.damage.load_source_files()

    def test_01_load_source_files(self):
        tables = [i[0] for i in self.damage.db.con.execute('show tables').fetchall()]
        self.assertIn("deck", tables)
        self.assertIn("axis", tables)
        self.assertIn("support", tables)
        self.assertIn("ascending", tables)
        self.assertIn("descending", tables)

    def test_02_build_point_geometry(self):
        self.pipeline.build_point_geometry()
        attributes = self.pipeline.get_attributes("ascending")
        self.assertIn("geom", attributes)
        attributes = self.pipeline.get_attributes("descending")
        self.assertIn("geom", attributes)

    def test_03_build_process_tables(self):
        self.pipeline.build_process_tables("EPSG:3857")
        attributes = self.pipeline.get_attributes("proc_deck")
        self.assertIn("geom", attributes)

    def test_04_process_axis(self):
        self.pipeline.process_axis()
        attributes = self.pipeline.get_attributes("proc_axis")
        self.assertIn("length", attributes)
        self.assertIn("azimuth", attributes)

    def test_05_process_deck(self):
        self.pipeline.process_deck(buffer_distance=100.0)
        attributes = self.pipeline.get_attributes("proc_deck")
        self.assertIn("span_count", attributes)
        self.assertIn("buffer", attributes)

    def test_06_relate_deck_axis(self):
        self.pipeline.relate_deck_axis()
        attributes = self.pipeline.get_attributes("proc_deck")
        self.assertIn("orientation", attributes)
        attributes = self.pipeline.get_attributes("proc_axis")
        self.assertIn("rdeck", attributes)

    def test_07_create_sectors(self):        
        self.pipeline.create_sectors()
        attributes = self.pipeline.get_attributes("sectors")
        self.assertIn("sector_tag", attributes)
        self.assertIn("center", attributes)

    def test_08_relate_deck_pspoints(self):
        self.pipeline.relate_deck_pspoints()
        attributes = self.pipeline.get_attributes("proc_ascending")
        self.assertIn("rdeck", attributes)
        self.assertIn("rsector", attributes)
        attributes = self.pipeline.get_attributes("proc_descending")
        self.assertIn("rdeck", attributes)
        self.assertIn("rsector", attributes)
    
    def test_09_relate_axis_pspoints(self):
        self.pipeline.relate_axis_pspoints()
        attributes = self.pipeline.get_attributes("proc_ascending")
        self.assertIn("ndist_axis", attributes)
        self.assertIn("proj_axis", attributes)
        attributes = self.pipeline.get_attributes("proc_descending")
        self.assertIn("ndist_axis", attributes)
        self.assertIn("proj_axis", attributes)
    
    def test_10_deck_edge_control(self):
        self.pipeline.deck_edge_control(buffer_distance=100.0)
        attributes = self.pipeline.get_attributes("proc_deck")
        self.assertIn("edge_check", attributes)

    def test_11_init_result_table(self):
        self.pipeline.init_result_table()
        attributes_ew = self.pipeline.get_attributes("result_ew")
        attributes_ns = self.pipeline.get_attributes("result_ns")
        self.assertIn("rdeck", attributes_ew)
        self.assertIn("rdeck", attributes_ns)

    def test_12_get_ns_bridge_uid(self):
        uids = self.pipeline.get_ns_bridge_uid()
        self.assertIsInstance(uids, list)
        self.assertGreaterEqual(len(uids), 0)
    
    def test_13_get_ew_bridge_uid(self):
        uids = self.pipeline.get_ew_bridge_uid()
        self.assertIsInstance(uids, list)
        self.assertGreaterEqual(len(uids), 0)

class TestPipeline02(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        deck = Deck(source_file='./examples/toy_data/deck.shp')
        axis = Axis(source_file='./examples/toy_data/axis.shp')
        support = Support(source_file='./examples/toy_data/support.shp', source_projection="EPSG:3857")
        ascending = Ascending(
            source_file='./examples/toy_data/ascending_data.csv',
            lat_field="latitude",
            lon_field="longitude",
            orbit_azimuth=45.0,
            incidence_angle=30.0
        )
        descending = Descending(
            source_file='./examples/toy_data/descending_data.csv',
            lat_field="latitude",
            lon_field="lonngitude",
            orbit_azimuth=135.0,
            incidence_angle=60.0
        )
        bridgedamage = BridgeDamage(deck=deck, axis=axis, support=support, ascending=ascending, descending=descending)
        cls.damage = DamageAssessment(deck=deck, axis=axis, support=support, ascending=ascending, descending=descending)
        cls.damage.db.setup()
        cls.pipeline = DBPipeline(bridgedamage, cls.damage.db.con)
        cls.damage.load_source_files()


    def test_build_point_geometry_valueerror(self):
        with self.assertRaises(ValueError):
            self.pipeline.build_point_geometry()

class TestDBQueries(unittest.TestCase):
    @classmethod
    
    def setUpClass(cls):
        cls.queries = DBQueries()
    
    def test_deck_geometry(self):
        query = self.queries.deck_geometry(deckuid=1, deck_table="proc_deck")
        expected_query = "SELECT ST_AsWKB(geom) FROM proc_deck WHERE uid = 1"
        self.assertEqual(query, expected_query)
    
    def test_buffer_geometry(self):
        query = self.queries.buffer_geometry(deckuid=1, table_name="proc_deck")
        expected_query = "SELECT ST_AsWKB(buffer) FROM proc_deck WHERE uid = 1"
        
        self.assertEqual(query, expected_query)
    
    def test_sector_geometry(self):
        query = self.queries.sector_geometry(deckuid=1)
        expected_query = "SELECT ST_AsWKB(geom) FROM sectors WHERE rdeck = 1"
        self.assertEqual(query, expected_query)
    
    def test_support_geometry(self):
        query = self.queries.support_geometry(deckuid=1, table_name="proc_support")
        expected_query = "SELECT ST_AsWKB(geom) FROM proc_support WHERE rdeck = 1"
        self.assertEqual(query, expected_query)
    
    def test_axis_geometry(self):
        query = self.queries.axis_geometry(deckuid=1, table_name="proc_axis")
        expected_query = "SELECT ST_AsWKB(geom) FROM proc_axis WHERE rdeck = 1"
        self.assertEqual(query, expected_query)
    
    def test_deck_edge(self):
        query = self.queries.deck_edge(deckuid=1, table_name="deck")
        expected_query = "SELECT ST_AsWKB(ST_StartPoint(deck_edge)) as st, ST_AsWKB(ST_EndPoint(deck_edge)) as ed, FROM proc_deck WHERE uid = 1"
        self.assertEqual(query, expected_query)
    
    def test_scatter_geometry(self):
        query = self.queries.scatter_geometry(deckuid=1, table_name="proc_ascending")
        expected_query = "SELECT ST_X(geom) as x, ST_Y(geom) as y FROM proc_ascending WHERE rdeck = 1"
        self.assertEqual(query, expected_query)
    
    def test_projected_scatters(self):
        query = self.queries.projected_scatters(deckuid=1, table_name="proc_descending")
        expected_query = "SELECT ST_X(proj_axis) as x, ST_Y(proj_axis) as y FROM proc_descending WHERE rdeck = 1"
        self.assertEqual(query, expected_query)
    
    def test_buffer_edge(self):
        query = self.queries.buffer_edge(deckuid=1, axis_name="proc_axis", deck_name="proc_deck")
        expected_query = """
                SELECT
                    ST_Distance(ST_StartPoint(deck.buffer_edge), ST_StartPoint(axis.geom)) / axis.length as p1,
                    ST_Distance(ST_EndPoint(deck.buffer_edge), ST_StartPoint(axis.geom)) / axis.length as p2,
                FROM (SELECT * FROM proc_axis WHERE rdeck = 1) as axis
                JOIN proc_deck as deck
                ON axis.rdeck = deck.uid
                """
        self.assertEqual(query.strip(), expected_query.strip())
    
    def test_deck_edge_graph(self):
        query = self.queries.deck_edge_graph(deckuid=1, axis_name="proc_axis", deck_name="proc_deck")
        expected_query = """
                SELECT
                    ST_Distance(ST_StartPoint(deck.deck_edge), ST_StartPoint(axis.geom)) / axis.length as p1,
                    ST_Distance(ST_EndPoint(deck.deck_edge), ST_StartPoint(axis.geom)) / axis.length as p2,
                FROM (SELECT * FROM proc_axis WHERE rdeck = 1) as axis
                JOIN proc_deck as deck
                ON axis.rdeck = deck.uid
                """
        self.assertEqual(query.strip(), expected_query.strip())
    
    def test_scatter_graph(self):
        query = self.queries.scatter_graph(deckuid=1, table_name="proc_ascending", name_fields=["field1", "field2"])
        expected_query = """ 
                SELECT 
                    proc_scatter.ndist_axis as x,
                    scatter.field2 - scatter.field1  as y,
                FROM (SELECT * FROM proc_proc_ascending WHERE rdeck = 1) as proc_scatter
                JOIN proc_ascending as scatter
                ON proc_scatter.uid = scatter.uid
                """
        self.assertEqual(query.strip(), expected_query.strip())
    
    def test_support_graph(self):
        query = self.queries.support_graph(deckuid="1", axis_name="proc_axis", support_name="proc_support")
        expected_query = """
              SELECT *
              FROM (
                SELECT 
                    ST_Distance(ST_StartPoint(axis.geom), ST_Centroid(ST_Intersection(support.geom, axis.geom))) / axis.length as p1,
                FROM (SELECT * FROM proc_proc_support WHERE rdeck = 1) AS support
                JOIN proc_proc_axis AS axis
                ON support.rdeck = axis.rdeck
              )
              WHERE p1 IS NOT NULL
              """
        self.assertEqual(query.strip(), expected_query.strip())


if __name__ == '__main__':
    unittest.main()