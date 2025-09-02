import os
import sys
sys.path.append('./src')
import unittest
from safebridge.pipeline import DBPipeline, DBQueries
from safebridge.data import BridgeDamage, Deck, Axis, Support, Ascending, Descending
from safebridge.damage_assessment import DamageAssessment
from duckdb import DuckDBPyConnection
from duckdb.duckdb import ParserException
import logging
logging.disable(logging.INFO)
logging.disable(logging.ERROR)

class TestDamageAssessment(unittest.TestCase):
    def setUp(self):
        
        self.deck = Deck(source_file='./examples/toy_data/deck.shp')
        self.axis = Axis(source_file='./examples/toy_data/axis.shp')
        self.support = Support(source_file='./examples/toy_data/support.shp', source_projection="EPSG:3857")
        self.ascending = Ascending(
            source_file='./examples/toy_data/ascending_data.csv',
            lat_field="lat",
            lon_field="lon",
            orbit_azimuth=45.0,
            incidence_angle=30.0,
        )
        self.descending = Descending(
            source_file='./examples/toy_data/descending_data.csv',
            lat_field="lat",
            lon_field="lon",
            orbit_azimuth=135.0,
            incidence_angle=60.0,
        )

    def tearDown(self):
        del self.deck, self.axis, self.support, self.ascending, self.descending
    
    @classmethod
    def setUpClass(cls):
        cls.damage = DamageAssessment(
            deck=Deck(source_file='./examples/toy_data/deck.shp'),
            axis=Axis(source_file='./examples/toy_data/axis.shp'),
            support=Support(source_file='./examples/toy_data/support.shp', source_projection="EPSG:3857"),
            ascending = Ascending(
                source_file='./examples/toy_data/ascending_data.csv',
                lat_field="lat",
                lon_field="lon",
                orbit_azimuth=45.0,
                incidence_angle=30.0,
            ),
            descending = Descending(
                source_file='./examples/toy_data/descending_data.csv',
                lat_field="lat",
                lon_field="lon",
                orbit_azimuth=135.0,
                incidence_angle=60.0,
            )
        )

    @classmethod
    def tearDownClass(cls):
        del cls.damage
        
    def test_nonvalid_source_file(self):
        self.deck.source_file = []
        with self.assertRaises(ValueError):
            DamageAssessment(
                deck=self.deck,
                axis=self.axis,
                support=self.support,
                ascending=self.ascending,
                descending=self.descending,
            )

        with self.assertRaises(ValueError):
            self.deck.source_file = ""
            DamageAssessment(
                deck=self.deck,
                axis=self.axis,
                support=self.support,
                ascending=self.ascending,
                descending=self.descending,
            )

    def test_nonvalid_projection(self):
        self.deck.source_projection = 1
        with self.assertRaises(ValueError):
            DamageAssessment(
                deck=self.deck,
                axis=self.axis,
                support=self.support,
                ascending=self.ascending,
                descending=self.descending,
            )
        self.deck.source_projection = ""
        with self.assertRaises(ValueError):
            DamageAssessment(
                deck=self.deck,
                axis=self.axis,
                support=self.support,
                ascending=self.ascending,
                descending=self.descending,
            )
        
    def test_nonvalid_tablename(self):
        self.deck.table_name = " "
        with self.assertRaises(ValueError):
            DamageAssessment(
                deck=self.deck,
                axis=self.axis,
                support=self.support,
                ascending=self.ascending,
                descending=self.descending,
            )
        self.deck.table_name = []
        with self.assertRaises(ValueError):
            DamageAssessment(
                deck=self.deck,
                axis=self.axis,
                support=self.support,
                ascending=self.ascending,
                descending=self.descending,
            )
    
    def test_00_with_unit(self):
        self.ascending.unit = "dm"
        with self.assertRaises(ValueError):
            DamageAssessment(
                deck=self.deck,
                axis=self.axis,
                support=self.support,
                ascending=self.ascending,
                descending=self.descending,
            )
        for i, j in [("m",1),("cm",1/100),("mm", 1/1000)]:
            self.ascending.unit = i
            self.descending.unit = i
            
            damage = DamageAssessment(
                deck=self.deck,
                axis=self.axis,
                support=self.support,
                ascending=self.ascending,
                descending=self.descending
            )
            self.assertEqual(damage.damage.ascending.scaling_factor, j)
            self.assertEqual(damage.damage.descending.scaling_factor, j)
        
    

    def test_load_source_files(self):
        self.damage.load_source_files()
        self.assertIsInstance(self.damage.db.con, DuckDBPyConnection)

    def test_prepocess(self):
        with self.assertRaises(AssertionError):
            # self.damage.load_source_files()
            self.damage.preprocess(computational_projection=1, buffer_distance=1)

        with self.assertRaises(AssertionError):
            # self.damage.load_source_files()
            self.damage.preprocess(computational_projection="EPSG:3947", buffer_distance="1")

        with self.assertRaises(AssertionError):
            # self.damage.load_source_files()
            self.damage.preprocess(computational_projection="EPSG:3947", buffer_distance=-1)
        
        self.damage.preprocess(computational_projection="EPSG:28992", buffer_distance=10)
        
        tables = [i[0] for i in self.damage.db.con.execute("show tables").fetchall()]
        self.assertTrue('result_ns' in tables)
        self.assertTrue('result_ew' in tables)

    

class TestDamageAssessmentFilter(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        
        cls.damage = DamageAssessment(
            deck=Deck(source_file='./examples/toy_data/deck.shp'),
            axis=Axis(source_file='./examples/toy_data/axis.shp'),
            support=Support(source_file='./examples/toy_data/support.shp'),
            ascending = Ascending(
                source_file='./examples/toy_data/ascending_data.csv',
                lat_field="lat",
                lon_field="lon",
                orbit_azimuth=45.0,
                incidence_angle=30.0,
            ),
            descending = Descending(
                source_file='./examples/toy_data/descending_data.csv',
                lat_field="lat",
                lon_field="lon",
                orbit_azimuth=135.0,
                incidence_angle=60.0,
            )
        )
        cls.damage.load_source_files()
        cls.damage.preprocess(computational_projection="EPSG:28992", buffer_distance=10)

    @classmethod
    def tearDownClass(cls):
        del cls.damage

    def test_filter_single_condition(self):
        condition = ('lat', '>', 52.0)
        self.damage.filter(self.damage.damage.ascending, condition)
        result = self.damage.db.con.execute("SELECT COUNT(*) FROM proc_ascending").fetchone()[0]
        self.assertGreaterEqual(result, 0)

    def test_filter_multiple_conditions_and_logic(self):
        conditions = [('lat', '>', 52.0), ('lon', '<', 5.0), ('lon', '>', 4.0), ('lat', '<', 80.0)]
        self.damage.filter(self.damage.damage.descending, conditions, logic="AND")
        result = self.damage.db.con.execute("SELECT COUNT(*) FROM proc_ascending").fetchone()[0]
        self.assertGreaterEqual(result, 0)

    def test_filter_invalid_condition_format(self):
        condition = ('lat', '>',)
        with self.assertRaises(ValueError):
            self.damage.filter(self.damage.damage.ascending, condition)

    def test_filter_invalid_logic(self):
        conditions = [('lat', '>', 52.0), ('lon', '<', 5.0)]
        with self.assertRaises(ValueError):
            self.damage.filter(self.damage.damage.ascending, conditions, logic="INVALID")

    def test_filter_invalid_column(self):
        condition = ('invalid_column', '=', 'value')
        with self.assertRaises(ValueError):
            self.damage.filter(self.damage.damage.ascending, condition)

    def test_filter_invalid_operator(self):
        condition = ('lat', 'INVALID_OPERATOR', 52.0)
        with self.assertRaises(ValueError):
            self.damage.filter(self.damage.damage.ascending, condition)

    def test_filter_in_operator(self):
        condition = ('lat', 'IN', [52.0, 53.0])
        self.damage.filter(self.damage.damage.ascending, condition)
        result = self.damage.db.con.execute("SELECT COUNT(*) FROM proc_ascending").fetchone()[0]
        self.assertGreaterEqual(result, 0)

    def test_filter_in_operator_invalid_value(self):
        condition = ('lat', 'IN', 'invalid_value')
        with self.assertRaises(ValueError):
            self.damage.filter(self.damage.damage.ascending, condition)
    
    def test_filter_in_operator_valid_value(self):
        condition = ('lat', '=', None)
        self.damage.filter(self.damage.damage.ascending, condition)
        result = self.damage.db.con.execute("SELECT COUNT(*) FROM proc_ascending").fetchone()[0]
        self.assertGreaterEqual(result, 0)
    
    def test_filter_in_operator_invalid_value3(self):
        condition = ('lat', '=', "This isn't a valid operator")
        with self.assertRaises(ParserException):
            self.damage.filter(self.damage.damage.ascending, condition)
    
class TestDamageAssessment_da(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        
        cls.damage = DamageAssessment(
            deck=Deck(source_file='./examples/toy_data/deck.shp'),
            axis=Axis(source_file='./examples/toy_data/axis.shp'),
            support=Support(source_file='./examples/toy_data/support.shp', source_projection="EPSG:3857"),
            ascending = Ascending(
                source_file='./examples/toy_data/ascending_data.csv',
                lat_field="lat",
                lon_field="lon",
                orbit_azimuth=45.0,
                incidence_angle=30.0,
            ),
            descending = Descending(
                source_file='./examples/toy_data/descending_data.csv',
                lat_field="lat",
                lon_field="lon",
                orbit_azimuth=135.0,
                incidence_angle=60.0,
            )
        )
        cls.damage.load_source_files()
        cls.damage.preprocess(computational_projection="EPSG:28992", buffer_distance=10)
    
    @classmethod
    def tearDownClass(cls):
        del cls.damage   
    
    def test_assesstment_method(self):
        self.damage.assess_damage()
        result = self.damage.db.con.execute("SELECT COUNT(*) FROM result_ns").fetchone()[0]
        self.assertGreater(result, 0)
    
    def test_report_generation(self):
        self.damage.assess_damage()
        
        with self.assertRaises(ValueError):
            self.damage.generate_report('some_col_name')
        
        self.damage.generate_report()
        pdf_files = [f for f in os.listdir("safebridgeDB") if f.endswith(".pdf")]
        self.assertGreaterEqual(len(pdf_files), 1)
    
    def test_output_parquet(self): 
        self.damage.export_results(filetype="parquet")
        parquet_files = [f for f in os.listdir("safebridgeDB") if f.endswith(".parquet")]
        self.assertGreater(len(parquet_files), 1)
    
    def test_output_shapefile(self):
        self.damage.export_results(filetype="shapefile")
        shapefile_files = [f for f in os.listdir("safebridgeDB") if f.endswith(".shp")]
        self.assertGreater(len(shapefile_files), 1)

    def test_output_notsupportedfile(self):
        with self.assertRaises(ValueError):
            self.damage.export_results(filetype="unsupported_filetype")


class TestDamageAssessment_damage(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.damage = DamageAssessment(
            deck=Deck(source_file='./examples/toy_data/deck.shp'),
            axis=Axis(source_file='./examples/toy_data/axis.shp'),
            support=Support(source_file='./examples/toy_data/support.shp', source_projection="EPSG:3857"),
            ascending = Ascending(
                source_file='./examples/toy_data/ascending_data.csv',
                lat_field="lat",
                lon_field="lon",
                orbit_azimuth=45.0,
                incidence_angle=30.0,
            ),
            descending = Descending(
                source_file='./examples/toy_data/descending_data.csv',
                lat_field="lat",
                lon_field="lon",
                orbit_azimuth=135.0,
                incidence_angle=60.0,
            )
        )
        cls.damage.load_source_files()
        cls.damage.preprocess(computational_projection="EPSG:28992", buffer_distance=10)
    
    @classmethod
    def tearDownClass(cls):
        del cls.damage
    
    def test_duckdb_file_connection(self):
        dbpath = self.damage.db._db_path 
        self.damage.connect_duckdb_file(db_path = self.damage.db._db_path)
        tables = [i[0] for i in self.damage.db.con.execute("show tables").fetchall()]
        self.assertEqual('deck' in tables, True)


class TestDamageAssessment_unit(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.damage = DamageAssessment(
            deck=Deck(source_file='./examples/toy_data/deck.shp'),
            axis=Axis(source_file='./examples/toy_data/axis.shp'),
            support=Support(source_file='./examples/toy_data/support.shp', source_projection="EPSG:3857"),
            ascending = Ascending(
                source_file='./examples/toy_data/ascending_data.csv',
                lat_field="lat",
                lon_field="lon",
                orbit_azimuth=45.0,
                incidence_angle=30.0,
                unit="m"
            ),
            descending = Descending(
                source_file='./examples/toy_data/descending_data.csv',
                lat_field="lat",
                lon_field="lon",
                orbit_azimuth=135.0,
                incidence_angle=60.0,
                unit="cm"
            )
        )
        cls.damage.load_source_files()
        cls.damage.preprocess(computational_projection="EPSG:28992", buffer_distance=10)
        
    
    @classmethod
    def tearDownClass(cls):
        del cls.damage
    
    def test_unit_mismatch(self):
        with self.assertRaises(ValueError):
            self.damage.assess_damage()
        

if __name__ == '__main__':
    unittest.main()