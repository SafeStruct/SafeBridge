import os
import unittest
from duckdb import DuckDBPyConnection
from safebridge.database import DataBase
import logging
logging.disable(logging.INFO)
logging.disable(logging.ERROR)

test_file            = "examples/toy_data/ascending_data.csv"
test_shape           = "examples/toy_data/deck.shp"
invalid_file         = "examples/toy_data/invalid_data.txt"
non_existent_file    = "examples/toy_data/non_existent.csv"

class TestDataBase(unittest.TestCase):
    
    def setUp(self):
        self.db = DataBase()
    
    def test_basic_setup(self):
        self.assertTrue(self.db.con is None)
    
    def test_setup_method(self):
        self.db.setup()
        self.assertIsInstance(self.db.con, DuckDBPyConnection)
        self.assertTrue(os.path.exists(self.db._db_path))
        self.assertTrue(self.db._db_path.endswith('.duckdb'))
        self.assertTrue(self.db._db_path.startswith('safebridgeDB/'))
    
        
    def test_load_file(self):
        self.db.setup()
        

        # Test loading a non-existent file
        with self.assertRaises(FileNotFoundError):
            self.db.load_file(non_existent_file, "test_table")
        

        # Test loading a file with invalid format
        with self.assertRaises(ValueError):
            self.db.load_file(invalid_file, "test_table")
        # Test loading a valid file and invalid table name
        with self.assertRaises(ValueError):
            self.db.load_file(invalid_file, 1)
        # test with whitespace as table name
        with self.assertRaises(ValueError):
            self.db.load_file(invalid_file, " ")
        # test with empty string as table name
        with self.assertRaises(ValueError):
            self.db.load_file(invalid_file, "")
        
        
        # # Test loading a file with an empty table name
        with self.assertRaises(ValueError):
            self.db.load_file(test_file, "")
        
        # Test loading valid file
        self.db.load_file(test_shape, "shape_table")
        self.assertTrue("shape_table" in self.db.con.sql("SHOW TABLES").fetchnumpy()['name'].tolist())

    # test loading recently created duckdbfile
    def test_connect_duckdbfile(self):
        self.db.setup()
        self.db.load_file(test_file, "test_table")

        # # reconnecting
        self.db.connect_duckdbfile(self.db._db_path)

        self.assertTrue("test_table" in self.db.con.sql("SHOW TABLES").fetchnumpy()['name'].tolist())

        with self.assertRaises(FileNotFoundError):
            self.db.connect_duckdbfile("./test/fixtures/non_existent.duckdb")



if __name__ == '__main__':
    unittest.main()