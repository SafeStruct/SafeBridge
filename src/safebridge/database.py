import os
import duckdb
from datetime import datetime

class DataBase:
    """ DataBase is a class that provides an interface for managing a DuckDB database. 
    It includes methods for setting up the database, loading files into tables, and connecting to existing DuckDB files. The class also supports spatial extensions for geospatial data processing.
    
    Methods
    -------
    __init__():
        Initializes the DataBase class and sets up the initial state.
    setup():
        Sets up the DuckDB connection and loads the spatial extension. Creates a new DuckDB database file in a specified directory.
    init_db_dir() -> str:
        Initializes the database directory and creates a new DuckDB database file with a unique timestamp-based name.
    load_file(source_file: str, table_name: str):
        Loads a file into the DuckDB database. Supports CSV and Shapefile formats. Creates a new table, adds a UID column, and generates unique IDs.
    connect_duckdbfile(duckdb_file: str):
        Connects to an existing DuckDB database file. Closes the current connection and establishes a new one to the specified file.
    
    Attributes
    ----------
    con : duckdb.DuckDBPyConnection or None
        The DuckDB connection object. Initially set to None.
    _db_path : str
        The path to the DuckDB database file. Set during setup or connection to an existing file.
    """
    
    def __init__(self):
        """Initialize the DataBase class."""
        
        self.con = None

    def setup(self):
        """Set up the DuckDB connection and load the spatial extension.

        This method initializes the database directory, creates a new DuckDB database file, and establishes a connection to it. It also loads the spatial extension if available.

        Raises
        ------
            RuntimeError: If the spatial extension fails to load.
        """

        self._db_path = self.init_db_dir()
        self.con = duckdb.connect(self._db_path)
        # Load the spatial extension if available
        self.con.load_extension("spatial")  

    def init_db_dir(self) -> str:
        """ Initialize the database directory and create a new DuckDB database file.

        This method creates a directory for the DuckDB database if it does not exist.
        It generates a new DuckDB database file with a timestamp to ensure uniqueness.
        
        Returns
        -------
            str: The full path to the newly created DuckDB database file.
        
        Raises
        ------
            RuntimeError: If the spatial extension fails to load.
        
        """

        db_folder = "safebridgeDB"
        os.makedirs(db_folder, exist_ok=True) # if the folder does not exist, it will be created if exist it will not raise an error
        
        # Unique name assingment in here 
        fname = os.path.join(db_folder, 
                             f"{datetime.now().strftime('%Y%m%d%H%M')}.duckdb"
                             )
        # return the full path to the database file in case it will be used
        return fname
    
    def load_file(self, source_file: str, table_name: str):
        """ Load a file into the DuckDB database.  
        
        This method checks the file extension to determine the appropriate loading method. It supports CSV and Shapefile formats. For CSV files, it uses `read_csv_auto`, and for Shapefiles, it uses `ST_Read`. It creates a new table with the specified name. If the table already exists, it will be dropped and recreated. This method also creates a sequence for the table to generate unique IDs and adds a UID column to the table. If the file does not exist, it raises a `FileNotFoundError`. If the file format is unsupported, it raises a `ValueError`.
        
        Arguments
        ---------
        source_file : str
            The path to the file to load.
        table_name : str
            The name of the table to create or append to.
        
        Raises
        ------
        FileNotFoundError: 
            If the specified file does not exist.
        ValueError
            If the file format is unsupported.
        """

        # Validate inputs
        for check in [source_file, table_name]:
            if not isinstance(check, str):
                raise ValueError(f"{check} must be a string.")
            if check == "":
                raise ValueError(f"{check} must not be empty.")
            if check.isspace():
                raise ValueError(f"{check} must not be an empty or whitespace-only string.")
        
        
        if source_file.endswith('.csv'):
            load_method = "read_csv_auto"
        elif source_file.endswith('.shp'):
            load_method = "ST_Read"
        else:
            raise ValueError("Unsupported file format. Only CSV and Shapefile are supported.")
        
        if not os.path.exists(source_file):
            raise FileNotFoundError(f"File {source_file} does not exist.")
        
        self.con.execute(f"""
                         CREATE OR REPLACE TABLE {table_name} AS SELECT * FROM {load_method}('{source_file}');
                         CREATE OR REPLACE SEQUENCE {table_name}_id;
                         ALTER TABLE {table_name} ADD COLUMN uid INTEGER DEFAULT nextval('{table_name}_id');
                         """
                         )

    def connect_duckdbfile(self, duckdb_file: str):
        """Connect to an existing DuckDB database file.

        This method allows you to connect to a DuckDB database file that already exists.
        It will close the current connection and establish a new one to the specified file.

        Arguments
        ---------
        duckdb_file : str
            The path to the DuckDB database file to connect to.
        Raises
        ------
        FileNotFoundError: 
            If the specified DuckDB file does not exist.
        ValueError:
            If the duckdb_file is not a valid string or is empty.
        RuntimeError:
            If the spatial extension fails to load.
        """
        if not os.path.exists(duckdb_file):
            raise FileNotFoundError(f"Database file {duckdb_file} does not exist.")
        
        self._db_path = os.path.splitext(duckdb_file)[0]
        self.con = duckdb.connect(duckdb_file)
        self.con.load_extension("spatial")