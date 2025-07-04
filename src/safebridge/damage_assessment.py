import re
import time
import warnings
import dateutil.parser as dsparser

from .data import Deck, Axis, Support, Ascending, Descending, BridgeDamage
from .pipeline import DBPipeline, DBQueries
from .database import DataBase
from .solvers import NS_Solver, EW_Solver
from .plotter import Plotter

from numpy import ndarray, array
from shapely.wkb import loads as wkbloads
from typing import Union
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore")

class DamageAssessment:
    """
    The DamageAssessment class is responsible for assessing the damage to a bridge based on the provided data.
    It initializes with the deck, axis, support, and persistent scatter data for both ascending and descending orbits. 
    It provides methods to process the provided data, damage assessement and generate reports. This class is designed to work with the SafeBridge data model, which includes various data classes for bridge components.
    
    Attributes
    ----------
    damage : BridgeDamage
        An instance of the BridgeDamage class containing the deck, axis, support, and persistent scatter data.
    
    db : DataBase
        An instance of the DataBase class for managing database connections and operations.
    query : DBQueries
        An instance of the DBQueries class for executing SQL queries.
    _plotter : Plotter
        An instance of the Plotter class for visualizing the data and results.
    dbpipeline : DBPipeline
        An instance of the DBPipeline class for processing the data through various steps.
    _buf_size : float
        The buffer distance used for processing geometries.
    
    Methods
    -------
    connect_duckdb_file(db_path: str)
        Connects to a DuckDB file at the specified path.
    load_source_files()
        Loads the source files for deck, axis, support, ascending, and descending data into the database.
    preprocess(computational_projection: str, buffer_distance: float)
        Preprocesses the data for damage assessment by building geometries, reprojecting them, and relating them.
    filter(safebridge_data: Union[Ascending, Descending, Deck, Axis, Support], condition: Union[tuple, list[tuple]], logic: str = "AND")
        Filters the data in the specified table based on the provided conditions and logic.
    assess_damage()

    """
    def __init__(self, deck:Deck, axis:Axis, support:Support, ascending = Ascending, descending = Descending):
        """
        Initialize the DamageAssessment with deck, axis, support and persistent scatter data.

        Arguments
        ----------
        deck : Deck
            The deck data.
        axis : Axis
            The axis data.
        support : Support
            The support data.
        ascending : Ascending
            The persistent scatter data for ascending orbit.
        descending : Descending
            The persistent scatter data for descending orbit.
        """
        
        for obj in [ascending, descending]:
            if obj.unit == "mm":
                obj.scaling_factor = 1/1000
            elif obj.unit == "cm":
                obj.scaling_factor = 1/100
            elif obj.unit == "m":
                obj.scaling_factor = 1.0
            else:
                raise ValueError(f"Invalid scaling factor for {obj} data. Use 'mm', 'cm', or 'm'.")
            
            # rest of the type and value checks in here for ascending and descending points
            
        for obj in [deck, axis, support, ascending, descending]:
            if not isinstance(obj.source_file, str):
                raise ValueError(f"The source_file for {obj.table_name} must be a string.")
            if not isinstance(obj.source_projection, str):
                raise ValueError(f"The source_projection for {obj.table_name} must be a string.")
            if obj.table_name in ["", " "]:
                raise ValueError(f"The table_name for {obj} must not be empty or contain only spaces.")
        
        self.damage = BridgeDamage(deck, axis, support, ascending, descending)
        self.db = DataBase()
        self.query = DBQueries()
        self._plotter = Plotter()
    
    def connect_duckdb_file(self, db_path: str):
        """ Connects to a DuckDB file using the db attribute of the class.
        This method establishes a connection to the DuckDB file at the specified path.
                
        Arguments
        ----------
        db_path : str 
            The path to the DuckDB file.
        """
        self.db.connect_duckdbfile(db_path)

    def load_source_files(self):
        """
        Loads the source files for deck, axis, support, ascending, and descending data.
        """
        # before the loading the new files to db initialize the db connection with a new DuckDB file
        if self.db.con is None:
            self.db.setup()
            
        for i in self.damage.__dataclass_fields__.keys():
            obj = self.damage.__getattribute__(i)
            st = time.time()
            self.db.load_file(obj.source_file, obj.table_name)
            print(f"{obj.table_name} dataset has been loaded to db from {obj.source_file} in {time.time() - st:.2f} seconds.")

    def preprocess(self, computational_projection: str, buffer_distance: float):
        """ Preprocess the data for damage assessment.

        Arguments
        ----------
        computational_projection : str
            The coordinate reference system for computations.
        buffer_distance : float
            The distance to buffer geometries.
        """
        assert isinstance(computational_projection, str), "Computational projection must be a string representing the EPSG code."
        assert isinstance(buffer_distance, (int, float)), "Buffer distance must be a numeric value."
        
        if buffer_distance <= 0:
            raise ValueError("Buffer distance must be greater than 0.")
        if computational_projection is None:
            raise ValueError("Computational projection must be specified.")
        
        self._buf_size = buffer_distance
        # Initialize the DBPipeline with the damage data and database connection
        self.dbpipeline = DBPipeline(self.damage, self.db.con)

        # 1-building point geometries for the ascending and descending constalliations
        self.dbpipeline.build_point_geometry()
        # 2-creating the proc_{table_name} tables with the reprojection of geometries to computational projection
        self.dbpipeline.build_process_tables(computational_projection)    
        # 3-reodering the axis vertices, calculating the length and azimuth
        self.dbpipeline.process_axis()
        # 4-calculating the deck span counts, generating the deck buffer geometries
        self.dbpipeline.process_deck(buffer_distance)
        # 5-relating the axis and the deck geometries, calculating the deck length
        self.dbpipeline.relate_deck_axis()
        # 6-creating sector geometries from the deck geometries
        self.dbpipeline.create_sectors()
        # 7-deck and sector assignment to the points
        self.dbpipeline.relate_deck_pspoints()
        # 8 project points on tho the axis line record the projected point as proj_axis and claculate the normalized distance along the axis line stating from the start point of the axis line as ndist_axis field
        self.dbpipeline.relate_axis_pspoints()
        # 9 check if there is at least one projected point for both orbital orientation within the radius of buffer_distance/2 at the both edge of the deck geometry
        self.dbpipeline.deck_edge_control(buffer_distance)
        self.dbpipeline.init_result_table()
    
    def filter(self, 
               safebridge_data: Union[Ascending, Descending, Deck, Axis, Support], 
               condition: Union[tuple, list[tuple]], 
               logic:str = "AND"
               ) -> None:
        """ Filter the data based on specified criteria.

        This method allows you to filter the data in the specified table based on the provided conditions.
        The conditions can be a single tuple or a list of tuples, where each tuple contains the column name, operator, and value to filter by. The method will apply the filter to table named `proc_{safebridge_data.table_name}`. IF tyhe user wants to start from scractch for the filter one needs to rerun the preprocess steps. 

        Arguments
        ----------
        safebridge_data : Union[Ascending, Descending, Deck, Axis, Support]
            The data to filter. Should be single instance of one of the data classes.
        condition : tuple
            A tuple containing the column name, opearator and value to filter by (column, operator, value)
        logic : str
            The logical operator to use for combining conditions. Defaults to "AND". Can be "AND" or "OR".
        Raises
        ------
            ValueError: If the condition format is invalid or if the logic is not "AND" or "OR".
        

        Example:
        filter(safebridge_data, ('column_name', '=', 'value'), logic='AND')
        filter(safebridge_data, [('column_name1', '=', 'value1'), ('column_name2', '>', 10)], logic='OR')
        filter(safebridge_data, [('column_name1', '=', 'value1'), ('column_name2', '>', 10), ('column_name3', 'LIKE', '%pattern%')], logic='AND')    
        """ 
        
        if isinstance(condition, tuple) and len(condition) == 3:
            conditions = [condition]
        elif isinstance(condition, list) and all(isinstance(c, tuple) and len(c) == 3 for c in condition):
            conditions = condition
        else:
            raise ValueError("Invalid condition format. Use (col, op, val) or list of them.")

        logic = logic.upper()
        if logic not in {"AND", "OR"}:
            raise ValueError("Logic must be 'AND' or 'OR'.")

        # valid_columns = get_column_names(safebridge_data.table_name, self.db.con)
        valid_columns = self.dbpipeline.get_attributes(safebridge_data.table_name)
        allowed_operators = {"=", "!=", "<", "<=", ">", ">=", "LIKE", "IN"} #"ILIKE"

        clause_parts = []
        ## inner method
        def _quote(value):
            """
            Quote the value for SQL query.
            
            Args:
                value: The value to quote.
            
            Returns:
                str: The quoted value.
            """
            if value is None:
                return "NULL"
            elif isinstance(value, (int, float)):
                return str(value)
            elif isinstance(value, (list, tuple)):
                return "(" + ", ".join(_quote(v) for v in value) + ")"
            else:
                escape = str(value).replace("'", "''")
                return f"{escape}"  # Basic SQL string escape    

        for col, op, val in conditions:
            if col not in valid_columns:
                raise ValueError(f"The {col} column does not exist in the {safebridge_data.table_name} table.")
            if op.upper() not in allowed_operators:
                raise ValueError(f"THE {op} operator is not allowed. Use one of {allowed_operators}.")
            if op.upper() == "IN":
                if not isinstance(val, (list, tuple)):
                    raise ValueError(f"'IN' operator requires a list or tuple value")
                clause_parts.append(f"{col} IN {_quote(val)}")
            else:
                clause_parts.append(f"{col} {op} {_quote(val)}")

        where_sql = f" {logic} ".join(clause_parts)
        final_query =  f"SELECT uid FROM {safebridge_data.table_name} WHERE {where_sql}"
        self.db.con.execute(f"""
            CREATE OR REPLACE TABLE proc_{safebridge_data.table_name} AS
            SELECT * FROM proc_{safebridge_data.table_name} WHERE uid IN ({final_query})
        """)
        
    def assess_damage(self):
        """ Assess damage based on the processed data.

        This method evaluates the damage to the bridge based on the ascending and descending data. It uses the full time range of the both ascending and descending datasets to determine the extent of damage. 
        It will assess all available data for `NS` oriented bridges for `EW` oriented ones it will use the overlapping time period of the ascending and descending data. The method will perform the necessary calculations to determine the extent of damage and will store the results in the database with a table called `result`.
        
        """
        
        timeOverlapInfo = self._get_timeoverlap()
        ns_decks = self.dbpipeline.get_ns_bridge_uid()

        st = time.time()
        for deckUid in ns_decks:
            ns_solver = NS_Solver(self._ns_solver_data(deckUid, timeOverlapInfo))

            self.db.con.execute(f"""
                INSERT INTO result (rdeck,
                                orient,
                                tilt_asc, 
                                defl_asc, 
                                tilt_dsc, 
                                defl_dsc, 
                                ns_quadratic_asc_x, 
                                ns_quadratic_asc_y,
                                ns_quadratic_dsc_x,
                                ns_quadratic_dsc_y, 
                                ns_analytical_asc_y,
                                ns_analytical_dsc_y) 
                                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) """, (
                    deckUid,
                    "NS",
                    ns_solver.quadratic_tilt('ascending') * self.damage.ascending.scaling_factor,
                    ns_solver.quadratic_deflection('ascending'),
                    ns_solver.quadratic_tilt('descending') * self.damage.descending.scaling_factor,
                    ns_solver.quadratic_deflection('descending'),
                    ns_solver._quadratic_x('ascending'),
                    ns_solver._quadratic_y('ascending'),
                    ns_solver._quadratic_x('descending'),
                    ns_solver._quadratic_y('descending'),
                    ns_solver.analytical_curve('ascending'),
                    ns_solver.analytical_curve('descending')
                )
            )
        print(f"NS solver completed in {time.time() - st:.2f} seconds.")
            
        ew_solver = EW_Solver(
            timeOverlapInfo,
            self.damage.ascending.incidence_angle,
            self.damage.descending.incidence_angle,
            self.damage.ascending.orbit_azimuth,
            self.damage.descending.orbit_azimuth,
            )
        
        ew_decks = self.dbpipeline.get_ew_bridge_uid()
        st1 = time.time()
        for deckUid in ew_decks:
            # QUERY
            sectors = self.db.con.sql(f""" SELECT uid, sector_tag, ndist FROM sectors WHERE rdeck = {deckUid} ORDER BY ndist ASC""").fetchnumpy()
            
            # QUERY
            deck_orientation_angle = self.db.con.sql(f"SELECT azimuth FROM proc_{self.damage.axis.table_name} WHERE rdeck = '{deckUid}'").fetchone()[0]
            deck_length = self.db.con.sql(f"SELECT deck_length FROM proc_{self.damage.deck.table_name} WHERE uid = {deckUid}").fetchone()[0]
            
            
            dataStore = dict()
            
            for indx, uid in enumerate(sectors['uid']):
                asc_ts = self._sector_mean_ts(uid, 'ascending' , timeOverlapInfo['ascending']['name'],  self.damage.ascending.scaling_factor)
                dsc_ts = self._sector_mean_ts(uid, 'descending', timeOverlapInfo['descending']['name'], self.damage.descending.scaling_factor)

                if (sectors['sector_tag'][indx] == 'N' or sectors['sector_tag'][indx] == 'S') and (asc_ts[0] is None or dsc_ts[0] is None):
                    break
                
                average_ts = ew_solver.average_ts( ascending_ts  = asc_ts, descending_ts = dsc_ts)
                long, vert = ew_solver.los_long_vert_displacement( average_ts['ascending'], average_ts['descending'], deck_orientation_angle)
                dataStore[sectors['sector_tag'][indx]] = dict(uid = uid, long = long[-1], vert = vert[-1])
            
            
            tilt = ew_solver.get_tilt(dataStore, deck_length)    
            deflection = ew_solver.get_deflection(dataStore, sectors['ndist'], deck_length)
            
            self.db.con.execute(f"""INSERT INTO result (rdeck, orient, tilt, defl) VALUES (?, ?, ?, ?)""", (deckUid, "EW", tilt * self.damage.ascending.scaling_factor, deflection))
        
        print(f"EW solver completed in {time.time() - st1:.2f} seconds.")
            
        
    #TODO: JOIN THE result and proc_{self.damage.deck.table_name} tables to write out the results

    def _sector_mean_ts(self, uid:str, orbit:str, name_fields:list[str], scaling_factor:float = 1.0):
        """
        Calculate the mean time series for a given sector and orbit.
        

        Args:
            uid (str): The unique identifier for the sector.
            orbit (str): The orbital orientation ('ascending' or 'descending').
            name_fields (list[str]): The list of field names to calculate the mean for.
        Returns:
            tuple: A tuple containing the unique identifier and the mean values for the specified fields.
        """
        selectStatement = ",".join([f"MEAN({i} - {name_fields[0]})*{scaling_factor} AS {i}" for i in name_fields])
        pointUIDs =f"SELECT uid FROM proc_{orbit} WHERE rsector == '{uid}'"

        return self.db.con.sql(f"SELECT {selectStatement} FROM {orbit} WHERE uid IN ({pointUIDs})").fetchall()[0]
    
    def _get_timeoverlap(self) -> dict:
        """
        Get the time overlap information between ascending and descending data.

        Returns:
            dict: A dictionary containing the start and end dates for both ascending and descending data.
        """
        # ascName, ascDate = extract_dates(get_column_names(self.damage.ascending.table_name, self.db.con))
        ascName, ascDate = self._extract_dates(self.dbpipeline.get_attributes(self.damage.ascending.table_name))
        # dscName, dscDate = extract_dates(get_column_names(self.damage.descending.table_name, self.db.con))
        dscName, dscDate = self._extract_dates(self.dbpipeline.get_attributes(self.damage.descending.table_name))

        lmin = max(ascDate.min(), dscDate.min())
        lmax = min(ascDate.max(), dscDate.max())

        lmax_bounder = "ascending" if (ascDate == lmax).any() else "descending"
        
        rmax = ascDate[ascDate >= lmax][0] if lmax_bounder == "descending" else dscDate[dscDate >= lmax][0] 
        return dict(
            rmin = lmin, rmax = rmax,
            ascending = dict(name = ascName, date = ascDate),
            descending = dict(name = dscName, date = dscDate),
        )
        
    def _ns_solver_data(self, deckUid:int, timeOverlapInfo:dict) -> dict:
        """
        Prepare the data for the NS solver.
        This method retrieves the necessary data for the NS solver based on the provided deck UID and time overlap information.
        Arguments
        ----------
            deckUid (int): The unique identifier for the deck.
            timeOverlapInfo (dict): A dictionary containing the time overlap information for ascending and descending data.
        Returns
        -------
            dict: A dictionary containing the deck, ascending, and descending data for the NS solver.
        """
        deck = self.db.con.sql(f"""
                SELECT span_count, deck_length
                FROM proc_{self.damage.deck.table_name}
                WHERE uid = {deckUid}
            """).fetchnumpy()

        asc = self.db.con.sql(f"""
            SELECT 
                first.ndist_axis as ndist,
                second.{timeOverlapInfo['ascending']['name'][-1]} - second.{timeOverlapInfo['ascending']['name'][0]}  as disp,
            FROM (SELECT uid, ndist_axis FROM proc_{self.damage.ascending.table_name} WHERE rdeck = {deckUid}) as first
            JOIN {self.damage.ascending.table_name} as second
            ON first.uid = second.uid
            ORDER BY first.ndist_axis ASC
        """).fetchnumpy()

        dsc = self.db.con.sql(f"""
            SELECT 
                first.ndist_axis as ndist,
                second.{timeOverlapInfo['descending']['name'][-1]} - second.{timeOverlapInfo['descending']['name'][0]}  as disp,
            FROM (SELECT uid, ndist_axis FROM proc_{self.damage.descending.table_name} WHERE rdeck = {deckUid}) as first
            JOIN {self.damage.descending.table_name} as second
            ON first.uid = second.uid
            ORDER BY first.ndist_axis ASC
        """).fetchnumpy()

        return dict(deck = deck, ascending = asc, descending = dsc)

    def generate_report(self, based_on:str=None):
        """ Generate a PDF report of the damage assessment results.
        This method generates a PDF report containing the damage assessment results for each deck in the database.
        The report will include plots for each deck based on the specified column name.
        
        Arguments
        ----------
        based_on (str): The column name to base the report on. This should be a valid column in the damage deck table.
        
        Raises
        ------
            ValueError: If the specified column does not exist in the damage deck table.
        """
    
        # valid_columns = get_column_names(self.damage.deck.table_name, self.db.con)
        valid_columns = self.dbpipeline.get_attributes(self.damage.deck.table_name)
        if based_on not in valid_columns and based_on is not None:
            raise ValueError(f"The {based_on} column does not exist in the {self.damage.deck.table_name} table.")
        
        with PdfPages(self.db._db_path.split('.')[0] + '_report.pdf') as pdf:
            # TODO: ITERATE OVER THE DECK UIDS FOR BOTH NS AND EW 
            # and generate the plots for each deck
            for deckuid in self.db.con.sql(f"SELECT rdeck FROM result").fetchnumpy()['rdeck']:
                fig, ax = self._plot(deckuid, self._buf_size/2)
                pdf.savefig(fig)
                plt.close(fig)
                

    def export_results(self, bridge_object:Union[Deck, Axis, Support, Ascending, Descending],
                      output_file:str):
        # TODO:
        pass
        
    def _extract_dates(self, column_names: list[str]) -> tuple[ndarray, ndarray]:
        """ Extracts date fields from a list of column names and returns name fields and date fields as numpy arrays.

        Arguments:
        ----------
            column_names (list[str]): A list of column names to search for date patterns.
        
        Returns:
        -------
            tuple[ndarray, ndarray]: A tuple containing two numpy arrays:
                - nameFields: Array of column names that contain date patterns.
                - dateFields: Array of dates extracted from the column names.

        """
        nameFields = []    
        dateFields = []
        for date in column_names:
            dateSearch = re.search(r'\d{8}', date)
            if dateSearch:
                pdate = dsparser.parse(dateSearch.group(), fuzzy=True).date()
                dateFields.append(pdate)
                nameFields.append(date)
        return array(nameFields), array(dateFields)

    def _plot(self, deckuid: int, buf_dist:float):
        """ Plot the damage assessment results for a specific deck UID.
        This method retrieves the necessary data for the specified deck UID and plots the damage assessment results using the Plotter class.
        
        Arguments
        ----------
        deckuid (int): The unique identifier for the deck.
        buf_dist (float): The buffer distance used for processing geometries.
        
        Returns
        -------
        fig: matplotlib.figure.Figure
            The figure object containing the plotted damage assessment results.
        """
        timeoverlapInfo = self._get_timeoverlap()

        deck = wkbloads(self.db.con.sql(self.query.deck_geometry(deckuid, f"proc_{self.damage.deck.table_name}")).fetchone()[0])
        axis = wkbloads(self.db.con.sql(self.query.axis_geometry(deckuid, f"proc_{self.damage.axis.table_name}")).fetchone()[0])
        support = wkbloads(self.db.con.sql(self.query.support_geometry(deckuid, f"proc_{self.damage.support.table_name}")).fetchall())
        sectors = wkbloads(self.db.con.sql(self.query.sector_geometry(deckuid)).fetchall())
        deck_edge = wkbloads(self.db.con.sql(self.query.deck_edge(deckuid, f"{self.damage.deck.table_name}")).fetchall()[0])
        asc_points = self.db.con.sql(self.query.scatter_geometry(deckuid, f'proc_{self.damage.ascending.table_name}')).fetchnumpy()
        asc_proj = self.db.con.sql(self.query.projected_scatters(deckuid, f'proc_{self.damage.ascending.table_name}')).fetchnumpy()
        dsc_points = self.db.con.sql(self.query.scatter_geometry(deckuid, f'proc_{self.damage.descending.table_name}')).fetchnumpy()
        dsc_proj = self.db.con.sql(self.query.projected_scatters(deckuid, f'proc_{self.damage.descending.table_name}')).fetchnumpy()
        buf_edges = self.db.con.sql(self.query.buffer_edge(deckuid, f"proc_{self.damage.axis.table_name}", f"proc_{self.damage.deck.table_name}")).fetchone()
        deck_edge_graph = self.db.con.sql(self.query.deck_edge_graph(deckuid, f"proc_{self.damage.axis.table_name}", f"proc_{self.damage.deck.table_name}")).fetchone()
        asc_geom_graph = self.db.con.sql(self.query.scatter_graph(deckuid, self.damage.ascending.table_name, timeoverlapInfo['ascending']['name'])).fetchnumpy()
        dsc_geom_graph = self.db.con.sql(self.query.scatter_graph(deckuid, self.damage.descending.table_name, timeoverlapInfo['descending']['name'])).fetchnumpy()
        support_graph = self.db.con.sql(self.query.support_graph(deckuid, self.damage.axis.table_name, self.damage.support.table_name)).fetchnumpy()
        deck_orientation = self.db.con.sql(f"SELECT orientation FROM proc_{self.damage.deck.table_name} WHERE uid = {deckuid}").fetchone()[0]
        ascending_quad_solution = self.db.con.sql(f"""SELECT ns_quadratic_asc_x as x, ns_quadratic_asc_y as y FROM result WHERE rdeck = {deckuid}""").fetchnumpy()
        descending_quad_solution = self.db.con.sql(f"""SELECT ns_quadratic_dsc_x as x, ns_quadratic_dsc_y as y FROM result WHERE rdeck = {deckuid}""").fetchnumpy()
        ascending_analytic_solution = self.db.con.sql(f"SELECT ns_analytical_asc_y FROM result WHERE rdeck = {deckuid}").fetchall()[0][0]
        descending_analytic_solution = self.db.con.sql(f"SELECT ns_analytical_dsc_y FROM result WHERE rdeck = {deckuid}").fetchall()[0][0]
        ascending_tilt_deflection = self.db.con.sql(f"SELECT tilt_asc, defl_asc FROM result WHERE rdeck = {deckuid}").fetchone()
        descending_tilt_deflection = self.db.con.sql(f"SELECT tilt_dsc, defl_dsc FROM result WHERE rdeck = {deckuid}").fetchone()
        ew_tilt_deflection = self.db.con.sql(f"SELECT tilt, defl FROM result WHERE rdeck = {deckuid}").fetchone()
        
        self._plotter.plot(
            deck_geom = deck,
            sector_geom = sectors,
            axis_geom = axis,
            support_geom = support,
            deck_edges = deck_edge,
            ascending_geom = asc_points,
            descending_geom = dsc_points,
            buf_dist = buf_dist,
            projected_ascending=asc_proj,
            projected_descending=dsc_proj,
            buffer_edges = buf_edges,
            deck_edge_graph = deck_edge_graph,
            ascending_geom_graph = asc_geom_graph,
            descending_geom_graph = dsc_geom_graph,
            support_graph = support_graph,
            deck_orientation = deck_orientation,
            ascending_quad_solution=ascending_quad_solution,
            descending_quad_solution=descending_quad_solution,
            ascending_analytic_solution=ascending_analytic_solution,
            descending_analytic_solution=descending_analytic_solution,
            ascending_tilt_deflection = ascending_tilt_deflection,
            descending_tilt_deflection = descending_tilt_deflection,
            ew_tilt_deflection = ew_tilt_deflection
        )
        
        
        self._plotter.postprocess(name_tag=deckuid)

        return self._plotter.get_figure()