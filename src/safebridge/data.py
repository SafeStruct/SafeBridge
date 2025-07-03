from dataclasses import dataclass
from typing import Literal

@dataclass
class BaseData:
    """ BaseData is a dataclass that represents the metadata for a data source.

    Attributes
    ----------
    source_file : str
        The path to the source file containing the data.
    table_name : str
        The name of the table associated with the data.
    source_projection : str
        The spatial reference system of the source data. Default is `"EPSG:4326"`.
    """
    source_file : str
    table_name : str
    source_projection : str = "EPSG:4326"

    
@dataclass 
class Deck(BaseData):
    """ Deck class representing the `'deck'` table in the database.
    
    Attributes
    ----------
    source_file : str
        The path to the source file containing the data.
    table_name : str
        The name of the table associated with the data. Default is `"deck"`.
    source_projection : str
        The spatial reference system of the source data. Default is `"EPSG:4326"`.
    """
    table_name : str = "deck"

@dataclass
class Axis(BaseData):
    """ Axis class representing the `'axis'` table in the database.
    
    Attributes
    ----------
    source_file : str
        The path to the source file containing the data.
    table_name : str
        The name of the table associated with the data. Default is `"axis"`.
    source_projection : str
        The spatial reference system of the source data. Default is `"EPSG:4326"`.
    """
    table_name : str = "axis"

@dataclass
class Support(BaseData):
    """ Support class representing the `'support'` table in the database.

    Attributes
    ----------
    source_file : str
        The path to the source file containing the data.
    table_name : str
        The name of the table associated with the data. Default is `"support"`.
    source_projection : str
        The spatial reference system of the source data. Default is `"EPSG:4326"`.
    """
    table_name : str = "support"

@dataclass
class Ascending(BaseData):
    """ Ascending class representing the `'ascending'` table in the database.
    
    Attributes
    ----------
    source_file : str
        The path to the source file containing the data.
    unit : "mm", "cm", "m"
        The unit of measurement for the data, default is `"m"`.
    lat_field : str
        The name of the field containing latitude data.
    lon_field : str
        The name of the field containing longitude data.
    table_name : str
        The name of the table associated with the data. Deafault is `"ascending"`.
    orbit_azimuth : float
        The azimuth angle of the orbit in degrees.
    incidence_angle : float
        The incidence angle of the orbit in degrees.
    source_projection : str
        The spatial reference system of the source data, Default is `"EPSG:4326"`.
    scaling_factor : float
        A scaling factor for the data, Default is `1.0`.
    """
    table_name : str = "ascending"
    unit : Literal["mm", "cm", "m"] = "m"
    lat_field : str = None
    lon_field : str = None
    orbit_azimuth : float = None
    incidence_angle : float = None
    scaling_factor : float = None
    

@dataclass
class Descending(Ascending):
    """ Descending class representing the `'descending'` table in the database.

    Attributes
    ----------
    source_file : str
        The path to the source file containing the data.
    unit : "mm", "cm", "m"
        The unit of measurement for the data, default is `"m"`.
    lat_field : str
        The name of the field containing latitude data.
    lon_field : str
        The name of the field containing longitude data.
    table_name : str
        The name of the table associated with the data. Default is `"descending"`.
    orbit_azimuth : float
        The azimuth angle of the orbit in degrees.
    incidence_angle : float
        The incidence angle of the orbit in degrees.
    source_projection : str
        The spatial reference system of the source data, defaulting to `"EPSG:4326"`.
    scaling_factor : float
        A scaling factor for the data, defaulting to `1.0`.
    """
    table_name : str = "descending"


@dataclass
class BridgeDamage:
    """ BridgeDamage is a dataclass that encapsulates the metadata for bridge damage data.

    Attributes
    ----------
    deck : Deck
        An instance of the Deck class representing the deck data.
    axis : Axis
        An instance of the Axis class representing the axis data.
    support : Support
        An instance of the Support class representing the support data.
    ascending : Ascending
        An instance of the Ascending class representing ascending orbit data.
    descending : Descending
        An instance of the Descending class representing descending orbit data.
    """
    deck: Deck
    axis: Axis
    support: Support
    ascending : Ascending
    descending : Descending
