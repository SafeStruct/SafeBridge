# imporitng the libs
from safebridge.damage_assessment import DamageAssessment
from safebridge.data import Deck, Axis, Support, Ascending, Descending

# creating the damage assessment object
damage_assessment = DamageAssessment(
    deck = Deck(source_file = "./toy_data/deck.shp"), 
    axis = Axis( source_file = "./toy_data/axis.shp" ), 
    support = Support(source_file = "./toy_data/support.shp"), 
    ascending = Ascending(
        source_file = "./toy_data/ascending_data.csv",
        unit = "mm",
        lat_field = "lat",
        lon_field = "lon",
        orbit_azimuth = 348.66,
        incidence_angle = 31.1,
    ),
    descending=Descending(
        source_file = "./toy_data/descending_data.csv",
        unit = "mm",
        lat_field = "lat",
        lon_field = "lon",
        orbit_azimuth = 190.72,
        incidence_angle = 35.4,
    )
)

# loading source files
damage_assessment.load_source_files()

# preprocessing the provided data
damage_assessment.preprocess(computational_projection="EPSG:28992", buffer_distance = 6)

# running the damage assessment
damage_assessment.assess_damage()

# generating the report
damage_assessment.generate_report()
