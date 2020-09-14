CrossSection - a python tool for creating ready to use cross-section data from high resolution elevation data for hydrologic modelling 
Florian Krebs (1,2), Philipp Kraft (1), Martin Bach (1),  Lutz Breuer (2)

(1) Research Centre for Biosystems, Land Use and Nutrition, Chair of Landscape, Water and Biogeochemical Cycles, Justus Liebig University Giessen, Heinrich-Buff-Ring 26, 35392 Giessen, Germany <br>
(2) knoell Germany GmbH, Konrad-Zuse-Ring 25, 68163 Mannheim, Germany

Correspondence to: Florian Krebs (fkrebs@knoell.com)

# Website (under construction)
Webiste: https://flokrebs.github.io/CrossSections/

# install with PIP
pip install crossections

# import tool
import crossections as cs

# run tools....

# example yaml file

Project: Example

Stream_shapefile: C:\Geodata\streams.shp
DTM: C:\Geodata\dtm.tif
Structures_shapefile: c:\Geodata\roads_lines.shp

Buffer_value: 5
Aggregation_field: MOD_ID
strahler_field: strahler

Line distance: 20
Cross-section width: 30

Depth extrapolation: True
strahler_depth: { 1: 0.1, 2: 0.2, 3: 0.3, 4: 0.4, 5: 0.5, 6: 0.6, 7: 0.7, 8: 0.8}
default_depth: .3
Point spacing: .5

inner_limit: 15
outer_limit: 25
Angle_sum: 90
Angle_max: 50
Angle_negative: 

Negative sum allowed: 1
Bottom Angle limit: 
Bottom height lim: 0.4
Filter bs max: 10
Filter depth min: .5
Depth method: lowest

Create_lwst_pt_shape: True
Plotting: True 
