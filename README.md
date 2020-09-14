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

Project: Example<br>

Stream_shapefile: C:\Geodata\streams.shp<br>
DTM: C:\Geodata\dtm.tif<br>
Structures_shapefile: c:\Geodata\roads_lines.shp<br>

Buffer_value: 5<br>
Aggregation_field: MOD_ID<br>
strahler_field: strahler<br>

Line distance: 20<br>
Cross-section width: 30<br>

Depth extrapolation: True<br>
strahler_depth: { 1: 0.1, 2: 0.2, 3: 0.3, 4: 0.4, 5: 0.5, 6: 0.6, 7: 0.7, 8: 0.8}<br>
default_depth: .3<br>
Point spacing: .5<br>

inner_limit: 15<br>
outer_limit: 25<br>
Angle_sum: 90<br>
Angle_max: 50<br>
Angle_negative: <br>

Negative sum allowed: 1<br>
Bottom Angle limit: <br>
Bottom height lim: 0.4<br>
Filter bs max: 10<br>
Filter depth min: .5<br>
Depth method: lowest<br>

Create_lwst_pt_shape: True<br>
Plotting: True <br>
