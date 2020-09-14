# -*- coding: utf-8 -*-
'''
Copyright (c) 2020 by Florian Krebs
This file is part of CrossSections Tool.

:author: Florian Krebs

:paper: Krebs et. al (2020) - CrossSection - a python tool for creating ready
 to use cross-section data from high resolution elevation data for hydrologic
 modelling. In prep.

This package enables the creation and analysis of stream cross sections.

:dependencies: - Numpy >1.8
               - Pandas >0.13
               - GeoPandas 1.5
               - Matplotlib >1.4
               - yaml 5.1.2
               - rasterio 1.0.21
               - shapely 1.6.4.post1
               
Please cite our paper, if you are using CrossSections.

If you want to participate:
Feel free to fork your repository and submit changes!

https://github.com/FloKrebs/CrossSections
'''

from . import Parameters  # Parameter class
from . import Geodata_package  # creates the cross-sections
from . import Analysis  # analyses the cross-sections
from . import Geodata_results  # creates resulting geodata
from . import Depth  # applies Depth extrapolation
from . import Plotting  # plots profiles
from . import Results  # filters the results
from . import main  # run runs with the complete workflow

print('CrossSections v1.1.0')
__version__ = '1.1.0'