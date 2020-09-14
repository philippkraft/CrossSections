# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 15:15:12 2019
CrossSection tool v01
@author: Florian Krebs
"""
from .Analysis import Analysis as Analysis
from .Geodata_package import Geodata_package as Geodata_package
from .Geodata_results import Geodata_results as Geodata_results
from .Depth import Depth as Depth
from .Plotting import Plotting as Plotting
from .Results import Results as Results
from .Parameters import Parameters as Parameters
    
def run(yaml_file):
    '''
    runs the complete workflow
    '''
    # imports all modules
    parameters = Parameters(yaml_file)
    # runs the Geodata package processing tool to create cross sections    
    Geodata_package(parameters)
    # Analyses the cross sections with the defined rules in the parameter file
    Analysis(parameters, parameters.CrossSections_csv,parameters.Xresults)
    # Extracts the lowest point from the point shapes
    Results(parameters, parameters.Xresults,parameters.Xresults_filter, '')
    # Extracts the lowest point from the point shapes
    Geodata_results(parameters)
    # runs the Depth extrapolation module
    Depth(parameters)
    # Plots the profiles
    Plotting(parameters)
    print('finished ...')