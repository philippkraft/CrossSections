# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 09:14:58 2019
parameters class to import parameters from .yaml file

@author: Florian Krebs


Example yaml file
----------------------------------------------------------
Project: Cross_sections_Example

Stream_shapefile: C:\stream.shp
DTM: c:\dtm.tif   
Structures_shapefile: c:\structures.shp

Buffer_value: 20
Aggregation_field: key
strahler_field: str

Line distance: 20
Cross-section width: 30
strahler_depth: {1: 0.1, 2: 0.2, 3: 0.3, 4: 0.4, 5: 0.5, 6: 0.6, 7: 0.7, 8: 0.8}
default_depth: 0.1
Point spacing: 1

inner_limit: 15
outer_limit: 30
Angle_sum: 90
Angle_max: 50
Angle_negative: 

Negative sum allowed: 1
Bottom Angle limit: 
Bottom height lim: 0.15
Filter bs max: 10
Filter depth min: .5
Depth method: lowest

Create_lwst_pt_shape: True
Plotting: True
--------------------------------------------------------------
"""
import yaml
from os import path

class Parameters:
    def __init__(self, config_file):
        self.Project = None
        self.project_path = None
        self.results_path = None
        self.yaml_file = config_file
        self.Stream_shapefile = None
        self.dtm = None
        self.structures_shp = None
        self.Buffer_value = None
        self.Aggregation_field = None
        self.strahler_field = None
        self.distance = None
        self.width = None
        self.ticks = None
        self.inner_limit = None
        self.outer_limit = None
        self.Angle_sum = None
        self.Angle_max = None
        self.angle_bot_lim  = None
        self.bottom_height_lim = None
        self.Angle_negative = None
        self.cr8_lpt_shp = None
        self.plotting = None
        self.depth_method = None
        self.n_sum_allowed = None
        self.XLines = None
        self.lowest_points = None
        self.XPoints = None
        self.Xresults = None
        self.CrossSections = None    
        self.Xresults = None
        self.Xresults_filter = None
        self.Xresults_dcorr = None
        self.Xresults_f_shp = None
        self.filter_bs_max = None
        self.filter_depth_min = None
        self.strahler_depth = None
        self.default_depth = None
        self.plotting_path_dc = None
        self.Intersect_pts = None
        self.depth_extr = False
        self.__set_parameters__()

    def __set_parameters__(self):
        '''
        Sets the parameters from .yaml file
        @type Project: <str>
        @param Project: Name of the Project
        @type project_path: <str>
        @param project_path: path to the project
        @type results_path: <str>
        @param results_path: result path
        @type Stream_shapefile: <str>
        @param Stream_shapefile: path to Stream geometries used as input
        (all file types compatible with geopandas read function
        https://geopandas.org/io.html)
        @type dtm: <str>
        @param dtm: path to dtm file (all file types compatible with 
        rasterio read function, https://rasterio.readthedocs.io)
        @type structures_shp: <str>
        @param structures_shp: path to structure geometries used as input
        (all file types compatible with geopandas read function
        https://geopandas.org/io.html)
        @type Buffer_value: <float>
        @param Buffer_value: value used for buffering of structures (lines)
        @type Aggregation_field: <str>
        @param Aggregation_field: field name used for aggregation
        @type strahler_field: <str>
        @param strahler_field: field name with strahler code
        @type distance: <float>
        @param distance: distance tick spacing of cross section creation
        @type width: <float>
        @param width: width of cross sections
        @type ticks: <float>
        @param ticks: distance of points along cross section.
        @type inner_limit: <int>
        @param inner_limit: defines search radius from the middle point to be
        used for the identification of the lowest point
        @type outer_limit: <int>
        @param outer_limit: defines the search limit of the profile extend 
        around the lowest point
        @type Angle_sum: <float>
        @param Angle_sum: limit of aggregated angles as construction condition
        (if sum is reached, the last point is considered cs limit)
        @type Angle_max: <float>
        @param Angle_max: limit of angle size as construction condition
        (if angle is reached, the last point is considered cs limit)
        @type angle_bot_lim: <float>
        @param angle_bot_lim: limit of angle size as construction condition for
        the bottom width. (if angle is reached, the last point is considered
        bw limit)
        @type bottom_height_lim: <float>
        @param bottom_height_lim: height limit as construction condition for
        the bottom width. (if height difference is reached, the last point is 
        considered bw limit)
        @type Angle_negative: <float>
        @param Angle_negative: limit of negative angle as construction condition 
        for the cs. (if negative angle is reached, the last point is considered
        cs limit)
        @type cr8_lpt_shp: <boolean>
        @param cr8_lpt_shp: creates point geometries with the lowest point
        @type plotting: <boolean>
        @param plotting: Enables/disables plotting function
        @type depth_method: <boolean>
        @param depth_method: Defines depth calculation method. 
        Options: "lowest" or "mean". lowest: lowest Z value is used. mean:
            mean Z value of cs points is used.
        @type n_sum_allowed: <int>
        @param n_sum_allowed: Defines how often the angle sum can be negative.
        This can be useful if a profile is noisy with slopy characters.
        Default: 1
        @type XLines: <str>
        @param XLines: Path to result file (CrossLines.gpkg)
        @type lowest_points: <str>
        @param lowest_points: Path to lowest point file ('lwst_pnt.gpkg')
        @type XPoints: <str>
        @param XPoints: cross section points file ('Points.gpkg')
        @type Xresults: <str>
        @param Xresults: cross section result file ('results.csv')
        @type Xresults_filter: <str>
        @param Xresults_filters: cross section result filter file
                                ('results_filter.csv')
        @type Xresults_f_shp: <str>
        @param Xresults_f_shp: cross section result filter file geometries
                                ('results_filter.gpkg')
        @type strahler_cs_width: <dict>
        @param strahler_cs_width: limits the cs analysis to a respective strahler
        order extent (e.g. {1 : 30, 2 : 40 ...})
        @type filter_bs_max: <float>
        @param filter_bs_max: value to filter the results by max bankslope
        (higher bankslope values are dropped)
        @type filter_depth_min: <float>
        @param filter_depth_min: value to filter the results by min depth
        (lower depth values are dropped)
        @type strahler_depth: <float>
        @param strahler_depth: value to create artificial depth by strahler
        @type default_depth: <float>
        @param default_depth: default value to create artificial depth
        '''
        with open(path.join(self.yaml_file)) as f:
            data = yaml.load(f, Loader=yaml.Loader)
            self.Project = data['Project']
            self.Stream_shapefile = data['Stream_shapefile']
            self.dtm = data['DTM']
            self.structures_shp = data['Structures_shapefile']
            self.Buffer_value = data['Buffer_value']
            self.strahler_field = data['strahler_field']
            self.Aggregation_field = data['Aggregation_field']
            self.distance = data['Line distance']
            self.width = data['Cross-section width']
            self.depth_extr = bool(data['Depth extrapolation'])
            self.strahler_depth = data['strahler_depth']
            self.default_depth = data['default_depth']
            self.ticks = data['Point spacing']
            self.inner_limit = data['inner_limit']
            self.outer_limit = data['outer_limit']
            self.Angle_sum = data['Angle_sum']
            self.Angle_max = data['Angle_max']
            self.Angle_negative = data['Angle_negative']
            self.angle_bot_lim = data['Bottom Angle limit']
            self.cr8_lpt_shp = data['Create_lwst_pt_shape']
            self.plotting = data['Plotting']
            self.depth_method = data['Depth method']
            self.n_sum_allowed = data['Negative sum allowed']
            self.filter_bs_max = data['Filter bs max']
            self.filter_depth_min = data['Filter depth min']
            self.bottom_height_lim = data['Bottom height lim']        
        self.project_path = path.dirname(path.abspath(self.yaml_file))
        self.results_path = self.project_path + '/' + self.Project
        self.XLines = path.join(self.results_path, 'CrossLines.gpkg')
        self.lowest_points = path.join(self.results_path, 'lwst_pnt.gpkg')
        self.XPoints= path.join(self.results_path, 'Points.gpkg')
        self.Xresults = path.join(self.results_path, 'results.csv')
        self.CrossSections_csv = path.join(self.results_path, 'CrossSections.csv')     
        self.Xresults_filter = path.join(self.results_path, 'results_filter.csv')
        self.Xresults_dcorr = path.join(self.results_path, 'results_dext.csv')
        self.Xresults_dcorr_filter = path.join(self.results_path, 'results_dext_filter.csv')
        self.Xresults_f_shp = path.join(self.results_path, 'results_filter.gpkg')
        self.plotting_path = path.join(self.results_path, 'Plot')
        self.l_pts_dict = path.join(self.results_path, 'lst_pts.csv')
        self.Intersect_pts = path.join(self.results_path, 'Intersect_pts.gpkg')