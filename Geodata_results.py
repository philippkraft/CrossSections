# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 09:48:16 2019
Filters the cross-section points for the lowest point
@author: Florian Krebs
"""

import geopandas as gpd
import pandas as pd
import csv


class Geodata_results:
    def __init__(self, parameters):
        # gets parameters from parameter class
        if type(parameters) != str:
            self.parameters = parameters
        elif type(parameters) == str:
            from Parameters import Parameters
            self.parameters = Parameters(parameters)
        self.main()

    def main(self):
        '''
        reads the filtered results and creates point geometries
        @type filter: Pandas DataFrame
        @param filter: DataFrame with filtered anaylsis results
        @type filter_IDs: <list>
        @param filter_IDs: list of profiles with regular and filtered
        cross section analysis results
        @type l_pts_dict: <dict>
        @param l_pts_dict: dictionary containing the information which
                            cs point is the lowest in the profile
        @type res: Geopandas GeoDataFrame
        @param res: Lowest points with geometries and attributes
        '''
        # reads filter results
        self.filter = pd.read_csv(self.parameters.Xresults_filter)
        # gets filter IDs
        self.filter_IDs = self.filter.ID.unique()
        # reads the lowest points for cs as dict
        l_pts_dict = self.get_dict()
        # extract lowest points by dict
        res = self.extract_lowest_point(l_pts_dict)
        # checks if lowest point shapefile should be created
        if self.parameters.cr8_lpt_shp is not None:
            res.to_file(self.parameters.lowest_points, driver="GPKG")
        # filters the geometries by filter ID
        res = res[res.ID.isin(self.filter_IDs)]
        # merges attributes to result file
        self.res = self.merge_analysis(res)
        # exports results, moves geometry col to end
        cols_at_end = ['geometry']
        self.res = self.res[[c for c in self.res if c not in cols_at_end] +
                            [c for c in cols_at_end if c in self.res]]
        self.res.to_file(self.parameters.Xresults_f_shp, driver="GPKG")

    def get_dict(self):
        '''
        reads the dictionary with the infromation on lowest poit IDs
        @type l_pts_dict: <dict>
        @param l_pts_dict: dictionary containing the information which
                            cs point is the lowest in the profile
        '''
        with open(self.parameters.l_pts_dict) as csv_file:
            reader = csv.reader(csv_file)
            l_pts_dict = dict(reader)
        return l_pts_dict

    def extract_lowest_point(self, l_pts_dict):
        '''
        extracts the lowest point from the cs profile points
        @type l_pts_dict: <dict>
        @param l_pts_dict: dictionary containing the information which
                            cs point is the lowest in the profile
        @type res: Geopandas GeoDataFrame
        @param res: Lowest points with geometries
        '''
        print('Extract lowest point')
        cs_pnts = gpd.read_file(self.parameters.XPoints)
        self.cs_pnts = cs_pnts
        crs = cs_pnts.crs
        new = []
        for index, row in cs_pnts.iterrows():
            try:
                PT = int(float(l_pts_dict[str(row.ID)]))
                if row.X == PT:
                    new.append(row)
            except:
                print('error in %s' % row.ID)
        df = pd.DataFrame(new)
        res = gpd.GeoDataFrame(df, crs=crs, geometry=df.geometry)
        return res

    def merge_analysis(self, fil):
        '''
        merges the analysis data to the point geometries
        @type res: Geopandas GeoDataFrame
        @param res: Lowest points with geometries and attributes
        '''
        print('merge with analysis data')
        new = []
        geo = []
        for index, row in fil.iterrows():            
            df =  self.filter[self.filter.ID == row.ID]
            new.append([r for i, r in df.iterrows()][0])
            geo.append(row.geometry)
        res = gpd.GeoDataFrame(new, geometry=geo, crs=fil.crs)
        return res