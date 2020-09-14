# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 14:16:48 2020

@author: Florian Krebs
"""
import pandas as pd
import numpy as np
import math
from .Analysis import Analysis as Analysis
from .Results import Results as Results


class Depth:
    def __init__(self, parameters):
        # get Parameters from Parameters class
        if type(parameters) != str:
            self.parameters = parameters
        elif type(parameters) == str:
            from Parameters import Parameters
            self.parameters = Parameters(parameters)
        self.cd_method = None
        self.__call__()

    def __call__(self):
        if self.parameters.depth_extr is True:
            print('Depth extrapolation')
            self.df = pd.read_csv(self.parameters.Xresults_filter)
            self.height_cols = [col for col in self.df.columns if 'X' in col]
            # checks if artifical depth should be added
            self.check_calculate_depth()
            # checks method
            if self.cd_method is not None:
                print('Method:', self.cd_method)
                res_depth_ext = self.iterate_list()
                res_depth_ext = Analysis(self.parameters, res_depth_ext,
                                   self.parameters.Xresults_dcorr)
                Results(self.parameters, self.parameters.Xresults_dcorr,
                        self.parameters.Xresults_dcorr_filter, 'dc')
            else:
                print('skipped.. no method defined')
            print('...done')

    def iterate_list(self):
        '''
        Iterrates the profiles to add Depth extrapolation. Uses rows of the
        results. Creates dataframe with height points for each ID wich are 
        replaced in create_new_cs function.
        @type self.corr_heights: pd.DataFrame
        @param self.corr_heights: Holds ID and Point(X..) Z correction values
        @type res_depth_ext: pd.DataFrame
        @param res_depth_ext: DataFrame with extrapolated cross sections 
        '''
        new = []
        h_corr = []
        for index, row in self.df.iterrows():
            row, heights = self.add_calculated_depth(row)
            new.append(row)
            if heights is not False:
                h_corr.append(heights)
        # concats single DataFrames       
        results = pd.DataFrame(new)
        # exports the Dataframe with Depth interpolation attributes
        results.to_csv(self.parameters.Xresults_filter, index=False)
        self.corr_heights = pd.concat(h_corr)
        res_depth_ext = self.create_new_cs(results)
        return res_depth_ext

    def create_new_cs(self, results):
        '''
        applies the new calculated depths to the cross section data. Selects
        the cross section points by ID from the corr_heights dataframe and
        replaces the row of the height point (X..) with the corresponding
        correction height.
        @type self.corr_heights: pd.DataFrame
        @param self.corr_heights: Holds ID and Point(X..) Z correction values
        @type res_depth_ext: pd.DataFrame
        @param res_depth_ext: DataFrame with extrapolated cross sections
        '''
        new = []
        for index, row in results.iterrows():
            # check if corrected height exist
            if row.ID in self.corr_heights.ID.unique():
                corr = self.corr_heights[self.corr_heights.ID == row.ID]
                for i, r in corr.iterrows():
                    row[r.point] = r.Z
            new.append(row)
        res_depth_ext = pd.DataFrame(new)
        return res_depth_ext

    def check_calculate_depth(self):
        '''
        Checks if the parameters default_depth or strahler_depth are defined.
        If defined a depth extension is calculated.
        @type cd_method: <str> or <boolean>
        @param cd_method: checks the type of depth calculation method
        '''      
        if 'art_depth' in self.df.columns:
            self.cd_method = 'art_depth'
        elif self.parameters.strahler_depth is not None:           
            self.strahler_dict = self.parameters.strahler_depth
            self.cd_method = 'strahler'
        elif self.parameters.default_depth is not None:
            self.cd_method = 'default'
        else:
            self.cd_method = 'infinite'
            
    def add_calculated_depth(self, row):
        '''
        Adds a calculated cross section depth extension. Especially for small
        profiles the definition of a fixed depth can lead to unplausile results
        (crossing of lines). Hence the construction is first checked by a
        triangle calculation (following common mathematical rules). By the
        construction of the triangle the height from side c (bottom width) to
        point C is derived. This is done by calculating the angles alpha
        (left side), beta (right side) togehter with c (Bottom width).
        The height (h) from c to C is calculated. If the height
        is lower than the parametrized depth value (strahler, default), the
        point C is defined as the lower depth extension limit.
        @type self.Xpts: <int>
        @param self.Xpts: Number ticks
        '''
        # get alpha of both sides
        # left
        a_l = row.depth
        b_l = row.width_l-row.bottom_width_l
        gamma = math.radians(90)
        alpha_l = math.atan2(a_l, b_l)
        row['alpha_l'] = np.degrees(alpha_l)
        # right
        a_r = row.depth
        b_r = row.width_r-row.bottom_width_r
        alpha_r = math.atan2(a_r, b_r)
        row['alpha_r'] = np.degrees(alpha_r)
        # construct new triangle down the bottomline
        c = row.bottom_width
        # alpha is calculated by alpha_l and alpha_r
        alpha = alpha_r
        beta = alpha_l
        gamma = math.radians(180)-alpha-beta
        row['alpha_c'] = np.degrees(alpha)
        row['beta_c'] = np.degrees(beta)
        row['gamma_c'] = np.degrees(gamma)
        # calculate b (law of sines)
        b = c*np.sin(beta) / np.sin(gamma)
        # calculate a (law of cosines)
        a = np.sqrt((b**2-2*b*c*np.cos(alpha)+c**2))
        # calculate the height h from c to C
        hc = a*np.sin(beta)
        row['hc'] = hc
        # gets the area of the triangle
        row['Area'] = (c*hc)/2
        # construct helper triangle to get c
        xc = self.c_pos_bottom(hc, beta)
        row, heights = self.check_hc(row, alpha, beta, c, xc, hc)
        return row, heights

    def check_hc(self, row, alpha, beta, c, xc, hc):
        '''
        Checks if the h value is lower then the maximum depth value.
        If the value is higher an additional bottom level can be added,
        as the bottom dimensions are big enough.
        Otherwise this would lead to unplausible results.
        '''
        if self.cd_method == 'default':
            depth = self.parameters.default_depth
        elif self.cd_method == 'strahler':
            try:
                depth = float(self.strahler_dict[
                        row[self.parameters.strahler_field]])
            except:
                print('Check Strahler ID %s value in parameter strahler_depth'
                      % row[self.parameters.strahler_field])
        # art_depth is always dominant
        # needs to be here because every row migh have a different depth value
        elif self.cd_method == 'art_depth':
            depth = row['art_depth']
        elif self.cd_method == 'infinite':
            depth = 999
        if hc == 0:
            row['bw_l_ext'] = np.NaN
            row['bw_r_ext'] = np.NaN
            row['depth_ext'] = np.NaN
            row['depth_corr'] = np.NaN
            heights = False
        elif np.isnan(hc) == True or np.isnan(xc) == True:
            row['bw_l_ext'] = np.NaN
            row['bw_r_ext'] = np.NaN
            row['depth_ext'] = np.NaN
            row['depth_corr'] = np.NaN
            heights = False
        elif hc <= depth:
            row['depth_corr'] = 1
            row['bw_l_ext'] = xc
            row['bw_r_ext'] = -99
            row['depth_ext'] = hc
            heights = self.correct_heights(alpha, beta, c, xc, hc, row, depth)
        else:
            left_c = self.c_pos_bottom(depth, beta)
            right_c = self.c_pos_bottom(depth, alpha)
            row['bw_l_ext'] = left_c
            row['bw_r_ext'] = right_c
            row['depth_ext'] = depth
            row['depth_corr'] = 1
            heights = self.correct_heights(alpha, beta, left_c,
                                           xc, hc, row, depth)
        return row, heights

    def correct_heights(self, alpha, beta, c, xc, hc, row, depth_limit):
        '''
        corrects the height by tick frequency.
        the triangle is divided in left / right side of the lowest point
        as the tool iterates from left to right, the position value needs
        to get back to zero. E.g bw = 1.5, left 1, right 0.5.
        freq count needs to count back from 1 to 0.5
        '''
        Z_low = row.Z_low
        freq = self.parameters.ticks
        # gets the first point of the bottom line
        pt = row.lowest_PT-(row.bottom_width_l/freq)
        # checks if alpha or beta should be used
        freq_count = 0
        # needed for right side
        opp_count = 0
        new = []
        # iterates until c is filled up with ticks
        if row.bottom_width == freq:
            # this part is used if the bottom width is equal to the
            # tick frequency - it lowers the two lowest cs points
            depth = row['depth_ext']
            new.append((row.ID, 0, "X%s" % int(pt), depth, Z_low-depth))
            new.append((row.ID, freq, "X%s" % int(pt+1), depth, Z_low-depth))
        else:
            while freq_count < row.bottom_width:
                # checks side, if alpha or beta is used
                if freq_count <= xc:
                    # gets the depth for the c point on bottom
                    depth = self.z_pos(beta, freq_count)
                    if depth > depth_limit:
                        depth = depth_limit
                    # raises the counter
                    freq_count = freq_count + freq                   
                elif freq_count > xc:
                    if opp_count == 0:
                        opp_count = row.bottom_width-freq_count+freq
                    # gets the depth for the c point on bottom  
                    depth = self.z_pos(alpha, opp_count)
                    if depth > depth_limit:
                        depth = depth_limit
                    # raises the counter
                    freq_count = freq_count + freq
                    # gets next position of the triangle
                    opp_count = opp_count - freq
                new.append((row.ID, freq_count, "X%s" % int(pt),
                            depth, Z_low-depth))
                # raises to next point
                pt = pt+1
        # writes the results to a DataFrame
        res = pd.DataFrame(new, columns=['ID', 'position', 'point',
                                         'depth', 'Z'])
        del new, freq, freq_count, opp_count, pt
        return res

    def z_pos(self, beta, c):
        '''
        constructs a triangle with beta, c = tick, and gamma = 90 to get
        length of b
        '''
        # set alpha to 90
        gamma = math.radians(90)
        b = c*np.sin(beta) / np.sin(gamma)
        return b

    def c_pos_bottom(self, depth, beta):
        '''
        constructs a triangle with alpha = 90, beta, and depth to get
        length of c
        '''
        # set alpha to 90
        alpha = math.radians(90)
        b = depth
        a = b*np.sin(alpha) / np.sin(beta)
        gamma = math.radians(180)-alpha-beta
        # get position of C (x axis)
        xc = np.sqrt((a**2-2*a*b*np.cos(gamma)+b**2))
        return xc