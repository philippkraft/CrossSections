# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 12:24:30 2019
Class for Analysing the created cross-section (cs) data
@author: Florian Krebs
"""

import pandas as pd
import numpy as np
import csv


class Analysis():
    def __init__(self, parameters, inp, out):
        # get Parameters from Parameters class or path
        if type(parameters) != str:
            self.parameters = parameters
        elif type(parameters) == str:
            from Parameters import Parameters
            self.parameters = Parameters(parameters)
        self.strahler_ids = np.array([])
        self.prints = False
        self.df_input = inp
        self.df_out = out
        self.__call__()

    def __call__(self):
        '''
        calls the main script to create a pandas DataFrame containing the
        cross sections and analsis results.
        Writes the lowest cs point ID to a dictionary
        Finally filters the values if a filter condition is set in the
        parameter class.
        Aggregates the filtered results if a aggregation field is set.
        @type self.df: pd.DataFrame
        @param self.df: List of cs points
        @type l_pts_dict: <dict>
        @param l_pts_dict: dictionary containing the information which
                            cs point is the lowest in the profile
        '''
        if isinstance(self.df_input, pd.DataFrame):
            self.df = self.df_input
        else:
            self.df = pd.read_csv(self.df_input)
        # runs the main analysis module
        self.res = self.main()
        # creates a dict with lowest profile point
        l_pts_dict = self.get_lwst_ID()
        # writes the dict to a csv output file
        self.write_dict(l_pts_dict)
        # writes result file
        self.res.to_csv(self.df_out, index=False)

    def main(self):
        '''
        reads the cs geodata created from the Geodata_package class and
        creates a pandas DataFrame. Checks if a strahler or aggregation field
        is set and adds the field to the result columns, respectively.
        Checks if artificial depth should be added. Selects all height
        columns ('X'), calculates the number of points and starts iterating
        the cs.
        @type self.df: pd.DataFrame
        @param self.df: List of cs points
        @type self.height_cols: <int>
        @param self.height_cols: Number of height columns in df
        @type self.n_points: <int>
        @param self.n_points: Number of height points from self.height_cols
        @type self.middle_point: <int>
        @param self.middle_point: cs middle point (rounded)
        @type self.columns: Pandas dataframe columns
        @param self.columns: Pandas dataframe columns
        '''
        print(self.df_input)
        # sets result columns
        self.columns = ['ID', 'intsct', 'geometry']
        # checks if parameter columns are added
        if self.parameters.strahler_field is not None:
            self.columns.insert(2, self.parameters.strahler_field)
            self.strahler_ids = self.df[
                    self.parameters.strahler_field].unique()
        if self.parameters.Aggregation_field is not None:
                self.columns.insert(4, self.parameters.Aggregation_field)
        if 'art_depth' in self.df.columns:
            self.columns = self.columns+['art_depth']
        # identifies height columns
        self.height_cols = [col for col in self.df.columns if 'X' in col]
        # gets number of height points
        self.n_points = len(self.height_cols)
        # gets middle point
        self.middle_point = round(self.n_points/2)
        # creates list of columns for column filter
        self.columns = self.columns+self.height_cols
        # starts iterrating through DataFrame to analyse the cs
        self.res = self.__iterrate_list()
        return self.res

    def write_dict(self, l_pts_dict):
        '''
        writes a dictionary to csv file
        @type l_pts_dict:	<dict>
        @param l_pts_dict: 	dictionary containing the information which
        cs point is the lowest in the profile
        '''
        with open(self.parameters.l_pts_dict, 'w', newline="") as csv_file:
            writer = csv.writer(csv_file)
            for key, value in l_pts_dict.items():
                writer.writerow([key, value])

    def get_lwst_ID(self):
        '''
        gets the lowest cs ID
        '''
        l_pts_dict = dict(zip(self.res['ID'].values, self.res.lowest_PT))
        return l_pts_dict

    def __iterrate_list(self):
        '''
        Iterrates the cs list and analyses the cs according to the defined
        parameters.
        - Creates a numpy array of the cs
        - Applies a inner and outer limit (search radius) if set in the
        parameter class
        - searches the lowest cs point
        - searches for the highest points in left and right direction from the
        lowest point
        - starts analysing the cs in the direction of the lower heighest point
        - writes the results to a result pandas DataFrame
        @type self.Xpts: <int>
        @param self.Xpts: Number ticks
        @type self.Ypts: <int>
        @param self.Ypts: Number of height points from self.height_cols
        @type self.middle_point: <int>
        @param self.middle_point: cs middle point (rounded)
        @type self.lwst_pnt: np.array
        @param self.lwst_pnt: lowest cs point
        @type self.res: pd.DataFrame
        @param self.lwst_pnt: Dataframe with results
        '''
        print('Analysing %s cross sections...' % str(len(self.df)))
        # new list to join results
        new = []
        # iterrate through list of cross sections
        for index, row in self.df.iterrows():
            # gets x directions
            self.Xpts = [i for i in range(0, self.n_points)]
            # gets y heights
            self.Ypts = [float(row[i]) for i in self.height_cols]
            # creates a numpy array cs
            points = []
            for i in self.Xpts:
                xy = [i, self.Ypts[i]]
                points.append(xy)
                del xy
            # create numpy array of points
            cs = np.array(points)
            # idetifies lowest point in cs and sets cs limit
            self.lwst_pnt = self.process_inner_cs(cs)
            # sets the cs outer limit and identifies highest point in
            # left and right direction
            hp_left, hp_right, cs = self.process_cs_limit(cs)            
            # identify lowest lwst_nbr direction (left [0] or right [1])
            # this is used to get first side of itertion
            values = [hp_left[1]]+[hp_right[1]]
            lwst_value = min(values)
            first_direction = values.index(lwst_value)
            # writes lowest point and Z value of lowest point to row
            row['lowest_PT'] = self.lwst_pnt[0]
            row['Z_low'] = self.lwst_pnt[1]
            if first_direction == 0:
                row = self.iterrate_cs(row, cs, lwst_value, -1, 1)
            elif first_direction == 1:
                row = self.iterrate_cs(row, cs, lwst_value, 1, -1)
            del cs, values, lwst_value, first_direction
            new.append(row)
        self.res = pd.DataFrame(new)
        self.res.to_csv(self.parameters.Xresults, index=False)
        return self.res

    def iterrate_cs(self, row, cs, lwst_value, v1, v2):
        '''
        iterates the cross section starting from the lowest point in left or
        right direction. Continues searching for corresponding neighbor in the
        opposite cs direction as long as no error in cs analysis occoures.
        - Analysis in first cs direction, returns side attributes (1)
        - Analysis in opposite direction (2) with (1) attributes
        - Re-Analyses first direction if (2) has changed the limiting
        attributes
        - calculates bankslope, depth and bottom width in corresponding modules
        - adds a plotting file path
        - return the analysis results as pandas row (series) in order to append
        to final results
        @type row: pandas.core.series.Series
        @param row: series to store results
        @type cs: np.array
        @param cs: cross section values
        @type lwst_value: <int>
        @param lwst_value: lowest point number in cs
        @type v1: <int>
        @param v1: cs direction
        @type v2: <int>
        @param v2: opposite cs direction
        '''
        if self.prints is True:
            print('iter:', 1)
        if v1 == -1:
            # attributes for left direction parametrisation
            btw_side = 'bottom_width_l'
            btw_side2 = 'bottom_width_r'
            w_side = 'width_l'
            w_side2 = 'width_r'
            bs_side1 = 'bankslope_l'
            bs_side2 = 'bankslope_r'
        elif v1 == 1:
            # attributes for right direction parametrisation
            btw_side = 'bottom_width_r'
            btw_side2 = 'bottom_width_l'
            w_side = 'width_r'
            w_side2 = 'width_l'
            bs_side1 = 'bankslope_r'
            bs_side2 = 'bankslope_l'
        # runs the analysis in the ininitial direction
        row, count1, nbr1, cs_error, height1, bw1 = self.an_dir(
                                            row, v1, cs, lwst_value)
        # if regular cs limit is found (cs_error = False), the tool iterates
        # in the opposite cs direction to find corresponding neighbor
        if cs_error is not True:
            row, cs_error = self.calc_bankslope(row, bs_side1, cs,
                                                nbr1, bw1, v1)
            row[btw_side] = bw1
            row[w_side] = count1*self.parameters.ticks
            if self.prints is True:
                print('iter:', 2)
            row, count2, nbr2, cs_error, height2, bw2 = self.an_dir(row, v2,
                                                                    cs,
                                                                    height1)
        # if regular cs limit is found (cs_error = False), the tool iterates
        # again in the opposite cs direction to find corresponding neighbor
        if cs_error is not True:
            row, cs_error = self.calc_bankslope(row, bs_side2,
                                                cs, nbr2, bw2, v2)
            row['bs_mean'] = float(np.mean([row[bs_side1], row[bs_side2]]))
            # identifies the direction
            if v2 == -1:
                self.calc_depth(row, height1, height2, nbr2, nbr1, cs, v2)
            if v2 == 1:
                self.calc_depth(row, height1, height2, nbr1, nbr2, cs, v2)
            row[btw_side2] = bw2
            row['bottom_width'] = bw1 + bw2
            row[w_side2] = count2*self.parameters.ticks
            row['width'] = (count1 + count2)*self.parameters.ticks
        # if regular cs limit is found (cs_error = False), the tool iterates
        # again in the opposite cs direction to find corresponding neighbor
        if cs_error is not True:
            if height2 < height1:
                if self.prints is True:
                    print('iter:', 3)
                row, count3, nbr3, cs_error, height3, bw3 = self.an_dir(
                        row, v1, cs, height2)
                if cs_error is not True:
                    row['width'] = (count2 + count3)*self.parameters.ticks
                    row[w_side] = count3*self.parameters.ticks
                    row, cs_error = self.calc_bankslope(row, 'bankslope_l',
                                                        cs, nbr3, bw3, v1)
                    row['bs_mean'] = float(np.mean([row['bankslope_l'],
                                                   row['bankslope_r']]))
                    row[btw_side] = bw3
                    row['bottom_width'] = bw3 + bw2
                if v1 == -1:
                    self.calc_depth(row, height1, height2, nbr2, nbr3, cs, v1)
                if v1 == 1:
                    self.calc_depth(row, height1, height2, nbr3, nbr2, cs, v1)
#        print('ID:', row.ID, row.reason_l, row.reason_r, round(row.lim_z_l, 2),
#              round(row.lim_z_r, 2), round(row['width'], 2))
        # A relative plotting path to the cross section plot is added
        row = self.add_plotting_path(row)
        del cs_error
        return row

    def add_plotting_path(self, row):
        '''
        Adds a path to plotting file
        @type row: pandas.core.series.Series
        @param row: series to store results
        '''
        row['plot'] = '%s\%s.png' % (self.parameters.plotting_path, row.ID)
        return row

    def calc_depth(self, row, height1, height2, nbr1, nbr2, cs, v):
        '''
        calculates the profile depth
        options:lowest - by lowest cs point (recommended)
                mean - by mean cs heights within cs (l +1  / r -1 )
        @type row: pandas.core.series.Series
        @param row: series to store results
        @type height1: <float>
        @param height1: height of point 1
        @type height2: <float>
        @param height2: height of point 2
        @type nbr1: <int>
        @param nbr1: next (inside) cs point to limit in 1st direction
        @type nbr2: <int>
        @param nbr2: next (inside) cs point to limit in 2nd direction
        @type cs: np.array
        @param cs: cross section values
        '''
        if v == -1:
            # attributes for left direction parametrisation
            depth_side1 = 'depth_r'
            depth_side2 = 'depth_l'
        elif v == 1:
            # attributes for right direction parametrisation
            depth_side1 = 'depth_l'
            depth_side2 = 'depth_r'
        if self.parameters.depth_method == 'lowest':
            row[depth_side1] = height1-row['Z_low']
            row[depth_side2] = height2-row['Z_low']
            row['depth'] = min(height1, height2)-row['Z_low']
        if self.parameters.depth_method == 'mean':
            nbr1, nbr2 = nbr1+1, nbr2-1
            p_lim = cs[nbr1:nbr2]
            mean_val = np.mean([i[1] for i in p_lim])
            row[depth_side1] = height1-row['Z_low']
            row[depth_side2] = height2-row['Z_low']
            row['depth'] = min(height1, height2) - mean_val
        if row['depth'] < 0:
            row['depth'] = np.NaN
            row['reason_l'] = 'Error: negative depth'
            row['reason_r'] = 'Error: negative depth'

    def calc_bankslope(self, row, bs_side, cs,
                       nbr, bw, v1):
        '''
        Calculates the inverse bankslope for each side
        @type row: pandas.core.series.Series
        @param row: series to store results
        @type bs_side: <str>
        @param bs_side: name of field for bs results e.g. 'bankslope_l'
        @type cs: np.array
        @param cs: cross section values
        @type nbr: <int>
        @param nbr: limiting cs point
        @type bw: <float>
        @param bw: bottom width
        @type v1: <int>
        @param v1: cs direction
        '''
        if v1 == -1:
            length = self.lwst_pos-nbr-bw
        elif v1 == 1:
            length = nbr-self.lwst_pos-bw
        bs = length / (cs[nbr, 1] - cs[self.lwst_pos, 1])
        if bs <= 0:
            self.write_error_row(row, v1)
            cs_error = True
        else:
            cs_error = False
            row[bs_side] = bs
        return row, cs_error

    def process_inner_cs(self, cs):
        '''
        Limits the inner cs calculation to a distance in left and right
        direction. The inner_limit is the point spacing from the middle-point
        of the cs. The reurn is the lowest point in the inner cs
        @type cs: np.Array
        @param cs: cross section points
        '''
        if self.parameters.inner_limit is None:
            left_lim = 1
            right_lim = self.n_points-1
        elif self.parameters.inner_limit >= self.n_points:
            left_lim = 1
            right_lim = self.n_points-1
        else:
            # identify lowest point between cs limit
            left_lim = round(
                    self.middle_point - self.parameters.inner_limit / 2 - 1)
            # cannot be smaller than the cs limit +1
            if left_lim <= 0:
                left_lim = 1
            right_lim = round(self.middle_point + self.parameters.inner_limit /
                              2 - self.n_points)
            # cannot be bigger than the cs limit-1
            if right_lim >= self.n_points:
                right_lim = self.n_points - 1
        self.lwst_pnt = self.get_lwst_point(cs, left_lim, right_lim)
        return self.lwst_pnt

    def process_cs_limit(self, cs):
        '''
        Limits the cs calculation to a distance in left and right
        direction from the lowest point. If no outer limit is set the cs limit
        is the highest point in the respective direction (left / right).
        If a outer limit is set the module searches the highest point within
        the outer limit boundary.
        @type cs: np.Array
        @param cs: cross section points
        @param self.n_points: number of ticks
        @type LINES: string
        @type DIST_LINE: float
        @param hp_left: <int>
        @type hp_left: Position of left cs limit
        @param hp_right: <int>
        @type hp_right: Position of right cs limit
        @param left_lim: <int>
        @type left_limt: Position of left limit
        @param right_lim: <int>
        @type right_limt: Position of right limit
        '''
        ref_point = int(self.lwst_pnt[0])
        outer_limit = self.parameters.outer_limit
        if outer_limit is None:
            left_lim = 0
            right_lim = self.n_points
        elif outer_limit >= self.n_points:
            left_lim = 0
            right_lim = self.n_points
        else:
            # creates index position for left limit (int)
            left_lim = int(ref_point - (outer_limit / 2)) - 1
            # cannot be smaller than the cs limit
            if left_lim < 0:
                left_lim = 0
            # creates index position for right limit (int)
            right_lim = int(ref_point + (outer_limit / 2)) + 1
            # cannot be bigger than the cs limit
            if right_lim > self.n_points:
                right_lim = self.n_points
        # slices the profile to the outer limit
        cs = cs[left_lim:right_lim]
        # derives highest points for left and right direction
        hp_left = self.get_hgst_point(cs, 0, ref_point-left_lim)
        hp_right = self.get_hgst_point(cs, ref_point-left_lim, -1)
        # gets the position of the lowset point in the new index
        self.lwst_pos = ref_point-left_lim
        # makes the left_prof_lim available
        self.left_prof_lim = left_lim
        return hp_left, hp_right, cs

    def get_lwst_point(self, cs, left_lim, right_lim):
        '''
        Identifies the lowest cs point within the left and right limit.
        @type cs: np.Array
        @param cs: cross section points
        @param left_lim: <int>
        @type left_limt: Position of left limit
        @param right_lim: <int>
        @type right_limt: Position of right limit
        '''
        cs_lp = cs[left_lim:right_lim]
        min_val = cs_lp.min(axis=0)[1]
        position = np.argmin(cs_lp, axis=0)[1]+left_lim
        self.lwst_pnt = np.array([int(position), min_val])
        return self.lwst_pnt

    def get_hgst_point(self, cs, left_lim, right_lim):
        '''
        Identifies the highest cs point within the left and right limit
        @type cs: np.Array
        @param cs: cross section points
        @param left_lim: <int>
        @type left_limt: Position of left limit
        @param right_lim: <int>
        @type right_limt: Position of right limit
        '''
        cs_lp = cs[left_lim:right_lim]
        max_val = cs_lp.max(axis=0)[1]
        position = np.argmax(cs_lp, axis=0)[1]+left_lim
        hgst = np.array([int(position), max_val])
        del cs_lp, max_val, position
        return hgst

    def an_dir(self, row, dr, cs, highest_point):
        '''
        Analysis the cross section direction.
        Loops through one cs side (left or right from the lowset point)
        and constructs the cs until a condition is fulfilled or the cs point
        limit is reached
        @param cs: cross section points
        @type cs: np.Array
        @param row: series to store results
        @type row: pandas.core.series.Series
        @param highest_point: Z value of highest point
        @type highest_point: <float>
        @param bott_add: bottom width of direction
        @type bott_add: <float>
        @param count: calculates the number of loops till cs construction
                    finished this is used to calculate the profile width
        @type count: <int>
        @param bott_add: bottom width of direction
        @type bott_add: <float>
        @param bot_fin: if True the bottom width calculation is done
        @type bot_fin: <boolean>
        @param cs_error: if True the cs could not be created
        @type cs_error: <boolean>
        @param angle_sum: total sum (aggregated differences) of angles
        @type angle_sum: <float>
        @param angle_add: single angle calculated from position to neighbor
        @type angle_add: <float>
        @param angle_former: value of last angle calculation
        @type angle_former: <float>
        @param angle_check: difference of angle to former angle
        @type angle_check: <float>
        @param neg_angle: used to define if a negative angle is calculated
        @type neg_angle: <float>
        @param count_neg: counts the number of negative angles (used for
                            self.parameters.n_sum_allowed)
        @type count_neg: <int>
        '''
        # sets the starting position to the lowest point
        position = self.lwst_pos
        if dr == 1:
            side_lim_z = 'lim_z_r'
            limit_position = 'lim_pt_r'
            reason = 'reason_r'
            # index position right is different to left
            n_positions = len(cs[position:])-1
        elif dr == -1:
            side_lim_z = 'lim_z_l'
            limit_position = 'lim_pt_l'
            reason = 'reason_l'
            n_positions = len(cs[:position])
        # sets the first neighbor in defined (dr) direction
        nbr = position+1*dr
        # gets Z value of lowest point
        lwst_Z = cs[position][1]
        # sets counter to 1
        count = 0
        # defines the condition parameters for comparison
        bott_add, angle_sum, angle_add, neg_angle, angle_former = 0, 0, 0, 0, 0
        # helper parameters for condtions
        bot_fin, cs_error = False, False  
        for num in range(1, n_positions+1):
            # iterrates through remaining points
            if count != 0:
                position = position + 1*dr
                nbr = nbr + 1*dr
            count = count + 1
            if self.prints is True:
                print(side_lim_z , position, nbr, num, n_positions,
                      cs[nbr][1], highest_point)
            # if position is already the last point, the loop stops
            if num == n_positions:
                if self.prints is True:
                    print(6)
                cs_error = True
                row[limit_position] = nbr
                row[reason] = 'Error: cs limit reached'
                break
            # calculates the angle of current position and neighbor
            angle_add = self.calc_sin(cs, nbr, position)
            if bot_fin is False:
                # checks / extents the bottom width
                bott_add, bot_fin = self.calc_bottom_width(cs, nbr, lwst_Z,
                                                           bott_add, bot_fin,
                                                           angle_add)
            # adds to the negative angle sum if angle is negative
            if angle_add < 0:
                neg_angle = neg_angle+angle_add
            # Counts the occourrance of negative angles while looping
            count_neg = 0
            # checks angle_sum < negative angle (check for slope)
            if angle_sum < neg_angle*-1:
                count_neg = count_neg + 1
            # checks if angle is positive
            elif angle_add > 0:
                # calculates the difference to former angle
                angle_check = angle_add-angle_former
#                if angle_check > 0:
                    # add diff angle1-angle2 to angle sum
                angle_sum = angle_sum+angle_check
                angle_former = angle_add
            if count_neg == self.parameters.n_sum_allowed:
                # checks if count_neg is within the condition
                row[limit_position] = (nbr - 1 * dr)+self.left_prof_lim
                row[side_lim_z] = cs[nbr - 1 * dr][1]
                row[reason] = 'angle_sum < neg_angle*-1'
                count = count-1
                break
            if self.parameters.Angle_negative is not None:
                # checks if condition is defined
                if neg_angle < self.parameters.Angle_negative:
                    # checks if neg_angle is within the condition
                    if self.prints is True:
                        print(2)
                    if dr == -1:
                        nbr = nbr+2
                    elif dr == 1:
                        nbr = nbr-2
                    row[limit_position] = nbr+self.left_prof_lim
                    row[side_lim_z] = cs[nbr][1]
                    row[reason] = 'angle < %s' % str(
                            self.parameters.Angle_negative)
                    count = count-1
                    cs_error = True
                    break
            elif self.parameters.Angle_sum is not None:
                # checks if condition is defined
                if angle_sum >= self.parameters.Angle_sum:
                    # checks if angle_sum is within the condition
                    if self.prints is True:
                        print(3)
                    row[limit_position] = nbr + self.left_prof_lim
                    row[side_lim_z] = cs[nbr][1]
                    row[reason] = 'angle sum > %s' % str(
                            self.parameters.Angle_sum)
                    break
            elif self.parameters.Angle_max is not None:
                # checks if condition is defined
                if angle_add >= self.parameters.Angle_max:
                    # checks if angle_add is within the condition
                    if self.prints is True:
                        print(4)
                    if dr == 1:
                        nbr += -1
                    elif dr == -1:
                        nbr += 1
                    row[limit_position] = nbr + self.left_prof_lim
                    row[side_lim_z] = cs[nbr][1]
                    row[reason] = 'angle max > %s' % str(
                            self.parameters.Angle_max)
                    count = count-1
                    break
            if cs_error is True:
                break
            elif cs[nbr][1] >= highest_point:
                # checks if height is bigger or equal to highest_point
                if self.prints is True:
                    print(5)
                row[side_lim_z] = cs[nbr][1]
                row[limit_position] = nbr + self.left_prof_lim
                row[reason] = 'height %s' % side_lim_z
                break
        del neg_angle, angle_sum, angle_add, angle_former
        if cs_error is True:
            # writes dummy NaN to result file
            row, height = self.write_error_row(row, dr)
        else:
            # writes the height limit value to row
            height = row[side_lim_z]
        if self.prints is True:
            print(side_lim_z, row[side_lim_z], row[limit_position])
        return row, count, nbr, cs_error, height, bott_add

    def calc_bottom_width(self, cs, nbr, lwst_Z, bott_add, bot_fin, angle_add):
        '''
        checks bottom width if bot_fin == False and / or parameter
        conditions are met.
        Part of module an_dir
        Returns bott_add, bot_fin
        @param cs: cross section points
        @type cs: np.Array
        @param nbr: <int>
        @type nbr: neighbor point in cs profile
        @param lwst_Z: <float>
        @type lwst_Z: Z value of lowest point
        @param lwst_Z: <float>
        @type lwst_Z: Z value of lowest point
        @param bott_add: bottom width of direction
        @type bott_add: <float>
        @param bot_fin: if True the bottom width calculation is done
        @type bot_fin: <boolean>
        @param angle_add: single angle calculated from position to neighbor
        @type angle_add: <float>
        '''
        if self.parameters.bottom_height_lim is not None:
            # checks if parameter bottom_height_lim is set
            if cs[nbr][1] - lwst_Z <= self.parameters.bottom_height_lim:
                bott_add = bott_add + self.parameters.ticks
                bot_fin = False
            else:
                if self.parameters.angle_bot_lim is not None:
                    # checks if parameter angle_bot_lim is set
                    if angle_add <= self.parameters.angle_bot_lim:
                        bott_add = bott_add + self.parameters.ticks
                        bot_fin = False
                else:
                    bot_fin = True
        else:
            if self.parameters.angle_bot_lim is not None:
                # checks if parameter angle_bot_lim is set
                if angle_add <= self.parameters.angle_bot_lim:
                    bott_add = bott_add + self.parameters.ticks
                    bot_fin = False
            else:
                bott_add = bott_add + self.parameters.ticks
                bot_fin = True
        return bott_add, bot_fin

    def write_error_row(self, row, dr):
        '''
        writes np.NaN to current cs row if error in analysis occoured
        Used in module an_dir
        @type row: pandas.core.series.Series
        @param row: series to store results
        '''
        row['lim_z_l'] = np.NaN
        row['lim_z_r'] = np.NaN
        row['lim_pt_l'] = np.NaN
        row['lim_pt_r'] = np.NaN
        row['bankslope_l'] = np.NaN
        row['bankslope_r'] = np.NaN
        row['depth'] = np.NaN
        row['width'] = np.NaN
        row['bs_mean'] = np.NaN
        if dr == -1:
            row['reason_r'] = 'Error: opposite side'
        elif dr == 1:
            row['reason_l'] = 'Error: opposite side'
        height = np.NaN
        return row, height

    def calc_sin(self, cs, point1, point2):
        '''
        calculates the sinus from two points and returns the angle
        '''
        a = cs[point1, 1]-cs[point2, 1]
        c = np.sqrt(a**2+self.parameters.ticks**2)
        angle = np.degrees(np.sin(a/c))
        return angle