# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 13:24:28 2019
Class Plot for the creation of Cross-Section Plots
@author: Florian Krebs
"""

import pandas as pd
import matplotlib.pyplot as plt
from os import path, makedirs
import numpy as np


class Plotting:
    def __init__(self, parameters):
        # get Parameters from Parameters class
        if type(parameters) != str:
            self.parameters = parameters
        elif type(parameters) == str:
            from Parameters import Parameters
            self.parameters = Parameters(parameters)
        self.main()

    def main(self):
        # creates the plotting folder
        self.create_folder(self.parameters.plotting_path)
        # reads results filter file
        self.Xresults = pd.read_csv(self.parameters.Xresults_filter)
        # get depth extrapolation if processed
        if self.parameters.depth_extr is True:
            self.dext = pd.read_csv(self.parameters.Xresults_dcorr)         
        # sets list of height columns
        self.height_cols = [col for col in self.Xresults.columns if 'X' in col]
        # checks is a inner limit is set and gets the coordinates
        if self.parameters.inner_limit != None:
           self.il_coords = self.inner_limit_coords()
        # plots all profiles in list
        self.plot_all_ID()

    def plot_all_ID(self):
        '''
        Plots all Profile IDs
        @type self.Xresults: pd.DataFrame
        @param self.Xresults: Pandas DataFrame including profile data
        '''
        for cs in self.Xresults.ID.unique():
            print('Plotting cs:', cs)
                    # creates the destination path of the plot
            self.file = path.join(
                    self.parameters.plotting_path,
                    '%s.png' % cs)
            if self.parameters.depth_extr is True:
                self.cs_depth_ext(cs)
            self.__select_cs(cs)

    def __select_cs(self, cs_ID):
        '''
        Selects the cs ID from the Pandas dataframe and plots x,y plots
        of the cs. Joins additional plotting in
        @type self.Xresults:   pd.DataFrame
        @param self.Xresults: Pandas DataFrame including cs data
        @type self.height_cols: <list>
        @type self.height_cols: List of height columns
        '''
        # filters the cs from the list
        cs = self.Xresults.loc[self.Xresults.ID == cs_ID]
        # creates a list of x_ticks for plotting
        x_ticks = np.arange(0, len(self.height_cols)*self.parameters.ticks,
                            self.parameters.ticks)
        # creates x and y points
        x = [i for i in x_ticks]
        y = [float(cs[i]) for i in self.height_cols]
        # gets the left and right limiting point coordinates
        l_lim = (cs.lim_pt_l.values[0]*self.parameters.ticks,
                 cs.lim_z_l.values[0])
        r_lim = (cs.lim_pt_r.values[0]*self.parameters.ticks,
                 cs.lim_z_r.values[0])
        # gets the lowest point coordinates
        lwst = (cs.lowest_PT.values[0]*self.parameters.ticks,
                cs.Z_low.values[0])       
        # gets attributes and schematic plots to add to 
        attrib, p_plot, c_plot = self.set_attributes(cs, lwst)
        # plots the cross section results
        self.plot(x, y, cs_ID, x_ticks, l_lim, r_lim,
                  lwst, attrib, p_plot, c_plot)

    def set_attributes(self, cs, lwst):
        '''
        sets the attributes for additional plotting informations
        @type attrib:   <list>
        @param attrib: List containing plot informations
        @type p_plot: <list>
        @type p_plot: List containing schematic cs points
        @type c_plot: <list>
        @type c_plot: List containing potential depth construction
        @type ol_coords: <list>
        @type ol_coords: List containing position of outer limit
        '''
        attrib = (cs.reason_l.values[0], 
                  cs.reason_r.values[0],
                  round(cs.depth.values[0],2),
                  round(cs.width.values[0],2),
                  round(cs.bankslope_l.values[0],2),
                  round(cs.bankslope_r.values[0],2),
                  cs.bottom_width.values[0])
        # checks if a strahler field is set
        if self.parameters.strahler_field != None:
            strahler = cs[self.parameters.strahler_field].values[0]
            if type(strahler) == float:
                strahler = int(strahler)
            attrib = attrib + (strahler,)
        else:
            attrib = attrib+(0,)
        # constructs the schematic cs for plotting
        p_plot = self.construct_cs(cs, lwst)
        # checks if depth construction was applied
        if self.parameters.depth_extr is True:
            attrib = attrib+(round(cs.depth_ext.values[0],2),)
            c_plot = self.construct_depth(cs, lwst)
        else:
            attrib = attrib+(0,)
            c_plot = None
        # checks if a outer limit is set
        if self.parameters.outer_limit != None:
           self.ol_coords = self.outer_limit_coords(cs)
        elif self.parameters.outer_limit == None:
           self.ol_coords = [0, len(cs)]        
        return attrib, p_plot, c_plot

    def inner_limit_coords(self):
        '''
        gets the inner limit to show in plot as lines
        '''
        left = ((len(self.height_cols)/2)-(
                self.parameters.inner_limit/2))*self.parameters.ticks
        right = ((len(self.height_cols)/2)+(
                self.parameters.inner_limit/2))*self.parameters.ticks
        il_coords = [left, right]
        return il_coords

    def outer_limit_coords(self, cs):
        '''
        gets the outer limit to show in plot as lines
        '''
        left = (cs.lowest_PT.values[0]-(
                self.parameters.outer_limit/2))*self.parameters.ticks
        right = (cs.lowest_PT.values[0]+(
                self.parameters.outer_limit/2))*self.parameters.ticks
        ol_coords = [left, right]
        return ol_coords

    def construct_cs(self, cs, lwst_pnt):
        '''
        constructs the schematic cs to show in plot
        '''
        # gets all parameters from cs
        depth = cs.depth.values[0]
        width_l = cs.width_l.values[0]
        width_r = cs.width_r.values[0]       
        bott_l = cs.bottom_width_l.values[0]
        bott_r = cs.bottom_width_r.values[0]
        z_ref = cs.Z_low.values[0]
        # assings the points
        p1 = (lwst_pnt[0]-bott_l, z_ref)
        p2 = (lwst_pnt[0]+bott_r, z_ref)
        p3 = (lwst_pnt[0]+width_r, z_ref+depth)
        p4 = (lwst_pnt[0]-width_l, z_ref+depth)
        p_plot = (p1, p2, p3, p4, p1)
        return p_plot

    def construct_depth(self, cs, lwst_pnt):
        '''
        (optional)
        constructs a potential depth cs calculated in depth_extr
        ''' 
        # gets all parameters from cs
        depth = cs.depth_ext.values[0]
        b_width_l = cs.bottom_width_l.values[0]
        b_width_r = cs.bottom_width_r.values[0]
        bott_l = cs.bw_l_ext.values[0]
        bott_r = cs.bw_r_ext.values[0]
        z_ref = cs.Z_low.values[0]
        p1 = (lwst_pnt[0]+b_width_r, z_ref)
        p2 = (lwst_pnt[0]-b_width_l, z_ref)
        # checks if trapezoid or triangle
        if bott_r != -99:
            p3 = (lwst_pnt[0]-b_width_l+bott_l, z_ref-depth)
            p4 = (lwst_pnt[0]+b_width_r-bott_r, z_ref-depth)
            c_plot = (p1, p2, p3, p4, p1)
        else:
            p3 = (lwst_pnt[0]-b_width_l+bott_l, z_ref-depth)
            c_plot = (p1, p2, p3, p1)
        return c_plot
    
    def cs_depth_ext(self, cs_ID):
        '''
        Selects the cs ID from the Pandas dataframe and plots x,y plots
        of the cs. Joins additional plotting in
        @type self.Xresults:   pd.DataFrame
        @param self.Xresults: Pandas DataFrame including cs data
        @type self.height_cols: <list>
        @type self.height_cols: List of height columns
        '''
        # filters the cs from the list
        
        cs = self.dext.loc[self.dext.ID == cs_ID]
        # creates a list of x_ticks for plotting
        x_ticks = np.arange(0, len(self.height_cols)*self.parameters.ticks,
                            self.parameters.ticks)
        # creates x and y points
        x = [i for i in x_ticks]
        y = [float(cs[i]) for i in self.height_cols]
        self.dext_cs = [x, y]
    
    def plot(self, x, y, cs_ID, x_ticks, l_lim, r_lim,
             lwst, attrib, p_plot, c_plot):
        '''
        Plots x,y points of  and adds all attributes
        '''
        # creates a figure
        fig, ax = plt.subplots(figsize=(12, 6))
        plt.ioff()
        ax.grid(True)
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='gray', linestyle='dashed', label='cross section')
        if self.parameters.depth_extr is True:
            ax.plot([x for x in self.dext_cs[0]],
                    [y for y in self.dext_cs[1]], '-', linestyle='dashed',
                    label='depth_ext', color='black')           
        # plots the regular cs
        ax.plot(x, y, '-')
        # shows the lowest point
        ax.scatter(lwst[0], lwst[1], marker='o', label='lowest point',
                   color='green')
        # shows the limiting points on both sides
        ax.scatter(l_lim[0], l_lim[1], marker='o', label='limiting point',
                   color='red')
        ax.scatter(r_lim[0], r_lim[1], marker='o', color='red')
        # plots the shematic plot
        ax.plot([i[0] for i in p_plot], [i[1] for i in p_plot], '-',
                color='gray')
        # (optional) plots the potential depth
        if c_plot is not None:
            ax.plot([i[0] for i in c_plot], [i[1] for i in c_plot], '--',
                    color='red', label='Potential depth')
        # shows the inner limit (if set)
        if self.parameters.inner_limit != None:
            plt.axvline(x=self.il_coords[0], color='orange',
                        linestyle='dashed',
                        label='Radius Limit')
            plt.axvline(x=self.il_coords[1], color='orange',
                        linestyle='dashed')
        # shows the outer limit (if set)
        if self.parameters.outer_limit != None:
            plt.axvline(x=self.ol_coords[0], color='green',
                        linestyle='dashed',
                        label='search boundary')
            plt.axvline(x=self.ol_coords[1], color='green',
                        linestyle='dashed')
        # adds labels
        plt.ylabel('Height [m]')
        plt.xlabel('Length [m]')
        # sets axis x limit
        plt.xlim(self.parameters.ticks*-1, len(
                self.height_cols)*self.parameters.ticks+self.parameters.ticks)
        ax.legend(loc=0, shadow=True, fontsize='small')
        # defines box properties
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        # constructs string sequence used for text box
        reason = "reason: [L] %s \n              [R] %s" % (attrib[0],
                          attrib[1])
        stats =  "d: %s (%s) w: %s bw: %s" % (attrib[2], attrib[8],
                        attrib[3],attrib[6])
        bs = "bankslope [L] %s [R] %s" % (attrib[4], attrib[5])
        # merges textstring
        if self.parameters.strahler_field != None:
            textstr = '%s \n %s \n %s   s: %s' % (reason, stats, bs, attrib[7])
        else:
            textstr = '%s \n %s \n %s' % (reason, stats, bs)
        # plots box
        plt.text(0.01, 0.02, textstr, fontsize=9, 
                 transform=plt.gcf().transFigure, color='black', bbox=props)
        # adds title
        plt.title('Cross section: %s' % cs_ID)
        plt.subplots_adjust(bottom=.2)
        # saves figure
        plt.savefig(self.file)
        plt.clf()
        plt.close()
        del fig, ax, x, y, l_lim, r_lim, lwst, attrib, p_plot, c_plot

    def create_folder(self, directory):
        '''
        creates folder in target directory
        '''
        try:
            if not path.exists(directory):
                makedirs(directory)
        except OSError:
            print('Error: Creating directory. ' + directory)    