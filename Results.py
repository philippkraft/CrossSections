# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 14:48:16 2020

@author: Florian Krebs
"""
import pandas as pd
import os

class Results:    
    def __init__(self, parameters, df_input, df_out, ext):
        # get Parameters from Parameters class
        if type(parameters) != str:
            self.parameters = parameters
        elif type(parameters) == str:
            from Parameters import Parameters
            self.parameters = Parameters(parameters)
        self.df_input = df_input
        self.df_out = df_out
        self.ext = ext
        self.res = None
        self.errors = ['Error: opposite side', 'Error: cs limit reached',
                       'angle_sum < neg_angle*-1']
        self.__call__()
    
    def __call__(self):
        self.res = pd.read_csv(self.df_input)
        self.__aggregate_reason(self.res)
        print(self.agg_results_r)
        self.res = self.filter_error()
        # aggregates the raw results (no filter)
        self.aggregate(self.res)
        filter_ = self.create_filter()
        self.aggregate(filter_)
        print(self.agg_results)
        return self.agg_results

    
    def create_filter(self):
        '''
        Filters the results by max bankslope, minimum depth and 
        intersection value if a respective parameter is set or calculated.
        
        @type self.res: Dataframe with results
        @param self.res: pd.DataFrame
        '''
        if self.parameters.filter_bs_max is not None:
            res = self.res[
                    self.res.bs_mean <= self.parameters.filter_bs_max]
        else:
            res = self.res
        if self.parameters.filter_depth_min is not None:    
            res = res[
                    res.depth >= self.parameters.filter_depth_min]
        else:
            res = self.res             
        l_intsct = res[res.intsct != 1]
        l_intsct.to_csv(self.df_out, index=False)
        return l_intsct
    
    
    def filter_error(self):
        res_no_err = self.res[~self.res['reason_l'].isin(self.errors)]
        res_no_err = res_no_err[~res_no_err['reason_r'].isin(self.errors)]
        return res_no_err
    
    def aggregate(self, res):
        '''
        looks for the appropriate aggregation function according to the 
        defined parameters in the parameter class
        
        @type self.res: Dataframe with results
        @param self.res: pd.DataFrame
        '''       
        if self.parameters.strahler_field is not None:
            self.__aggregate_by_strahler(res)
            if self.parameters.Aggregation_field is not None:
                self.__aggregate_by_field_strahler(res)                
        elif self.parameters.Aggregation_field is not None:
                self.__aggregate_by_field(res)

    def __aggregate_reason(self, res):
        '''
        Aggregates the results by cs limiting reason
        
        @type self.res: Dataframe with results
        @param self.res: pd.DataFrame
        '''
        df = res[res['intsct'] != 1]
        self.agg_results_r = df.groupby(
               ['reason_l', 'reason_r']
            ).agg(
                {    'reason_l': 'count',
                     'reason_r': 'count'
                }
            )
        self.agg_results_r.to_csv(os.path.join(self.parameters.results_path,
                                             'results_agg_%s.csv' % self.ext))

    def __aggregate_by_field(self, res):
        '''
        Aggregates the results on the defined aggregation_field and 
        cs limiting reason
        
        @type self.res: Dataframe with results
        @param self.res: pd.DataFrame
        '''
        df = res[res['intsct'] != 1]
        self.agg_results = df.groupby(
               [self.parameters.Aggregation_field]
            ).agg(
                {
                     'width': "mean",  # median value of width
                     'depth': 'mean',  # mean value of depth
                     'bs_mean': 'mean',  # median value of bankslope
                     'bottom_width' : 'mean',
                     'reason_l': 'count',
                     'reason_r': 'count'
                }
            )
        self.agg_results.to_csv(os.path.join(self.parameters.results_path,
                                             'results_agg_%s_%s.csv' %
                                             (self.parameters.Aggregation_field,
                                              self.ext)))

    def __aggregate_by_strahler(self, res):
        '''
        Aggregates the results on the defined strahler_field by 
        cs limiting reason
        
        @type self.res: Dataframe with results
        @param self.res: pd.DataFrame
        '''
        df = res[res['intsct'] != 1]
        self.agg_results = df.groupby(
               [self.parameters.strahler_field]
            ).agg(
                {
                     'width': "mean",  # median value of width
                     'depth': 'mean',  # mean value of depth
                     'bs_mean': 'mean',  # median value of bankslope
                     'bottom_width' : 'mean',
                     'reason_l': 'count',
                     'reason_r': 'count'
                }
            )
        self.agg_results.to_csv(os.path.join(self.parameters.results_path,
                                             'results_agg_%s_%s.csv'
                                             % (self.parameters.strahler_field,
                                                self.ext)))

    def __aggregate_by_field_strahler(self, res):
        '''
        Aggregates the results on the defined aggregation_field and strahler 
        field by cs limiting reason       
        @type self.res: Dataframe with results
        @param self.res: pd.DataFrame
        '''
        df = res[res['intsct'] != 1]
        self.agg_results = df.groupby(
               [self.parameters.Aggregation_field,
                self.parameters.strahler_field]
            ).agg(
                {
                     'width': "mean",  # median value of width
                     'depth': 'mean',  # mean value of depth
                     'bs_mean': 'mean',  # median value of bankslope
                     'bottom_width' : 'mean',
                     'reason_l': 'count',
                     'reason_r': 'count'
                }
            )
        self.agg_results.to_csv(os.path.join(self.parameters.results_path,
                                             'results_agg_%s_%s_%s.csv' %
                                             (self.parameters.strahler_field,
                                              self.parameters.Aggregation_field,
                                              self.ext)))

if __name__ == "__main__":
#    from Parameters import Parameters
    config_file = r"C:\Users\Scuff\OneDrive\phd\Topics\CrossSections\Example.yaml"
    file = r"C:\Users\Scuff\OneDrive\phd\Topics\CrossSections\Example\results_depth_corr.csv"
#    Parameters = Parameters(config_file)
    x = Results(config_file, file, 'dc')