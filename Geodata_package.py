# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 11:01:35 2019

@author: Florian Krebs
"""

from shapely.geometry import LineString, Point
import shapely.ops as ops
import math
import geopandas as gpd
import pandas as pd
from shapely.ops import split, snap
import rasterio as rio
import csv
import os
import numpy as np


class Geodata_package:
    def __init__(self, parameters):
        self.parameters = parameters
        self.point_gdf = None
        self.point = None
        self.cross_sections = None
        self.crs = None
        self.stream_gpd = None
        self.fields_add = None
        self.__call__()

    def __call__(self):
        '''
        runs the geodata creation class
        - creates a results folder
        - reads the stream geometries
        - gets and checks crs of stream geometries and dtm
        - creates the cross sections (ID)
        - creates points on the cross sections (ID as identifier) by tick
        - adds Z information to points
        - writes a result file
        @param self.crs: <string>
        @type self.crs: Reference coordinate system of stream geometries
        @param self.steam_gpd: <gpd.GeodataFrame>
        @type self.stream_gpd: Stream geometries with attributes
        @param self.split_line: <gpd.GeodataFrame>
        @type self.split_line: Stream splitted by parameters.distance
        @param self.cross_sections: <gpd.GeodataFrame>
        @type self.cross_sectionse: Cross sections
        '''
        # creates a results folder
        self.create_folder(self.parameters.results_path)
        # reads the stream geometries
        self.stream_gpd = gpd.read_file(self.parameters.Stream_shapefile)
        # gets and checks crs of stream geometries and dtm
        self.crs = self.stream_gpd.crs
        # checks if crs match
        self.check_crs()
        print('Create Crossections')
        # creates split lines
        self.split_line = self.create_split_lines()
        # creates cross sections on split lines
        self.cross_sections = self.create_CrossSections()
        self.create_points()
        print('CS done')

    def check_crs(self):
        '''
        checks the input crs dtm
        @param dtm_open: <rasterio raster>
        @type self.cross_sectionse: Cross sections
        @param dtm_crs: <EPSG>
        @type dtm_crs: coordinate system of dtm
        '''
        dtm_open = rio.open(self.parameters.dtm)
        dtm_crs = dtm_open.crs.data
        dtm_open.close
        if dtm_crs != self.crs:
            print('Warning!!! CRS of dtm and shapefile do not match')
        else:
            print('CRS check successful')

    def split_MLS_dist(self, MultiLineString, distance):
        '''
        Splits the MultiLineString stream geometries into defined distance
        segments used for cs creation. Creates points along the line to split
        @param MultiLineString: shapely MultiLineString
        @type MultiLineStrings: Stream geometries MultiLineString
        @param distance: <float>
        @type distance: distance used to spolit the lines
        @param current_dist: <float>
        @type current_dist: sum of distances (aggregates to total length)
        @param segments: <list>
        @type segments: list containing the new created line segments
        @param current_dist: <float>
        @type current_dist: distance added by each distance to sum up length
        '''
        # list to hold the line segments
        segments = []
        for line in MultiLineString:
            list_points = []
            # set the current distance to place the point
            current_dist = distance
            # get the total length of the line
            line_length = line.length
            # append the starting coordinate to the list
            list_points.append(Point(list(line.coords)[0]))
            # while the current cumulative distance is less than the
            # total length of the line
            while current_dist < line_length:
                # use interpolate and increase the current distance
                list_points.append(line.interpolate(current_dist))
                current_dist += distance
            # append end coordinate to the list
            list_points.append(Point(list(line.coords)[-1]))
            points = ops.cascaded_union(list_points)
            # snaps the points to line and splits the line (tolerance 0.01)
            line_snap = snap(line, points, 0.01)
            line_split = split(line_snap, points)
            for i in line_split:
                segments.append(i)
        return segments

    def split_LS_dist(self, line, distance):
        '''
        Splits a single LineString into defined distance segments
        Used for cross section creation from stream shapefile and from the 
        created cross section LineString in the module create points.
        @type distance: distance used to spolit the lines
        @param current_dist: <float>
        @type current_dist: sum of distances (aggregates to total length)
        @param line_split: <list>
        @type line_split: list containing the new created line segments
        @param current_dist: <float>
        @type current_dist: distance added by each distance to sum up length
        @param position: <int>
        @type position: Position of point within the cross section
        @param points_lst: <list>
        @type points_lst: holding the cross section point geometries
        '''
        # list to hold all the point coords
        list_points = []
        # set the current distance to place the poin
        current_dist = distance
        # get the total length of the line
        line_length = line.length
        # append the starting coordinate to the list
        list_points.append(Point(list(line.coords)[0]))
        # while the current cumulative distance is less than the total
        # length of the line
        while current_dist < line_length:
            # use interpolate and increase the current distance
            list_points.append(line.interpolate(current_dist))
            current_dist += distance
        # append end coordinate to the list
        list_points.append(Point(list(line.coords)[-1]))
        # uninions all points
        points = ops.cascaded_union(list_points)
        # snaps the lines to points
        line_snap = snap(line, points, 0.01)
        # splits the line by snapped points
        line_split = split(line_snap, points)
        # this is used for the creation of the final points in the module
        # create_points.
        points_lst = []
        count = 0
        position = []
        for i in points:
            points_lst.append(i)
            position.append(count)
            count = count + 1
        del count
        return line_split, points_lst, position

    def create_split_lines(self):
        '''
        Creates split lines from stream geometries. Uses the distance
        parameter. Checks the geometry first if stream shapees are a 
        Multilinestring or Linestring for processing purposes. 
        @param gdf: geopandas GeoDataFrame
        @type gdf: Stream geometries gpd. self.stream_gpd
        @param distance: <float>
        @type distance: distance used to split the lines. parametes.distance
        @param line_gpd: <gpd.GeoDataFrame>
        @type line_gpd: Splitted lines
        '''
        distance = self.parameters.distance
        gdf = self.stream_gpd
        print('Create split lines with %s m ticks' % distance)
        # checks the type of geometry
        if gdf.geometry[0].geom_type == 'MultiLineString':
            geoms = []
            for i, r in gdf.iterrows():
                # creates a list of LineStrings
                [geoms.append(i) for i in r.geometry]
        elif gdf.geometry[0].geom_type == 'LineString':
            geoms = []
            for i, r in gdf.iterrows():
                # creates a list of LineStrings
                geoms.append(r.geometry)
        # merges the single geometries to one feature
        merged_line = ops.linemerge(geoms)
        # checks if a multi or single linestring is created
        if merged_line.geom_type == 'MultiLineString':
            line_split_d = self.split_MLS_dist(merged_line, distance)
        elif merged_line.geom_type == 'LineString':
            line_split_d, points_lst, position = self.split_LS_dist(
                    merged_line, distance)
        # creates GeoPandasDataframe
        line_gpd = gpd.GeoDataFrame(crs=self.crs,
                                    geometry=[i for i in line_split_d])
        return line_gpd

    def __join_attributes(self, left_df):
        '''
        joins cross-sections and stream attributes
        used in module join_intsct_attributes
        @param left_df: geopandas GeoDataFrame
        @type left_df: geodatafram used to join to cross sections
        '''
        self.fields_add = []
        if self.parameters.strahler_field is not None:
            self.fields_add.append(self.parameters.strahler_field)
        if self.parameters.Aggregation_field is not None:
            self.fields_add.append(self.parameters.Aggregation_field)
        if 'art_depth' in self.stream_gpd.dtypes:
            self.fields_add.append('art_depth')
        if self.fields_add != []:
            right_df = self.stream_gpd[self.fields_add+['geometry']]
            left_df = gpd.sjoin(left_df, right_df, how='inner',
                                op='intersects')
            left_df = left_df.drop('index_right', axis=1)
            return left_df
        else:
            return left_df

    def join_intsct_attributes(self, cross_sections, structures):
        '''
        Joins the attributes from cross section and structures to the
        cross sections
        joins cross-sections and create_structures
        @param cross_sections: geopandas GeoDataFrame
        @type cross_sections: cross sections
        @param structures: geopandas GeoDataFrame
        @type structures: structures
        '''
        left_df = cross_sections
        left_df = left_df.to_crs(crs=self.crs)
        right_df = structures
        right_df = right_df.to_crs(crs=self.crs)
        right_df['intsct'] = 1
        right_df = right_df[['intsct', 'geometry']]
        join = gpd.sjoin(left_df, right_df, how='left', op='intersects')
        join['intsct'] = join['intsct_right']
        join = join.drop(columns=['intsct_left', 'intsct_right',
                                  'index_right'])
        return join

    def create_CrossSections(self):
        '''
        creates the cross sections along the splitted lines. gets the first
        and last point of the segmented line. The angle of both points
        is calculated and used as bearing angle to create 90° / -90° points for
        left and right side of the line using the width/2 for each side.
        @param geoline: pandas DataFrame
        @type geoline: dataframe containing geometries and angles
        @param cross_sections: geopandas GeoDataFrame
        @type cross_sections: cross sections
        @param structures: geopandas GeoDataFrame
        @type structures: structures
        '''
        cs_line = []
        cs_points = []
        # iterrates the splitted line geometries
        for index, row in self.split_line.iterrows():
            # gets first and last point of line
            firstPoint, scndPoint = self.get_points(row.geometry)
            # sets width from parameter
            width = self.parameters.width
            # calculates the angle between the first and last point
            angle1 = self.getAngle(firstPoint, scndPoint)
            # gets the first cross section point
            p1 = self.getPoint1(firstPoint, angle1, width/2)
            # calculates the angle of the first point and new created point
            angle2 = self.getAngle(p1, firstPoint)
            # uses the first cross section point, width and angle to create
            # the opposite point
            p2 = self.getPoint2(p1, angle2, width)
            # creates the cross section LineString using both calculated points
            line_cross = LineString([p1, p2])
            # append all attributes to list
            cs_line.append((line_cross, angle1, angle2))
            cs_points.append(p1)
            cs_points.append(p2)
        # create a dataframe from list
        geoline = pd.DataFrame(cs_line, columns=['geometry',
                                                 'angle1',
                                                 'angle2'])
        # create a Geopandas dataframe from dataframe
        cross_sections = gpd.GeoDataFrame(geoline, crs=self.crs,
                                          geometry=geoline.geometry)
        # defines a unique cross section ID
        cross_sections['ID'] = cross_sections.index
        # sets the intersect value to 0
        cross_sections['intsct'] = 0
        # checks if intersections should be caclulated
        if self.parameters.structures_shp is not None:
            print('run intersection analysis')
            # runs the intersection analysis
            structures = self.create_structures()
            # joins intersection information to cross sections
            cross_sections = self.join_intsct_attributes(cross_sections,
                                                         structures)
        # joins fields from the stream shape to cross sections
        cross_sections = self.__join_attributes(cross_sections)
        cross_sections = cross_sections.drop_duplicates(
                                            subset='geometry', keep='first')
        # writes cross sections to Geopackage
        cross_sections.to_file(self.parameters.XLines, driver='GPKG')
        return cross_sections

    def getAngle(self, pt1, pt2):
        '''
        calculates the angle of two points and returns the tangens
        @param x: coordinate
        @type x: x coordinate
        @param y: coordinate
        @type y: y coordinate
        '''
        x_diff = pt2.x - pt1.x
        y_diff = pt2.y - pt1.y
        return math.degrees(math.atan2(y_diff, x_diff))

    def getPoint1(self, pt, bearing, dist):
        '''
        creates the first point of the cross section using the angle
        and distance from the first / second point relation
        @param x: coordinate
        @type x: x coordinate
        @param y: coordinate
        @type y: y coordinate
        @param bearing: angle
        @type bearing: angle of first and second point relation
        '''
        angle = bearing + 90
        bearing = math.radians(angle)
        x = pt.x + dist * math.cos(bearing)
        y = pt.y + dist * math.sin(bearing)
        return Point(x, y)

    def getPoint2(self, pt, bearing, dist):
        '''
        creates the second point of the cross section using the angle
        and distance from the first / second point relation
        @param x: coordinate
        @type x: x coordinate
        @param y: coordinate
        @type y: y coordinate
        @param bearing: angle
        @type bearing: angle of first and second point relation
        '''
        bearing = math.radians(bearing)
        x = pt.x + dist * math.cos(bearing)
        y = pt.y + dist * math.sin(bearing)
        return Point(x, y)

    def get_points(self, line):
        '''
        gets first and second point of a line
        @param line: shapely LineString
        @type line: line segment used for cs creation
        @param firstPoint: shapely point
        @type firstPoint: first point
        @param scndPoint: shapely point
        @type scndPoint: first point
        '''
        firstPoint = list(line.coords)[0]
        scndPoint = list(line.coords)[1]
        return Point(firstPoint), Point(scndPoint)
    
    def create_structures(self):
        '''
        prepares shapefiles for intersections calculation to be joined to
        the cross-section shapefile. Line intersections create points. Points
        are buffered. Polygons are buffered, respectively or used without 
        buffer.
        @param self.intsct_gdf: <GeoDataFrame>
        @type self.intsct_gdf: Contains geometries for intersection
        '''
        print(' - create_structures with %s m buffer' % 
              self.parameters.Buffer_value)
        self.intsct_gdf = gpd.read_file(self.parameters.structures_shp)
        self.intsct_gdf = self.intsct_gdf.to_crs(self.crs)
        geom = self.intsct_gdf.geom_type.unique()[0]
        if geom == 'LineString':
            print('Intersection with Lines')
            points = gpd.GeoDataFrame(crs=self.crs, geometry=list(
                    self.stream_gpd.unary_union.intersection(
                            self.intsct_gdf.unary_union)))
            if self.parameters.Buffer_value is not None:
                points['geometry'] = points.geometry.buffer(
                        self.parameters.Buffer_value)
                self.intsct_gdf = points
            elif self.value is None:
                points['geometry'] = points.geometry.buffer(1)
                self.intsct_gdf = points
        elif geom == 'Polygon':
            print('Intersection with Polygons')
            if self.parameters.Buffer_value is not None:
                self.intsct_gdf['geometry'] = self.intsct_gdf.geometry.buffer(
                        self.parameters.Buffer_value)
        elif geom == 'Point':
            if self.parameters.Buffer_value is not None:
                self.intsct_gdf['geometry'] = self.intsct_gdf.geometry.buffer(
                        self.parameters.Buffer_value)
            elif self.parameters.Buffer_value is None:
                self.intsct_gdf['geometry'] = self.intsct_gdf.geometry.buffer(1)
        self.intsct_gdf.to_file(self.parameters.Intersect_pts, driver='GPKG')
        return self.intsct_gdf

    def tiny_window(self, dtm, x, y):
        '''
        creates a window to read the x, y Z value of the dtm
        @param dtm: raster
        @type dtm: dtm used to grab Z values
        '''
        r, c = dtm.index(x, y)
        return ((r, r+1), (c, c+1))

    def get_columns(self):
        '''
        get columns and create all optional columns if assigned
        '''
        # get number of columns for X
        cols = math.ceil(self.parameters.width/self.parameters.ticks)+1
        # set column names
        col_name = ['ID']+[str(
                'X%s' % i) for i in range(0, cols)]+['intsct', 'geometry']
        # create a list with additional columns to include in the results
        include = []
        if self.parameters.strahler_field is not None:
            col_name = col_name + [self.parameters.strahler_field]
            include = include + [self.parameters.strahler_field]
        if self.parameters.Aggregation_field is not None:
            col_name = col_name + [self.parameters.Aggregation_field]
            include = include + [self.parameters.Aggregation_field]
        if 'art_depth' in self.stream_gpd.dtypes:
            col_name = col_name + ['art_depth']
            include = include + ['art_depth']
        return cols, col_name, include

    def create_points(self):
        '''
        Creates points on crossection lines and adds Z infromation to point.
        Finally creates a csv file with ID, X tick, z, intersect and
        geometry.
        @param dtm: raster
        @type dtm: dtm used to grab Z values
        '''
        print('Creating points with %s m ticks' % self.parameters.ticks)
        # opens the dtm in read mode
        dtm = rio.open(self.parameters.dtm)
        # creates a list of points to fill up
        pts = []
        # creates columns
        cols, col_name, include, = self.get_columns()
        # opens the csv to write the data and columns to
        with open(self.parameters.CrossSections_csv, 'w',
                  newline='') as csv_file:
            writer = csv.writer(csv_file, delimiter=',', quotechar='"',
                                quoting=csv.QUOTE_MINIMAL)
            writer.writerow(col_name)
            # iterrate through the cross section data
            for i, r in self.cross_sections.iterrows():
                # splits the cross section line py point tick
                line_split, points_lst, position = self.split_LS_dist(
                        r.geometry, self.parameters.ticks)
                print('Profile:', r.ID)
                # creates a list to store the z values for profile ID
                z_list = []
                for num in position[0:cols]:
                    # reads the point z value by a tiny window. This is very
                    # useful as large dtm files don't need to be completely
                    # loaded
                    z = dtm.read(window=self.tiny_window(
                            dtm, points_lst[num].x, points_lst[num].y))
                    try:
                        if z.item(0) > 0: 
                            # try to assign z value
                            z_list = z_list + [z.item(0)]
                        else:
                            # if a dtm nodata value of e.g -99 is set
                            z_list = z_list + [np.NaN]
                    except:
                        print('error: out of DTM extent')
                    # creates a list of x and z, intersect
                    var = [i, num, z.item(0), r.intsct, points_lst[num]]
                    # append to points list
                    pts.append(tuple(list(var) + [r[x] for x in include]))
                # write all informations as row to the csv file
                writer.writerow([i] + [s for s in z_list] + [r.intsct] +
                                [str(points_lst[num])] +
                                [r[x] for x in include])
                del line_split, points_lst, position
        # closes the dtm
        dtm.close()
        # closes the csv file
        csv_file.close()
        # creates cols for geopandas geodataframe
        colx = ['ID', 'X', 'z', 'intsct', 'geometry']+[str(x) for x in include]
        # writes the points to a Pandas dataframe
        geopoint = pd.DataFrame(pts, columns=colx)
        # converts the Pandas dataframe to Geopandas dataframe
        geopoint = gpd.GeoDataFrame(geopoint, crs=self.crs,
                                    geometry=geopoint.geometry)
        # Writes the points to Geopackage
        geopoint.to_file(self.parameters.XPoints, driver='GPKG')

    def create_folder(self, directory):
        '''
        creates folder in target directory
        '''
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            print('Error: Creating directory. ' + directory)