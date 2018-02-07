#import sys
#workdir = "/home/sven/Documents/github/shyft-gis"
#sys.path.append(workdir)
#import layer_settings as settings

import os
import qgis 
from qgis.core import QgsMapLayerRegistry
from qgis.analysis import QgsZonalStatistics
from qgis.core import QgsVectorLayer
import process_layers as proc
from shutil import copy

from netCDF4 import Dataset






class catchment_layers:
    def __init__(self, name, ptp, ptr, plf, pcn, rlf, dff):
        self.__name__ = name
        self._path_to_polygon = ptp
        self._path_to_raster = ptr
        self.polygon_layer_files = plf  #grid files
        self.polygon_centroid_names = pcn #centroid

        self.polygon_layer_files = [self._path_to_polygon + f for f in self.polygon_layer_files]
        #self.polygon_centroid_files = [self._path_to_polygon + f for f in self.polygon_centroid_files]
        self.raster_layer_files = rlf # reservoirs only zero at the moment
        self.raster_layer_files = [self._path_to_raster + f for f in self.raster_layer_files]
        self.fraction_attribute_names = ["forest_",
                           "lake_",
                           "glacier_",
                           "reservoir_"]
        
        # DEM specific layers
        self._path_to_dem = ptd
        self.dem_feature_files = dff
        self.dem_feature_files = [self._path_to_dem + f for f in self.dem_feature_files]
        self.dem_feature_attribute_names = ["elevation_", "slope_", "aspect_"]
        
        # coordinate system specifications
        self.epsg_code = "EPSG:32633"
        self.grid_mapping_name = "transverse_mercator"
        self.proj4 = "+proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


    def copy_files(self, path_interm):
    	'''
    	to generate a work copy of original files
    	(avoid changes at orig files)

    	path_interm: path where to copy files
    	'''
    	for fn in self.polygon_layer_files:
    		copy(fn, path_interm)
            copy(fn[:-4]+".dbf", path_interm)
            copy(fn[:-4]+".prj", path_interm)
            copy(fn[:-4]+".qpj", path_interm) 
            copy(fn[:-4]+".shx", path_interm)

    	self.polygon_layer_files = [DATA_INT + f for f in self.polygon_layer_files]


    def calculate_means(self):
    # polygon_layer_files: list of absolute path polygon layer files
    # raster_layer_files: list of absolute path raster layer files
    # attribute_names: name prefix of newly calculated attributes (means from raster to polygon)
    
    #for landcover fractions
    for poly in self.polygon_layer_files:
        for rast, attr_name in zip(self.raster_layer_files, self.fraction_attribute_names):
            assert os.path.isfile(poly), "File {} does not exist.".format(poly)
            assert os.path.isfile(rast), "File {} does not exist.".format(rast)
            poly_object = QgsVectorLayer(poly, 'zonepolygons', 'ogr')
            zoneStat = QgsZonalStatistics(poly_object, rast, attr_name, 1, QgsZonalStatistics.Mean)
            check = zoneStat.calculateStatistics(None)
            assert check == 0, "zoneStat.calculateStatistics(None) returned non-zero value... check input layers!"

    #for dem characteristics (elev, slope, aspect)
    for poly in self.polygon_layer_files:
        for rast, attr_name in zip(self.raster_layer_files, self.dem_feature_attribute_names):
            assert os.path.isfile(poly), "File {} does not exist.".format(poly)
            assert os.path.isfile(rast), "File {} does not exist.".format(rast)
            poly_object = QgsVectorLayer(poly, 'zonepolygons', 'ogr')
            zoneStat = QgsZonalStatistics(poly_object, rast, attr_name, 1, QgsZonalStatistics.Mean)
            check = zoneStat.calculateStatistics(None)
            assert check == 0, "zoneStat.calculateStatistics(None) returned non-zero value... check input layers!"

    refresh_layers() # TODO: refresh does not work somehow... in qgis: reload layer and delete old one as workaround




def setup_finse():
    ptp = DATA_PATH+'shyft_catchments/' # path to polygon
	ptr = DATA_PATH+'shyft_catchments/' # path to raster
	ptd = DATA_PATH+'dem/'

	plf = ["12.166_shyft_gridded/12.166_gridded.shp", 
				"50.2_shyft_gridded/50.2_gridded.shp", 
				"50.5_shyft_gridded/50.5_gridded.shp", 
				"50.8_shyft_gridded/50.8_gridded.shp", 
				"50.10_shyft_gridded/50.10_shyft_gridded.shp", 
				"50.11_shyft_gridded/50.11_shyft_gridded.shp", 
				"50.13_shyft_gridded/50.13_gridded.shp", 
				"50.38_shyft_gridded/50.38_gridded.shp"] # polygon layer files

	pcn = ["12.166_centroids", "50.2_centroids", "50.5_centroids", "50.8_centroids", "50.10_centroids", "50.11_centroids", "50.13_centroids", "50.38_centroids"]

	rlf = ["g100_clc12_V18_5a/corine_forest_fraction.tif", 
                "g100_clc12_V18_5a/corine_lake_fraction.tif",
                "g100_clc12_V18_5a/corine_glacier_fraction.tif",
                "g100_clc12_V18_5a/corine_reservoir_fraction.tif"] # raster layer files
	dff = ["dem_finse.tif", "slope_finse.tif", "aspect_finse.tif"] # dem feature files
	finse = catchment_layers("finse", ptp, ptr, plf, pcn, rlf, dff)
	return finse




def atnsjoen_fractions():
    # calculate fractions for each subcatchment-features
    atnsjoen = atnsjoen_layer()
    proc.calculate_means(atnsjoen.polygon_layer_files, atnsjoen.raster_layer_files, atnsjoen.attribute_names)
    
def atnsjoen_dem_features():
    # calculate DEM-feature (elevation, slope, aspect) for subcatchment features
    atnsjoen = atnsjoen_layers()
    proc.calculate_means(atnsjoen.polygon_layer_files, atnsjoen.dem_feature_files, atnsojoen.dem_feature_attribute_names)

def get_nr_cells(catchment):
    nr_cells = 0
    for poly in catchment.polygon_layer_files:
        poly_object = QgsVectorLayer(poly, 'zonepolygons', 'ogr')
        nr_cells += len([f for f in poly_object.getFeatures()])
    return nr_cells

def write_atnsjoen_cell_data():
    atnsjoen = atnsjoen_layers()
    create_cell_data_files(atnsjoen)










def create_cell_data_files(catchment, nc_path='/home/felixm/00_PhD/enki/enki_data/qgis_scripts/cell_data'):
    nr_cells = get_nr_cells(catchment)

    with Dataset('/'.join([nc_path,catchment.__name__+'_cell_data.nc']), 'w', format='NETCDF4') as dset:
        cell = dset.createDimension('cell', nr_cells) # only dimension in nr of cells
        
        # write fractions
        fractions = ["forest-fraction", "reservoir-fraction", "lake-fraction", "glacier-fraction"]
        layer_attributes = ["forest_mean", "reservoir_mean", "lake_mean", "glacier_mean"]
        layer_attributes = [attr[0:10] for attr in layer_attributes] # shapefile attribute names limited to 10 characters
        for frac,layer_attr in zip(fractions, layer_attributes):
            fraction_var = dset.createVariable(frac,'f8',('cell',))
            values = []
            for poly in catchment.polygon_layer_files:
                for value in proc.get_attribute_list(poly,layer_attr):
                    assert value>=0 and value <=1, "Value in {} is too large or small: {}".format(layer_attr, value)
                    values.append(value)
            assert len(values) == nr_cells, "Number of fractions and cell number not the same!"
            fraction_var[:] = values
            fraction_var.units = '-'
            fraction_var.grid_mapping = 'crs'
            fraction_var.coordinates = 'y x z'

        # write x,y coordinates
        coordinates = ["x","y"]
        x_values = []
        y_values = []
        for poly in catchment.polygon_centroid_names:
            proc.get_geometry_list(poly)
            x_list, y_list = proc.get_geometry_list(poly)
            for x,y in zip(x_list,y_list):
                x_values.append(x)
                y_values.append(y)
        x_var = dset.createVariable('x','f8',('cell',))
        x_var.units = 'm'
        x_var.axis = 'X'
        x_var.standard_name = 'projection_x_coordinate'
        assert len(x_values) == nr_cells, "Number of x-values and cell number not the same!"
        x_var[:] = x_values

        y_var = dset.createVariable('y','f8',('cell',))
        y_var.units = 'm'
        y_var.axis = 'Y'
        y_var.standard_name = 'projection_y_coordinate'
        assert len(y_values) == nr_cells, "Number of y-values and cell number not the same!"
        y_var[:] = y_values

        # write elevation
        z_var = dset.createVariable('z','f8',('cell',))
        values = []
        for poly in catchment.polygon_layer_files:
            for value in proc.get_attribute_list(poly, 'elevation_mean'[0:10]):
                assert value>=0 and value<9000, "Value in {} is too large or small: {}".format(layer_attr, value)
                values.append(value)

        assert len(values) == nr_cells, "Number of elevations and cell number not the same!"
        z_var[:] = values
        z_var.units = 'm'
        z_var.long_name = "height above mean sea level"
        z_var.axis = 'Z'
        z_var.standard_name = 'height'

        # write area
        area_var = dset.createVariable('area','f8',('cell',))
        values = []
        for poly in catchment.polygon_layer_files:
            for value in proc.get_attribute_list(poly, 'area'):
                assert value>=0, "Value in {} is too large or small: {}".format(layer_attr, value)
                values.append(value)
        assert len(values) == nr_cells, "Number of areas and cell number not the same!"
        area_var[:] = values
        area_var.units = 'm^2'
        area_var.grid_mapping = 'crs'
        area_var.coordinates = 'y x z'
        
        # write catchment id
        id_var = dset.createVariable('catchment_id','i4',('cell',))
        values = []
        for poly in catchment.polygon_layer_files:
            for value in proc.get_attribute_list(poly, 'catch_id'):
                values.append(value)
        assert len(values) == nr_cells, "Number of areas and cell number not the same!"
        id_var[:] = values
        id_var.units = '-'
        id_var.grid_mapping = 'crs'
        id_var.coordinates = 'y x z'
       
        # write crs
        crs_dim = dset.createDimension('crs_dim', None) # only dimension in nr of cells
        crs_var = dset.createVariable('crs','i4',('crs_dim',))
        crs_var.epsg_code = catchment.epsg_code
        crs_var.grid_mapping_name = catchment.grid_mapping_name
        crs_var.proj4 = catchment.proj4 

        # write aspect
        aspect_var = dset.createVariable('aspect','f8',('cell',))
        values = []
        for poly in catchment.polygon_layer_files:
            for value in proc.get_attribute_list(poly, 'aspect_mean'[0:10]):
                assert value>=0 and value <= 360, "Value in {} is too large or small: {}".format("aspect", value)
                values.append(value)
        assert len(values) == nr_cells, "Number of areas and cell number not the same!"
        aspect_var[:] = values
        aspect_var.units = 'degree'
        aspect_var.grid_mapping = 'crs'
        aspect_var.coordinates = 'y x z'
        aspect_var.standard_name = "apect angle of area"
        aspect_var.long_name = "apect angle of area; 0 is north; 90 is east"
        
        # write slope
        slope_var = dset.createVariable('slope','f8',('cell',))
        values = []
        for poly in catchment.polygon_layer_files:
            for value in proc.get_attribute_list(poly, 'slope_mean'[0:10]):
                assert value>=0 and value <= 90, "Value in {} is too large or small: {}".format("slope", value)
                values.append(value)
        assert len(values) == nr_cells, "Number of areas and cell number not the same!"
        slope_var[:] = values
        slope_var.units = 'degree'
        slope_var.grid_mapping = 'crs'
        slope_var.coordinates = 'y x z'
        slope_var.standard_name = "slope angle of area"


if __name__ == "__main__":

    HOME= os.environ['HOME']
    DATA_PATH= HOME+'/Documents/finse/GIS_data/'
    DATA_INT=DATA_PATH+'shyft_catchments/int'

	finse = setup_finse()
    finse.copy_files(DATA_INT) # create intermediate files and point them to polygon_layer_files
    finse.calculate_means()