# first run:
# import sys
# workdir = "/home/felixm/00_PhD/enki/enki_data/qgis_scripts"
# sys.path.append(workdir)
# import layer_settings as settings

import qgis 
from qgis.core import QgsMapLayerRegistry
from qgis.analysis import QgsZonalStatistics
from qgis.core import QgsVectorLayer
import process_layers as proc

from netCDF4 import Dataset

class viksvatn_layers:
    def __init__(self):
        self.__name__ = "viksvatn"
        self._path_to_polygon = "/home/felixm/00_PhD/enki/enki_data/norway/viksvatn/catchments/subcatch/"
        self._path_to_raster = "/home/felixm/00_PhD/enki/enki_data/norway/general_rasters/landcover/outfiles/"
        self.polygon_layer_files = ["subcatch_1_32N.shp",
                          "subcatch_2_32N.shp",
                          "subcatch_3_32N.shp",
                          "subcatch_4_32N.shp",
                          "subcatch_5_32N.shp",
                          "subcatch_6_32N.shp"]
        self.polygon_centroid_names = ["subcatch_1_centroid_32N",
                          "subcatch_2_centroid_32N",
                          "subcatch_3_centroid_32N",
                          "subcatch_4_centroid_32N",
                          "subcatch_5_centroid_32N",
                          "subcatch_6_centroid_32N"]
        self.polygon_layer_files = [self._path_to_polygon + f for f in self.polygon_layer_files]
        #self.polygon_centroid_files = [self._path_to_polygon + f for f in self.polygon_centroid_files]
        self.raster_layer_files = ["corine_forest_mask_32N.tif",
                        "corine_lake_mask_32N.tif",
                        "corine_glacier_mask_32N.tif",
                        "corine_reservoir_mask_32N.tif"] # reservoirs only zero at the moment
        self.raster_layer_files = [self._path_to_raster + f for f in self.raster_layer_files]
        self.fraction_attribute_names = ["forest_",
                           "lake_",
                           "glacier_",
                           "reservoir_"]
        
        # DEM specific layers
        self.dem_feature_files = ["/home/felixm/00_PhD/enki/enki_data/norway/viksvatn/dem/viksvatn_mosaic_from_saga.tif",
                                  "/home/felixm/00_PhD/enki/enki_data/norway/viksvatn/dem/outfiles/viksvatn_slope.tif",
                                  "/home/felixm/00_PhD/enki/enki_data/norway/viksvatn/dem/outfiles/viksvatn_aspect.tif"]
        self.dem_feature_attribute_names = ["elevation_", "slope_", "aspect_"]
        
        # coordinate system specifications
        self.epsg_code = "EPSG:32632"
        self.grid_mapping_name = "transverse_mercator"
        self.proj4 = "+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

class atnsjoen_layers:
    def __init__(self):
        self.__name__ = "atnsjoen"
        self._path_to_polygon = "/home/felixm/00_PhD/enki/enki_data/norway/atnasjo/catchments/subcatch/"
        self._path_to_raster = "/home/felixm/00_PhD/enki/enki_data/norway/general_rasters/landcover/outfiles/"
        self.polygon_layer_files = ["atnsjoen_subcatch_1_33N.shp"]
        self.polygon_centroid_names = ["atnsjoen_subcatch_1_centroid_33N"]

        self.polygon_layer_files = [self._path_to_polygon + f for f in self.polygon_layer_files]
        #self.polygon_centroid_files = [self._path_to_polygon + f for f in self.polygon_centroid_files]
        self.raster_layer_files = ["corine_forest_mask_33N.tif",
                        "corine_lake_mask_33N.tif",
                        "corine_glacier_mask_33N.tif",
                        "corine_reservoir_mask_33N.tif"] # reservoirs only zero at the moment
        self.raster_layer_files = [self._path_to_raster + f for f in self.raster_layer_files]
        self.fraction_attribute_names = ["forest_",
                           "lake_",
                           "glacier_",
                           "reservoir_"]
        
        # DEM specific layers
        self.dem_feature_files = ["/home/felixm/00_PhD/enki/enki_data/norway/atnasjo/dem/atnsjo_dem.tif",
                                  "/home/felixm/00_PhD/enki/enki_data/norway/atnasjo/dem/outfile/atnsjo_slope.tif",
                                  "/home/felixm/00_PhD/enki/enki_data/norway/atnasjo/dem/outfile/atnsjo_aspect.tif"]
        self.dem_feature_attribute_names = ["elevation_", "slope_", "aspect_"]
        
        # coordinate system specifications
        self.epsg_code = "EPSG:32633"
        self.grid_mapping_name = "transverse_mercator"
        self.proj4 = "+proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

def viksvatn_fractions():
    # calculate fractions for each subcatchment-features
    viksvatn = viksvatn_layer()
    proc.calculate_means(viksvatn.polygon_layer_files, viksvatn.raster_layer_files, viksvatn.attribute_names)

def atnsjoen_fractions():
    # calculate fractions for each subcatchment-features
    atnsjoen = atnsjoen_layer()
    proc.calculate_means(atnsjoen.polygon_layer_files, atnsjoen.raster_layer_files, atnsjoen.attribute_names)
    
def viksvatn_dem_features():
    # calculate DEM-feature (elevation, slope, aspect) for subcatchment features
    viksvatn = viksvatn_layers()
    proc.calculate_means(viksvatn.polygon_layer_files, viksvatn.dem_feature_files, viksvatn.dem_feature_attribute_names)

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


def write_viksvatn_cell_data():
    viksvatn = viksvatn_layers()
    create_cell_data_files(viksvatn)

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
