import os
from shutil import copy
from netCDF4 import Dataset

# QGIS imports
from qgis.core import QgsMapLayerRegistry
from qgis.analysis import QgsZonalStatistics
from qgis.core import QgsVectorLayer

class catchment_layers:
    def __init__(self, name, ptp, ptd, ptr, plf, pcn, rlf, dff):
        self.__name__ = name
        self._path_to_polygon = ptp
        self._path_to_raster = ptr
        self.polygon_layer_files = plf  #grid files
        self.polygon_centroid_names = pcn #centroid
        self.polygon_layer_files_tmp = None # for temporary copy of shape files

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

        # cells
        self.nr_cells = self.get_nr_cells()

    def get_nr_cells(self):
        nr_cells = 0
        for poly in self.polygon_layer_files:
            poly_object = QgsVectorLayer(poly, 'zonepolygons', 'ogr')
            nr_cells += len([f for f in poly_object.getFeatures()])
        return nr_cells

    def get_layer(self, name):
        layer = None
        for lyr in QgsMapLayerRegistry.instance().mapLayers().values():
            if lyr.name() == name:
                layer = lyr
                break
        return layer

    def get_attribute_list(self, layer_path, attr_name):
        print("Accessing {}, searching for attribute {}.".format(layer_path,attr_name))
        layer = QgsVectorLayer(layer_path , 'zonepolygons', 'ogr')
        attr_lst = []
        idx = layer.fieldNameIndex(attr_name)
        assert idx >= 0, "Attribute field not know. FieldNameIdx = {}".format(idx)
        for f in layer.getFeatures():
            attr_lst.append(f.attributes()[idx])
        return attr_lst
    
    def get_geometry_list(self, layer_name):
        layer = self.get_layer(layer_name)
        x_lst,y_lst = [],[]
        for f in layer.getFeatures():
            x_lst.append(f.geometry().geometry().x())
            y_lst.append(f.geometry().geometry().y())
        return x_lst, y_lst

    def copy_files(self, path_tmp):
        """
    	to generate a work copy of original files
    	(avoid changes at orig files)

    	path_tmp: path where to copy files
    	"""
        for fn in self.polygon_layer_files:
            copy(fn, path_tmp)
            copy(fn[:-4]+".dbf", path_tmp)
            copy(fn[:-4]+".prj", path_tmp)
            copy(fn[:-4]+".qpj", path_tmp)
            copy(fn[:-4]+".shx", path_tmp)
        self.polygon_layer_files_tmp = [path_tmp +'/'+ f.split('/')[-1] for f in self.polygon_layer_files]

    def _calculate_means(self, polyg_layer_files, raster_layer_files, attribute_names):
    # polygon_layer_files: list of absolute path polygon layer files
    # raster_layer_files: list of absolute path raster layer files
    # attribute_names: name prefix of newly calculated attributes (means from raster to polygon)
        for poly in polyg_layer_files:
            for rast, attr_name in zip(raster_layer_files, attribute_names):
                assert os.path.isfile(poly), "File {} does not exist.".format(poly)
                assert os.path.isfile(rast), "File {} does not exist.".format(rast)
                poly_object = QgsVectorLayer(poly, 'zonepolygons', 'ogr')
                zoneStat = QgsZonalStatistics(poly_object, rast, attr_name, 1, QgsZonalStatistics.Mean)
                check = zoneStat.calculateStatistics(None)
                assert check == 0, "zoneStat.calculateStatistics(None) returned non-zero value... check input layers!"
    
    def calculate_landcover_attributes(self):
        print("Calculating landcover franctions for catchment {}".format(self.__name__))
        self._calculate_means(self.polygon_layer_files_tmp, self.raster_layer_files, self.fraction_attribute_names)
    
    def  calculate_topography_attributes(self):
        print("Calculating topographic attributes for catchment {}".format(self.__name__))
        self._calculate_means(self.polygon_layer_files_tmp, self.dem_feature_files, self.dem_feature_attribute_names)


    def create_cell_data_files(self, outfile):
        """
        catchment: instance of catchment_layers class
        outfile: *.nc file
        """
        with Dataset(outfile, 'w', format='NETCDF4') as dset:
            cell = dset.createDimension('cell', self.nr_cells) # only dimension in nr of cells
            
            # write fractions
            fractions = ["forest-fraction", "reservoir-fraction", "lake-fraction", "glacier-fraction"]
            layer_attributes = ["forest_mean", "reservoir_mean", "lake_mean", "glacier_mean"]
            layer_attributes = [attr[0:10] for attr in layer_attributes] # shapefile attribute names limited to 10 characters
            for frac,layer_attr in zip(fractions, layer_attributes):
                fraction_var = dset.createVariable(frac,'f8',('cell',))
                values = []
                for poly in self.polygon_layer_files_tmp:
                    for value in self.get_attribute_list(poly,layer_attr):
                        assert value>=0 and value <=1, "Value in {} is too large or small: {}".format(layer_attr, value)
                        values.append(value)
                assert len(values) == self.nr_cells, "Number of fractions and cell number not the same!"
                fraction_var[:] = values
                fraction_var.units = '-'
                fraction_var.grid_mapping = 'crs'
                fraction_var.coordinates = 'y x z'

            # write x,y coordinates
            coordinates = ["x","y"]
            x_values = []
            y_values = []
            for poly in self.polygon_centroid_names:
                x_list, y_list = self.get_geometry_list(poly)
                for x,y in zip(x_list,y_list):
                    x_values.append(x)
                    y_values.append(y)
            x_var = dset.createVariable('x','f8',('cell',))
            x_var.units = 'm'
            x_var.axis = 'X'
            x_var.standard_name = 'projection_x_coordinate'
            assert len(x_values) == self.nr_cells, "Number of x-values and cell number not the same!"
            x_var[:] = x_values

            y_var = dset.createVariable('y','f8',('cell',))
            y_var.units = 'm'
            y_var.axis = 'Y'
            y_var.standard_name = 'projection_y_coordinate'
            assert len(y_values) == self.nr_cells, "Number of y-values and cell number not the same!"
            y_var[:] = y_values

            # write elevation
            z_var = dset.createVariable('z','f8',('cell',))
            values = []
            for poly in self.polygon_layer_files_tmp:
                for value in self.get_attribute_list(poly, 'elevation_mean'[0:10]):
                    assert value>=0 and value<9000, "Value in {} is too large or small: {}".format(layer_attr, value)
                    values.append(value)

            assert len(values) == self.nr_cells, "Number of elevations and cell number not the same!"
            z_var[:] = values
            z_var.units = 'm'
            z_var.long_name = "height above mean sea level"
            z_var.axis = 'Z'
            z_var.standard_name = 'height'

            # write area
            area_var = dset.createVariable('area','f8',('cell',))
            values = []
            for poly in self.polygon_layer_files_tmp:
                for value in self.get_attribute_list(poly, 'area'):
                    assert value>=0, "Value in {} is too large or small: {}".format(layer_attr, value)
                    values.append(value)
            assert len(values) == self.nr_cells, "Number of areas and cell number not the same!"
            area_var[:] = values
            area_var.units = 'm^2'
            area_var.grid_mapping = 'crs'
            area_var.coordinates = 'y x z'
            
            # write catchment id
            id_var = dset.createVariable('catchment_id','i4',('cell',))
            values = []
            for poly in self.polygon_layer_files_tmp:
                for value in self.get_attribute_list(poly, 'catch_id'):
                    values.append(value)
            assert len(values) == self.nr_cells, "Number of areas and cell number not the same!"
            id_var[:] = values
            id_var.units = '-'
            id_var.grid_mapping = 'crs'
            id_var.coordinates = 'y x z'
           
            # write crs
            crs_dim = dset.createDimension('crs_dim', None) # only dimension in nr of cells
            crs_var = dset.createVariable('crs','i4',('crs_dim',))
            crs_var.epsg_code = self.epsg_code
            crs_var.grid_mapping_name = self.grid_mapping_name
            crs_var.proj4 = self.proj4 

            # write aspect
            aspect_var = dset.createVariable('aspect','f8',('cell',))
            values = []
            for poly in self.polygon_layer_files_tmp:
                for value in self.get_attribute_list(poly, 'aspect_mean'[0:10]):
                    assert value>=0 and value <= 360, "Value in {} is too large or small: {}".format("aspect", value)
                    values.append(value)
            assert len(values) == self.nr_cells, "Number of areas and cell number not the same!"
            aspect_var[:] = values
            aspect_var.units = 'degree'
            aspect_var.grid_mapping = 'crs'
            aspect_var.coordinates = 'y x z'
            aspect_var.standard_name = "apect angle of area"
            aspect_var.long_name = "apect angle of area; 0 is north; 90 is east"
            
            # write slope
            slope_var = dset.createVariable('slope','f8',('cell',))
            values = []
            for poly in self.polygon_layer_files_tmp:
                for value in self.get_attribute_list(poly, 'slope_mean'[0:10]):
                    assert value>=0 and value <= 90, "Value in {} is too large or small: {}".format("slope", value)
                    values.append(value)
            assert len(values) == self.nr_cells, "Number of areas and cell number not the same!"
            slope_var[:] = values
            slope_var.units = 'degree'
            slope_var.grid_mapping = 'crs'
            slope_var.coordinates = 'y x z'
            slope_var.standard_name = "slope angle of area"
