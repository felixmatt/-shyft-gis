# first run
# import sys
# workdir = "/home/felixm/00_PhD/enki/enki_data/qgis_scripts"
# sys.path.append(workdir)
# import process_layers as proc

import os

import qgis
from qgis.core import QgsMapLayerRegistry
from qgis.analysis import QgsZonalStatistics
from qgis.core import QgsVectorLayer
import layer_settings
def refresh_layers():
    layers = qgis.utils.iface.legendInterface().layers()
    for layer in layers:
        layer.setCacheImage(None)
        layer.triggerRepaint()

def get_layer(name):
    layer = None
    for lyr in QgsMapLayerRegistry.instance().mapLayers().values():
        if lyr.name() == name:
            layer = lyr
            break
    return layer

def get_attribute_list(layer_path, attr_name):
    layer = QgsVectorLayer(layer_path, 'zonepolygons', 'ogr')
    attr_lst = []
    idx = layer.fieldNameIndex(attr_name)
    assert idx >= 0, "Attribute field not know. FieldNameIdx = {}".format(idx)
    for f in layer.getFeatures():
        attr_lst.append(f.attributes()[idx])
    return attr_lst

def get_geometry_list(layer_name):
    layer = get_layer(layer_name)
    x_lst,y_lst = [],[]
    for f in layer.getFeatures():
        x_lst.append(f.geometry().geometry().x())
        y_lst.append(f.geometry().geometry().y())
    return x_lst, y_lst


def calculate_means(polygon_layer_files, raster_layer_files, attribute_names):
    # polygon_layer_files: list of absolute path polygon layer files
    # raster_layer_files: list of absolute path raster layer files
    # attribute_names: name prefix of newly calculated attributes (means from raster to polygon)
    for poly in polygon_layer_files:
        for rast, attr_name in zip(raster_layer_files, attribute_names):
            assert os.path.isfile(poly), "File {} does not exist.".format(poly)
            assert os.path.isfile(rast), "File {} does not exist.".format(rast)
            poly_object = QgsVectorLayer(poly, 'zonepolygons', 'ogr')
            zoneStat = QgsZonalStatistics(poly_object, rast, attr_name, 1, QgsZonalStatistics.Mean)
            check = zoneStat.calculateStatistics(None)
            assert check == 0, "zoneStat.calculateStatistics(None) returned non-zero value... check input layers!"
    refresh_layers() # TODO: refresh does not work somehow... in qgis: reload layer and delete old one as workaround
