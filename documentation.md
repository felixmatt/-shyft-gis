#Generatecell-data.nc for shyft

##DEM

- Load tifs

- merge them to one(default settings) Raster-Misselaneous-Merge

- clip to your area(not necessary but to reduce data volume)

- calulate slope andaspect from your dem (Raster-Analysis-DEM)

  - for slope chooseMode ‘Slope’
  - for aspect chooseMode “Aspect’

  when calculatingaspect tick “Return 0 for flat (instead of -9999)”

- save the output

10m DEM from Norwayis available here: <https://kartkatalog.geonorge.no/metadata/kartverket/dtm-10-terrengmodell-utm33/dddbb667-1303-4ac5-8640-7ec04c0e3918>



##CORINE Land Cover

- load corine layer

- 67800, 105500,6697000, 6749000

- for forest fraction

  ```
  ("corine_clip@1"=23)/23 OR ("corine_clip@1"=24)/24OR ("corine_clip@1"=25)/25 OR (("corine_clip@1"<23)AND ("corine_clip@1">25))*0
  ```

- for lake fraction

  ```
  ("corine_clip@1"=41)/41 OR (("corine_clip@1"<41)AND ("corine_clip@1">41))*0
  ```

- for glacier fraction

  ```
  ("corine_clip@1"=34)/34 OR (("corine_clip@1"<34)AND ("corine_clip@1">34))*0
  ```

- for reservoir fraction

  ```
  "corine_clip@1"*0
  ```

- set UTM33(EPSG:32633)

- set min to 0 and maxto 1 to see the correct results!

The complete CorineLand Cover map is available here:<https://land.copernicus.eu/pan-european/corine-land-cover>



##CATCHMENTS

Tools to derive catchment boundaries from outlet positions: <http://nevina.nve.no/>

- load originalcatchments
- by using the toolsVector-Geoprocessing-Difference and Clip, avoid overlapping catchment areas
  - result for Finse:
- create grid layer that covers the extension of your study area (Vector-ResearchTools-Vector Grid)
  - Set X and Y spacingto 1000 as we want a 1x1km grid cell size.
- Save the layer
- use toVector-Geoprocessing Tools-Clip tool to get the grid on your catchments. 
  - Input layer is the grid, clip layer is the catchment layer
- save the layer



- add fields(attribute table – edit mode – add field)

  - fields:

    - catch_id (int,Length 5)

    - area (double,Length 10, precision 2)

      to calculate area,select field, type ‘$area’ and click ‘Update all’.

- Calculate centroids for gridded catchments. The centroid of every polygon…
  (Vector - Geometry Tools - Polygon Vector)

- Save the output

  ​

##PYTHON CONSOLE

- create ‘/int’ directory for temp-files
- Open Python Console(Plugins – Python Console)
- load the file ‘layersettings.py’

### Explanations to‘layer_settings.py’

- *workdir*: path to process_layer.py
- *ptp, ptr, ptd*: here all the same directory ('project')
  COMMENT: Could be cleaned up a bit...
- *plf*: list of all your shape files representing your gridded shapefiles
- *pcn*: list of the centroids layers as they are named in qgis
- *rlf*: list of your land cover fraction raster layer (lake, forest, glacier, reservoir)
- *dff*: list of your topography characteristic raster layers (dem, slope, aspect)
- *DATA_PATH*: path to the 'project' directory
- *DATA_INT*: path to the directory for the intermediate files
- *OUT*: path to the' directory where the output should be saved



