# shyft-gis

This repository provides py-scripts and documentation for creating a Shyft cell data file 
(netcdf format) that can be read using the 
[CFRegionModelRepository](https://github.com/statkraft/shyft/blob/master/shyft/repository/netcdf/cf_region_model_repository.py) 
repository. The scripts rely on a working Qgis 2.x install (Qgis 3.x not tested).

An exmaple is proveded in finse_qgis.tar.gz

To run the example (tested for Ubuntu 14.04 LTS and Qgis 2.18):

- clone the repository to $HOME (`cd $HOME; git clone https://github.com/felixmatt/shyft-gis.git`)
- untar finse_qgis.tar.gz (e.g. `tar -xvzf finse_qgis.tar.gz` in Ubuntu)
- open the Qgis project `project.qgs` located in directory `finse_qgis` using Qgis
- open the Qgis python console (`Plugins` -> `Python Console`)
- open the `layer_settings.py` file using the editor provided in the Qgis python console
- run the script in the python console

If the script runs succesfully, a `cell_data.nc` is created in `$HOME/shyft-gis`.

In order to setup other regions, the Qgis project in combination with the [layer_settings.py](https://github.com/felixmatt/shyft-gis/blob/master/layer_settings.py)
gives an idea on which maps are required to run the scripts.
