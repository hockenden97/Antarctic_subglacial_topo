# Antarctic_subglacial_topo

This repository contains the code associated with the manuscript 'Complex mesoscale landscapes beneath Antarctica mapped from space'. The code contained within this repository allows the user to apply the Ice Flow Perturbation Analysis methodology detailed in that manuscript to ice-surface elevation and velocity datasets, and invert for the bed elevation and basal slipperiness. It also allows the user to reproduce the figures within that manuscript. 

### REMA_on_its.py
To apply the IFPA method, the ice surface data and velocity data need to be on the same coordinate grid. This script can be used to interpolate the REMA ice surface data (100 m resolution) onto the velocity grid (450 m resolution for the ITSLIVE/MEaSURES grid.

### file_processing.py and inversion_module_v3.py (and inversion module_v2.py) 
These files contain various functions which apply the Ice Flow Perturbation Analysis methodology over a specified domain. Inversion_module_v3.py shouldn't produce error messages due to dividing by 0 (unlike inversion_module_v2.py), but the map in the manuscript was produced with inversion_module_v2.py. 

### mpi_inversion.py 
This script allows the application of the Ice Flow Perturbation Analysis methodology over multiple domains, and uses mpi to process those domains in parallel to reduce the overall processing time. 

### Stitching.py 
This script joins together the multiple outputs produced by running mpi_inversion.py into one single bed topography .nc file. 

### metrics.py and Antarctic_IFPA_Metrics.ipynb and Antarctic_IFPA_Metrics_bedmap3.ipynb
This script was used to calculate various metrics for the IFPA bed topography, as detailed in the methods section of the manuscript. The notebook Antarctic_IFPA_Metrics.ipynb can also be used for this. 

### Antarctic_FIGURES.ipynb
This Jupyter notebook contains the code required to reproduce all the figures in the manuscript. This notebook can either be used with the datasets provided for these figures, or to reproduce these datasets from the original data, depending on user preference. Currently, the full IFPA map has not been provided in the Zenodo repository which accompanies this notebook as it is under embargo.

### Antarctic_IFPA_worked_example.ipynb
This notebook works through several examples of applying the IFPA code to regions in Antarctica, explaining the steps in the IFPA code.  

