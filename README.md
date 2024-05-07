# Antarctic_subglacial_topo

This repository contains the code associated with the manuscript 'Antarctic subglacial topography mapped from space reveals complex mesoscale landscape dynamics'. The code contained within this repository allows the user to apply the Ice Flow Perturbation Analysis methodology detailed in that manuscript to ice-surface elevation and velocity datasets, and invert for the bed elevation and basal slipperiness. It also allows the user to reproduce the figures within that manuscript. 

### file_processing.py and inversion_module_v2.py 
These files contain various functions which apply the Ice Flow Perturbation Analysis methodology over a specified domain.

### mpi_inversion.py 
This script allows the application of the Ice Flow Perturbation Analysis methodology over multiple domains, and uses mpi to process those domains in parallel to reduce the overall processing time. 

### Stitching.py 
This script joins together the multiple outputs produced by running mpi_inversion.py into one single bed topography .nc file. 

### metrics.py 
This script was used to calculate various metrics for the IFPA bed topography, as detailed in the methods section of the manuscript. 

### Antarctic_subglacial_topography.ipynb
This Jupyter notebook contains the code required to reproduce figure 1 - A3 from the manuscript. 

### Figure_A1abc.py
This python script contains the code required to reproduce figure A1 panels a and b from the manuscript. (Provided as a separate script as the datasets to plot the bed topography for the entire continent are large and load slowly in notebook form.). 

### Map_plots.py
This python script contains the code required to reproduce figures A4 - A39 from the manuscript - comparative figures of the new IFPA bed topography, and the bed topography from Bedmachine Antarctica (Morlighem et al. 2020). 

Morlighem, M. et al. Deep glacial troughs and stabilizing ridges unveiled beneath the margins of the Antarctic ice sheet. Nature Geoscience (2020). Dataset available at DOI: 10.5067/FPSU0V1MWUB6. 

