With the code and scripts in this directory one can reproduce the numerical experiments and plots discussed in the publication:

Riha, S., Peliz, A. (2013): A 2-Layer Primitive Equation Model of an Idealized Strait of Gibraltar Connected to an Eastern Basin. Ocean Dyanamics.

The configuration files for the Hallberg Isopycnal Model (HIM) are in the 'run*' directories. You need to download HIM before you can run the experiments. Search for the file 'HIM1.10_beta2_code.tar.gz' on the web. You will also need the following python-packages: numpy, scipy, matplotlib and netCDF4.

Runs 148-151 are the BASIC expermients.  
runs 178-181: DENSE  
runs 182-185: LIGHT  
runs 170-173: CAPE  
The other runs are for plotting the spinup phase (same as above but more frequent output).


The scripts for plotting the figures and printing the tables are located in the 'postprocessing' folder. To plot Fig. 1, you also have to get the file 'ETOPO1_Ice_g_gdal.grd.gz', unpack it and move it to the 'postprocessing' folder. The file can be downloaded here: http://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gdal.grd.gz

Apparently git cannot store empty folders. You have to create these yourself:  
>mkdir runXYZ/saves  
>mkdir runXYZ/data  
>mkdir postprocess/figures  
>mkdir postprocess/tables_tex  

Before compiling HIM, you have to edit the runXZY/Makefiles to suit your computing environment and:  
*) point SRCDIR to the HIM source code  
*) change INPUTDIR in runXZY/init.h  

'>python prep_run.py' will build the input files for HIM.  
'>./make.sh' compiles HIM.  
'>./run.sh' runs the experiments.  


Note that the the instructions for compiling and running the model assume that you use 'mpirun' (openmpi). If this is not the case, then according changes have to be made in the 'run*/init.h', 'runXZY/Makefiles' and './run.sh' files. We used a single multi-core AMD Opteron Processor 3280, and one 3-year experiment (i.e. one of the experiments in the runXZY directories) took about half a day to complete.


To make the figures after all expermients are completed, move into the 'postprocessing' directory and type '>python fig2.py' etc. The .pdfs are written to the 'figures' directory. There is also a script for printing out the numbers of the tables in the publication in a latex booktabs environment.


Stefan Riha   hoitaus@gmail.com

[![githalytics.com alpha](https://cruel-carlota.gopagoda.com/9da48a37900dc43d54671cb680f0faa8 "githalytics.com")](http://githalytics.com/poidl/article_gib1)
