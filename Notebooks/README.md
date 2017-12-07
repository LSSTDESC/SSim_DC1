## DC1 Data Access Notebooks

These notebooks contain demonstration material convering how to access and use the DC1 datasets. 
To run these notebooks at NERSC, follow the instructions for setting up [Jupyter-dev](https://github.com/LSSTDESC/Monitor/blob/master/doc/jupyter-dev.md). 
You will need a more modern version (as of 2017-12-04) of a Python3 kernel, which we provide through [lsst_kernel.sh](../scripts/lsst-kernel.sh).

These notebooks include instructions for accessing the DC1 data directly using the Butler and indirectly through HDF5 files extracted from the database and read with [Pandas](https://pandas.pydata.org/).

The ["Dask" subdirectory](./Dask) contains examples of how to use [Dask](https://dask.pydata.org/en/latest/) to do out-of-memory operations on large HDF5 files and build smaller Pandas dataframes that you can use in memory.  You can also download the HDF5 files to your own machine.

While we've done our best to validate that these notebooks work on a single system (NERSC) at a specific point in time (December 2017), they may not "just work" for you without modification. They are written by various group members and downloaded from SSim meeting presentations. For this reason, you may sometimes have to modify the expected input file locations or obtain the needed input files etc.  However, the information in these files will give you the background you need to write your own analysis and validation programs. Several of the oldest and most obsolete notebooks have been archived in the ["Deprecated" subdirectory](./Deprecated).

## Description

* Dask
** [DC1_Dask_Access.ipynb](./Dask/DC1_Dask_Access.ipynb) - Introduction to indirect data access from HDF5 files using Dask.
** [DC1_Dask_Validation.ipynb](./Dask/DC1_Dask_Validation.ipynb) - Some data validation examples using the Dask access interface.
* [DC1_Butler_Access.ipynb](./DC1_Butler_Access.ipynb) - Introduction to direct Butler access.
* [DESC-SSim_Flags.ipynb](./DESC-SSim_Flags.ipynb) - "Best practices" for selecting objects using flags
* [DESC-SSim_Patch_Geometry.ipynb](./DESC-SSim_Patch_Geometry.ipynb) - Survey geometry plotting examples using SkyMap
* [PSF_test.ipynb](./PSF_test.ipynb) - Some simple PSF testing examples.