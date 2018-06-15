# Contour advection calculations for Mars polar vortex

## Directory structure

* `cas/` : code for running contour advection (CA)
  * `cas_output/` : output from CA calculations goes here (contents is ignored by git)
  * `input contours/` : initial contours for CA calculations, calculated by `python/make_contour.py`
  * `winds/` : input winds for CA script
  * `cas` : this is the executable for the CA calculation, which can be compiled using gfortran with the following command `gfortran -o cas cas.f uvw.f`. CA code is comprised of the following (see `cas/README` for more details):
    * `cas.f`	: bulk of the CA code
    * `uvw.f`	: subroutines to initialize wind grids and to read in wind data
    * `cas.h`	: "include" common blocks
  * `output.cas` : 'default' output file for CA output, then copied to a more descriptive filename in `cas_output/`
  * `output.log` : logfile for CA output
  * `pv.in` : 'default' input contour file (copied over from `input_contours/`)
  * `run_multi_contour.py` : This script runs the CA code over a series of contours from `input_contours/`, and using winds from `../model_output/`, then puts the output in `cas_output/`. It is currently a python script but could be easily changed to a shell script.
  * `run_params_will.in` : input for the `cas` executable, which is run with the following command: `./cas < run_params_will.in`. The only thing you need to change here is the number on the final line which is the start day for the calculation.

* `model_output/` : Binary model output, netCDF files made from model output, and wind files for input into CA code. Too large for github, but can dowload here: https://www.dropbox.com/s/1lhhv52l1cc56z6/model_output.zip?dl=0 (about 2GB).

* `plots/` : currently shows results of CA calculations for two runs with same wave 2 topography but one with zero relaxation and one with relaxation giving an annular vortex.

* `python/` : python scripts for making inputs for CA code and plotting outputs
  * `interpolate_winds.py` : produces the wind input files for the CA code (found in `../model_output/winds_*`)
  * `make_contour.py` : makes input contours for CA code
  * `make_netCDF.py` : produces a netCDF file of q, u, and v from swbob binary output
  * `plot_cas_multicontour.py`: plots results of CA code (as in `../plots/` directory)
  * `plot_pv.py` : produces quick plots of PV for checking output of swbob runs
  * `read_BOB.py` : reads swbob binary output and converts to a (lat,lon) array

* `swbob/` : code for running the shallow water model. Can add more info here if needed.


## swbob run naming convention

A typical file
