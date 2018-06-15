# Contour advection calculations for Mars polar vortex

## Getting started

First try compiling the contour advection code. In directory `cas`, run `gfortran -o cas cas.f uvw.f`. Now see if it will run: try `./cas < run_params_will.in`. This should produce a file `output.cas`

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
  * `make_contour.py` : makes input contours for CA code. Note the length of the CA calculation is determined by the number in the contour file header (currently set to 50 time steps (i.e. 25 days)).
  * `make_netCDF.py` : produces a netCDF file of q, u, and v from swbob binary output
  * `plot_cas_multicontour.py`: plots results of CA code (as in `../plots/` directory)
  * `plot_pv.py` : produces quick plots of PV for checking output of swbob runs
  * `read_BOB.py` : reads swbob binary output and converts to a (lat,lon) array

* `swbob/` : code for running the shallow water model. Can add more info here if needed.


## swbob run naming convention

A typical filename from swbob output (and other files derived from this output) might look like the following:

`res_test57.-70.-nu4-urlx-kt0.0.c-0020sat200.0.T85`

Here's what it all means:

`res_test` : name of the type of run. I called it this when I was testing different resolutions, should probably update to something more descriptive.
`57.-70.` : inner and outer latitudes of the initial annular vortex
`nu4-urlx` : doesn't change among runs
`kt0.0` : relaxation parameter. In this case, 0 relaxation. A relaxation parameter 4.0 would mean a timescale of 0.25 days (i.e. kt = 1/t_relax).
`c-0020`: angular frequency of the tropography
`sat200.0`: maximum topography amplitude
`T85`: spectral resolution
