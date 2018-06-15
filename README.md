# Contour advection calculations for Mars polar vortex

## Directory structure

* `cas` : code for running contour advection (CA)
  * `cas_output/` : output from CA calculations goes here (contents is ignored by git)
  * `input contours/` : initial contours for CA calculations, calculated by `python/make_contour.py`
  * `winds/` : input winds for CA script
  * `cas` : this is the executable for the CA calculation, which can be compiled using gfortran with the following command `gfortran -o cas cas.f uvw.f`
  * `cas.f`	: bulk of the CA code
  * `uvw.f`	: subroutines to initialize wind grids and to read in wind data
  * `cas.h`	: "include" common blocks
