#!/usr/bin/env bash
# LSST Jupyter-dev kernel for SSim DC1 notebook examples.
# For installation and setup instructions, see:
# https://github.com/LSSTDESC/Monitor/blob/master/doc/jupyter-dev.md
#
# Remember that this file needs to be executable:
# chmod ug+rx lsst-kernel.sh

INST_DIR=/global/common/software/lsst/cori-haswell-gcc/stack/
source $INST_DIR/setup_w_2017_46_py3_gcc6.sh
setup lsst_sims
setup afw

exec python -m ipykernel $@
