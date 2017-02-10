#!/bin/bash

ulimit -c ${CORE_LIMIT:-1000} # Limit core dump
set -e # exit on error

# Workaround for EUPS trying to write to home directory
export HOME=`pwd`

# Workaround for low level libraries such as OpenBLAS allocating many threads
export OMP_NUM_THREADS=1 

# Find script to run
export SCRIPT=${SCRIPT_LOCATION}/${PIPELINE_PROCESS:-$1}

# Set up Twinkles environment and invoke process specific script
# If DM already set up then no need to do it again
if [ -v PILOT_SET_ENV ]
then
set -xe; export SHELLOPTS; source ${SCRIPT} 
else
source ${DM_DIR}/${DM_SETUP}; set -xe; export SHELLOPTS; source ${SCRIPT}
fi
