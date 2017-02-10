#!/bin/bash -l
#SBATCH -p regular   #Submit to the regular 'partition'
#SBATCH -N 1         #Use 1 node
#SBATCH -t 48:00:00  #Set up to 48 hour time limit
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use Haswell nodes
#SBATCH --output=/global/homes/d/desc/jobcontrol/logs/dc1-%j.out

# Set up CD1 environment
cd /global/homes/d/desc/dc1_workflow/DM
source setupStack-dc1.sh
export PILOT_SET_ENV=true 

export CLASSPATH=~desc/jobcontrol/org-srs-jobcontrol-2.1.1-SNAPSHOT-jar-with-dependencies.jar
unset LS_COLORS
export P2_SENDMAIL=/global/homes/b/bvan/bsub/bridge.bash

module load java
java org.srs.jobcontrol.pilot.JobControlPilot -L SCRATCH -C haswell -p dc1-dm "$@"
