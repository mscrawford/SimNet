#!/bin/bash

#$ -N PPA
#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -S /bin/bash
#$ -l h_rt=10:00:00
#$ -l h_vmem=20G
#$ -pe openmpi-orte 20

ml load foss/2018b
ml load R/3.6.0

export R_PROFILE=$EBROOTR/lib64/R/library/snow/RMPISNOWprofile
export NSLOTS=20

APP="Rscript --no-save $HOME/PPA.R"
mpirun --display-allocation --verbose -np $NSLOTS $APP

