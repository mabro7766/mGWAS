#!/bin/bash

echo "Running htu_gwas script ..."

#echo "Current job ID : $SGE_TASK_ID "

#$ -t 2
  # -t 2-4 this is an array job job with tasks from 2 to 4 indicating that the GWAS will be performed on trait in col 2, 3, 4 of Y

module load gridengine
module load R/x86_64/3.2.2

R --max-ppsize=500000 --no-save --args $SGE_TASK_ID < htu_gwas_vesto_v3_noplots.R > htu_gwas_noplots.log
