#!/bin/bash

echo "append ID ..."

# module load gridengine
module load R/x86_64/3.2.2

R --max-ppsize=500000 --no-save < append.R > appendR.log
