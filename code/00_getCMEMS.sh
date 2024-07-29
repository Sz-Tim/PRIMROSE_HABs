#!/bin/bash

source ~/miniforge3/bin/activate cmems_env_py3.9
Rscript code/00_getCMEMS.R
source ~/miniforge3/bin/deactivate