#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-40:00
#SBATCH -p pleiades
#SBATCH --mem=8000
#SBATCH -o contour_%j.out
#SBATCH -e contour_%j.err

source run_contourplot.sh