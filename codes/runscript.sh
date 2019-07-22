#!/bin/bash
#SBATCH -n 64               # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-10:10          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p bigmem   # Partition to submit to
#SBATCH --mem=500000          # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o myoutput_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid

python border_check.py