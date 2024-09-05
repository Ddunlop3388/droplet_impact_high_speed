#!/bin/bash

#SBATCH -J1Node_cpu-mem=7G
#SBATCH --account=gts-ssuresh313
#SBATCH -N1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=7G
#SBATCH -t360
#SBATCH -qembers
#SBATCH -oReport-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wdunlop3@gatech.edu
cd $SLURM_SUBMIT_DIR

C99='mpicc -std=c99' /storage/home/hcoda1/8/wdunlop3/basilisk/src/qcc -Wall -O2 -D_MPI=1 -grid=octree droplet_impact.c -o droplet_impact -L/storage/home/hcoda1/8/wdunlop3/basilisk/src/gl -lglutils -lfb_tiny -lm


mpirun -np 1 ./droplet_impact 2.28e-3 1.0 8
