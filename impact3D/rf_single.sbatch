#!/bin/bash

#SBATCH -JDropletImpactExample3D
#SBATCH --account=gts-ssuresh313
#SBATCH -N1 --ntasks-per-node=4
#SBATCH --mem-per-cpu=2G
#SBATCH -t360
#SBATCH -qembers
#SBATCH -oReport-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wdunlop3@gatech.edu
cd $SLURM_SUBMIT_DIR

/storage/home/hcoda1/8/wdunlop3/basilisk/src/qcc -Wall -O2 -grid=octree droplet_impact.c -o droplet_impact -L/storage/home/hcoda1/8/wdunlop3/basilisk/src/gl -lglutils -lfb_tiny -lm


./droplet_impact 1e-3 0.5 8
