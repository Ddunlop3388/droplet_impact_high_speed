#!/bin/bash

#SBATCH -JDropImpact
#SBATCH --account=gts-ssuresh313
#SBATCH -N1 --ntasks-per-node=8
#SBATCH --mem-per-cpu=7G
#SBATCH -t360
#SBATCH -qembers
#SBATCH -oReport-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wdunlop3@gatech.edu
cd $SLURM_SUBMIT_DIR

/storage/home/hcoda1/8/wdunlop3/basilisk/src/qcc -Wall -O2 -grid=octree droplet_impact.c -o droplet_impact -L/storage/home/hcoda1/8/wdunlop3/basilisk/src/gl -lglutils -lfb_tiny -lm


./droplet_impact 2.28e-3 1.0 8

rm droplet_impact

/storage/home/hcoda1/8/wdunlop3/basilisk/src/qcc -Wall -O2 -grid=octree droplet_impact.c -o droplet_impact -L/storage/home/hcoda1/8/wdunlop3/basilisk/src/gl -lglutils -lfb_tiny -lm

./droplet_impact 2.28e-3 1.0 9

rm droplet_impact

