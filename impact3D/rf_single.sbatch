#!/bin/bash

#SBATCH -JBurzynskiValidation
#SBATCH --account=gts-ssuresh313
#SBATCH -N1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=7G
#SBATCH -t360
#SBATCH -qembers
#SBATCH -oReport-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wdunlop3@gatech.edu
cd $SLURM_SUBMIT_DIR

uptime_seconds=$(cat /proc/uptime | awk '{print int($1)}')

/storage/home/hcoda1/8/wdunlop3/basilisk/src/qcc -Wall -O2 -grid=octree droplet_impact.c -o droplet_impact_${uptime_seconds} -L/storage/home/hcoda1/8/wdunlop3/basilisk/src/gl -lglutils -lfb_tiny -lm


./droplet_impact_${uptime_seconds} 3.7e-3 10.0 8

rm droplet_impact_${uptime_seconds}
