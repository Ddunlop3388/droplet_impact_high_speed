#!/bin/bash

#SBATCH -JBursynskiValidation
#SBATCH --account=gts-ssuresh313
#SBATCH -N1 --ntasks-per-node=10
#SBATCH --mem-per-cpu=16GB
#SBATCH -t360
#SBATCH -qinferno
#SBATCH -oReport-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wdunlop3@gatech.edu
cd $SLURM_SUBMIT_DIR

uptime_seconds=$(cat /proc/uptime | awk '{print int($1)}')
numtasks=10

C99='mpicc -std99' /storage/home/hcoda1/8/wdunlop3/basilisk/src/qcc -Wall -O2 -grid=octree -D_MPI=$numtasks droplet_impact_clsvof.c -o droplet_impact_${uptime_seconds} -L/storage/home/hcoda1/8/wdunlop3/basilisk/src/gl -lglutils -lfb_tiny -lm -L/usr/local/pace-apps/manual/packages/mvapich/2.3.6-ouywal -lmpi


srun --ntasks=$numtasks ./droplet_impact_${uptime_seconds} 2.28e-3 1.0 9

rm droplet_impact_${uptime_seconds}
