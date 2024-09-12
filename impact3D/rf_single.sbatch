#!/bin/bash

#SBATCH -JMoitaValidation
#SBATCH --account=gts-ssuresh313
#SBATCH -N1 --ntasks-per-node=8
#SBATCH --mem-per-cpu=16GB
#SBATCH -t360
#SBATCH -qinferno
#SBATCH -oReport-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wdunlop3@gatech.edu
cd $SLURM_SUBMIT_DIR

uptime_seconds=$(cat /proc/uptime | awk '{print int($1)}')
numtasks=8

C99='mpicc -std99' /storage/home/hcoda1/8/wdunlop3/basilisk/src/qcc -Wall -O2 -grid=octree -D_MPI=$numtasks droplet_impact.c -o droplet_impact_${uptime_seconds} -L/storage/home/hcoda1/8/wdunlop3/basilisk/src/gl -lglutils -lfb_tiny -lm -L/usr/local/pace-apps/manual/packages/mvapich/2.3.6-ouywal -lmpi


srun --ntasks=$numtasks ./droplet_impact_${uptime_seconds} 3.214e-3 3.3225 8

rm droplet_impact_${uptime_seconds}
