#!/bin/bash

#SBATCH -JDropletImpactExample
#SBATCH --account=gts-ssuresh313
#SBATCH -N1 --ntasks-per-node=4
#SBATCH --mem-per-cpu=2G
#SBATCH -t1440
#SBATCH -qinferno
#SBATCH -oReport-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wdunlop3@gatech.edu
cd $SLURM_SUBMIT_DIR

/storage/home/hcoda1/8/wdunlop3/basilisk/src/qcc droplet_impact.c -o droplet_impact -L/storage/home/hcoda1/8/wdunlop3/basilisk/src/gl -lglutils -lfb_tiny -lm 

u0=20

while [ $u0 -le 250 ]
do
        echo $u0

        dia_number=1

        drop_dia=$(awk "BEGIN {print $dia_number * 10^-6}")
        ./droplet_impact $drop_dia $u0
        echo $drop_dia
        dia_number=10

        while [ $dia_number -le 100 ]
        do
                drop_dia=$(awk "BEGIN {print $dia_number * 10^-6}")

                echo $drop_dia
                ./droplet_impact $drop_dia $u0
                dia_number=$((dia_number+10))
        done
        u0=$((u0+20))
done

./droplet_impact 1e-4 1.0
