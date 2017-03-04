#!/bin/sh

# Connect through SSH (-Y : graphic)
# ssh -Y laura.domine@hopper.calcul.labos.polytechnique.fr
# Send a file to hopper
# scp utils.cpp laura.domine@hopper.calcul.labos.polytechnique.fr:~
# Get a file from hopper (-r : recursive)
# scp -rp laura.domine@hopper.calcul.labos.polytechnique.fr:~/displit/histograms/ .
# Launch interactive run
# srun --ntasks=1 ./main 1
# Clear directory
# rm `ls mpi* | grep -v "_FINAL$"`

#SBATCH --job-name=displitGenerate
##SBATCH --time=120  # max 120 minutes
#SBATCH --account=small # 96 coeurs timelimit 10 jours
#SBATCH --mail-user=laura.domine@polytechnique.edu --mail-type=ALL
## Pour reserver 1 coeur
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
##SBATCH --partition=interactive


module add root/6.06.08
module add openmpi/gcc/64/1.8.4

/usr/bin/time -o displitGenerate.log mpiexec ./main generate-mpi parameters
# mpiexec --mca orte_base_help_aggregate 0 -np $SLURM_NTASKS ./main generate