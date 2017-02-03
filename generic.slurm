#!/bin/sh

# Connect through SSH
# ssh -Y laura.domine@hopper.calcul.labos.polytechnique.fr
# Send a file to hopper
# scp utils.cpp laura.domine@hopper.calcul.labos.polytechnique.fr:~
# Get a file from hopper
# scp laura.domine@hopper.calcul.labos.polytechnique.fr:~/displit/tree.root .
# Launch interactive run
# srun --ntasks=1 ./main 1

#SBATCH --job-name=test
#SBATCH --time=120  # max 120 minutes
#SBATCH --account=medium # 96 coeurs timelimit 10 jours
#SBATCH --mail-user=laura.domine@polytechnique.edu --mail-type=ALL
## Pour reserver 1 coeur
#SBATCH --ntasks=1
#SBATCH --partition=interactive   // On force l'utilisation de la queue interactive.
 
 
module add root
 
./main 1