#!/bin/bash
#SBATCH --nodes=1
#SBATCH -p tornado-k40
#SBATCH -t 00-02:00:00
#SBATCH -J Ren_cgm
#SBATSH -o ./result.txt
#SBATSH -e ./error.txt
if [ -f /etc/profile.d/modules-basis.sh ]; 
then
source /etc/profile.d/modules-basis.sh
fi
module purge
module load compiler/gcc/9

gcc shi.c -lpthread -o shi
./shi
