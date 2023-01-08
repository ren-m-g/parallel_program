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
module add python/3.9
module add mpi/openmpi/4.0.1/gcc/9


mpiexec -n 8 python3 pmpi.py