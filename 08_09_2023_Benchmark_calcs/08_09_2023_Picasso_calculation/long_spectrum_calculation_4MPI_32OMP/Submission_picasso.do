#!/bin/bash

# The name to show in queue lists for this job:
#SBATCH -J Malakhov

# Number of desired cores:
#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=32


# Amount of RAM needed for this job:
#SBATCH --mem=100gb

# The time the job will be running, 10 hours:
#SBATCH --time=48:00:00
#SBATCH --constraint=sr
##SBATCH --nodelist=sr[045-049,126-127,140-147]
##SBATCH --nodelist=sr[046-047]





num_mpi_tasks=$SLURM_NTASKS
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd `pwd`


mpirun -np $num_mpi_tasks ../../../../mpicbwe.x input.txt 1> output.log 2> output.err




exit 0




