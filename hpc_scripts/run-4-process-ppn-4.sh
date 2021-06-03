#!/bin/bash

#PBS -l select=1:ncpus=4:pmem=1gb

#PBS -l walltime=0:35:00

#PBS -q short_cpuQ

module load mpich-3.2
# module load openmpi-4.0.4
echo "Number of processes: 4"
echo "Process per node: 1"
mpiexec -n 4 -ppn 1 /home/matteo.bortolon/ipfp-hpc-cluster/build/distributedIPFP /home/matteo.bortolon/data/aggregate_visit_matrix.txt /home/matteo.bortolon/data/poi_marginals_2020_03_02.txt /home/matteo.bortolon/data/cbg_marginals_2020_03_02.txt /home/matteo.bortolon/output