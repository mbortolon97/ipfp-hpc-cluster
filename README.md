# ipfp-hpc-cluster
This repository contain a implementation of the IPFP that work on large HPC cluster

# Compile from source

To compile from source execute:
```
cd ipfp-hpc-cluster
mkdir build
cd build
cmake -DMPI_CXX_COMPILER=/usr/local/bin/mpicxx -DMPI_C_COMPILER=/usr/local/bin/mpicc -DMPIEXEC_EXECUTABLE=/usr/local/bin/mpiexec.hydra ..
make -j $(nproc)
```