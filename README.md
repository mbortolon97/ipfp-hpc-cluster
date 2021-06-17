# [HPC Project] High performance IPFP

This repository contain an implementation of the **Iterative Proportional Fitting Procedure** (IPFP) which can be run on large HPC cluster.
This implementation is thought to speed-up the the IPFP step used by [Albani et al.](https://www.nature.com/articles/s41598-021-88281-w#citeas), where they use this algorithm to forecast how many people will frequent a **point of interest** (eg. restaurants, shops, ...) in a given hour of the day.

### Useful Links

- Resources : [Google Drive folder](https://drive.google.com/drive/folders/1l0M-Kd9Mqq0Dy5ns5Z3sh17xlMcemMrd?usp=sharing)

### How to run the code

To compile from source execute:

```
cd ipfp-hpc-cluster
mkdir build
cd build
cmake -DMPI_CXX_COMPILER=/usr/local/bin/mpicxx -DMPI_C_COMPILER=/usr/local/bin/mpicc -DMPIEXEC_EXECUTABLE=/usr/local/bin/mpiexec.hydra ..
make -j $(nproc)
```
