# CSCI7551 - Parallel and Distributed Processing 
- __Author__: Javier Pastorino, University Of Colorado, Denver
- __Date__: September 2016
# Final Project

## File Description
- pfilter_seq.cpp   :    Particle Filters Sequential Algorithm
- pfilter_omp.cpp   :    Particle Filters OpenMp Algorithm
- pfilter_cuda.cu   :    Particle Filters CUDA Algorithm
- pfilter_mixed.cu  :    Particle Filters CUDA & OpenMP Algorithm

## Compiling Code On Hydra

### Sequential
g++  pfilter_seq.cpp  -o pfilter_seq

### OpenMP
g++ -O -fopenmp  pfilter_omp.cpp  -o pfilter_omp

### CUDA
nvcc pfilter_cuda.cu -o pfilter_cuda

### Mixed 
nvcc -Xcompiler -fopenmp  pfilter_mixed.cu -o pfilter_mixed

## Running Code On Hydra

### Sequential
bpsh 8 ./pfilter_seq 1000 100 0 | tee -a  ./sequential.out

### OpenMP
bpsh 9 ./pfilter_omp 1000 100 0 | tee -a  ./omp.out

### CUDA
bpsh 12 ./pfilter_cuda 1000 100 0 | tee -a  ./cuda.out

### Mixed 
bpsh 13 ./pfilter_mixed 1000 100 0 | tee -a  ./mixed.out



