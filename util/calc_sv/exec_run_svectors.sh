#! /bin/sh

source activate bth-wrf_sv

# Compile svectors.exe
gfortran -o svectors.exe -fcheck=bounds -freal-4-real-8 svectors.f lancz4p.f lancz4.f atx93_klm.f90 -mcmodel=large
# Run svectors.exe
./svectors.exe

