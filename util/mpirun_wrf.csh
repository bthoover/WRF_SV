#! /bin/csh

set rundir=${1}
set nproc=${2}

cd ${rundir}
mpirun -np ${nproc} ./wrf.exe

