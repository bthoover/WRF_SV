#! /bin/sh

source activate bth-wrf_sv

root_dir=/Users/hoover/WRF_SV

curr_dir=${root_dir}/workspace
util_dir=${root_dir}/util
wrfplus_repo_dir=/Users/morgan/WRFPLUSV3
namelist_repo_dir=${root_dir}/namelists
work_dir=${curr_dir}/post_process
sv_repo_dir=${curr_dir}/tlm_adj_sv
nonlin_dir=${curr_dir}/nonlin_trj

set yyyy_beg=2020
set mm_beg=03
set dd_beg=06
set hh_beg=12
set fcst_len_hrs=24
set update_bdy_hrs=1
set wrf_tstep_sec=60

echo "moving from ${curr_dir} to ${work_dir}..."

mkdir -p ${work_dir}
cd ${work_dir}
#
# read_vectors
#
cp ${sv_repo_dir}/vector.dat .
cp ${util_dir}/post_process/read_vectors.f .
gfortran -o read_vectors.exe read_vectors.f
./read_vectors.exe
# Return
echo "running, returning from ${work_dir} to ${curr_dir}..."
cd ${curr_dir}

