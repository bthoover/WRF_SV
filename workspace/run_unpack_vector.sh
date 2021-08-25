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

yyyy_beg=2020
mm_beg=03
dd_beg=06
hh_beg=12
fcst_len_hrs=24
update_bdy_hrs=1
wrf_tstep_sec=60
J_val=6

echo "moving from ${curr_dir} to ${work_dir}..."

mkdir -p ${work_dir}
cd ${work_dir}
#
# unpack_vectors
#
date_beg=${yyyy_beg}-${mm_beg}-${dd_beg}_${hh_beg}:00:00
cp ${nonlin_dir}/wrfout_d01_${date_beg} .
cp ${util_dir}/post_process/unpack_vectors.py .

echo ${work_dir}/
echo ${fcst_len_hrs}
echo ${yyyy_beg}
echo ${mm_beg}
echo ${dd_beg}
echo ${hh_beg}
echo ${update_bdy_hrs}
echo ${wrf_tstep_sec}
echo ${J_val}


python3 unpack_vectors.py << EOF
${work_dir}/
${fcst_len_hrs}
${yyyy_beg}
${mm_beg}
${dd_beg}
${hh_beg}
${update_bdy_hrs}
${wrf_tstep_sec}
${J_val}
EOF
# Return
echo "running, returning from ${work_dir} to ${curr_dir}..."
cd ${curr_dir}

