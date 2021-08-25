#! /bin/sh

source activate bth-wrf_sv

root_dir=/Users/hoover/WRF_SV

curr_dir=${root_dir}/workspace
util_dir=${root_dir}/util
work_dir=${curr_dir}/sv_pert
nonlin_dir=${curr_dir}/nonlin_trj
pproc_dir=${curr_dir}/post_process

yyyy_beg=2020
mm_beg=03
dd_beg=06
hh_beg=12
max_u=3.0
max_v=3.0
max_t=1.5
date_beg=${yyyy_beg}-${mm_beg}-${dd_beg}_${hh_beg}:00:00

echo "moving from ${curr_dir} to ${work_dir}..."

mkdir -p ${work_dir}
cd ${work_dir}
cp ${nonlin_dir}/wrfout_d01_${date_beg} .
cp ${pproc_dir}/init_pert_d01_J1_${date_beg} .
cp ${util_dir}/compute_sv_pert.py .
# Run python program
python3 compute_sv_pert.py << EOF
${work_dir}/
${yyyy_beg}
${mm_beg}
${dd_beg}
${hh_beg}
${max_u}
${max_v}
${max_t}
EOF
# Return
echo "returning from ${work_dir} to ${curr_dir}..."
cd ${curr_dir}

