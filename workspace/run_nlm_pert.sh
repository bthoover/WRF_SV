#! /bin/sh

source activate bth-wrf_sv

root_dir=/Users/hoover/WRF_SV
curr_dir=${root_dir}/workspace
work_base_dir=${curr_dir}/tlm_nlm_pert
nlm_work_dir=${work_base_dir}/em_real
tlm_work_dir=${work_base_dir}/em_tlm

echo "moving from ${curr_dir} to ${nlm_work_dir}..."
cd ${nlm_work_dir}
# Add perturbation to NLM initialization file
echo "adding perturbation to NLM initialization"
mv wrfinput_d01 wrfinput_d01_uptd
python3 perturb_wrfinit_nlm.py << EOF
wrfinput_d01_uptd
wrfinput_d01_ptd
SV_pert.nc
EOF
# Run TLM
ln -s wrfinput_d01_ptd wrfinput_d01
./mpirun_wrf.csh ${nlm_work_dir} 12

# Return
echo "finished, returning to ${curr_dir}"
cd ${curr_dir}

