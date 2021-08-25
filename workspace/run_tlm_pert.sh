#! /bin/sh

source activate bth-wrf_sv

root_dir=/Users/hoover/WRF_SV
curr_dir=${root_dir}/workspace
work_base_dir=${curr_dir}/tlm_nlm_pert
nlm_work_dir=${work_base_dir}/em_real
tlm_work_dir=${work_base_dir}/em_tlm

echo "moving from ${curr_dir} to ${tlm_work_dir}..."
cd ${tlm_work_dir}
# Add perturbation to TLM initialization file
echo "adding perturbation to TLM initialization"
python3 perturb_wrfinit_tlm.py << EOF
init_pert_d01
SV_pert.nc
EOF
# Run TLM
./mpirun_wrf.csh ${tlm_work_dir} 12

# Return
echo "finished, returning to ${curr_dir}"
cd ${curr_dir}

