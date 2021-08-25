#!/bin/tcsh

set root_dir=/Users/hoover/WRF_SV

set curr_dir=${root_dir}/workspace
set util_dir=${root_dir}/util
set wrfplus_repo_dir=/Users/morgan/WRFPLUSV3
set namelist_repo_dir=${root_dir}/namelists
set work_base_dir=${curr_dir}/tlm_nlm_pert
set nlm_work_dir=${work_base_dir}/em_real
set tlm_work_dir=${work_base_dir}/em_tlm
set init_repo_dir=${root_dir}/ICs
set sv_pert_dir=${curr_dir}/sv_pert
set wrfplus_trj_dir=${curr_dir}/nonlin_trj

set yyyy_beg=2020
set mm_beg=03
set dd_beg=06
set hh_beg=12
set date_beg=${yyyy_beg}-${mm_beg}-${dd_beg}_${hh_beg}:00:00

echo "moving from ${curr_dir} to ${work_base_dir}..."

mkdir -p ${work_base_dir}
cd ${work_base_dir}

#
# Stage NLM
#
echo "staging NLM in ${nlm_work_dir}..."
mkdir -p ${nlm_work_dir}
# Copy WRF model-files
ln -sf ${wrfplus_repo_dir}/main/wrf.exe ${nlm_work_dir}/.
ln -sf ${wrfplus_repo_dir}/run/RRTM_DATA_DBL ${nlm_work_dir}/RRTM_DATA
ln -sf ${wrfplus_repo_dir}/run/RRTMG_LW_DATA_DBL ${nlm_work_dir}/RRTMG_LW_DATA
ln -sf ${wrfplus_repo_dir}/run/RRTMG_SW_DATA_DBL ${nlm_work_dir}/RRTMG_SW_DATA
ln -sf ${wrfplus_repo_dir}/run/SOILPARM.TBL ${nlm_work_dir}/.  
ln -sf ${wrfplus_repo_dir}/run/VEGPARM.TBL ${nlm_work_dir}/.  
ln -sf ${wrfplus_repo_dir}/run/GENPARM.TBL ${nlm_work_dir}/.  
ln -sf ${wrfplus_repo_dir}/run/LANDUSE.TBL ${nlm_work_dir}/.
# Soft-link namelist for NLM
ln -sf ${namelist_repo_dir}/namelist.input.nlm ${nlm_work_dir}/namelist.input
# Copy input files and run-script
cp ${init_repo_dir}/wrfinput_d01_uptd ${nlm_work_dir}/wrfinput_d01
cp ${init_repo_dir}/wrfbdy_d01 ${nlm_work_dir}/.
cp ${util_dir}/mpirun_wrf.csh ${nlm_work_dir}/.
cp ${util_dir}/plus.io_config ${nlm_work_dir}/.
# Copy IC perturbation program
cp ${util_dir}/perturb_wrfinit_nlm.py ${nlm_work_dir}/.
# Copy SV perturbation file
cp ${sv_pert_dir}/SV_pert.nc ${nlm_work_dir}/.
#
# Stage TLM
#
echo "staging TLM in ${tlm_work_dir}..."
mkdir -p ${tlm_work_dir}
# Copy WRF model-files
ln -sf ${wrfplus_repo_dir}/main/wrf.exe ${tlm_work_dir}/.
ln -sf ${wrfplus_repo_dir}/run/RRTM_DATA_DBL ${tlm_work_dir}/RRTM_DATA
ln -sf ${wrfplus_repo_dir}/run/RRTMG_LW_DATA_DBL ${tlm_work_dir}/RRTMG_LW_DATA
ln -sf ${wrfplus_repo_dir}/run/RRTMG_SW_DATA_DBL ${tlm_work_dir}/RRTMG_SW_DATA
ln -sf ${wrfplus_repo_dir}/run/SOILPARM.TBL ${tlm_work_dir}/.  
ln -sf ${wrfplus_repo_dir}/run/VEGPARM.TBL ${tlm_work_dir}/.  
ln -sf ${wrfplus_repo_dir}/run/GENPARM.TBL ${tlm_work_dir}/.  
ln -sf ${wrfplus_repo_dir}/run/LANDUSE.TBL ${tlm_work_dir}/.
# Soft-link namelist for TLM
ln -sf ${namelist_repo_dir}/namelist.input.tlm ${tlm_work_dir}/namelist.input
# Copy (unperturbed) NLM trajectory files
cp ${wrfplus_trj_dir}/wrfout_d01_* ${tlm_work_dir}/.
ln -sf ${wrfplus_trj_dir}/auxhist6_d01_* ${tlm_work_dir}/.
mv ${tlm_work_dir}/wrfout_d01_${date_beg} ${tlm_work_dir}/wrfout_init_saved
cp ${tlm_work_dir}/wrfout_init_saved ${tlm_work_dir}/init_pert_d01
ln -s ${tlm_work_dir}/init_pert_d01 ${tlm_work_dir}/wrfout_d01_${date_beg}
# Copy initial condition files and run-script, io_config file
cp ${init_repo_dir}/wrfinput_d01_uptd ${tlm_work_dir}/wrfinput_d01
cp ${init_repo_dir}/wrfbdy_d01 ${tlm_work_dir}/.
cp ${util_dir}/mpirun_wrf.csh ${tlm_work_dir}/.
cp ${util_dir}/plus.io_config ${tlm_work_dir}/.
# Copy IC perturbation program
cp ${util_dir}/perturb_wrfinit_tlm.py ${tlm_work_dir}/.
# Copy SV perturbation file
cp ${sv_pert_dir}/SV_pert.nc ${tlm_work_dir}/.

# Return
echo "finished, returning to ${curr_dir}"
cd ${curr_dir}

