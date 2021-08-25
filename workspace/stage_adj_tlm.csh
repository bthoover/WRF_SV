#!/bin/tcsh

set root_dir=/Users/hoover/WRF_SV

set curr_dir=${root_dir}/workspace
set util_dir=${root_dir}/util
set wrfplus_repo_dir=/Users/morgan/WRFPLUSV3
set namelist_repo_dir=${root_dir}/namelists
set work_base_dir=${curr_dir}/tlm_adj_sv
set adj_work_dir=${work_base_dir}/em_adj
set tlm_work_dir=${work_base_dir}/em_tlm
set init_repo_dir=${root_dir}/ICs
set wrfplus_trj_dir=${curr_dir}/nonlin_trj

set yyyy_beg=2020
set mm_beg=03
set dd_beg=06
set hh_beg=12
set fcst_len_hrs=24
set yyyy_end=2020
set mm_end=03
set dd_end=07
set hh_end=12
set lpo_kmin=1
set lpo_kmax=41
set lpo_jmin=63
set lpo_jmax=83
set lpo_imin=135
set lpo_imax=155
set update_bdy_hrs=1
set wrf_tstep_sec=60
set wrf_root_dir=${work_base_dir}

set date_beg=${yyyy_beg}-${mm_beg}-${dd_beg}_${hh_beg}:00:00
set date_end=${yyyy_end}-${mm_end}-${dd_end}_${hh_end}:00:00


echo "moving from ${curr_dir} to ${work_base_dir}..."

mkdir -p ${work_base_dir}
cd ${work_base_dir}

echo "staging SV in ${work_base_dir}..."

cp ${util_dir}/calc_sv/svectors.f .
cp ${util_dir}/calc_sv/lancz4p.f .
cp ${util_dir}/calc_sv/lancz4.f .
cp ${util_dir}/calc_sv/atx93_klm.f90 .
cp ${util_dir}/calc_sv/exec_run_svectors.sh .
cp ${util_dir}/calc_sv/govern_wrf_tlm_adj.py .

echo "${yyyy_beg}" > governor_inputs.txt
echo "${mm_beg}" >> governor_inputs.txt
echo "${dd_beg}" >> governor_inputs.txt
echo "${hh_beg}" >> governor_inputs.txt
echo "${fcst_len_hrs}" >> governor_inputs.txt
echo "${lpo_kmin}" >> governor_inputs.txt
echo "${lpo_kmax}" >> governor_inputs.txt
echo "${lpo_jmin}" >> governor_inputs.txt
echo "${lpo_jmax}" >> governor_inputs.txt
echo "${lpo_imin}" >> governor_inputs.txt
echo "${lpo_imax}" >> governor_inputs.txt
echo "${update_bdy_hrs}" >> governor_inputs.txt
echo "${wrf_tstep_sec}" >> governor_inputs.txt
echo "${wrf_root_dir}" >> governor_inputs.txt

echo "staging ADJ in ${adj_work_dir}..."
mkdir -p ${adj_work_dir}

ln -sf ${wrfplus_repo_dir}/main/wrf.exe ${adj_work_dir}/.
ln -sf ${wrfplus_repo_dir}/run/RRTM_DATA_DBL ${adj_work_dir}/RRTM_DATA
ln -sf ${wrfplus_repo_dir}/run/RRTMG_LW_DATA_DBL ${adj_work_dir}/RRTMG_LW_DATA
ln -sf ${wrfplus_repo_dir}/run/RRTMG_SW_DATA_DBL ${adj_work_dir}/RRTMG_SW_DATA
ln -sf ${wrfplus_repo_dir}/run/SOILPARM.TBL ${adj_work_dir}/.  
ln -sf ${wrfplus_repo_dir}/run/VEGPARM.TBL ${adj_work_dir}/.  
ln -sf ${wrfplus_repo_dir}/run/GENPARM.TBL ${adj_work_dir}/.  
ln -sf ${wrfplus_repo_dir}/run/LANDUSE.TBL ${adj_work_dir}/.

ln -sf ${namelist_repo_dir}/namelist.input.adj ${adj_work_dir}/namelist.input

cp ${wrfplus_trj_dir}/wrfout_d01_* ${adj_work_dir}/.
ln -sf ${wrfplus_trj_dir}/auxhist6_d01_* ${adj_work_dir}/.

cp ${init_repo_dir}/wrfinput_d01_uptd ${adj_work_dir}/wrfinput_d01
cp ${init_repo_dir}/wrfbdy_d01 ${adj_work_dir}/.
cp ${util_dir}/mpirun_wrf.csh ${adj_work_dir}/.
cp ${util_dir}/plus.io_config ${adj_work_dir}/.

cp ${adj_work_dir}/wrfout_d01_${date_end} ${adj_work_dir}/fcst.nc
mv ${adj_work_dir}/wrfout_d01_${date_end} ${adj_work_dir}/wrfout_final_saved
ln -sf ${adj_work_dir}/wrfinput_d01 ${adj_work_dir}/xref.nc
mv ${adj_work_dir}/fcst.nc ${adj_work_dir}/final_sens_d01_${date_end}
ln -sf ${adj_work_dir}/final_sens_d01_${date_end} ${adj_work_dir}/wrfout_d01_${date_end}
ln -sf ${adj_work_dir}/final_sens_d01_${date_end} ${adj_work_dir}/final_sens_d01

echo "staging TLM in ${tlm_work_dir}..."
mkdir -p ${tlm_work_dir}

ln -sf ${wrfplus_repo_dir}/main/wrf.exe ${tlm_work_dir}/.
ln -sf ${wrfplus_repo_dir}/run/RRTM_DATA_DBL ${tlm_work_dir}/RRTM_DATA
ln -sf ${wrfplus_repo_dir}/run/RRTMG_LW_DATA_DBL ${tlm_work_dir}/RRTMG_LW_DATA
ln -sf ${wrfplus_repo_dir}/run/RRTMG_SW_DATA_DBL ${tlm_work_dir}/RRTMG_SW_DATA
ln -sf ${wrfplus_repo_dir}/run/SOILPARM.TBL ${tlm_work_dir}/.  
ln -sf ${wrfplus_repo_dir}/run/VEGPARM.TBL ${tlm_work_dir}/.  
ln -sf ${wrfplus_repo_dir}/run/GENPARM.TBL ${tlm_work_dir}/.  
ln -sf ${wrfplus_repo_dir}/run/LANDUSE.TBL ${tlm_work_dir}/.

ln -sf ${namelist_repo_dir}/namelist.input.tlm ${tlm_work_dir}/namelist.input

cp ${wrfplus_trj_dir}/wrfout_d01_* ${tlm_work_dir}/.
ln -sf ${wrfplus_trj_dir}/auxhist6_d01_* ${tlm_work_dir}/.
mv ${tlm_work_dir}/wrfout_d01_${date_beg} ${tlm_work_dir}/wrfout_init_saved
cp ${tlm_work_dir}/wrfout_init_saved ${tlm_work_dir}/init_pert_d01
ln -s ${tlm_work_dir}/init_pert_d01 ${tlm_work_dir}/wrfout_d01_${date_beg}

cp ${init_repo_dir}/wrfinput_d01_uptd ${tlm_work_dir}/wrfinput_d01
cp ${init_repo_dir}/wrfbdy_d01 ${tlm_work_dir}/.
cp ${util_dir}/mpirun_wrf.csh ${tlm_work_dir}/.
cp ${util_dir}/plus.io_config ${tlm_work_dir}/.

cd ${curr_dir}

