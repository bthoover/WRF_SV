#!/bin/tcsh

set root_dir=/Users/hoover/WRF_SV

set curr_dir=${root_dir}/workspace
set util_dir=${root_dir}/util
set wrfplus_repo_dir=/Users/morgan/WRFPLUSV3
set namelist_repo_dir=${root_dir}/namelists
set work_dir=${curr_dir}/nonlin_trj
set init_repo_dir=${root_dir}/ICs

echo "moving from ${curr_dir} to ${work_dir}..."

mkdir -p ${work_dir}
cd ${work_dir}
ln -sf ${wrfplus_repo_dir}/main/wrf.exe .
ln -sf ${wrfplus_repo_dir}/run/RRTM_DATA_DBL RRTM_DATA
ln -sf ${wrfplus_repo_dir}/run/RRTMG_LW_DATA_DBL RRTMG_LW_DATA
ln -sf ${wrfplus_repo_dir}/run/RRTMG_SW_DATA_DBL RRTMG_SW_DATA
ln -sf ${wrfplus_repo_dir}/run/SOILPARM.TBL .  
ln -sf ${wrfplus_repo_dir}/run/VEGPARM.TBL .  
ln -sf ${wrfplus_repo_dir}/run/GENPARM.TBL .  
ln -sf ${wrfplus_repo_dir}/run/LANDUSE.TBL .

ln -sf ${namelist_repo_dir}/namelist.input.trj namelist.input
cp ${init_repo_dir}/wrfinput_d01_uptd wrfinput_d01
cp ${init_repo_dir}/wrfbdy_d01 .
cp ${util_dir}/plus.io_config .
mpirun -np 12 ./wrf.exe

echo "running, returning from ${work_dir} to ${curr_dir}..."
cd ${curr_dir}

