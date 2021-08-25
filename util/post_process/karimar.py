#!/usr/bin/env python
#import atx93_klm

import numpy as np
import scipy as sp
import scipy.linalg
import netCDF4
from netCDF4 import Dataset
from scipy.io import FortranFile
from wrf import getvar, extract_times
import os
import subprocess
from subprocess import PIPE
import glob
import datetime as dt
import scipy.optimize
import shutil 
import numpy as np

print('Post Processing: STARTS')

#Input initial time (EDIT)
duration = 24  #duration of the OTI
st = dt.datetime(2020,11,10,0) #year,month,day, hour (in this order) (EDIT)

boundary_condition_update = 6
time_step_seconds = 180
dt = dt.timedelta(hours=duration)
et = st+dt
st_date_string = f'{st.year}-{str(st.month).zfill(2)}-{str(st.day).zfill(2)}_{str(st.hour).zfill(2)}:{str(st.minute).zfill(2)}:{str(st.second).zfill(2)}'
et_date_string = f'{et.year}-{str(et.month).zfill(2)}-{str(et.day).zfill(2)}_{str(et.hour).zfill(2)}:{str(et.minute).zfill(2)}:{str(st.second).zfill(2)}'

print(st_date_string)
print(et_date_string)


"""POST PROCESS  vector.dat"""
print("Post Process vector.dat")
    
"""This code reads one of the eigenvectors the output (units 30+j, j = 1, N) of program read_vector.f . The user must specify which eigenvector to post-process to create an initial condition file."""

J = 1   # (EDIT)

path ='/tornado/home1/class/fall18/ledesmamaldo/'
f    = FortranFile(path+'wrfdata_eta_24/fort.3'+str(J), 'r' )
XIC  = f.read_reals( dtype='float64' )
f.close
print('data read, vector of length:', XIC.size)

path_fwd_initial     = path+'WRFPLUSV3/test/em_fwd_sv/wrfout_d01_'+st_date_string
path_post_processing = path+'wrfdata_eta_24/' 
nc_tlm_path          = path+'/WRFPLUSV3/test/em_tlm_postprocessing/' 

shutil.copy(path_fwd_initial,path_post_processing) # copy wrfout_d01_initialtime file to post processing directory

os.chdir(path_post_processing) # cd to post processing directory

nc_fwd_fname = path_post_processing+'wrfout_d01_'+st_date_string
ncf_fwd      = Dataset(nc_fwd_fname, 'r') #To open wrfout_d01_initialtime file


# Get variables for size, not value
u = getvar(ncf_fwd,"U",meta=False)
v = getvar(ncf_fwd,"V",meta=False)
t = getvar(ncf_fwd,"T",meta=False)

nz=t.shape[0]
ny=u.shape[1]
nx=v.shape[2]
print(nz,nx,ny)

U = np.zeros([nz,ny,nx+1])
V = np.zeros([nz,ny+1,nx])
T = np.zeros([nz,ny,nx])

us=0
uf=nz*ny*(nx+1)

vs=uf
vf=vs+nz*(ny+1)*nx

ts = vf
tf=ts+nz*ny*nx

#########
u_wght = 0.5
v_wght = 0.5
t_wght = 0.5*(9.8*9.8)/(300.e-2)**2.
print(u_wght, v_wght, t_wght, ' <<<<< WEIGHTINGS')
#########

G_U = np.reshape(XIC[us:uf], (nz, ny, nx+1))/np.sqrt(u_wght,dtype=np.float64)
G_V = np.reshape(XIC[vs:vf], (nz, ny+1, nx))/np.sqrt(v_wght,dtype=np.float64)
G_T = np.reshape(XIC[ts:tf], (nz, ny, nx))/np.sqrt(t_wght,dtype=np.float64)

uu = G_U*np.sqrt(u_wght,dtype=np.float64)
vv = G_V*np.sqrt(v_wght,dtype=np.float64)
tt = G_T*np.sqrt(t_wght,dtype=np.float64)


#Norm

energy = np.sum(uu*uu,dtype=np.float64)+np.sum(vv*vv,dtype=np.float64)+np.sum(tt*tt,dtype=np.float64)

print(' for J = ', J, 'energy = ', energy)
if(energy!=0.0):
    G_U = G_U/np.sqrt(energy,dtype=np.float64)
    G_V = G_V/np.sqrt(energy,dtype=np.float64)
    G_T = G_T/np.sqrt(energy,dtype=np.float64)
else:
    print('SV', J, 'is a null vector. Check it')

####
print('normalized_energy = ', np.sum(G_U*G_U,dtype=np.float64)+np.sum(G_V*G_V,dtype=np.float64)+np.sum(G_T*G_T,dtype=np.float64))


os.chdir(path_post_processing) # cd  to post processing directory

#Rename wrfout_d01_initialtime to init_pert_d01
os.rename('wrfout_d01_'+st_date_string, 'init_pert_d01')

nc_fname = path_post_processing+'init_pert_d01'
ncf_tlm_new = Dataset(nc_fname, 'a')  #To open init_pert_d01 file

ncf_tlm_new.variables['G_U'][0] = ncf_tlm_new.variables['G_U'][0] + G_U[:]
ncf_tlm_new.variables['G_V'][0] = ncf_tlm_new.variables['G_V'][0] + G_V[:]
ncf_tlm_new.variables['G_T'][0] = ncf_tlm_new.variables['G_T'][0] + G_T[:]

ncf_tlm_new.close()

# Copy init_pert_d01 file to the em_tlm directory
shutil.copy(path_post_processing+'init_pert_d01',nc_tlm_path) 

#Rename init_pert_d01 to init_pert_d01J_d01 file in the post_processing directory 
os.rename('init_pert_d01', 'init_pert_00'+str(J)+'_d01')
          

print("Python: Run TLM post process")
os.chdir(nc_tlm_path)  # cd to em_tlm directory
subprocess.run('./run_wrf_tlm.csh') # run tlm
print("Python: TLM DONE post process")

"""POST PROCESS OUTPUT NEW ENERGY U V THETA"""
print("Post Process output new energy")

# Check energy of (first) leading SV at t = 0

nc_fname = path_post_processing+'init_pert_00'+str(J)+'_d01'  
ncf = Dataset(nc_fname, 'r') #To open init_pert_d01J file

G_U = getvar(ncf,"G_U",meta=False)
G_V = getvar(ncf,"G_V",meta=False)
G_T = getvar(ncf,"G_T",meta=False)

ncf.close()

nz=G_T.shape[0]
ny=G_U.shape[1]
nx=G_V.shape[2]
print(nz,nx,ny)
#

us=0
uf=nz*ny*(nx+1)

vs=uf
vf=vs+nz*(ny+1)*nx

ts = vf
tf=ts+nz*ny*nx


u_wght = 0.5
v_wght = 0.5
g = 9.8
theta_ref = 300
N = .01
t_wght = 0.5*(g*g)/(theta_ref*N)**2.

G_U = G_U*np.sqrt(u_wght,dtype=np.float64)
G_V = G_V*np.sqrt(v_wght,dtype=np.float64)
G_T = G_T*np.sqrt(t_wght,dtype=np.float64)


energy_init = np.sum(G_U*G_U)+np.sum(G_V*G_V)+np.sum(G_T*G_T)

# Copy tlout_d01_finaltime file to post processing directory
shutil.copy(nc_tlm_path+'tlout_d01_'+et_date_string,path_post_processing)  
os.chdir(path_post_processing) # cd  to post processing directory

#Rename tlout_d01_finaltime file to final_pert_00J_d01
os.rename('tlout_d01_'+et_date_string, 'final_pert_00'+str(J)+'_d01')

nc_fname = path_post_processing+'final_pert_00'+str(J)+'_d01'  
ncfi = Dataset(nc_fname, 'r') #To open final_pert_00J file

#Local projector Operator
LPO_u=np.zeros([nz,ny,nx+1])
LPO_u[1:-1,57:100,45:105]=1.0

LPO_v=np.zeros([nz,ny+1,nx])
LPO_v[1:-1,57:100,45:105]=1.0

LPO_T=np.zeros([nz,ny,nx])
LPO_T[1:-1,57:100,45:105]=1.0

u_final = getvar(ncfi,"G_U",meta=False)*np.sqrt(u_wght,dtype=np.float64)*LPO_u[:,:,:]
v_final = getvar(ncfi,"G_V",meta=False)*np.sqrt(v_wght,dtype=np.float64)*LPO_v[:,:,:]
t_final = getvar(ncfi,"G_T",meta=False)*np.sqrt(t_wght,dtype=np.float64)*LPO_T[:,:,:]


energy_final=np.sum(u_final*u_final,dtype=np.float64)+np.sum(v_final*v_final,dtype=np.float64)+np.sum(t_final*t_final,dtype=np.float64)

print('For SV, initial energy = ,', energy_init, 'and final energy = ', energy_final)

"""   PERTRUB NLM IC    """
print("Perturb NLM IC")

path_fwd_wrfinput      = path+'WRFPLUSV3/test/em_fwd_sv/wrfinput_d01'
path_fwd_posprocessing = path+'WRFPLUSV3/test/em_fwd_posprocessing/'

##### POSITIVE Pertubation ######

shutil.copy(path_fwd_wrfinput,path_fwd_posprocessing) # copy wrfinput_d01 to em_fwd post processing

os.chdir(path_fwd_posprocessing) # cd to em_fwd post processing directory
nc_fname = path_fwd_posprocessing+'wrfinput_d01'   
ncf_initial = Dataset(nc_fname, "a", format="NETCDF4")  #to open the wrfinput_d01 file

nc_fname = path_post_processing+'init_pert_00'+str(J)+'_d01'
ncf_initial_p = Dataset(nc_fname,'r')

G_U = ncf_initial_p.variables['G_U'][0, :, :, :]     # u-wind
G_V = ncf_initial_p.variables['G_V'][0, :, :, :]     # v-wind
G_T = ncf_initial_p.variables['G_T'][0, :, :, :]     # theta

print('max G_U is: ', np.amax(G_U))
print('max G_V is: ', np.amax(G_V))
print('max G_T is: ', np.amax(G_T))

#Add perturbation

ncf_initial.variables['U'][0,:] = ncf_initial.variables['U'][0,:] + G_U*50
ncf_initial.variables['V'][0,:] = ncf_initial.variables['V'][0,:] + G_V*50
ncf_initial.variables['T'][0,:] = ncf_initial.variables['T'][0,:] + G_T*50

ncf_initial.close()

print("Python: Run FWD post process +pert")
os.chdir(path_fwd_posprocessing) # cd to em_fwd post processing directory
subprocess.run('./run_wrf_fwd.csh') # To run fwd
print("Python: FWD DONE post process +pert")

shutil.copy(path_fwd_posprocessing+'wrfout_d01_'+st_date_string,path_post_processing) # copy wrfout_d01_initialtime to post processing directory
os.chdir(path_post_processing) # cd  to post processing directory
os.rename('wrfout_d01_'+st_date_string, 'wrfout_d01'+str(J)) # rename wrfout_d01_initialtime to wrfout_d01J


##### NEGATIVE Perturbation ######

shutil.copy(path_fwd_wrfinput,path_fwd_posprocessing)  # copy wrfinput_d01 to em_fwd post processing

os.chdir(path_fwd_posprocessing) # cd to em_fwd post processing directory
nc_fname = path_fwd_posprocessing+'wrfinput_d01'  
ncf_initial = Dataset(nc_fname, "a", format="NETCDF4")  #to open the wrfinput_d01 file

nc_fname = path_post_processing+'init_pert_00'+str(J)+'_d01'
ncf_initial_p = Dataset(nc_fname,'r')

G_U = ncf_initial_p.variables['G_U'][0, :, :, :]     # u-wind
G_V = ncf_initial_p.variables['G_V'][0, :, :, :]     # v-wind
G_T = ncf_initial_p.variables['G_T'][0, :, :, :]     # theta

print('max G_U is: ', np.amax(G_U))
print('max G_V is: ', np.amax(G_V))
print('max G_T is: ', np.amax(G_T))

#Add perturbation

ncf_initial.variables['U'][0,:] = ncf_initial.variables['U'][0,:] - G_U*50
ncf_initial.variables['V'][0,:] = ncf_initial.variables['V'][0,:] - G_V*50
ncf_initial.variables['T'][0,:] = ncf_initial.variables['T'][0,:] - G_T*50

ncf_initial.close()

print("Python: Run FWD post process -pert")
os.chdir(path_fwd_posprocessing) # cd to em_fwd post processing directory
subprocess.run('./run_wrf_fwd.csh') #to run fwd
print("Python: FWD DONE post process -pert")

shutil.copy(path_fwd_posprocessing+'wrfout_d01_'+st_date_string,path_post_processing) # copy wrfout_d01_initialtime to post processing directory
os.chdir(path_post_processing) # cd  to post processing directory
os.rename('wrfout_d01_'+st_date_string, 'wrfout_d01'+str(J)+'_neg') # rename wrfout_d01_initialtime to wrfout_d01J_neg

print('Post Processing: DONE')