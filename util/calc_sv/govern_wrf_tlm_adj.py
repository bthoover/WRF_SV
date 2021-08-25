######################################################################################
#
# PYTHON 3 PROGRAM
#
#
#
#
######################################################################################
#
# Import necessary modules.routines
#
from scipy.io import FortranFile #.................................................... FORTRAN unformatted i/o module
import os #........................................................................... operating system module
import subprocess #................................................................... subprocess management module
from subprocess import PIPE #......................................................... pipeline to standard stream routine
import numpy as np #.................................................................. array module
import glob #......................................................................... global functions module
from netCDF4 import Dataset #......................................................... netCDF i/o routine
from wrf import getvar, extract_times #............................................... WRF variable/time extraction routines
import datetime as dt #............................................................... datetime module
import scipy as sp #.................................................................. scientific functions module
import scipy.optimize #............................................................... function optimization routine
import scipy.linalg #................................................................. linear algebra routine
import shutil #....................................................................... shell utility module
#
######################################################################################
#
# Get user inputs from text file (governor_inputs.txt, generated before running)
# 
inp_hdl = open('governor_inputs.txt') #............................................... file-handler for input text file (generated before running)
inp_list = [line.rstrip('\n') for line in inp_hdl] #.................................. list of inputs from file-handler (newline character stripped)
beg_yyyy       = inp_list[0] #........................................................ beginning year (string, YYYY)
beg_mm         = inp_list[1] #........................................................ beginning year (string, MM)
beg_dd         = inp_list[2] #........................................................ beginning year (string, DD)
beg_hh         = inp_list[3] #........................................................ beginning year (string, HH)
fcst_len_hrs   = inp_list[4] #........................................................ forecast duration (string, hours)
lpo_kmin       = inp_list[5] #........................................................ local projection operator (LPO): Minimum k (string)
lpo_kmax       = inp_list[6] #........................................................ local projection operator (LPO): Maximum k (string)
lpo_jmin       = inp_list[7] #........................................................ local projection operator (LPO): Minimum j (string)
lpo_jmax       = inp_list[8] #........................................................ local projection operator (LPO): Maximum j (string)
lpo_imin       = inp_list[9] #........................................................ local projection operator (LPO): Minimum i (string)
lpo_imax       = inp_list[10] #....................................................... local projection operator (LPO): Maximum i (string)
update_bdy_hrs = inp_list[11] #....................................................... update interval of boundary conditions (string, hours)
wrf_tstep_sec  = inp_list[12] #....................................................... WRF model time-step (string, seconds)
wrf_root_dir   = inp_list[13] #....................................................... WRF model root-directory (contains em_adj, em_tlm, fort.16)
inp_hdl.close()
#
######################################################################################
#
# Initialize some variables
#
duration = float(fcst_len_hrs) #...................................................... forecast duration (float, hours)
boundary_condition_update = float(update_bdy_hrs) #................................... update interval of boundary conditions (float, hours)
time_step_seconds = float(wrf_tstep_sec) #............................................ WRF model time-step(float, seconds)
st = dt.datetime.strptime(beg_yyyy+beg_mm+beg_dd+beg_hh,'%Y%m%d%H') #................. starting datetime
et = st + dt.timedelta(hours=duration) #.............................................. ending (forecast verification) datetime
st_date_string = dt.datetime.strftime(st,'%Y-%m-%d_%H:%M:%S') #....................... starting date-string (WRF-file formatted)
et_date_string = dt.datetime.strftime(et,'%Y-%m-%d_%H:%M:%S') #....................... ending date-string (WRF-file formatted)
path = wrf_root_dir
nc_tlm_path_fwd = path+'/em_tlm/wrfout_d01_'+st_date_string #......................... full-path to WRF initial condition file (residing in TLM directory)
nc_fwd = Dataset(nc_tlm_path_fwd, 'r') #.............................................. WRF initial condition file-handle
#
######################################################################################
#
# Set up initial perturbation variables, dimensions, constants
#
u = getvar(nc_fwd,"G_U",meta=False) #................................................. perturbation zonal flow (initialized to default vals) 
v = getvar(nc_fwd,"G_V",meta=False) #................................................. perturbation merid flow (initialized to default vals)
T = getvar(nc_fwd,"G_T",meta=False) #................................................. perturbation temperature (initialized to default vals)
#
nz = T.shape[0] #..................................................................... number of levels
ns = u.shape[1] #..................................................................... number of north-south points
ew = v.shape[2] #..................................................................... number of east-west points
print('BTH:','nz=',nz,'ns=',ns,'ew=',ew)
#
"""Goal: K= (E^^(-1/2))P*EP(E^^(-1/2)) where P=TLM, P*=ADJ, E = weighted matrix"""
"""E is a matrix of weights derivable from the analytic expression for the quantity used to \
 introduce the norm. A diagonal matrix."""
#
# constants
g = 9.8                                      # gravity
cp = 1006.                                   # J/kg K Specific heat capacity
N = 0.01                                     # reference value for the Brunt-Vaisaila frequency
theta_ref =300.                              # theta reference
w = 0.5                                      # weighting for u and v
wT = 0.5*((g**2)/(((theta_ref)**2)*(N**2)))  # weighting for theta
#
""" To Remove Lateral Boundaries """
# LPO limits
kmin = int(lpo_kmin) - 1                     # JPO minimum k (integer)
kmax = int(lpo_kmax)                         # JPO maximum k (integer)
jmin = int(lpo_jmin) - 1                     # JPO minimum j (integer)
jmax = int(lpo_jmax)                         # JPO maximum j (integer)
imin = int(lpo_imin) - 1                     # JPO minimum i (integer)
imax = int(lpo_imax)                         # JPO maximum i (integer)
# Define LPO
LPO_u=np.zeros([nz,ns,ew+1]) #........................................................ LPO for u-field (initialized to zero)
LPO_u[kmin:kmax,jmin:jmax,imin:imax]=1.0
LPO_v=np.zeros([nz,ns+1,ew]) #........................................................ LPO for v-field (initialized to zero)
LPO_v[kmin:kmax,jmin:jmax,imin:imax]=1.0
LPO_T=np.zeros([nz,ns,ew]) #.......................................................... LPO for T-field (initialized to zero)
LPO_T[kmin:kmax,jmin:jmax,imin:imax]=1.0
#
######################################################################################
#
# Initialize WRF TLM
#
# STEP 1: Read in IC from atx93 call
f16 = FortranFile('fort.16', 'r' ) #.................................................. file-handle for fort.16
XIC = f16.read_reals( dtype='float64') #.............................................. initial condition state vector
print(len(XIC))
# STEP 2: Reshape XIC, multiply the state vector with matrix  E^^(-1/2) and remove 
#         lateral boundaries
x1 = XIC[0:nz*ns*u.shape[2]].reshape([nz,ns,ew+1]) #.................................. u-field component of state vector, reshaped to [nz,ns,ew+1]
x1 = x1 * (w**(-0.5))
x2 = XIC[nz*ns*u.shape[2]:nz*ns*u.shape[2]+ nz*ew*v.shape[1]].reshape([nz,ns+1,ew]) #. v-field component of state vector, reshaped to [nz,ns+1,ew]
x2 = x2 * (w**(-0.5))
x3 = XIC[nz*ns*u.shape[2]+ nz*ew*v.shape[1]:len(XIC)].reshape([nz,ns,ew]) #........... T-field component of state vector, reshaped to [nz,ns,ew]
x3 = x3 * (wT**(-0.5))
# STEP 3: Create the initial condition file (init_pert_d01) for the WRF_TLM
nc_tlm_path = path+'/em_tlm' #........................................................ path to em_tlm directory
# copy path_fwd_initial (file) to nc_tlm_path (subdirectory)
#shutil.copy(path_fwd_initial,nc_tlm_path)
# move to nc_tlm_path
os.chdir(nc_tlm_path)
# rename initial wrfout file to init_pert_d01 
#os.rename('wrfout_d01_'+st_date_string, 'init_pert_d01')
nc_tlm = nc_tlm_path+'/init_pert_d01' #............................................... full-path to init_pert_d01 file
ncf_tlm = Dataset(nc_tlm, 'a') #...................................................... file-handle for init_pert_d01 file
# write to G_* variables
ncf_tlm.variables['G_U'][0]= x1[:]
ncf_tlm.variables['G_V'][0]= x2[:]
ncf_tlm.variables['G_T'][0]= x3[:]
ncf_tlm.close()
#
######################################################################################
#
# Run TLM, produce input for ADJ
#
# STEP 1: Run TLM (integrate the resulting vector on the TLM from to to t)
print("Python: Run TLM")
os.chdir(path +'/em_tlm') 
pathto_TLM = path +'/em_tlm'
subprocess.run(['./mpirun_wrf.csh',path+'/em_tlm','12'])
print("Python: TLM DONE")
# STEP 2: Multiply the resulting vector with matrix E 
nc_tlm_final = path+'/em_tlm/tlout_d01_'+et_date_string #............................. file name of TLM output at final time
ncf_tlm_final = Dataset(nc_tlm_final, 'r') #......................................... file-handle for TLM output at final time
# Apply weighting and local projection operator to TLM u,v,T
u_tlm = getvar(ncf_tlm_final,"G_U",meta=False)*(w)*LPO_u[:,:,:] #.................... u-field from TLM, multiplied by E (with LPO)
v_tlm = getvar(ncf_tlm_final,"G_V",meta=False)*(w)*LPO_v[:,:,:] #.................... v-field from TLM, multiplied by E (with LPO)
T_tlm = getvar(ncf_tlm_final,"G_T",meta=False)*(wT)*LPO_T[:,:,:] #................... T-field from TLM, multiplied by E (with LPO)
#STEP 3: Create the initial condition file (final_sens_d01) to run the WRF_ADJ
#path_fwd_final=path+'/em_fwd/wrfout_d01_'+et_date_string #............................ full-path file name for WRF output at final time
#nc_adj_path = path+'/em_adj' #....................................................... path to WRF adjoint directory
# copy path_fwd_final (file) to nc_adj_path (subdirectory)
#shutil.copy(path_fwd_final,nc_adj_path) 
# move to em_ad subdirectory
os.chdir(path+'/em_adj') 
# rename to final_sens_d01
#os.rename('wrfout_d01_'+et_date_string, 'final_sens_d01') 
nc_adj_final =path+'/em_adj/final_sens_d01_'+et_date_string #........................ full-path file name for final_sens_d01
ncf_adj_final = Dataset(nc_adj_final, 'a') #......................................... file-handle for final_sens_d01
# write weighted TLM variables to G_* variables in final_sens_d01
ncf_adj_final.variables['G_U'][0]= u_tlm[:] 
ncf_adj_final.variables['G_V'][0]= v_tlm[:]
ncf_adj_final.variables['G_T'][0]= T_tlm[:]
ncf_adj_final.close()
#
######################################################################################
#
# Run ADJ, generate output for fort.17
#
# STEP 1: Run ADJ (integrate the resulting vector on the ADJ backward from t to to)
print("Python: Run ADJ")
os.chdir(path+'/em_adj') 
subprocess.run(['./mpirun_wrf.csh',path+'/em_adj','12'])
print("Python: ADJ DONE")
# STEP 2: Multiply the resulting vector with matrix E^^(-1/2) and remove lateral 
#         boundaries
nc_adj_initial=path+'/em_adj/gradient_wrfplus_d01_'+st_date_string #................... full-path file name for gradient_wrfplus_d01 at initial time
nc_adj = Dataset(nc_adj_initial, 'r') #............................................... file-handle for gradient_wrfplus_d01
u_adj= getvar(nc_adj,"A_U",meta=False) #.............................................. u-field from ADJ
v_adj= getvar(nc_adj,"A_V",meta=False) #.............................................. v-field from ADJ
T_adj= getvar(nc_adj,"A_T",meta=False) #.............................................. T-field from ADJ
# Apply weighting
xy_u=u_adj[:,:,:]*(w**(-0.5)) #....................................................... weighted u-field
xy_v=v_adj[:,:,:]*(w**(-0.5)) #....................................................... weighted v-field
xy_T=T_adj[:,:,:]*(wT**(-0.5)) #...................................................... weighted T-field
# STEP 3: Concatenation of variables in 1D array
XIC2=np.concatenate([xy_u,xy_v,xy_T], axis=None) #.................................... state vector for output
print(len(XIC2))
#
######################################################################################
#
# Left/Right testing and output
# 
""" SV Testing """
#XIC=x'
#XIC2=(E^^(-1/2))L*EL(E^^(-1/2))x' where L=TLM, L*=ADJ
#A=E^^(1/2)LE^^(-1/2)x'
# Extract fields from TLM for left portion of left/right test
a1=getvar(ncf_tlm_final,"G_U",meta=False)*(w**(0.5))*LPO_u[:,:,:] #.................... u-field from TLM (weighted, LPO)
a2=getvar(ncf_tlm_final,"G_V",meta=False)*(w**(0.5))*LPO_v[:,:,:] #.................... v-field from TLM (weighted, LPO)
a3=getvar(ncf_tlm_final,"G_T",meta=False)*(wT**(0.5))*LPO_T[:,:,:] #................... T-field from TLM (weighted, LPO)
# concatenate into a single state vector array
A=np.concatenate([a1,a2,a3], axis=None) #............................................. state vector for 'left'
# compute left and right
right=np.dot(XIC2,XIC) #.............................................................. 'right'
left=np.dot(A,A) #.................................................................... 'left'
# report left/right values
print('right values:',right)
print('left values:',left)
# if right vales = leftf values, it means it is symmetric
# Create fort.17 file and save XIC2 vector
# Write XIC2 in fort.17 
# move to wrf_sv subdirectory
#os.chdir('/Users/karimar/wrf_sv') 
f17 = FortranFile(path+'/fort.17', 'w') #............................................. file-handle for fort.17
f17.write_record(np.array(XIC2, dtype=np.float64))
f17.close()
#
print("Python:DONE")
#
######################################################################################
#
# DONE
#
