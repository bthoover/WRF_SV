#import atx93_klm
########################################################################
#
# PYTHON 3 PROGRAM
#
# This program will unpack individual fort.## files representing 
# different SVs, produced by running read_vectors.exe operating on 
# vector.dat.
#
#########################################################################
# Load necessary modules
#
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
#
#########################################################################
# Collect user inputs
#
postproc_dir = input() #................................................. post-processing directory (full-path, ends in '/')
dur_hh = input() #....................................................... duration of optimal time interval (OTI)
ini_yyyy = input() #..................................................... initialization year (yyyy)
ini_mm = input() #....................................................... initialization month (mm)
ini_dd = input() #....................................................... initialization day (dd)
ini_hh = input() #....................................................... initialization hour (hh)
update_bdy_hrs = input() #............................................... update interval of boundary conditions (hrs)
wrf_tstep_sec = input() #................................................ WRF model timestep (seconds)
J_num = input() #........................................................ integer-number of singular vector to unpack
print('Post Processing: STARTS')
# Digitize user inputs (input as strings), compute some more variables
duration = int(dur_hh) #................................................. duration of the OTI (integer)
st = dt.datetime.strptime(ini_yyyy+ini_mm+ini_dd+ini_hh,'%Y%m%d%H') #.... initialization datetime
boundary_condition_update = int(update_bdy_hrs) #........................ update interval of boundary conditions (integer, hours)
time_step_seconds = int(wrf_tstep_sec) #................................. WRF model timestep (integer, seconds)
oti = dt.timedelta(hours=duration) #..................................... timedelta of OTI
et = st+oti #............................................................ ending datetime
st_date_string = dt.datetime.strftime(st,'%Y-%m-%d_%H:%M:%S') #.......... starting date-string (WRF-file formatted)
et_date_string = dt.datetime.strftime(et,'%Y-%m-%d_%H:%M:%S') #.......... ending date-string (WRF-file formatted)
print(st_date_string)
print(et_date_string)
J = int(J_num) #......................................................... integer-number of singular vector to unpack (integer)
#
#########################################################################
#
# Post-process vector J from vector.dat
#
"""POST PROCESS  vector.dat"""
print("Post Process vector.dat")   
"""This code reads one of the eigenvectors the output (units 30+j, j = 1, N) of program read_vector.f . The user must specify which eigenvector to post-process to create an initial condition file."""
# Read vector from SV-J file
path = postproc_dir #.................................................... path to directory containing necessary files
f    = FortranFile(path+'fort.3'+str(J), 'r' ) #......................... fortran file for SV-J
XIC  = f.read_reals( dtype='float64' ) #................................. vector from f
f.close
print('data read, vector of length:', XIC.size)
# Extract dimension-sizes from WRF output file for initialization time
nc_fwd_fname = path+'wrfout_d01_'+st_date_string #....................... WRF file for initialization date
ncf_fwd      = Dataset(nc_fwd_fname, 'r') #.............................. WRF file-handle
nc_dims = ncf_fwd.dimensions #........................................... WRF dimensions library
#
nz = nc_dims['bottom_top'].size #........................................ WRF (unstaggered) vertical dimension
ny = nc_dims['south_north'].size #....................................... WRF (unstaggered) merid dimension
nx = nc_dims['west_east'].size #......................................... WRF (unstaggered) zonal dimension
print(nz,nx,ny)
ncf_fwd.close()
# Initialize SV-J grids [nz,ny,nx]
U = np.zeros([nz,ny,nx+1]) #............................................. SV-J zonal flow (W-E staggered) [nz,ny,nx+1]
V = np.zeros([nz,ny+1,nx]) #............................................. SV-J merid flow (S-N staggered) [nz,ny+1,nx]
T = np.zeros([nz,ny,nx]) #............................................... SV-J temperature (unstaggered) [nz,ny,nx]
# Define start/finish indices of each variable in XIC
us=0 #................................................................... SV-J start-index of zonal flow
uf=nz*ny*(nx+1) #........................................................ SV-J finish-index of zonal flow
vs=uf #.................................................................. SV-J start-index of merid flow
vf=vs+nz*(ny+1)*nx #..................................................... SV-J finish-index of merid flow
ts = vf #................................................................ SV-J start-index of temperature
tf=ts+nz*ny*nx #......................................................... SV-J finish-index of temperature
# Define SV weighting coefficients
u_wght = 0.5 #........................................................... zonal flow SV weighting coefficient
v_wght = 0.5 #........................................................... merid flow SV weighting coefficient
t_wght = 0.5*(9.8*9.8)/(300.e-2)**2. #................................... temperature SV weighting coefficient
print(u_wght, v_wght, t_wght, ' <<<<< WEIGHTINGS')
# Extract and reshape SV-J U, V, T fields (non-normalized)
uu = np.reshape(XIC[us:uf],(nz,ny,nx+1)) #............................... SV-J U field (non-normalized)
vv = np.reshape(XIC[vs:vf],(nz,ny+1,nx)) #............................... SV-J V field (non-normalized)
tt = np.reshape(XIC[ts:tf],(nz,ny,nx)) #................................. SV-J T field (non-normalized)
# Create a field normalized by the square-root of the weighing coeff.
G_U = uu/np.sqrt(u_wght,dtype=np.float64) #.............................. SV-J U field (normalized)
G_V = vv/np.sqrt(v_wght,dtype=np.float64) #.............................. SV-J V field (normalized)
G_T = tt/np.sqrt(t_wght,dtype=np.float64) #.............................. SV-J T field (normalized)
#########################################################################
#
# Compute energy and normalized energy
#
energy = ( #............................................................. total energy of SV-J (non-normalized)
           np.sum(uu*uu,dtype=np.float64) + 
           np.sum(vv*vv,dtype=np.float64) + 
           np.sum(tt*tt,dtype=np.float64)
         )
# For a non-null vector, normalize fields by square-root of energy
print(' for J = ', J, 'energy = ', energy)
if(energy!=0.0):
    G_U = G_U/np.sqrt(energy,dtype=np.float64)
    G_V = G_V/np.sqrt(energy,dtype=np.float64)
    G_T = G_T/np.sqrt(energy,dtype=np.float64)
else:
    print('SV', J, 'is a null vector. Check it')
normalized_energy = ( #.................................................. total energy of SV-J (normalized)
                      np.sum(G_U*G_U,dtype=np.float64) + 
                      np.sum(G_V*G_V,dtype=np.float64) + 
                      np.sum(G_T*G_T,dtype=np.float64)
                    )
print('normalized_energy = ', normalized_energy)
#
#########################################################################
#
# Write normalized SV-J initial-fields to file
#
# Rename wrfout_d01_<initial-time> to init_pert_d01_J#_<initial-time>
SV_J_fname = 'init_pert_d01_J'+str(J)+'_'+st_date_string #................... filename for SV-J
os.rename('wrfout_d01_'+st_date_string,SV_J_fname)
# Open SV_J_fname file
nc_fname = path+SV_J_fname #............................................ SV-J file-name
ncf_SV_J = Dataset(nc_fname, 'a') #..................................... SV-J file-handle
# Write normalized fields to file
ncf_SV_J.variables['G_U'][0] = G_U[:]
ncf_SV_J.variables['G_V'][0] = G_V[:]
ncf_SV_J.variables['G_T'][0] = G_T[:]
ncf_SV_J.close()
#
#########################################################################
#
# END PROGRAM
#
