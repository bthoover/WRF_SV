######################################################################
#
# PYTHON 3 PROGRAM
#
# Collects SVs (u,v,T) from post-processing and computes a scaled
# linear combination of all SVs (evenly weighted), such that the
# maximum perturbation on each level will (at most) equal a provided
# threshold value. 
#
# Outputs perturbation to a separate file.
#
######################################################################
#
# Import necessary modules
#
import numpy as np #.................................................. array module
from netCDF4 import Dataset #......................................... netCDF module
import wrf #.......................................................... WRF module
from scipy.interpolate import interp1d #.............................. interpolation module
from glob import glob #............................................... global (search) function
#
######################################################################
#
# Take inputs from user
#
SV_filedir = input() #................................................ directory containing individual SV files (ends in '/')
ini_yyyy = input() #.................................................. initialization year
ini_mm = input() #.................................................... initialization month
ini_dd = input() #.................................................... initialization day
ini_hh = input() #.................................................... initialization hour
max_u_str = input() #................................................. zonal flow threshold  applied to calibration of SVs (string)
max_v_str = input() #................................................. merid flow threshold  applied to calibration of SVs (string)
max_t_str = input() #................................................. temperature threshold  applied to calibration of SVs (string)
# Compute initialization date-string
ini_datestr = ini_yyyy+'-'+ini_mm+'-'+ini_dd+'_'+ini_hh+':00:00' #.... initialization date-string
# Compute thresholds as float
max_u = float(max_u_str) #............................................ zonal flow threshold applied to calibration of SVs (float)
max_v = float(max_v_str) #............................................ merid flow threshold applied to calibration of SVs (float)
max_t = float(max_t_str) #............................................ temperature threshold applied to calibration of SVs (float)
#
######################################################################
#
# Collect individual SV information (u,v,T)
#
SV_file_list = sorted( #.............................................. sorted list of all SV files
                      glob(SV_filedir+'init_pert_d01_J*_'+ini_datestr)
                     )
n_SV = len(SV_file_list) #............................................ number of individual SVs
# Read first SV file, extract full model pressure and SVs
SV_file=SV_file_list[0] #............................................. first SV file
SV_hdl=Dataset(SV_file) #............................................. SV file-handle
pre=np.asarray(wrf.getvar(SV_hdl,'pressure')).squeeze() #............. WRF pressure [ny,nx] (hPa)
nz,ny,nx = np.shape(pre) #............................................ dimensions of WRF grid
# Initialize SV arrays
u_svs = np.nan*np.ones((nz,ny,nx,n_SV)) #............................. zonal flow SVs [nz,ny,nx,n_SV] (initialized to NaN)
u_stag_svs = np.nan*np.ones((nz,ny,nx+1,n_SV)) #...................... zonal flow SVs on staggered-grid [nz,ny,nx+1,n_SV] (initialized to NaN)
v_svs = np.nan*np.ones((nz,ny,nx,n_SV)) #............................. merid flow SVs [nz,ny,nx,n_SV] (initialized to NaN)
v_stag_svs = np.nan*np.ones((nz,ny+1,nx,n_SV)) #...................... merid flow SVs on staggered-grid [nz,ny+1,nx,n_SV] (initialized to NaN)
t_svs = np.nan*np.ones((nz,ny,nx,n_SV)) #............................. temperature SVs [nz,ny,nx,n_SV] (initialized to NaN)
# Extract SVs from first file
ux=np.asarray(SV_hdl.variables['G_U']).squeeze() #.................... zonal flow SV [nz,ny,nx+1]
vx=np.asarray(SV_hdl.variables['G_V']).squeeze() #.................... merid flow SV [nz,ny+1,nx]
tx=np.asarray(SV_hdl.variables['G_T']).squeeze() #.................... temperature SV [nz,ny,nx+1]
# Commit ux, vx to staggered SV arrays
u_stag_svs[:,:,:,0] = ux
v_stag_svs[:,:,:,0] = vx
# Destagger u, v
ux=wrf.destagger(ux,stagger_dim=2)
vx=wrf.destagger(vx,stagger_dim=1)
# Commit ux, vx, tx to SV arrays
u_svs[:,:,:,0] = ux
v_svs[:,:,:,0] = vx
t_svs[:,:,:,0] = tx
# Close file-handle and remove file from list
SV_hdl.close()
SV_file_list.remove(SV_file)
# Loop through remaining files, commit ux, vx, tx to SV arrays
for i in range(len(SV_file_list)):
    SV_file=SV_file_list[i]
    SV_hdl=Dataset(SV_file)
    ux=np.asarray(SV_hdl.variables['G_U']).squeeze()
    vx=np.asarray(SV_hdl.variables['G_V']).squeeze()
    tx=np.asarray(SV_hdl.variables['G_T']).squeeze()
    u_stag_svs[:,:,:,i+1] = ux
    v_stag_svs[:,:,:,i+1] = vx
    ux=wrf.destagger(ux,stagger_dim=2)
    vx=wrf.destagger(vx,stagger_dim=1)
    u_svs[:,:,:,i+1] = ux
    v_svs[:,:,:,i+1] = vx
    t_svs[:,:,:,i+1] = tx
    SV_hdl.close()
#
######################################################################
#
# Compute uncalibrated total-SV as linear combination of each SV
#
u_sv_uncal = np.sum(u_svs,axis=3).squeeze() #......................... uncalibrated zonal flow SV-perturbation [nz,ny,nx]
v_sv_uncal = np.sum(v_svs,axis=3).squeeze() #......................... uncalibrated merid flow SV-perturbation [nz,ny,nx]
t_sv_uncal = np.sum(t_svs,axis=3).squeeze() #......................... uncalibrated temperature SV-perturbation [nz,ny,nx]
u_stag_sv_uncal = np.sum(u_stag_svs,axis=3).squeeze() #............... uncalibrated zonal flow SV-perturbation on staggered grid [nz,ny,nx+1]
v_stag_sv_uncal = np.sum(v_stag_svs,axis=3).squeeze() #............... uncalibrated merid flow SV-perturbation on staggered grid [nz,ny+1,nx]
# Calculate (absolute value) ratio of max threshold to uncalibrated
# SV perturbation
u_pe_rat = np.abs(max_u*np.power(u_sv_uncal,-1.)) #................... zonal flow thresh/sv ratio [nz,ny,nx]
v_pe_rat = np.abs(max_v*np.power(v_sv_uncal,-1.)) #................... merid flow thresh/sv ratio [nz,ny,nx]
t_pe_rat = np.abs(max_t*np.power(t_sv_uncal,-1.)) #................... temperature thresh/sv ratio [nz,ny,nx]
# Find minimum ratio to define SV_scalar
SV_scalar = min([ #................................................... scalar multiplier to calibrate SV perturbation
                  np.min(u_pe_rat) , 
                  np.min(v_pe_rat) , 
                  np.min(t_pe_rat)
                ])
#
######################################################################
#
# Compute calibrated total-SV by applying SV_scalar to uncalibrated
# total-SV (on staggered grids where appropriate)
u_stag_sv_cal = SV_scalar * u_stag_sv_uncal #......................... calibrated zonal flow SV-perturbation on staggered grid [nz,ny,nx+1]
v_stag_sv_cal = SV_scalar * v_stag_sv_uncal #......................... calibrated merid flow SV-perturbation on staggered grid [nz,ny+1,nx]
t_sv_cal = SV_scalar * t_sv_uncal #................................... calibrated temperature SV-perturbation [nz,ny,nx]
#
######################################################################
#
# Write calibrated total-SV to file
#
nc_out = Dataset( #...................................................... Dataset object for output
                  SV_filedir+'SV_pert.nc'  , # Dataset input: Output file name
                  "w"              , # Dataset input: Make file write-able
                  format="NETCDF4" , # Dataset input: Set output format to netCDF4
                )
# Dimensions
lat_vstag  = nc_out.createDimension( #................................... Output dimension
                                     "lat_vstag" , # nc_out.createDimension input: Dimension name
                                      None         # nc_out.createDimension input: Dimension size limit ("None" == unlimited)
                                   )
lat  = nc_out.createDimension( #......................................... Output dimension
                               "lat" , # nc_out.createDimension input: Dimension name
                                None   # nc_out.createDimension input: Dimension size limit ("None" == unlimited)
                             )
lon_ustag  = nc_out.createDimension( #................................... Output dimension
                                     "lon_ustag" , # nc_out.createDimension input: Dimension name
                                      None         # nc_out.createDimension input: Dimension size limit ("None" == unlimited)
                                   )
lon  = nc_out.createDimension( #......................................... Output dimension
                               "lon" , # nc_out.createDimension input: Dimension name
                               None    # nc_out.createDimension input: Dimension size limit ("None" == unlimited)
                             )
lev = nc_out.createDimension( #.......................................... Output dimension
                               "lev" , # nc_out.createDimension input: Dimension name
                               None    # nc_out.createDimension input: Dimension size limit ("None" == unlimited)
                             )
param = nc_out.createDimension( #........................................ Output dimension
                               "param" , # nc_out.createDimension input: Dimension name
                               None      # nc_out.createDimension input: Dimension size limit ("None" == unlimited)
                             )
# Variables
SV_PERT_U = nc_out.createVariable( #.................................... Output variable
                                    "SV_PERT_U"  , # nc_out.createVariable input: Variable name
                                    "f8"          , # nc_out.createVariable input: Variable format
                                    (
                                      "lev"       , # nc_out.createVariable input: Variable dimension
                                      "lat"       , # nc_out.createVariable input: Variable dimension
                                      "lon_ustag"   # nc_out.createVariable input: Variable dimension
                                    )
                                  )
SV_PERT_V = nc_out.createVariable( #.................................... Output variable
                                    "SV_PERT_V"  , # nc_out.createVariable input: Variable name
                                    "f8"          , # nc_out.createVariable input: Variable format
                                    (
                                      "lev"       , # nc_out.createVariable input: Variable dimension
                                      "lat_vstag" , # nc_out.createVariable input: Variable dimension
                                      "lon"         # nc_out.createVariable input: Variable dimension
                                    )
                                  )
SV_PERT_T = nc_out.createVariable( #.................................... Output variable
                                    "SV_PERT_T"  , # nc_out.createVariable input: Variable name
                                    "f8"          , # nc_out.createVariable input: Variable format
                                    (
                                      "lev"       , # nc_out.createVariable input: Variable dimension
                                      "lat"       , # nc_out.createVariable input: Variable dimension
                                      "lon"         # nc_out.createVariable input: Variable dimension
                                    )
                                  )
UMAX = nc_out.createVariable( #......................................... Output variable
                               "UMAX"   , # nc_out.createVariable input: Variable name
                               "f8"      , # nc_out.createVariable input: Variable format
                               (
                                 "param"   # nc_out.createVariable input: Variable dimension
                               )
                             )
VMAX = nc_out.createVariable( #......................................... Output variable
                               "VMAX"   , # nc_out.createVariable input: Variable name
                               "f8"      , # nc_out.createVariable input: Variable format
                               (
                                 "param"   # nc_out.createVariable input: Variable dimension
                               )
                             )
TMAX = nc_out.createVariable( #......................................... Output variable
                               "TMAX"   , # nc_out.createVariable input: Variable name
                               "f8"      , # nc_out.createVariable input: Variable format
                               (
                                 "param"   # nc_out.createVariable input: Variable dimension
                               )
                             )
# Fill netCDF arrays via slicing
SV_PERT_U[:,:,:] = u_stag_sv_cal
SV_PERT_V[:,:,:] = v_stag_sv_cal
SV_PERT_T[:,:,:] = t_sv_cal
UMAX[:] = max_u
VMAX[:] = max_v
TMAX[:] = max_t
# Close netCDF file
nc_out.close()
#
#########################################################################
#
# END
#
