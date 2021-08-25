#########################################################################
#
# PYTHON 3 PROGRAM
#
# Given a WRF-TLM input file (init_pert_d<domain>) and an optimal 
# perturbation file (u,v,t), this program will add the optimal 
# perturbation to the u, v, and t components to the G_U, G_V, and G_T 
# grids and write the updated grids to the WRF-TLM input file.
#
# 05/03/21: Adding compensating moisture term at this point in the code.
#           This is being computed here, rather than at the point where
#           optimal perturbations are computed, because the scaling of
#           the compensating moisture term is not *perfectly* linear with
#           the scaling of the optimal (temperature) perturbation. It is
#           very *nearly* perfectly linear, but in testing there appear
#           to be a minority of gridpoints that deviate strongly from
#           the ratio (say, from a ratio of 2 when P_MAG=0.5). This is
#           probably just goofy stuff showing up around a ptd_q value of
#           zero. But, to be safe, we are computing the entire term here
#           instead of there.
#
#
#########################################################################
#
# Import necessary modules
#
import numpy as np #..................................................... array module
from netCDF4 import Dataset #............................................ netCDF i/o module
import wrf #............................................................. WRF i/o module
from shutil import copyfile #............................................ file copy module
#
#########################################################################
#
# Obtain inputs from user
#
tlm_input_file = input() #............................................... old (pre-perturbed) WRF input file (full-path)
opt_pert_file = input() #................................................ file containing optimal perturbations
#
########################################################################
#
# Open wrfinput and optimal perturbation files and extract u, v, t fields
#
wrfinput_hdl = Dataset(tlm_input_file,'r+') #............................ WRF input file-handle (opened as appendable read)
u0 = wrfinput_hdl.variables['G_U'] #..................................... unperturbed zonal flow [ntim=1,nlev,nlat,nlon]
v0 = wrfinput_hdl.variables['G_V'] #..................................... unperturbed merid flow [ntim=1,nlev,nlat,nlon]
t0 = wrfinput_hdl.variables['G_T'] #..................................... unperturbed temperature [ntim=1,nlev,nlat,nlon]
opt_pert_hdl = Dataset(opt_pert_file) #.................................. optimal perturbation file-handle
up = opt_pert_hdl.variables['SV_PERT_U'] #............................... optimal zonal flow perturbation [nlev,nlat,nlon]
vp = opt_pert_hdl.variables['SV_PERT_V'] #............................... optimal merid flow perturbation [nlev,nlat,nlon]
tp = opt_pert_hdl.variables['SV_PERT_T'] #............................... optimal temperature perturbation [nlev,nlat,nlon]
#
#########################################################################
#
# Extract additional (or unit-specific) fields from wrfinput for 
# computing the compensating moisture term
#
init_p = np.asarray(wrf.getvar(
                               wrfinput_hdl ,
                               'pres'       ,
                               units='hPa'
                              )
                   ).squeeze() #. full model pressure (hPa) [lev,lat,lon]
init_t = np.asarray(wrf.getvar(
                               wrfinput_hdl ,
                               'temp'       ,
                               units='degC'
                              )
                   ).squeeze() # full model temperature (deg C) [lev,lat,lon]
q0 = wrfinput_hdl.variables['QVAPOR'] #.................................. unperturbed water vapor mixing ratio (kg/kg) [lev,lat,lon]
#
# Compute saturation vapor pressure, according to the equation provided here:
#
# https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf
#
es_exp = np.divide(7.5*init_t,237.3+init_t) #............................ exponent of saturation vapor pressure
init_es = 6.11*np.power(10.,es_exp) #.................................... full model saturation vapor pressure (hPa) [lev,lat,lon]
#
# Compute saturation vapor pressure after perturbation
#
es_exp = np.divide(7.5*init_t+tp,237.3+init_t+tp) #...................... exponent of saturation vapor pressure
init_es_p = 6.11*np.power(10.,es_exp) #.................................. full model (perturbed) saturation vapor pressure (hPa) [lev,lat,lon]
#
# Compute saturation mixing ratio, according to the equation provided here:
#
# https://www.weather.gov/media/epz/wxcalc/mixingRatio.pdf
#
# We are immediately computing in (kg/kg) by dividing the coefficient by
# 1000.0 here, i.e. 621.97 --> 0.62197
#
init_ws = 0.62197*np.divide(init_es,init_p-init_es) #.................... full model saturation mixing ratio (kg/kg) [lev,lat,lon]
#
# Compute saturation mixing ratio after perturbation
#
init_ws_p = 0.62197*np.divide(init_es_p,init_p-init_es_p) #.............. full model (perturbed) saturation mixing ratio (kg/kg) [lev,lat,lon]
#
# Compute the necessary mixing ratio to keep relative humidity constant:
#
# (q1/ws1) = (q2/ws2) --> q2 = q1(ws2/ws1)
#
q_p = q0 * np.divide(init_ws_p,init_ws) #................................ full model (perturbed) mixing ratio to maintain unperturbed rh (kg/kg) [lev,lat,lon]
#
# Compute initial mixing ratio perturbation to keep relative humidity
# constant
#
qp = q_p - q0 #.......................................................... NON-optimal mixing ratio perturbation to keep rh constant [lev,lat,lon]
#
#########################################################################
#
# Compute total field: <var>0 + <var>p
#
utot = np.asarray(u0) + np.asarray(up) #................................. perturbed zonal flow [ntim=1,nlev,nlat,nlon]
vtot = np.asarray(v0) + np.asarray(vp) #................................. perturbed merid flow [ntim=1,nlev,nlat,nlon]
ttot = np.asarray(t0) + np.asarray(tp) #................................. perturbed temperature [ntim=1,nlev,nlat,nlon]
# Suppressing compensating moisture term (for now)
#qtot = np.asarray(q0) + np.asarray(qp) #................................. perturbed water vapor mixing ratio [ntim=1,nlev,nlat,nlon]
#
#########################################################################
#
# Write total fields back to input fields via slicing
#
u0[:,:,:,:] = utot
v0[:,:,:,:] = vtot
t0[:,:,:,:] = ttot
# Suppressing compensating moisture term (for now)
#q0[:,:,:,:] = qtot
# Close file-handles
wrfinput_hdl.close()
opt_pert_hdl.close()
#
#########################################################################
#
# END
#
