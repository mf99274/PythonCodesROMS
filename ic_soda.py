import os
import time
import pyroms
import netCDF4 as nc
import numpy as np
from gridfill import fill
import mpl_toolkits.basemap as mp
from mpl_toolkits.basemap import Basemap
from netcdftime import num2date, date2num, datetime

# name of the simulation
simName = 'Baja'

# Select Year, Month, day, and hour of the input data
year = '1990'
month = '12'
day = '31'

time_units='days since 1900-01-01 00:00:00'

# Path to my input data
inputPath = '/media/mfocean/easystore/SODA_5days/'
os.chdir(inputPath)

# Path to coarse grid
gridPath = '/home/mfocean/Documents/2-PhD/98-models/2-ROMS/0-Grids/'
gridName = 'southernCCS_4km.nc'
gridData = nc.Dataset(gridPath+gridName)

ICPath = '/home/mfocean/Documents/2-PhD/98-models/2-ROMS/3-IC/'
ICname = 'IC_soda.nc'

kw = dict(eps=1e-4, relax=0.6, itermax=1e4, initzonal=False,
          cyclic=True, verbose=True)

#--------------- Don't need to change things after here ---------------#
print('Getting values from your coarse grid')
h = gridData.variables['h'][:]
cs_r = gridData.variables['Cs_r'][:]
cs_w=gridData.variables['Cs_w'][:] 
lon_rho = gridData.variables['lon_rho'][:]
lat_rho = gridData.variables['lat_rho'][:]
lon_u = gridData.variables['lon_u'][:]
lat_u = gridData.variables['lat_u'][:]
lon_v = gridData.variables['lon_v'][:]
lat_v = gridData.variables['lat_v'][:]
mask_rho = gridData.variables['mask_rho'][:]
mask_u = gridData.variables['mask_u'][:]
mask_v = gridData.variables['mask_v'][:]

# creating mask_u and mask_v in 3d to use later
mask_u_3d = np.ones((cs_r.shape[0],mask_u.shape[0],mask_u.shape[1]))
mask_u_3d[:,:,:] = mask_u[np.newaxis,:,:] == 1. 
mask_v_3d = np.ones((cs_r.shape[0],mask_v.shape[0],mask_v.shape[1]))
mask_v_3d[:,:,:] = mask_v[np.newaxis,:,:] == 1. 

print('Opening file based on the selection above...')
inputFileName ='soda3.3.2_'+year+'_'+month+'_'+day+'.nc'
inputData = nc.Dataset(inputFileName)

xt_ocean = inputData.variables['xt_ocean'][:]
yt_ocean = inputData.variables['yt_ocean'][:]
st_ocean = inputData.variables['st_ocean'][:]
xu_ocean = inputData.variables['xu_ocean'][:]
yu_ocean = inputData.variables['yu_ocean'][:]
timeInput = inputData.variables['time']

# getting the bounds so the data is smaller
lon_bound = np.logical_and(xt_ocean>=lon_rho.min()-0.5,xt_ocean<=lon_rho.max()+0.5)
lat_bound = np.logical_and(yt_ocean>=lat_rho.min()-0.5,yt_ocean<=lat_rho.max()+0.5)
lonu_bound = np.logical_and(xu_ocean>=lon_rho.min()-0.5,xu_ocean<=lon_rho.max()+0.5)
latu_bound = np.logical_and(yu_ocean>=lat_rho.min()-0.5,yu_ocean<=lat_rho.max()+0.5)
depth_bound = np.logical_and(st_ocean>=h.min()-5,st_ocean<=h.max()+100)

# new grid values
lonSoda = xt_ocean[lon_bound]
latSoda = yt_ocean[lat_bound]
lonuSoda = xu_ocean[lonu_bound]
latuSoda = yu_ocean[latu_bound]
depthSoda = st_ocean[depth_bound]

print('Loading data from input file')
ssh = inputData.variables['ssh']
salinity = inputData.variables['salt'][0,depth_bound,lat_bound,lon_bound]
temperature = inputData.variables['temp'][0,depth_bound,lat_bound,lon_bound]
u = inputData.variables['u'][0,depth_bound,latu_bound,lonu_bound]
v = inputData.variables['v'][0,depth_bound,latu_bound,lonu_bound]

print('Creating arrays for sea surface height')
dst_ssh = np.zeros((lon_rho.shape[0],lat_rho.shape[1]))

print('interpolating ssh')
ssh1,converged = fill(ssh[0,lat_bound,lon_bound], 0, 1, **kw)
ssh2=mp.interp(ssh1,lonSoda,latSoda,lon_rho,lat_rho,checkbounds=False,\
              masked=False,order=1)

ssh2 = ssh2 * mask_rho
dst_ssh = np.ma.MaskedArray(ssh2,mask_rho==0.)

dst_sal = np.zeros((cs_r.shape[0],lon_rho.shape[0],lat_rho.shape[1]))
dst_temp = np.ones((cs_r.shape[0],lon_rho.shape[0],lat_rho.shape[1]))
dst_u = np.ones((cs_r.shape[0],lon_u.shape[0],lat_u.shape[1]))
dst_v = np.ones((cs_r.shape[0],lon_v.shape[0],lat_v.shape[1]))

sal_arr = []
temp_arr = []
u_arr=[]
v_arr=[]

for k in range(depthSoda.shape[0]):
    "Modified SODA2ROMS from Trond Kristiansen"
    "Date of modification: Nov,13th 2015 "
    "Author: Matheus Fagundes"

    salinity1 = salinity[k,:,:]
    temperature1 = temperature[k,:,:]
    u1 = u[k,:,:]
    v1 = v[k,:,:]

    salty,converged = fill(salinity1,0,1,**kw)
    tempy,converged = fill(temperature1,0,1,**kw)
    uy,converged = fill(u1,0,1,**kw)
    vy,converged = fill(v1,0,1,**kw)

    salty = mp.interp(salty,lonSoda,latSoda,lon_rho,lat_rho,\
                 checkbounds=False,masked=False,order=1)
    tempy = mp.interp(tempy,lonSoda,latSoda,lon_rho,lat_rho,\
                 checkbounds=False,masked=False,order=1)
    uy = mp.interp(uy,lonuSoda,latuSoda,lon_u,lat_u,\
                 checkbounds=False,masked=False,order=1) 
    vy = mp.interp(vy,lonuSoda,latuSoda,lon_v,lat_v,\
                 checkbounds=False,masked=False,order=1)

    temp_arr.append(tempy)
    sal_arr.append(salty)
    u_arr.append(uy)
    v_arr.append(vy)

print('Converting list of temp and salinity in numpy array')
print('and interpolating vertically')
sal=np.array(sal_arr)*mask_rho
temp=np.array(temp_arr)*mask_rho

for x in range(lon_rho.shape[0]):
    for y in range(lat_rho.shape[1]):
        z_new = cs_r * h[x,y]
        dst_sal[:,x,y] = np.interp(z_new,-depthSoda[::-1],sal[::-1,x,y])
        dst_temp[:,x,y] = np.interp(z_new,-depthSoda[::-1],temp[::-1,x,y])
   
dst_sal = np.ma.masked_where(dst_sal==0.,dst_sal)
dst_temp = np.ma.masked_where(dst_temp==0.,dst_temp)

print('Converting list of u and v in numpy array')
print('and interpolating vertically')
uu = np.array(u_arr)*mask_u
vv = np.array(v_arr)*mask_v

for x in range(lon_u.shape[0]):
    for y in range(lat_u.shape[1]):
        z_new = cs_r * h[x,y]
        dst_u[:,x,y] = np.ma.masked_where(mask_u[x,y]==0.,np.interp(z_new,-depthSoda[::-1],uu[::-1,x,y]))

for x in range(lon_v.shape[0]):
    for y in range(lat_v.shape[1]):
        z_new = cs_r * h[x,y]
        dst_v[:,x,y] = np.interp(z_new,-depthSoda[::-1],vv[::-1,x,y])

dst_u1 = np.ma.masked_where(mask_u_3d==0.,dst_u)
dst_v1 = np.ma.masked_where(mask_v_3d==0.,dst_v)

######################### Computing Ubar and Vbar
dst_grd = pyroms.grid.get_ROMS_grid(simName)

######################### Ubar
z_u = 0.5 * (dst_grd.vgrid.z_w[0,:,:,:-1] + dst_grd.vgrid.z_w[0,:,:,1:])
dst_ubar = np.zeros((dst_u1.shape[1], dst_u1.shape[2]))

for i in range(dst_ubar.shape[1]):
    for j in range(dst_ubar.shape[0]):
        dst_ubar[j,i] = (dst_u1[:,j,i] * np.diff(z_u[:,j,i])).sum() / -z_u[0,j,i]

dst_ubar=np.ma.masked_invalid(dst_ubar)

######################### Vbar
z_v = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,:] + dst_grd.vgrid.z_w[0,:,1:,:])
dst_vbar = np.zeros((dst_v1.shape[1], dst_v1.shape[2]))

for i in range(dst_vbar.shape[1]):
    for j in range(dst_vbar.shape[0]):
        dst_vbar[j,i] = (dst_v1[:,j,i] * np.diff(z_v[:,j,i])).sum() / -z_v[0,j,i]

dst_vbar=np.ma.masked_invalid(dst_vbar)

######################### Converting time
tt = num2date(timeInput[0],units=timeInput.units,calendar=timeInput.calendar)
ocean_time = date2num(tt,units=time_units,calendar='julian')

######################### Creating netcdf ic file
fill_value = dst_sal.fill_value

# Create netCDF file

#dimensions:
#	xi_u = 904 ;
#	xi_v = 905 ;
#	xi_rho = 905 ;
#	eta_u = 323 ;
#	eta_v = 322 ;
#	eta_rho = 323 ;
#	s_rho = 30 ;
#	ocean_time = 1 ;
#	tracer = 2 ;
#variables:
#	double ocean_time(ocean_time) ;
#		ocean_time:long_name = "days since 1900-01-01 00:00:00" ;
#		ocean_time:units = "days" ;
#	double temp(ocean_time, s_rho, eta_rho, xi_rho) ;
#		temp:_FillValue = NaN ;
#		temp:long_name = "potential temperature" ;
#		temp:units = "Celsius" ;
#		temp:field = "temperature, scalar, series" ;
#	double sal(ocean_time, s_rho, eta_rho, xi_rho) ;
#		sal:_FillValue = NaN ;
#		sal:long_name = "salinity" ;
#		sal:units = "PSU" ;
#		sal:field = "salinity, scalar, series" ;
#	double zeta(ocean_time, eta_rho, xi_rho) ;
#		zeta:_FillValue = NaN ;
#		zeta:long_name = "free-surface" ;
#		zeta:units = "meter" ;
#		zeta:field = "free-surface, scalar, series" ;
#	double u(ocean_time, s_rho, eta_u, xi_u) ;
#		u:_FillValue = NaN ;
#		u:long_name = "3D u-momentum component" ;
#		u:units = "meter second-1" ;
#		u:field = "u-velocity, scalar, series" ;
#	double v(ocean_time, s_rho, eta_v, xi_v) ;
#		v:_FillValue = NaN ;
#		v:long_name = "3D u-momentum component" ;
#		v:units = "meter second-1" ;
#		v:field = "v-velocity, scalar, series" ;
#	double ubar(ocean_time, eta_u, xi_u) ;
#		ubar:_FillValue = NaN ;
#		ubar:long_name = "2D u-momentum component" ;
#		ubar:units = "meter second-1" ;
#		ubar:field = "ubar-velocity,, scalar, series" ;
#	double vbar(ocean_time, eta_v, xi_v) ;
#		vbar:_FillValue = NaN ;
#		vbar:long_name = "2D v-momentum component" ;
#		vbar:units = "meter second-1" ;
#		vbar:field = "vbar-velocity,, scalar, series" ;
#
#// global attributes:
#		:type = "Initial file ROMS" ;
#		:author = "Matheus Fagundes" ;
#		:history = "CreatedTue Mar  1 19:28:23 2016" ;

ncd = nc.Dataset(ICPath+ICname, 'w', format='NETCDF4')    

ncd.createDimension('xi_u', mask_u.shape[1])
ncd.createDimension('xi_v', mask_v.shape[1])
ncd.createDimension('xi_rho', mask_rho.shape[1])
ncd.createDimension('eta_u', mask_u.shape[0])
ncd.createDimension('eta_v', mask_v.shape[0])
ncd.createDimension('eta_rho', mask_rho.shape[0])
ncd.createDimension('s_rho', cs_r.shape[0])
ncd.createDimension('s_w', cs_w.shape[0])
ncd.createDimension('time', 1)
ncd.createDimension('tracer', 2)

# Ocean_time
ncd.createVariable('ocean_time', 'f8', ('time', ))
ncd.variables['ocean_time'].long_name = "time since initialization"
ncd.variables['ocean_time'].units = time_units
ncd.variables['ocean_time'][:] = ocean_time

# temperature
ncd.createVariable('temp', 'f8', ('time', 's_rho','eta_rho','xi_rho'),zlib=True)
ncd.variables['temp'].long_name = 'potential temperature'
ncd.variables['temp'].units = 'Celsius'
ncd.variables['temp'].field = 'temperature, scalar, series'
ncd.variables['temp'][:,:,:,:] = dst_temp

#Salinity
ncd.createVariable('salt', 'f8', ('time', 's_rho','eta_rho','xi_rho'),zlib=True)
ncd.variables['salt'].long_name = 'salinity'
ncd.variables['salt'].units = 'PSU'
ncd.variables['salt'].field = 'salinity, scalar, series'
ncd.variables['salt'][:,:,:,:] = dst_sal

#SSH
ncd.createVariable('zeta', 'f8', ('time','eta_rho','xi_rho'),zlib=True, fill_value=fill_value)
ncd.variables['zeta'].long_name = 'free-surface'
ncd.variables['zeta'].units = 'meter'
ncd.variables['zeta'].field = 'free-surface, scalar, series'
ncd.variables['zeta'][:,:,:] = dst_ssh

#Zonal velocity
ncd.createVariable('u', 'f8', ('time', 's_rho', 'eta_u', 'xi_u'),zlib=True, fill_value=fill_value)
ncd.variables['u'].long_name = '3D u-momentum component'
ncd.variables['u'].units = 'meter second-1'
ncd.variables['u'].field = 'u-velocity, scalar, series'
ncd.variables['u'][:,:,:,:] = dst_u1

#Meridional velocity
ncd.createVariable('v', 'f8', ('time', 's_rho', 'eta_v', 'xi_v'),zlib=True,fill_value=fill_value)
ncd.variables['v'].long_name = '3D u-momentum component'
ncd.variables['v'].units = 'meter second-1'
ncd.variables['v'].field = 'v-velocity, scalar, series'
ncd.variables['v'][:,:,:,:] = dst_v1

#UBAR
ncd.createVariable('ubar', 'f8', ('time', 'eta_u', 'xi_u'),zlib=True,fill_value=fill_value)
ncd.variables['ubar'].long_name = '2D u-momentum component'
ncd.variables['ubar'].units = 'meter second-1'
ncd.variables['ubar'].field = 'ubar-velocity,, scalar, series'
ncd.variables['ubar'][:,:,:] = dst_ubar

#VBAR
ncd.createVariable('vbar', 'f8', ('time', 'eta_v', 'xi_v'),zlib=True,fill_value=fill_value)
ncd.variables['vbar'].long_name = '2D v-momentum component'
ncd.variables['vbar'].units = 'meter second-1'
ncd.variables['vbar'].field = 'vbar-velocity,, scalar, series'
ncd.variables['vbar'][:,:,:] = dst_vbar

ncd.type='Initial file ROMS'
ncd.author='Cobialab'
ncd.history='Created '+time.ctime(time.time())
ncd.title='Mexico'

ncd.close()



