import os
import glob
import time
import pyroms
import numpy as np
import netCDF4 as nc
from gridfill import fill
import mpl_toolkits.basemap as mp
#from mpl_toolkits.basemap import Basemap
from netcdftime import num2date, date2num, datetime

# name of the simulation
simName = 'Baja'

# Select Year, Month, day, and hour to start bry data
year = '1990'
month = '12'
day = '31'

# Select Year, Month, day, and hour to finish bry data
year1 = '1995'
month1 = '11'
day1 = '05'

time_units='days since 1900-01-01 00:00:00'

# Path to my input data
inputPath = '/media/mfocean/easystore/SODA_5days/'
os.chdir(inputPath)

# Path to coarse grid
gridPath = '/home/mfocean/Documents/2-PhD/98-models/2-ROMS/0-Grids/'
gridName = 'southernCCS_4km.nc'
gridData = nc.Dataset(gridPath+gridName)

BCPath = '/home/mfocean/Documents/2-PhD/98-models/2-ROMS/4-BC/'
BCname = 'BC_soda.nc'

kw = dict(eps=1e-4, relax=0.6, itermax=1e4, initzonal=False,
          cyclic=True, verbose=True)

fill_value = 1e+20

# choosing the boundaries
north = 'yes'
south = 'yes'
west = 'yes'
east = 'no' #east boundary interpolation isn't set as of right now (09/19/2020)

#--------------- Don't need to change things after here ---------------#
################ Loading grid data
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

################# creating mask_u and mask_v in 3d to use later
mask_u_3d = np.ones((cs_r.shape[0],mask_u.shape[0],mask_u.shape[1]))
mask_u_3d[:,:,:] = mask_u[np.newaxis,:,:] == 1. 
mask_v_3d = np.ones((cs_r.shape[0],mask_v.shape[0],mask_v.shape[1]))
mask_v_3d[:,:,:] = mask_v[np.newaxis,:,:] == 1. 

############### 
dst_grd = pyroms.grid.get_ROMS_grid(simName) #Loading values from my grid

print('Getting the number of files to create the bdry condition')

listFiles = glob.glob("soda3.3.2_*.nc")
listFiles.sort()

print('Checking if the beginning of bry and the end are in the path')

for length in range(len(listFiles)):
    if listFiles[length] == 'soda3.3.2_'+year+'_'+month+'_'+day+'.nc':
        print('The file '+ 'soda3.3.2_'+year+'_'+month+'_'+day+'.nc'+' '+'is in the folder')
        idx_ini = length

for length in range(len(listFiles)):
    if listFiles[length] == 'soda3.3.2_'+year1+'_'+month1+'_'+day1+'.nc':
        print('The file '+ 'soda3.3.2_'+year1+'_'+month1+'_'+day1+'.nc'+' '+'is in the folder')
        idx_end = length

print('Opening file based on the selection above...')
inputFileName ='soda3.3.2_'+year+'_'+month+'_'+day+'.nc'
inputData = nc.Dataset(listFiles[idx_ini])

xt_ocean = inputData.variables['xt_ocean'][:]
yt_ocean = inputData.variables['yt_ocean'][:]
st_ocean = inputData.variables['st_ocean'][:]
xu_ocean = inputData.variables['xu_ocean'][:]
yu_ocean = inputData.variables['yu_ocean'][:]


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

print('Opening netcdf in the hardrive...')
############### Initializing netcdf file
ncd = nc.Dataset(BCPath+BCname, 'w', format='NETCDF4')

ncd.createDimension('xi_u', mask_u.shape[1])
ncd.createDimension('xi_v', mask_v.shape[1])
ncd.createDimension('xi_rho', mask_rho.shape[1])
ncd.createDimension('eta_u', mask_u.shape[0])
ncd.createDimension('eta_v', mask_v.shape[0])
ncd.createDimension('eta_rho', mask_rho.shape[0])
ncd.createDimension('s_rho', cs_r.shape[0])
ncd.createDimension('s_w', cs_w.shape[0])
ncd.createDimension('ocean_time', idx_end+1)
ncd.createDimension('tracer', 2)

print('Creating ocean_time variable')
ncd.createVariable('ocean_time', 'f8', ('ocean_time'))
ncd.variables['ocean_time'].long_name = 'open boundary conditions time'
ncd.variables['ocean_time'].units = time_units

if west == 'yes':
    ncd.createVariable('salt_west', 'f8', ('ocean_time', 's_rho', 'eta_rho'),zlib=True)
    ncd.variables['salt_west'].long_name = 'salinity western boundary condition'
    ncd.variables['salt_west'].units = 'PSU'
    ncd.variables['salt_west'].field = 'salt_west, scalar, series'
    ncd.variables['salt_west'].time='ocean_time'
    
    ncd.createVariable('temp_west', 'f8', ('ocean_time', 's_rho', 'eta_rho'),zlib=True)
    ncd.variables['temp_west'].long_name = 'potential temperature western boundary condition'
    ncd.variables['temp_west'].units = 'Celsius'
    ncd.variables['temp_west'].field = 'temp_west, scalar, series'
    ncd.variables['temp_west'].time='ocean_time'

    ncd.createVariable('zeta_west', 'f8', ('ocean_time', 'eta_rho'),zlib=True)
    ncd.variables['zeta_west'].long_name = 'free-surface western boundary condition'
    ncd.variables['zeta_west'].units = 'meter'
    ncd.variables['zeta_west'].field = 'ssh_west, scalar, series'
    ncd.variables['zeta_west'].time='ocean_time'

    ncd.createVariable('u_west', 'f8', ('ocean_time', 's_rho', 'eta_u'),zlib=True)
    ncd.variables['u_west'].long_name = '3D u-momentum western boundary condition'
    ncd.variables['u_west'].units = 'meter second-1'
    ncd.variables['u_west'].field = 'u_west, scalar, series'
    ncd.variables['u_west'].time='ocean_time'

    ncd.createVariable('v_west', 'f8', ('ocean_time', 's_rho', 'eta_v'),zlib=True)
    ncd.variables['v_west'].long_name = '3D u-momentum western boundary condition'
    ncd.variables['v_west'].units = 'meter second-1'
    ncd.variables['v_west'].field = 'v_west, scalar, series'
    ncd.variables['v_west'].time='ocean_time'

    ncd.createVariable('ubar_west', 'f8', ('ocean_time', 'eta_u'),zlib=True)
    ncd.variables['ubar_west'].long_name = '2D u-momentum western boundary condition'
    ncd.variables['ubar_west'].units = 'meter second-1'
    ncd.variables['ubar_west'].field = 'ubar_west, scalar, series'
    ncd.variables['ubar_west'].time='ocean_time'

    ncd.createVariable('vbar_west', 'f8', ('ocean_time', 'eta_v'),zlib=True)
    ncd.variables['vbar_west'].long_name = '2D u-momentum western boundary condition'
    ncd.variables['vbar_west'].units = 'meter second-1'
    ncd.variables['vbar_west'].field = 'vbar_west, scalar, series'
    ncd.variables['vbar_west'].time='ocean_time'

if east == 'yes':
    ncd.createVariable('salt_east', 'f8', ('ocean_time', 's_rho', 'eta_rho'),zlib=True)
    ncd.variables['salt_east'].long_name = 'salinity eastern boundary condition'
    ncd.variables['salt_east'].units = 'PSU'
    ncd.variables['salt_east'].field = 'salt_east, scalar, series'
    ncd.variables['salt_east'].time='ocean_time'

    ncd.createVariable('temp_east', 'f8', ('ocean_time', 's_rho', 'eta_rho'),zlib=True)
    ncd.variables['temp_east'].long_name = 'potential temperature easstern boundary condition'
    ncd.variables['temp_east'].units = 'Celsius'
    ncd.variables['temp_east'].field = 'temp_east, scalar, series'
    ncd.variables['temp_east'].time='ocean_time'

    ncd.createVariable('zeta_east', 'f8', ('ocean_time', 'eta_rho'),zlib=True)
    ncd.variables['zeta_east'].long_name = 'free-surface eastern boundary condition'
    ncd.variables['zeta_east'].units = 'meter'
    ncd.variables['zeta_east'].field = 'zeta_east, scalar, series'
    ncd.variables['zeta_east'].time='ocean_time'

    ncd.createVariable('u_east', 'f8', ('ocean_time', 's_rho', 'eta_u'),zlib=True)
    ncd.variables['u_east'].long_name = '3D u-momentum eastern boundary condition'
    ncd.variables['u_east'].units = 'meter second-1'
    ncd.variables['u_east'].field = 'u_east, scalar, series'
    ncd.variables['u_east'].time='ocean_time'

    ncd.createVariable('v_east', 'f8', ('ocean_time', 's_rho', 'eta_v'),zlib=True)
    ncd.variables['v_east'].long_name = '3D u-momentum eastern boundary condition'
    ncd.variables['v_east'].units = 'meter second-1'
    ncd.variables['v_east'].field = 'v_east, scalar, series'
    ncd.variables['v_east'].time='ocean_time'

    ncd.createVariable('ubar_east', 'f8', ('ocean_time', 'eta_u'),zlib=True)
    ncd.variables['ubar_east'].long_name = '2D u-momentum eastern boundary condition'
    ncd.variables['ubar_east'].units = 'meter second-1'
    ncd.variables['ubar_east'].field = 'ubar_east, scalar, series'
    ncd.variables['ubar_east'].time='ocean_time'

    ncd.createVariable('vbar_east', 'f8', ('ocean_time', 'eta_v'),zlib=True)
    ncd.variables['vbar_east'].long_name = '2D u-momentum eastern boundary condition'
    ncd.variables['vbar_east'].units = 'meter second-1'
    ncd.variables['vbar_east'].field = 'vbar_east, scalar, series'
    ncd.variables['vbar_east'].time='ocean_time'

if north == 'yes':
    ncd.createVariable('salt_north', 'f8', ('ocean_time', 's_rho', 'xi_rho'),zlib=True)
    ncd.variables['salt_north'].long_name = 'salinity northern boundary condition'
    ncd.variables['salt_north'].units = 'PSU'
    ncd.variables['salt_north'].field = 'salt_north, scalar, series'
    ncd.variables['salt_north'].time='ocean_time'

    ncd.createVariable('temp_north', 'f8', ('ocean_time', 's_rho', 'xi_rho'),zlib=True)
    ncd.variables['temp_north'].long_name = 'potential temperature northern boundary condition'
    ncd.variables['temp_north'].units = 'Celsius'
    ncd.variables['temp_north'].field = 'temp_north, scalar, series'
    ncd.variables['temp_north'].time='ocean_time'

    ncd.createVariable('zeta_north', 'f8', ('ocean_time', 'xi_rho'),zlib=True)
    ncd.variables['zeta_north'].long_name = 'free-surface northern boundary condition'
    ncd.variables['zeta_north'].units = 'meter'
    ncd.variables['zeta_north'].field = 'ssh_north, scalar, series'
    ncd.variables['zeta_north'].time='ocean_time'

    ncd.createVariable('u_north', 'f8', ('ocean_time', 's_rho', 'xi_u'),zlib=True)
    ncd.variables['u_north'].long_name = '3D u-momentum northern boundary condition'
    ncd.variables['u_north'].units = 'meter second-1'
    ncd.variables['u_north'].field = 'u_north, scalar, series'
    ncd.variables['u_north'].time='ocean_time'

    ncd.createVariable('v_north', 'f8', ('ocean_time', 's_rho', 'xi_v'),zlib=True)
    ncd.variables['v_north'].long_name = '3D u-momentum northern boundary condition'
    ncd.variables['v_north'].units = 'meter second-1'
    ncd.variables['v_north'].field = 'v_north, scalar, series'
    ncd.variables['v_north'].time='ocean_time'

    ncd.createVariable('ubar_north', 'f8', ('ocean_time', 'xi_u'),zlib=True)
    ncd.variables['ubar_north'].long_name = '2D u-momentum northern boundary condition'
    ncd.variables['ubar_north'].units = 'meter second-1'
    ncd.variables['ubar_north'].field = 'ubar_north, scalar, series'
    ncd.variables['ubar_north'].time='ocean_time'

    ncd.createVariable('vbar_north', 'f8', ('ocean_time', 'xi_v'),zlib=True)
    ncd.variables['vbar_north'].long_name = '2D u-momentum northern boundary condition'
    ncd.variables['vbar_north'].units = 'meter second-1'
    ncd.variables['vbar_north'].field = 'vbar_north, scalar, series'
    ncd.variables['vbar_north'].time='ocean_time'

if south == 'yes':
    ncd.createVariable('salt_south', 'f8', ('ocean_time', 's_rho', 'xi_rho'),zlib=True)
    ncd.variables['salt_south'].long_name = 'salinity southern boundary condition'
    ncd.variables['salt_south'].units = 'PSU'
    ncd.variables['salt_south'].field = 'salt_south, scalar, series'
    ncd.variables['salt_south'].time='ocean_time'

    ncd.createVariable('temp_south', 'f8', ('ocean_time', 's_rho', 'xi_rho'),zlib=True)
    ncd.variables['temp_south'].long_name = 'potential temperature southern boundary condition'
    ncd.variables['temp_south'].units = 'Celsius'
    ncd.variables['temp_south'].field = 'temp_south, scalar, series'
    ncd.variables['temp_south'].time='ocean_time'

    ncd.createVariable('zeta_south', 'f8', ('ocean_time','xi_rho'),zlib=True)
    ncd.variables['zeta_south'].long_name = 'free-surface southern boundary condition'
    ncd.variables['zeta_south'].units = 'meters'
    ncd.variables['zeta_south'].field = 'ssh_south, scalar, series'
    ncd.variables['zeta_south'].time='ocean_time'

    ncd.createVariable('u_south', 'f8', ('ocean_time', 's_rho', 'xi_u'),zlib=True)
    ncd.variables['u_south'].long_name = '3D u-momentum southern boundary condition'
    ncd.variables['u_south'].units = 'meter second-1'
    ncd.variables['u_south'].field = 'u_south, scalar, series'
    ncd.variables['u_south'].time='ocean_time'

    ncd.createVariable('v_south', 'f8', ('ocean_time', 's_rho', 'xi_v'),zlib=True)
    ncd.variables['v_south'].long_name = '3D u-momentum southern boundary condition'
    ncd.variables['v_south'].units = 'meter second-1'
    ncd.variables['v_south'].field = 'v_south, scalar, series'
    ncd.variables['v_south'].time='ocean_time'

    ncd.createVariable('ubar_south', 'f8', ('ocean_time', 'xi_u'), zlib=True)
    ncd.variables['ubar_south'].long_name = '2D u-momentum southern boundary condition'
    ncd.variables['ubar_south'].units = 'meter second-1'
    ncd.variables['ubar_south'].field = 'ubar_south, scalar, series'
    ncd.variables['ubar_south'].time='ocean_time'

    ncd.createVariable('vbar_south', 'f8', ('ocean_time', 'xi_v'),zlib=True)
    ncd.variables['vbar_south'].long_name = '2D u-momentum southern boundary condition'
    ncd.variables['vbar_south'].units = 'meter second-1'
    ncd.variables['vbar_south'].field = 'vbar_south, scalar, series'
    ncd.variables['vbar_south'].time='ocean_time'

print('netcdf created...')

########################## Starting interpolating
listIndex = idx_ini
while listIndex != idx_end + 1:

    dst_sal = np.zeros((cs_r.shape[0],lon_rho.shape[0],lat_rho.shape[1]))
    dst_temp = np.ones((cs_r.shape[0],lon_rho.shape[0],lat_rho.shape[1]))
    dst_u = np.ones((cs_r.shape[0],lon_u.shape[0],lat_u.shape[1]))
    dst_v = np.ones((cs_r.shape[0],lon_v.shape[0],lat_v.shape[1]))

    print('Opening file based on the selection above...')

    inputData = nc.Dataset(listFiles[listIndex])
    timeInput = inputData.variables['time']
    tt = num2date(timeInput[0],units=timeInput.units,calendar=timeInput.calendar)
    ocean_time = date2num(tt,units=time_units,calendar='julian')

    print('The units will be:'+' '+'hours since '+year+'-'+month+'-'+day)
    ncd.variables['ocean_time'][listIndex] = ocean_time

    print('Loading data from: '+listFiles[listIndex])
    ssh = inputData.variables['ssh']

    print('interpolating ssh')
    ssh1,converged = fill(ssh[0,lat_bound,lon_bound], 0, 1, **kw)
    ssh2=mp.interp(ssh1,lonSoda,latSoda,lon_rho,lat_rho,checkbounds=False,\
              masked=False,order=1)

    ssh2 = ssh2 * mask_rho
    dst_ssh = np.ma.MaskedArray(ssh2,mask_rho==0.)

    salinity = inputData.variables['salt'][0,depth_bound,lat_bound,lon_bound]
    temperature = inputData.variables['temp'][0,depth_bound,lat_bound,lon_bound]
    u = inputData.variables['u'][0,depth_bound,latu_bound,lonu_bound]
    v = inputData.variables['v'][0,depth_bound,latu_bound,lonu_bound]

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
    #dst_grd = pyroms.grid.get_ROMS_grid(simName)

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

    if north == 'yes':
        ncd.variables['zeta_north'][listIndex,:] = dst_ssh[-1,:]
        ncd.variables['salt_north'][listIndex,:,:] = dst_sal[:,-1,:]
        ncd.variables['temp_north'][listIndex,:,:] = dst_temp[:,-1,:]
        ncd.variables['u_north'][listIndex,:,:] = dst_u1[:,-1,:]
        ncd.variables['v_north'][listIndex,:,:] = dst_v1[:,-1,:]
        ncd.variables['ubar_north'][listIndex,:] = dst_ubar[-1,:]
        ncd.variables['vbar_north'][listIndex,:] = dst_vbar[-1,:]
    if south == 'yes':
        ncd.variables['zeta_south'][listIndex,:] = dst_ssh[0,:]
        ncd.variables['salt_south'][listIndex,:,:] = dst_sal[:,0,:]
        ncd.variables['temp_south'][listIndex,:,:] = dst_temp[:,0,:]
        ncd.variables['u_south'][listIndex,:,:] = dst_u1[:,0,:]
        ncd.variables['v_south'][listIndex,:,:] = dst_v1[:,0,:]
        ncd.variables['ubar_south'][listIndex,:] = dst_ubar[0,:]
        ncd.variables['vbar_south'][listIndex,:] = dst_vbar[0,:]
    if west == 'yes':
        ncd.variables['zeta_west'][listIndex,:]  = dst_ssh[:,0]
        ncd.variables['salt_west'][listIndex,:,:] = dst_sal[:,:,0]
        ncd.variables['temp_west'][listIndex,:,:] = dst_temp[:,:,0]
        ncd.variables['u_west'][listIndex,:,:] = dst_u1[:,:,0]
        ncd.variables['v_west'][listIndex,:,:] = dst_v1[:,:,0]
        ncd.variables['ubar_west'][listIndex,:] = dst_ubar[:,0]
        ncd.variables['vbar_west'][listIndex,:] = dst_vbar[:,0]
    if east == 'yes':
        ncd.variables['zeta_east'][listIndex,:]  = dst_ssh[:,-1]  
        ncd.variables['salt_east'][listIndex,:,:] = dst_sal[:,:,-1]
        ncd.variables['temp_east'][listIndex,:,:] = dst_temp[:,:,-1]
        ncd.variables['u_east'][listIndex,:,:] = dst_u1[:,:,-1]
        ncd.variables['v_east'][listIndex,:,:] = dst_v1[:,:,-1]
        ncd.variables['ubar_east'][listIndex,:] = dst_ubar[:,-1]
        ncd.variables['vbar_east'][listIndex,:] = dst_vbar[:,-1]

    listIndex += 1

ncd.close()
