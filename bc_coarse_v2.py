import os
import glob
import time
import pyroms
import numpy as np
import netCDF4 as nc
from gridfill import fill
import mpl_toolkits.basemap as mp
from mpl_toolkits.basemap import Basemap
from netcdftime import num2date, date2num, datetime

# name of the simulation
simName = 'Baja'

# Select Year, Month, day, and hour to start bry data
year = '2000'
month = '01'
day = '01'
hour = '00'

# Select Year, Month, day, and hour to finish bry data
year1 = '2000'
month1 = '12'
day1 = '30'
hour1 = '21'

time_units='days since 1900-01-01 00:00:00'

# Path to my input data
inputPath = '/media/mfocean/easystore/hycom/'+year
os.chdir(inputPath)

# Path to coarse grid
gridPath = '/home/mfocean/Documents/2-PhD/98-models/2-ROMS/9-simulations/new_domain/'
gridName = 'southernCCS1_4km.nc'
gridData = nc.Dataset(gridPath+gridName)

BCPath = '/home/mfocean/Documents/2-PhD/98-models/2-ROMS/9-simulations/new_domain/'
BCname = 'southernCCS1_BC_hycom'+str(year)+'.nc'

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
mask_rho = gridData.variables['mask_rho'][:]
cs_r = gridData.variables['Cs_r'][:]
cs_w=gridData.variables['Cs_w'][:] 
lon_rho = gridData.variables['lon_rho'][:]
lat_rho = gridData.variables['lat_rho'][:]
mask_u = gridData.variables['mask_u'][:]
mask_v = gridData.variables['mask_v'][:]

dst_grd = pyroms.grid.get_ROMS_grid(simName) #Loading values from my grid
z_u_n = 0.5 * (dst_grd.vgrid.z_w[0,:,-1,:-1] + dst_grd.vgrid.z_w[0,:,-1,1:])
z_u_s = 0.5 * (dst_grd.vgrid.z_w[0,:,0,:-1] + dst_grd.vgrid.z_w[0,:,0,1:])
z_u_w = 0.5 * (dst_grd.vgrid.z_w[0,:,:,0] + dst_grd.vgrid.z_w[0,:,:,1])
#    z_u_e = 0.5 * (dst_grd.vgrid.z_w[0,:,:,-1] + dst_grd.vgrid.z_w[0,:,:,-2])   
 
z_v_n = 0.5 * (dst_grd.vgrid.z_w[0,:,-1,:] + dst_grd.vgrid.z_w[0,:,-2,:])
z_v_s = 0.5 * (dst_grd.vgrid.z_w[0,:,0,:] + dst_grd.vgrid.z_w[0,:,1,:])
z_v_w = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,0] + dst_grd.vgrid.z_w[0,:,1:,0])

print('Getting the number of files to create the bdry condition')

listFiles = glob.glob("hycom_*.nc")
listFiles.sort()

print('Checking if the beginning of bry and the end are in the path')

for length in range(len(listFiles)):
    if listFiles[length] == 'hycom_'+year+month+day+'12_t0'+hour+'.nc':
        print('The file '+ 'hycom_'+year+month+day+'12_t0'+hour+'.nc'+' '+'is in the folder')
        idx_ini = length

for length in range(len(listFiles)):
    if listFiles[length] == 'hycom_'+year1+month1+day1+'12_t0'+hour1+'.nc':
        print('The file '+ 'hycom_'+year1+month1+day1+'12_t0'+hour1+'.nc'+' '+'is in the folder')
        idx_end = length

inputData = nc.Dataset(listFiles[0])
latInput = inputData.variables['lat'][:]
lonInput = inputData.variables['lon'][:]

if north == 'yes':
    lat_n = np.where(np.round(latInput,1)==np.round(lat_rho[-1,0],1))[0][0]
    lon_n,h1 = np.meshgrid(lon_rho[-1,:],cs_r)
    hh_n = np.zeros((cs_r.shape[0],lon_rho.shape[1]))
    for j_xi  in range(lon_rho.shape[1]):
            hh_n[:,j_xi] = cs_r*h[-1,j_xi]
    
if south == 'yes':
    lat_s = np.where(np.round(latInput,1)==np.round(lat_rho[0,0],1))[0][0]
    lon_s,h1 = np.meshgrid(lon_rho[0,:],cs_r)
    hh_s = np.zeros((cs_r.shape[0],lon_rho.shape[1]))
    for j_xi  in range(lon_rho.shape[1]):
            hh_s[:,j_xi] = cs_r*h[0,j_xi]

if west == 'yes':
    lon_w = np.where(np.round(lonInput,1)==np.round(lon_rho[0,0],1))[0][0]
    lat_w,h1 = np.meshgrid(lat_rho[:,0],cs_r)
    hh_w = np.zeros((cs_r.shape[0],lat_rho.shape[0]))
    for i_xi  in range(lat_rho.shape[0]):
            hh_w[:,i_xi] = cs_r*h[i_xi,0]

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

listIndex = idx_ini
while listIndex != idx_end + 1:

    print('Opening file based on the selection above...')

    inputData = nc.Dataset(listFiles[listIndex])
    
    depthInput = inputData.variables['depth'][:]
    latInput = inputData.variables['lat'][:]
    lonInput = inputData.variables['lon'][:]

    print('Loading time...')
    timeInput = inputData.variables['time']

    print('Converting time')
    print('Note to myself...')
    print('When using num2date function from netcdftime, the timestamp')
    print('has 12 hours more added since it is HOURS since 2000, and therefore, by subtracting 12')
    print('I get the correct timestamp')
    print('....')
    print('....')
    print('....')
    print('Subtracting 12 from the timeInput')

    timeValues = num2date(timeInput[0]-12.,timeInput.units,\
                          calendar=timeInput.calendar)

    print('Converting to my ocean time')
    oceanTime = datetime(timeValues.year, timeValues.month,\
                   timeValues.day, timeValues.hour, dayofwk=0)

    print('The units will be:'+' '+'hours since '+year+'-'+month+'-'+day+' '+hour+':00:00')
    ocean_time = date2num(oceanTime, units=time_units,calendar='julian')
    
    ncd.variables['ocean_time'][listIndex] = ocean_time

    print('Loading data from: '+listFiles[listIndex])
    salinity = inputData.variables['salinity']
    temperature = inputData.variables['water_temp']
    ssh = inputData.variables['surf_el']
    u = inputData.variables['water_u']
    v = inputData.variables['water_v']

    ssh0,converged = fill(ssh[0,:,:], 0, 1, **kw)
    ssh1=mp.interp(ssh0,lonInput,latInput,lon_rho,lat_rho,\
                  checkbounds=True,masked=True,order=1)

    ssh1 = ssh1 * mask_rho
    dst_ssh = np.ma.MaskedArray(ssh1,mask_rho==0.)

    if north == 'yes':
        ncd.variables['zeta_north'][listIndex,:] = dst_ssh[-1,:]
    if south == 'yes':
        ncd.variables['zeta_south'][listIndex,:] = dst_ssh[0,:]
    if west == 'yes':
        ncd.variables['zeta_west'][listIndex,:]  = dst_ssh[:,0]
    if east == 'yes':
        ncd.variables['zeta_east'][listIndex,:]  = dst_ssh[:,-1]  

    if west == 'yes':

        tempp = mp.interp(temperature[0,::-1,:,lon_w],latInput,\
             np.flipud(-depthInput),lat_w,hh_w,checkbounds=True,order=1)
        salty = mp.interp(salinity[0,::-1,:,lon_w],latInput,\
             np.flipud(-depthInput),lat_w,hh_w,checkbounds=True,order=1)
        vv = mp.interp(v[0,::-1,:,lon_w],latInput,\
             np.flipud(-depthInput),lat_w,hh_w,checkbounds=True,order=1)
        uu = mp.interp(u[0,::-1,:,lon_w],latInput,\
             np.flipud(-depthInput),lat_w,hh_w,checkbounds=True,order=1)

        uu_w = 0.5*(uu[:,:]+uu[:,:])
        vv_w =  0.5*(vv[:,:-1]+vv[:,1:])
        ncd.variables['salt_west'][listIndex,:,:] = salty
        ncd.variables['temp_west'][listIndex,:,:] = tempp
        ncd.variables['u_west'][listIndex,:,:] = uu_w
        ncd.variables['v_west'][listIndex,:,:] = vv_w
     
        # ubar and vbar for north/south boundaries
        print('Calculating ubar for north and south boundaries')
        ubar_west = np.zeros((mask_u.shape[0]))
        vbar_west = np.zeros((mask_v.shape[0]))

        for n in range(ubar_west.shape[0]):
            ubar_west[n] = (uu_w[:,n] * np.diff(z_u_w[:,n])).sum() / -z_u_w[0,n]

        for n in range(vbar_west.shape[0]):
            vbar_west[n] = (vv_w[:,n] * np.diff(z_v_w[:,n])).sum() / -z_v_w[0,n]

        vbar_west = np.ma.masked_invalid(vbar_west)
        ubar_west = np.ma.masked_invalid(ubar_west)

        ncd.variables['ubar_west'][listIndex,:] = ubar_west[:]
        ncd.variables['vbar_west'][listIndex,:] = vbar_west[:]

    if north == 'yes':
        
        tempp = mp.interp(temperature[0,::-1,lat_n,:],lonInput,\
             np.flipud(-depthInput),lon_n,hh_n,checkbounds=True,order=1)
        salty = mp.interp(salinity[0,::-1,lat_n,:],lonInput,\
             np.flipud(-depthInput),lon_n,hh_n,checkbounds=True,order=1)
        vv = mp.interp(v[0,::-1,lat_n,:],lonInput,\
             np.flipud(-depthInput),lon_n,hh_n,checkbounds=True,order=1)
        uu = mp.interp(u[0,::-1,lat_n,:],lonInput,\
             np.flipud(-depthInput),lon_n,hh_n,checkbounds=True,order=1)

        uu_n = 0.5*(uu[:,:-1]+uu[:,1:])
        vv_n =  0.5*(vv[:,:]+vv[:,:])
        ncd.variables['salt_north'][listIndex,:,:] = salty
        ncd.variables['temp_north'][listIndex,:,:] = tempp
        ncd.variables['u_north'][listIndex,:,:] = uu_n
        ncd.variables['v_north'][listIndex,:,:] = vv_n

        # ubar and vbar for north/south boundaries
        print('Calculating ubar for north and south boundaries')
        ubar_north = np.zeros((mask_u.shape[1]))
        vbar_north = np.zeros((mask_v.shape[1]))

        for n in range(ubar_north.shape[0]):
            ubar_north[n] = (uu_n[:,n] * np.diff(z_u_n[:,n])).sum() / -z_u_n[0,n]

        for n in range(vbar_north.shape[0]):
            vbar_north[n] = (vv_n[:,n] * np.diff(z_v_n[:,n])).sum() / -z_v_n[0,n]

        vbar_north = np.ma.masked_invalid(vbar_north)
        ubar_north = np.ma.masked_invalid(ubar_north)

        ncd.variables['ubar_north'][listIndex,:] = ubar_north[:]
        ncd.variables['vbar_north'][listIndex,:] = vbar_north[:]

    if south == 'yes':

        tempp = mp.interp(temperature[0,::-1,lat_s,:],lonInput,\
             np.flipud(-depthInput),lon_s,hh_s,checkbounds=True,order=1)
        salty = mp.interp(salinity[0,::-1,lat_s,:],lonInput,\
             np.flipud(-depthInput),lon_s,hh_s,checkbounds=True,order=1)
        vv = mp.interp(v[0,::-1,lat_s,:],lonInput,\
             np.flipud(-depthInput),lon_s,hh_s,checkbounds=True,order=1)
        uu = mp.interp(u[0,::-1,lat_s,:],lonInput,\
             np.flipud(-depthInput),lon_s,hh_s,checkbounds=True,order=1)

        uu_s = 0.5*(uu[:,:-1]+uu[:,1:])
        vv_s =  0.5*(vv[:,:]+vv[:,:])
        ncd.variables['salt_south'][listIndex,:,:] = salty
        ncd.variables['temp_south'][listIndex,:,:] = tempp
        ncd.variables['u_south'][listIndex,:,:] = uu_s
        ncd.variables['v_south'][listIndex,:,:] = vv_s

        # ubar and vbar for north/south boundaries
        print('Calculating ubar for north and south boundaries')
        ubar_south = np.zeros((mask_u.shape[1]))
        vbar_south = np.zeros((mask_v.shape[1]))

        for n in range(ubar_south.shape[0]):
            ubar_south[n] = (uu_s[:,n] * np.diff(z_u_s[:,n])).sum() / -z_u_s[0,n]

        for n in range(vbar_south.shape[0]):
            vbar_south[n] = (vv_s[:,n] * np.diff(z_v_s[:,n])).sum() / -z_v_s[0,n]

        vbar_south = np.ma.masked_invalid(vbar_south)
        ubar_south = np.ma.masked_invalid(ubar_south)

        ncd.variables['ubar_south'][listIndex,:] = ubar_south[:]
        ncd.variables['vbar_south'][listIndex,:] = vbar_south[:]

    listIndex += 1

ncd.close()

