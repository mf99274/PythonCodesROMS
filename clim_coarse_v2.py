import os
import time
import pyroms
import netCDF4 as nc
from gridfill import fill
import mpl_toolkits.basemap as mp
from mpl_toolkits.basemap import Basemap
from netcdftime import num2date, date2num, datetime

# files names
inputFileName = 'southernCCS_2000_2006_mon.nc'
climName = 'tsclim_nemo.nc'
gridName = 'southernCCS1_4km.nc'

# Path to my input data
inputPath = '/home/mfocean/Documents/2-PhD/98-models/2-ROMS/00-Inputs/bdr/'
gridPath = '/home/mfocean/Documents/2-PhD/98-models/2-ROMS/9-simulations/new_domain/'
climPath = '/home/mfocean/Documents/2-PhD/98-models/2-ROMS/5-Clim/'

# name of the simulation
simName = 'Baja'

# Select Year, Month, day, and hour of the input data
year = '2000'
month = '01'
day = '15'
hour = '00'

time_units='days since 1900-01-01 00:00:00'
fill_value = 1e+20

# fill values for extrapolation
kw = dict(eps=1e-4, relax=0.6, itermax=1e4, initzonal=False,
          cyclic=True, verbose=True)

#--------------- Don't need to change things after here ---------------#
inputData = nc.Dataset(inputPath+inputFileName)
depthInput = inputData.variables['depth'][:]
latInput = inputData.variables['latitude'][:]
lonInput = inputData.variables['longitude'][:]
timeInput = inputData.variables['time']

timeValues = num2date(timeInput[:],timeInput.units,calendar=timeInput.calendar) 

ocean_time = np.zeros(timeValues.shape[0])
for t in range(timeValues.shape[0]):
    oTime = timeValues[t]
    if t == 0:
        oceanTime = datetime(oTime.year, oTime.month, 1, int(hour), dayofwk=0)
        ocean_time[t] = date2num(oceanTime, units=time_units,calendar='julian')
    else:
        oceanTime = datetime(oTime.year, oTime.month, int(day), int(hour), dayofwk=0)
        ocean_time[t] = date2num(oceanTime, units=time_units,calendar='julian')
    

################ Loading grid data
gridData = nc.Dataset(gridPath+gridName)
h = gridData.variables['h'][:]
mask_rho = gridData.variables['mask_rho'][:]
cs_r = gridData.variables['Cs_r'][:]
cs_w=gridData.variables['Cs_w'][:] 
lon_rho = gridData.variables['lon_rho'][:]
lat_rho = gridData.variables['lat_rho'][:]
mask_u = gridData.variables['mask_u'][:]
mask_v = gridData.variables['mask_v'][:]

############### Initializing netcdf file
ncd = nc.Dataset(climPath+climName, 'w', format='NETCDF4')

ncd.createDimension('xi_u', mask_u.shape[1])
ncd.createDimension('xi_v', mask_v.shape[1])
ncd.createDimension('xi_rho', mask_rho.shape[1])
ncd.createDimension('eta_u', mask_u.shape[0])
ncd.createDimension('eta_v', mask_v.shape[0])
ncd.createDimension('eta_rho', mask_rho.shape[0])
ncd.createDimension('s_rho', cs_r.shape[0])
ncd.createDimension('s_w', cs_w.shape[0])
#ncd.createDimension('ocean_time',  ocean_time.shape[0])
ncd.createDimension('salt_time', ocean_time.shape[0])
ncd.createDimension('temp_time',ocean_time.shape[0])
ncd.createDimension('tracer', 2)

print('Creating lon_rho variable')
ncd.createVariable('lon_rho', 'f8', ('eta_rho','xi_rho'))
ncd.variables['lon_rho'].long_name = 'lon_rho'
ncd.variables['lon_rho'].units = 'degrees'
ncd.variables['lon_rho'].field = "xp, scalar, series"
ncd.variables['lon_rho'].FillValue = lon_rho.fill_value
ncd.variables['lon_rho'].missing_value = lon_rho.fill_value

print('Creating lat_rho variable')
ncd.createVariable('lat_rho', 'f8', ('eta_rho','xi_rho'))
ncd.variables['lat_rho'].long_name = 'lat_rho'
ncd.variables['lat_rho'].units = 'degrees'
ncd.variables['lat_rho'].field = "yp, scalar, series"
ncd.variables['lat_rho'].FillValue = lat_rho.fill_value
ncd.variables['lat_rho'].missing_value = lat_rho.fill_value

#print('Creating ocean_time variable')
#ncd.createVariable('ocean_time', 'f8', ('ocean_time'))
#ncd.variables['ocean_time'].long_name = 'open boundary conditions time'
#ncd.variables['ocean_time'].units = time_units
#ncd.variables['ocean_time'][:] = ocean_time[:]

ncd.createVariable('temp_time', 'f8', ('temp_time', ))
ncd.variables['temp_time'].long_name = 'temp_time'
ncd.variables['temp_time'].units = time_units
ncd.variables['temp_time'].field = "ocean_time, scalar, series"
ncd.variables['temp_time'][:] = ocean_time[:]

#print('Creating salt_time variable')
ncd.createVariable('salt_time', 'f8', ('salt_time', ))
ncd.variables['salt_time'].long_name = 'salt_time'
ncd.variables['salt_time'].units = time_units
ncd.variables['salt_time'].field = "salt_time, scalar, series"
ncd.variables['salt_time'][:] = ocean_time[:]

ncd.createVariable('salt', 'f8', ('salt_time', 's_rho','eta_rho','xi_rho'),zlib=True,fill_value=lon_rho.fill_value)
ncd.variables['salt'].long_name = 'salinity'
ncd.variables['salt'].units = 'PSU'
ncd.variables['salt'].field = 'salinity, scalar, series'

ncd.createVariable('temp', 'f8', ('temp_time', 's_rho','eta_rho','xi_rho'),zlib=True)
ncd.variables['temp'].long_name = 'potential temperature'
ncd.variables['temp'].units = 'Celsius'
ncd.variables['temp'].field = 'temperature, scalar, series'

listIndex = 0
while listIndex != timeValues.shape[0]:
 
    salinity = inputData.variables['so'][listIndex]
    temperature = inputData.variables['thetao'][listIndex]

    print('The file is:',str(listIndex))
    print('The time for the'+' '+str(listIndex)+' '+'is:', num2date(ocean_time[listIndex],units=time_units,calendar='julian'))

    temp_arr = []
    sal_arr = []
    u_arr=[]
    v_arr=[]

    for k in range(len(depthInput)):
        "Modified SODA2ROMS from Trond Kristiansen"
        "Date of modification: Nov,13th 2015 "
        "Author: Matheus Fagundes"
    
        kk=k+1
        print('Interpolating horizontally layer:', kk)

        salinity1 = salinity[k,:,:]
        salty,converged = fill(salinity1,0,1,**kw)
        salty = mp.interp(salty,lonInput,latInput,lon_rho,lat_rho,checkbounds=True,masked=True,order=1)
        #salty = mp.interp(salinity1,lonInput,latInput,lon_rho,\
        #         lat_rho,checkbounds=True,masked=True,order=1)        

        temperature1 = temperature[k,:,:]
        tempp,converged = fill(temperature1,0,1,**kw)
        tempp = mp.interp(tempp,lonInput,latInput,lon_rho,lat_rho,checkbounds=True,masked=True,order=1)
        #tempp = mp.interp(temperature1,lonInput,latInput,lon_rho,\
        #                lat_rho,checkbounds=True,masked=True,order=1) 

        sal_arr.append(salty)
        temp_arr.append(tempp)

    print('Converting list of temp and salinity in numpy array')
    sal=np.array(sal_arr)*mask_rho
    temp=np.array(temp_arr)*mask_rho

    print('Creating arrays for salinity and temperature')
    dst_sal = np.zeros((cs_r.shape[0],lon_rho.shape[0],lat_rho.shape[1]))
    dst_temp = np.ones((cs_r.shape[0],lon_rho.shape[0],lat_rho.shape[1]))

    print('Interpolating vertically temp and salinity')
    for x in range(lon_rho.shape[0]):
        for y in range(lat_rho.shape[1]):
            z_new = cs_r * h[x,y]
            dst_temp[:,x,y] = np.interp(z_new, flipud(-depthInput),\
                               np.transpose(temp[::-1,x,y]))
            dst_sal[:,x,y] = np.interp(z_new, flipud(-depthInput),\
                               np.transpose(sal[::-1,x,y]))

    ncd.variables['salt'][listIndex,:,:,:] = dst_sal
    ncd.variables['temp'][listIndex,:,:,:] = dst_temp

    listIndex += 1

ncd.close()

#######################################################################
#import netCDF4 as nc
#import netcdftime as nf

#oceanTime = nf.datetime(2000, 1, 1, 0, dayofwk=0)
#ocT = nf.date2num(oceanTime, units='days since 1900-01-01 00:00:00',calendar='julian')

#climPath = '/home/mfocean/Documents/2-PhD/98-models/2-ROMS/5-Clim/'
#climName = 'tsclim_nemo.nc'

#ncd = nc.Dataset(climPath+climName, 'a', format='NETCDF4')

#ncd.variables['ocean_time'][0] = ocT

#ncd.close()



