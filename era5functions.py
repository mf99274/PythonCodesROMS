########################################################################
# This code gets all the netcdf data from era5 and just converts       #
# the variables for ROMS readable variables, also, you don't need to   #
# tell the code what netcdf file is, it does by itself                 #
# 03/28/22, Matheus Fagundes                                           #
########################################################################
#import time
#import glob
import netCDF4 as nc
import netcdftime as nf

def savingVariables(dictionary,name,n4file,fillValue):
    """This function saves the variable and use era5 dictionary to created
       the netcdf 
       The dictionary needs to have the following in order:
       - name that ROMS can read
       - long name
       - units
     """

    if name == 'str' or name == 'ssr':
        n4file.createVariable(dictionary[name][0],'f8',('time','lat','lon'),zlib=True, fill_value=fillValue)
        n4file.variables[dictionary[name][0]].long_name = dictionary[name][1]
        n4file.variables[dictionary[name][0]].units = dictionary[name][2]
        n4file.variables[dictionary[name][0]].time = 'time'
        n4file.variables[dictionary[name][0]].positive_value = dictionary[name][4]
        n4file.variables[dictionary[name][0]].negative_value = dictionary[name][5]
        n4file.variables[dictionary[name][0]].coordinates = 'lon lat'
    else:
        n4file.createVariable(dictionary[name][0],'f8',('time','lat','lon'),zlib=True, fill_value=fillValue)
        n4file.variables[dictionary[name][0]].long_name = dictionary[name][1]
        n4file.variables[dictionary[name][0]].units = dictionary[name][2]
        n4file.variables[dictionary[name][0]].time = 'time'
        n4file.variables[dictionary[name][0]].coordinates = 'lon lat'
            

def convertingTime(netcdfFile,n4file,time_units='days since 1900-01-01 00:00:00',calendar='julian'):
    """time_units is the time_units you used to create BC for ROMS
       calendar is also the same calendar used to create BC for ROMS """

    aa = netcdfFile
    time = aa.variables['time']
    timeValues = nf.num2date(time[:],time.units,calendar=time.calendar)
    # converting back to number but using ROMS settings
    timeValues2 = nf.date2num(timeValues,units=time_units,calendar=calendar)
    
    n4file.createDimension('time',timeValues.shape[0])
    n4file.createVariable('time', 'f8', ('time'))
    n4file.variables['time'].long_name = 'atmospheric forcing time'
    n4file.variables['time'].units = time_units
    n4file.variables['time'][:] = timeValues2

def creatingLonLat(netcdfFile,n4file):
    """Changing the names on longitude and latitude """

    n4file.createDimension('lon', netcdfFile.variables['longitude'].shape[0])
    n4file.createVariable('lon', 'f8', ('lon'))
    n4file.variables['lon'].long_name = 'atmospheric forcing longitude'
    n4file.variables['lon'][:] = netcdfFile.variables['longitude'][:]

    n4file.createDimension('lat', netcdfFile.variables['latitude'].shape[0])
    n4file.createVariable('lat', 'f8', ('lat'))
    n4file.variables['lat'].long_name = 'atmospheric forcing latitude'
    n4file.variables['lat'][:] = netcdfFile.variables['latitude'][:]
    
def convertingVar(netcdfFile,name,dictionary,n4file):
    """This function converts atmospheric variables that 
       have different units from the original dataset.
        """
    if name == 'u10' or name == 'v10':
        n4file.variables[dictionary[name][0]][:,:,:] = netcdfFile.variables[name][:]

    elif name == 't2m':
        n4file.variables[dictionary[name][0]][:,:,:] = netcdfFile.variables[name][:] - 273.15

    elif name == 'sp':
        n4file.variables[dictionary[name][0]][:,:,:] = netcdfFile.variables[name][:]*0.01 

    elif name == 'tp' or name == 'e':
        rho_w = 1000 #kg/m3
        n4file.variables[dictionary[name][0]][:,:,:] = netcdfFile.variables[name][:]*(rho_w/(3600))

    elif name == 'tcc':
        n4file.variables[dictionary[name][0]][:,:,:] = netcdfFile.variables[name][:]




