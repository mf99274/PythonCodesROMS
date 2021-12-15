import time
import glob
import netCDF4 as nc
import netcdftime as nf
from gridfill import fill
import mpl_toolkits.basemap as mp

def netcdfFunction(dictionary,name,fill_value=-9999.):

    if name == 'str' or name == 'ssr':
        ncd.createVariable(dictionary[name][0],'f8',('time','lat','lon'),zlib=True, fill_value=fill_value)
        ncd.variables[dictionary[name][0]].long_name = dictionary[name][1]
        ncd.variables[dictionary[name][0]].units = dictionary[name][2]
        ncd.variables[dictionary[name][0]].time = 'time'
        ncd.variables[dictionary[name][0]].positive_value = dictionary[name][4]
        ncd.variables[dictionary[name][0]].negative_value = dictionary[name][5]
        ncd.variables[dictionary[name][0]].coordinates = 'lon lat'
    else:
        ncd.createVariable(dictionary[name][0],'f8',('time','lat','lon'),zlib=True, fill_value=fill_value)
        ncd.variables[dictionary[name][0]].long_name = dictionary[name][1]
        ncd.variables[dictionary[name][0]].units = dictionary[name][2]
        ncd.variables[dictionary[name][0]].time = 'time'
        ncd.variables[dictionary[name][0]].coordinates = 'lon lat'

def timeData(timeInput,year=2000):

    timeData = []
    timeIdx = []
    for t in range(timeInput[:].shape[0]):
        # calculating time
        timeValues = nf.num2date(timeInput[t],timeInput.units,calendar=timeInput.calendar)

        if timeValues.year == year:
            timeValues = nf.num2date(timeInput[t]-9.,timeInput.units,calendar=timeInput.calendar)
            print('The date is: ', timeValues)
            atmTime = nf.datetime(timeValues.year, timeValues.month, timeValues.day, timeValues.hour, dayofwk=0)
            timeData.append(nf.date2num(atmTime, units=time_units,calendar=calendar))

            # getting index
            ttIdx = np.where(timeInput==timeInput[t])
            timeIdx.append(ttIdx[0][0])
            
        else:
            continue
            
    return timeData,timeIdx

def convertFunction(dataValue,dictionary,times,timeIndex,name):
    """This function get the data from the any file and interpolate to the right 
       time and domain
       dataValue ->  the data based on the file
       dictionary -> the name of the variables for the data you got
       times -> just to make sure you got the right data time
       timeIndex -> the right indexes for the data
       name  -> based on the list of all variables 
       hrs -> the constant to divide the atmospheric fluxes
       """
    variable = []
    times = times
    timeIndex = timeIndex
    scale_factor = dataValue.scale_factor
    add_offset = dataValue.add_offset
    for t in range(len(timeIndex)):

        timeValues = nf.num2date(times[t],units=time_units,calendar=calendar)
        print('The date is: ', timeValues)
        
        if name == 't2m':
            variable.append(dataValue[timeIndex[t],::-1,:] - 273.15)
        
        elif name == 'ssrd':
            variable.append(dataValue[timeIndex[t],::-1,:]/(3600))

        elif name == 'strd':
            variable.append(dataValue[timeIndex[t],::-1,:]/(3600))

        elif name == 'str':
            variable.append(dataValue[timeIndex[t],::-1,:]/(3600))

        elif name == 'slhf':
            variable.append(dataValue[timeIndex[t],::-1,:]/(3600))

        elif name == 'sshf':
            variable.append(dataValue[timeIndex[t],::-1,:]/(3600))

        elif name == 'msl':
            variable.append(dataValue[timeIndex[t],::-1,:]*0.01)
 
        elif name == 'e':
            rho_w = 1000 #kg/m3
            variable.append(dataValue[timeIndex[t],::-1,:]*(rho_w/(3600)))

        elif name == 'tp':
            rho_w = 1000 #kg/m3
            variable.append(dataValue[timeIndex[t],::-1,:]*(rho_w/(3600)))

        else:
            variable.append(dataValue[timeIndex[t],::-1,:])

    ncd.variables[dictionary[name][0]][:,:,:] = np.asarray(variable) #* mask_rho


def airHumidity(dataValue1,dataValue2,times,timeIndex,name):
    """Calculating air humidity in percentage
       QAIR = 100 * (E/Es) 
       where:
       E  = 6.11 * 10.0 ** (7.5 * d2m / (237.7 + d2m))    vapor pressure (mb)
                                                         d2m in Celsius
       Es = 6.11 * 10.0 ** (7.5 * t2m / (237.7 + t2m))    saturation vapor
                                                         pressure (mb)
                                                         t2m in Celsius"""

    variable = []
    times = times
    timeIndex = timeIndex    

    for t in range(len(timeIndex)):

        timeValues = nf.num2date(times[t],units=time_units,calendar=calendar)
        print('The date is: ', timeValues) 
        t2m = dataValue1[timeIndex[t],::-1,:] - 273.15
        t2d = dataValue2[timeIndex[t],::-1,:] - 273.15
        E = 6.11*10.**(7.5*t2d/(237.7 + t2d))
        Es = 6.11*10.**(7.5*t2m/(237.7 + t2m))

        qair = 100 * (E/Es)
        variable.append(qair)

    ncd.variables['Qair'][:,:,:] = np.asarray(variable)

# Path to my input data
inputPath = '/home/mfocean/Documents/2-PhD/98-models/2-ROMS/10-atmosphericData/era5/'
listFiles = glob.glob(inputPath+'adaptor.mars.*.nc')

# setting time units like the model
time_units='days since 1900-01-01 00:00:00'
calendar = 'julian'

# selecting fill value
fill_value = -9999.

# dictionary for era5 atmospheric data
era5Names = {}
era5Names['u10'] = 'Uwind','surface u-wind component','meter second-1','time'
era5Names['v10'] = 'Vwind','surface v-wind component','meter second-1','time'
era5Names['msl'] = 'Pair','surface air pressure','millibar'
#era5Names['ssr'] = 'swrad','solar shortwave radiation','Watts meter-2','time','downward flux, heating','upward flux, cooling'
#era5Names['str'] = 'lwrad','net longwave radiation flux','Watts meter-2','time','downward flux, heating','upward flux, cooling'
era5Names['t2m'] = 'Tair','surface air temperature','Celsius'
era5Names['tcc'] = 'cloud','cloud fraction','nondimensional','time'
#era5Names['strd'] = 'lwrad_down','downwelling longwave radiation flux','Watts meter-2','time'
#era5Names['ssrd'] = 'swrad','solar shortwave radiation','Watts meter-2','time'
#era5Names['slhf'] = 'latent', 'net latent heat flux','watt meter-2'
#era5Names['sshf'] = 'sensible', 'net sensible heat flux', 'watt meter-2'
era5Names['e'] = 'evaporation','evaporation rate','kilogram meter-2 second-1'
era5Names['tp'] = 'rain','rain fall rate','kilogram meter-2 second-1'
era5Names['d2m'] = 'Qair','surface air relative humidity','percentage','time'

############################ year nc files
#year = np.arange(2000,2007,1)

yr = 2006
# atmopheric file
AtmPath = inputPath
AtmName = 'atm_exceptRad'+str(yr)+'_v1.nc'

listData = {}
listNames = []
for i in range(len(listFiles)):
    number = i +1
    aa = nc.Dataset(listFiles[i]).variables
    bb = list(aa.keys())
    listData['file'+str(i+1)] = nc.Dataset(listFiles[i]).variables

    for j in range(len(bb)):
        listNames.append(bb[j])

    listNames1 = list(set(listNames))

listData1 = list(listData.values())

# creating netcdf file
ncd = nc.Dataset(AtmPath+AtmName, 'w', format='NETCDF4')

for k in listNames1:
    k2 = list(listData.keys())
    for i in range(len(listData1)):
        if k == 'time':
            tt,ttIdx = timeData(listData[k2[i]][k],year=yr)
       
            # time
            ncd.createDimension('time',len(tt))
            ncd.createVariable('time', 'f8', ('time'))
            ncd.variables['time'].long_name = 'atmospheric forcing time'
            ncd.variables['time'].units = time_units
            ncd.variables['time'][:] = tt
            break
        elif k == 'longitude':
            lonInput = listData[k2[i]][k][:]

        elif k == 'latitude':
            latInput = listData[k2[i]][k][::-1]


ncd.createDimension('lon', lonInput.shape[0])
ncd.createVariable('lon', 'f8', ('lon'))
ncd.variables['lon'].long_name = 'atmospheric forcing longitude'
ncd.variables['lon'][:] = lonInput.data

ncd.createDimension('lat', latInput.shape[0])
ncd.createVariable('lat', 'f8', ('lat'))
ncd.variables['lat'].long_name = 'atmospheric forcing latitude'
ncd.variables['lat'][:] = latInput.data

for k in listNames1:
    k2 = list(listData.keys())
    for i in range(len(listData1)):
        if k in listData[k2[i]]:
            if k in era5Names:
                print(k)
                data = listData[k2[i]][k]                
                netcdfFunction(era5Names,k,fill_value=1e+20)
                print('netcdf var created')
                if k == 'd2m':
                    data1 = listData[k2[i]][k]  
                    data2 = listData['file9']['t2m']
                    airHumidity(data2,data1,tt,ttIdx,k)
                else:
                    convertFunction(data,era5Names,tt,ttIdx,k)                
                print('')

ncd.type='Atmophere'+' '+str(yr)
ncd.author='Cobialab'
ncd.history='Created '+time.ctime(time.time())
ncd.title='Mexico'

ncd.close()

