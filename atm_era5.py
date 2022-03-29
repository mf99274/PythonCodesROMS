import glob
import time
import netCDF4 as nc
import era5functions as era5

# Path to my input data
inputPath = '/home/mfocean/Documents/2-PhD/98-models/2-ROMS/00-Inputs/atmosphere/'
listFiles = glob.glob(inputPath+'adaptor.mars.*.nc')

# setting time units like the model
time_units='days since 1900-01-01 00:00:00'
calendar = 'julian'

# selecting fill value
fill_value = -9999.

# dictionary for era5 atmospheric data
#era5list = ['u10','v10','t2m','sp','e','tp','tcc']

era5Names = {}
era5Names['u10'] = 'Uwind','surface u-wind component','meter second-1'
era5Names['v10'] = 'Vwind','surface v-wind component','meter second-1'
era5Names['t2m'] = 'Tair','surface air temperature','Celsius'
era5Names['tp'] = 'rain','rain fall rate','kilogram meter-2 second-1'
era5Names['e'] = 'evaporation','evaporation rate','kilogram meter-2 second-1'
era5Names['sp'] = 'Pair','surface air pressure','millibar'
era5Names['tcc'] = 'cloud','cloud fraction','nondimensional','time'
#era5Names['ssr'] = 'swrad','solar shortwave radiation','Watts meter-2','time','downward flux, heating','upward flux, cooling'


era5list = ['e','tp','tcc']

for idx in range(len(era5list)):

    for i in range(len(listFiles)):
        aa = nc.Dataset(listFiles[i]).variables
        bb = list(aa.keys())
        for var in bb:
            if era5list[idx] != var:
                continue
            else:
                print(var) 
                data = nc.Dataset(listFiles[i])
                #val1 = data.variables[var]            
                ncd = nc.Dataset(inputPath+'atm_'+var+'.nc', 'w', format='NETCDF4')
                era5.convertingTime(nc.Dataset(listFiles[i]),ncd)
                era5.creatingLonLat(nc.Dataset(listFiles[i]),ncd)
                era5.savingVariables(era5Names,var,ncd,fillValue=fill_value)
                era5.convertingVar(nc.Dataset(listFiles[i]),var,era5Names,\
                ncd)

                ncd.type='Atmophere'+' '+var
                ncd.author='Cobialab'
                ncd.history='Created '+time.ctime(time.time())
                ncd.title='Mexico'

                ncd.close()


