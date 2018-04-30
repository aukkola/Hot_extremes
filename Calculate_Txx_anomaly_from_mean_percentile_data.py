# -*- coding: utf-8 -*-
"""
Created on Thu 3 August 2017

@author: annaukkola

"""


from netCDF4 import Dataset # to work with NetCDF files
import numpy as np
import glob
import sys
import os
import calendar
import datetime as dt
from cdo import *

cdo = Cdo()



### Set paths ###

root_path = "/srv/ccrc/data04/z3509830/Fluxnet_data/CMIP5_project/Processed_data_1950_2099/"



#Find continents to process
continents = os.listdir(root_path)


vars=["tasmax_month_mean"]
varname=["tasmax_tseries"]

day_ind_var=["hot_day_ind"]



#################
### Set years ###
#################

start_yr= 1990
end_yr  = 2010

#Number of days to extract prior to Txx day
timesteps = 90

no_yrs = len(range(start_yr, end_yr+1))

#Define missing values
miss_val = -9999


#No of days per month
no_days = calendar.mdays[1:13]

#indices for month start and end
month_end = np.cumsum(no_days)
month_start = month_end - no_days


##################
### Load files ###
##################


for k in range(len(continents)):

    #Find files for max day indices
    models = os.listdir(root_path + "/" + continents[k] + "/percentile_data_97.5/")


    #Loop through models
    for m in range(len(models)):

        #Progress
        print 'Processing', models[m]

        #Find model file for day of Txx
        maxday_file = glob.glob(root_path + "/" + continents[k] + "/percentile_data_97.5/" + models[m] + "/*.nc")

        #Read lon and lat to set up output array
        fh = Dataset(maxday_file[0], mode='r')

        #Get lat and lon
        try:
            lat = fh.variables['latitude'][:]
            lon = fh.variables['longitude'][:]
        except:
            try:
                lat = fh.variables['northing'][:]
                lon = fh.variables['easting'][:]
            except:
                lat = fh.variables['lat'][:]
                lon = fh.variables['lon'][:]


        #Get Tmax and hot day index
        tmax = fh.variables['tasmax_tseries'][:]
        maxday=fh.variables['hot_day_ind'][:]


	    #Loop through variables
        for v in range(len(vars)):

      	    #Find model file for monthly mean Txx
            var_file = glob.glob(root_path + "/" + continents[k] + "/" + vars[v] + 
                                 "/" + models[m] + "/*" + str(start_yr) + "-" + str(end_yr) +"*.nc")

            #Load data
            nc = Dataset(var_file[0], mode='r')
            tmax_mean = nc.variables['tasmax'][:]

            try:
                lat_mean = nc.variables['latitude'][:]
            except:
                try:
                    lat_mean = nc.variables['northing'][:]
                except:
                    lat_mean = nc.variables['lat'][:]


            #Flip if doesn't match Tmax data
            if lat[0] != lat_mean[0]:
                tmax_mean = tmax_mean[:,::-1,:]


	    	#Initialise output arrays (years, days, lat, lon)
            out_var=np.zeros((maxday.shape[0], len(lat), len(lon))) + miss_val

	    	#Loop through years
            for day in np.arange(0, maxday.shape[0]):


				#Mask
            	maxday[maxday.mask==True] = miss_val
                tmax[tmax.mask==True] = miss_val
                tmax_mean[tmax_mean.mask==True] = miss_val
                


            	### Calculate anomaly from monthly mean Txx ###

                for i in range(len(lat)):

                    for j in range(len(lon)):

                        #Extract data if cell not missing
                        if maxday[day,i,j] != miss_val:

                            #Max day
                            ind = maxday[day,i,j]

                            #Find month when Txx occurs (setting to 1990 here, slight error if a leap year)
                            date  = dt.datetime(1990, 1, 1) + dt.timedelta(ind - 1)
                            month = date.month

                            #Find corresponding mean monthly Txx
                            mean_txx = tmax_mean[month-1, i, j]

                            #Calculate anomaly from mean Txx
                            out_var[day, i, j] = tmax[day,i,j] - mean_txx


            ##############################
            ### Write result to NetCDF ###
            ##############################

            #Create output path and filename
            out_path = (root_path + "/" + continents[k] + "/Tmax_anomaly_percentile_data" + "/" + models[m])


            if not os.path.exists(out_path):
                os.makedirs(out_path)

            out_file = (out_path + "/" + models[m] + "_Tmax_anomaly_from_mean_" + 
                        str(start_yr) +  "_" + str(end_yr) + ".nc")


            # open a new netCDF file for writing.
            ncfile = Dataset(out_file,'w', format="NETCDF4_CLASSIC")
            # create the output data.

            # create the x, y and time dimensions
            ncfile.createDimension('lat', lat.shape[0])
            ncfile.createDimension('lon', lon.shape[0])
            ncfile.createDimension('day', maxday.shape[0])


            # create variables
            # first argument is name of variable, second is datatype, third is
            # a tuple with the names of dimensions.

            longitude = ncfile.createVariable("lon",  'f8', ('lon',))
            latitude  = ncfile.createVariable("lat",  'f8', ('lat',))
            day      = ncfile.createVariable("day", 'i4', ('day',))


            data = ncfile.createVariable(varname[v], 'f8',('day', 'lat','lon'), fill_value=miss_val)


            #Set variable information
            longitude.units = 'degrees_east'
            latitude.units = 'degrees_north'

            data.long_name = 'Tmax anomaly from 1990-2010 mean'


            # write data to variable.
            longitude[:]= lon
            latitude[:] = lat
            day[:]     = range(0, maxday.shape[0])


            data[:,:,:] = out_var

            # close the file.
            ncfile.close()
