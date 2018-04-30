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

from cdo import *

cdo = Cdo()



### Set paths ###

root_path = "/srv/ccrc/data04/z3509830/Fluxnet_data/CMIP5_project/Processed_data_1950_2099/"

PRECIP_DIR="/srv/ccrc/data45/z5025137/markus/pr"


#Find continents to process
continents = os.listdir(root_path)


vars=["tasmax_tseries", "hfls_tseries", "hfss_tseries"]
varname=["tasmax", "hfls", "hfss"]




#################
### Set years ###
#################

start_yr= 1950
end_yr  = 2099

#Number of days to extract prior to Txx day
timesteps = 90

no_yrs = len(range(start_yr, end_yr+1))

#Define missing values
miss_val = -9999


##################
### Load files ###
##################


for k in range(len(continents)):

    #Find files for max day indices
    models = os.listdir(root_path + "/" + continents[k] + "/" + "hfls")


    #Loop through models
    for m in range(len(models)):

        #Progress
        print 'Processing', models[m]

        #Find model file for day of Txx
        maxday_file = glob.glob(root_path + "/" + continents[k] + "/" + "MaxDay" + "/" + models[m] + "/*.nc")

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



	#Loop through variables
        for v in range(len(vars)):

      	    #Find model file for time series
            var_file = glob.glob(root_path + "/" + continents[k] + "/" + vars[v] + "/" + models[m] + "/*.nc")

	    	#Initialise output arrays (years, days, lat, lon)
	    	out_var=np.zeros((no_yrs, timesteps, len(lat), len(lon))) + miss_val

	    	#Loop through years
            for year in np.arange(start_yr, end_yr+1):

         	    #Extract year
            	maxday_yr = cdo.selyear(year, input=maxday_file[0], returnArray='DOY')
				var_yr = cdo.selyear(year, input=var_file[0], returnArray=varname[v])


				#Mask
            	maxday_yr[maxday_yr.mask==True] = miss_val
	        	var_yr[var_yr.mask==True] = miss_val


            	### Find lagged days ###

        		for i in range(len(lat)):

            	    for j in range(len(lon)):

                		#Extract data if cell not missing
                    	if maxday_yr[:,i,j] != miss_val and any(var_yr[:,i,j,] != miss_val):

                            ind = maxday_yr[:,i,j][0]
                            out_var[year-start_yr, :, i,j] = var_yr[np.arange(ind-timesteps, ind), i, j]






            ##############################
            ### Write result to NetCDF ###
            ##############################

            #Create output path and filename
            out_path = (root_path + "/" + continents[k] + "/" + varname[v] + "_lagged" + "/" + models[m])


			if not os.path.exists(out_path):
            	os.makedirs(out_path)

            out_file = (out_path + "/" + models[m] + "_" + varname[v] + "_lagged_" + str(timesteps) + "_days_1950_2099.nc")


            # open a new netCDF file for writing.
            ncfile = Dataset(out_file,'w', format="NETCDF4_CLASSIC")
            # create the output data.

            # create the x, y and time dimensions
            ncfile.createDimension('lat', lat.shape[0])
            ncfile.createDimension('lon', lon.shape[0])
            ncfile.createDimension('year', no_yrs)
            ncfile.createDimension('day', timesteps)


            # create variables
            # first argument is name of variable, second is datatype, third is
            # a tuple with the names of dimensions.

            longitude = ncfile.createVariable("lon",  'f8', ('lon',))
            latitude  = ncfile.createVariable("lat",  'f8', ('lat',))
            year      = ncfile.createVariable("year", 'i4', ('year',))
            day       = ncfile.createVariable("day",'i4', ('day',))


            data = ncfile.createVariable(varname[v], 'f8',('year', 'day', 'lat','lon'), fill_value=miss_val)


            #Set variable information
            longitude.units = 'degrees_east'
            latitude.units = 'degrees_north'

            data.long_name = (varname[v] + ' lagged from day Txx occurs')


            # write data to variable.
            longitude[:]= lon
            latitude[:] = lat
            year[:]     = range(start_yr, end_yr+1)
            day[:]      = range(1, timesteps+1)


            data[:,:,:,:] = out_var

            # close the file.
            ncfile.close()
