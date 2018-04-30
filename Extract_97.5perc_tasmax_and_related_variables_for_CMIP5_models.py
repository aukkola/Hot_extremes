# -*- coding: utf-8 -*-
"""
Created on Thu 3 August 2017

@author: annaukkola

"""


###--- Extract data for years 1990-2010 ---###


from netCDF4 import Dataset # to work with NetCDF files
import numpy as np
import glob
import sys
import os

from cdo import *

cdo = Cdo()



### Set paths ###

root_path = "/srv/ccrc/data04/z3509830/Fluxnet_data/CMIP5_project/Processed_data_1950_2099/"



#Source functions
lib_path  = root_path + '/../scripts/Python/functions'

sys.path.append(os.path.abspath(lib_path))
from find_days_above_percentile import *




###############################
### Set years and variables ###
###############################


### Set variables ###

### 1. Temperature variable (daily tasmax) ###

#Variable path
tasmax_var  = "tasmax_tseries"

#Variable name in netcdf file
tasmax_name = "tasmax"


### 2. Other variables ###

# saved on the day hot extreme occurs
# Variable paths
other_vars  = ["hfls_tseries", "hfss_tseries", "pr_tseries", "tasmax_tseries"]

#Variable names in Netcdf file
other_names = ["hfls", "hfss", "pr", "tasmax"]


### Set percentile ###

#code finds days above or equal to this percentile
percentile = 97.5


### Set years ###
start_yr= 1990
end_yr  = 2010


#Number of days to average prior to Txx day
#If don't want lagged variables, set to float('nan')
lags = [3,5,7,10]



#Other options
no_yrs = len(range(start_yr, end_yr+1))

#Define missing values
miss_val = -9999



### Get continents to process ###

#Find continents to process
continents = os.listdir(root_path)


#Sanity check
if len(continents) != 2:
    sys.exit("Not retrieving correct directories!")


##################
### Load files ###
##################

#Loop through continents
for k in range(len(continents)):


        #Find models
        models = os.listdir(root_path + "/" + continents[k] + "/" + tasmax_var)


        #Loop through models
        for m in range(len(models)):

            #Progress
            print 'Processing', models[m]

            #Find model file for day of Txx
            tasmax_file = glob.glob(root_path + "/" + continents[k] + "/" + tasmax_var + "/" + models[m] + "/*.nc")

            #Read lon and lat to set up output array
            fh = Dataset(tasmax_file[0], mode='r')

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



            ### Load other variables ###

            #Initialise empty character string for file names
            other_files = ["" for x in range(len(other_vars))]


            #Find files for other variables
            for v in range(len(other_vars)):
                
                #Find model file for time series
                other_files[v] = glob.glob(root_path + "/" + continents[k] + "/" +
                                           other_vars[v] + "/" + models[m] + 
                                           "/*regrid*.nc")[0]
    
    

            #Initialise output matrices for other variables
            #The function creates extra variables if using lags,
            #need to take this into account
            #Not a fan of this as code will fail if function changes...
            #But can't think of a better way right now
            
            #Set output variables, other vars first so can create lagged vars
            out_vars = other_vars
            
            #Then add extra lagged variables that code creates
            for v in range(len(out_vars)):
                for l in range(len(lags)):
                
                    lag_varname = out_vars[v] + '_lag' + str(lags[l])
                    out_vars = np.append(out_vars, lag_varname)
            
            
            #Add hot day index variable
            out_vars = np.append(out_vars, 'hot_day_ind')
            
            
            #Then add tasmax variable (changed code so already gets added above
            #as tasmax_tseries)
            #out_vars = np.append(out_vars, tasmax_name)


	    	#Initialise output arrays (days, days, lat, lon)
            #The number of days above the percentile will vary by grid cells
            #Set it to twice the number of days expected based on percentile to
            #be on the safe side 
            #e.g. if percentile 97.5, expect 365 * 0.25 days (~9 days per year)
            no_days = int(365 * ((100 - percentile) / 100) * 2) + 1
    
            #Initialise dictionary for output data
            out_data = {}
            
            for v in range(len(out_vars)):
                
                out_data[out_vars[v]] =  np.zeros((no_yrs * no_days, 
                                                   len(lat), len(lon))) + miss_val


           
	    	### Loop through years ###
            for year in np.arange(start_yr, end_yr+1):


                ### tasmax ###
                
         	    #Extract year
            	tasmax_yr = cdo.selyear(year, input=tasmax_file[0], returnArray=tasmax_name[0])

				#Mask
            	tasmax_yr[tasmax_yr.mask==True] = miss_val

                
                ### Other variables ###
                in_other_vars = {}
                
                #Get data and mask
                for v in range(len(other_files)):
                    
                    #Get data
                    var_data = cdo.selyear(year, input=other_files[v], 
                                           returnArray=other_names[v])
                    
                    #Mask 
                    var_data[var_data.mask == True] = miss_val
                        
                    #Add to dictionary
                    in_other_vars[other_vars[v]] = var_data



            	### Find tasmax days ###

                for i in range(len(lat)):

                    for j in range(len(lon)):

                		#Extract data if cell not missing
                    	if all(tasmax_yr[:,i,j] != miss_val) :


                            #Extract grid cell for other vars (needs to be a dictionary)                            
                            in_others = {}
                            
                            for v in in_other_vars:
                                in_others[v] = in_other_vars[v].data[:,i,j]
                            
                            
                            #Find days above percentile and associated variables
                            outs = find_days_above_percentile(mod_vec = tasmax_yr[:,i,j],
                                                              no_days = no_days,
                                                              percentile = percentile, 
                                                              other_vars = in_others,
                                                              lags = lags, 
                                                              miss_val = miss_val)


                            #Organise into matrix again
                            for v in out_data:
                                
                                #get start and end indices (year + 18 days)
                                ind_start = (year-start_yr) * no_days
                                ind_end   = ind_start + no_days
                                
                                #save to outputs
                                out_data[v][ind_start : ind_end, i, j] = outs[v]





            ##############################
            ### Write result to NetCDF ###
            ##############################

            #Create output path and filename
            out_path = (root_path + "/" + continents[k] + "/" + "percentile_data_" + 
                        str(percentile) + "/" + models[m])


            if not os.path.exists(out_path):
            	os.makedirs(out_path)

            out_file = (out_path + "/" + models[m] + "_data_for_days_above_percentile_" + 
                        str(percentile) + "_years_" + str(start_yr) + "_" + str(end_yr) + ".nc")


            # open a new netCDF file for writing.
            ncfile = Dataset(out_file,'w', format="NETCDF4_CLASSIC")
            # create the output data.

            # create the x, y and time dimensions
            ncfile.createDimension('lat', lat.shape[0])
            ncfile.createDimension('lon', lon.shape[0])
            ncfile.createDimension('days', out_data[v].shape[0])


            # create variables
            # first argument is name of variable, second is datatype, third is
            # a tuple with the names of dimensions.

            longitude = ncfile.createVariable("lon",  'f8', ('lon',))
            latitude  = ncfile.createVariable("lat",  'f8', ('lat',))
            day       = ncfile.createVariable("days",'i4', ('days',))


            nc_out_data = {}
            for v in out_data:
                nc_out_data[v] = ncfile.createVariable(v, 'f8',
                                ('days', 'lat','lon'), fill_value=miss_val)


            #Set variable information
            longitude.units = 'degrees_east'
            latitude.units = 'degrees_north'

            #data.long_name = (varname[v] + ' lagged from day Txx occurs')


            # write data to variable.
            longitude[:]= lon
            latitude[:] = lat
            day[:]      = range(1, no_yrs * no_days + 1)

            for v in nc_out_data:
                nc_out_data[v][:,:,:] = out_data[v]

            # close the file.
            ncfile.close()
            
            
            
            
            
            
            
            
            
            
