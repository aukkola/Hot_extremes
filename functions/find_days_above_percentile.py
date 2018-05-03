# -*- coding: utf-8 -*-
"""
Created on Fri May  6 14:54:23 2016

@author: annaukkola
"""


def find_days_above_percentile(mod_vec, no_days, percentile=97.5, other_vars={},
                               lags=float('nan'), miss_val=float('nan')):


    # source packages
    import sys
    import os
    import numpy as np


    ########################
    ### Find tasmax days ###
    ########################

    #Convert mod_vec into np.array in case not already, indexing might fail otherwise
    mod_vec = np.array(mod_vec)
    

    ### Calculate percentile threshold ###
    
    threshold = np.percentile(mod_vec, percentile)


    ### Find tasmax days above or equal to threshold ###
    
    #Weird problem using percentile = 100.0, hacky fix
    if percentile == 100.0:
        #if found several, pick first one
        hot_days = np.where(mod_vec == np.max(mod_vec))[0]
        if len(hot_days) > 1:
            hot_days = hot_days[0]
    else:
        hot_days = np.where(mod_vec >= threshold)[0]
    

    #Initialise as the correct size
    tasmax = np.zeros(no_days) + miss_val


    #Then find temperatures of hot days
    tasmax[0:len(hot_days)] = mod_vec[hot_days]


    #Collate outputs (add 1 to hot day indices for normal 1-based indexing)
    #Similarly make hot days the correct size for outputting
    hot_days_out = np.zeros(no_days) + miss_val 
    hot_days_out[0:len(hot_days)] = hot_days + 1 

    outs={'tasmax': tasmax, 'hot_day_ind': hot_days_out}


    #######################################################
    ### Find additional variables on the day and lagged ###
    #######################################################


    #If other variables to process
    if len(other_vars) > 0 :
        
        #Get variable names
        keys = other_vars.keys()
        
        #Initialise dictionary for outputs
        vars_on_the_day = {}
        
        
        #Loop through variables
        for k in other_vars:
                        
            #Convert dict to np.array or indexing might fail
            data = np.array(other_vars[k])
            
            #Get hot day values
            vars_on_the_day[k] = data[hot_days]
            
            #Make data the correct size (no_days) for outputting
            if(len(vars_on_the_day[k]) < no_days):
                
                vars_on_the_day[k] = np.append(vars_on_the_day[k], 
                                     np.zeros(no_days - len(vars_on_the_day[k])) + 
                                     miss_val)
                
                
            
            #Also calculate lagged variables
            if all(np.isnan(lags)) == False:
                
                
                #Then loop through lags
                for l in range(len(lags)):
                    
                    lag_varname = k + '_lag' + str(lags[l])

                    #Initialise as correct length
                    lag_data = np.zeros(len(hot_days)) + miss_val
                    
            
            
                    #Loop through hot days
                    for h in range(len(hot_days)):            
                    
                        #If can't count back, leave as missing
                        if (hot_days[h] - (lags[l] - 1)) < 0:
                            continue
                    
                        else:
                            ind = np.arange(hot_days[h] - (lags[l]-1), hot_days[h]+1)
                            lag_data[h] = np.mean(data[ind])
                        
            
                    #Add lag data variable to output dictionary
                    vars_on_the_day[lag_varname] = lag_data
                    
                    if(len(vars_on_the_day[lag_varname]) < no_days):
                            
                            vars_on_the_day[lag_varname] = np.append(vars_on_the_day[lag_varname], 
                                                 np.zeros(no_days - len(vars_on_the_day[lag_varname])) + 
                                                 miss_val)
                            

            
        #Add additional variables to outputs
        outs.update(vars_on_the_day)

                

    #Return outputs
    return outs;
