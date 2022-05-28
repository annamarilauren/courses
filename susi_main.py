# -*- coding: utf-8 -*-
"""
Created on Mon May 21 18:38:10 2018

@author: lauren
"""
import numpy as np
import pandas as pd
import datetime

from canopygrid import CanopyGrid
from mosslayer import MossLayer
from strip import StripHydrology, drain_depth_development

from susi_utils import  rew_drylimit

class Susi():
    def __init(self):
        pass

    def run_susi(self, forc, cpara, org_para, spara, start_yr, end_yr, hc, LAI): 
        
        print ('******** Susi-peatland simulator hydro (2022) c Annamari Laurén *********************')
        print ('           ')    
        print ('Initializing stand and site:') 
         
        dtc = cpara['dt']                                                          # canopy model timestep
    
        start_date = datetime.datetime(start_yr,1,1) 
        end_date=datetime.datetime(end_yr,12,31)
        length = (end_date - start_date).days +1                                   # simulation time in days
        n = spara['n']                                                             # number of columns along the strip        
        
        
        #********* Above ground hydrology initialization ***************
        cmask = np.ones(spara['n'])                                                # compute canopy and moss for each soil column (0, and n-1 are ditches)
        cstate = cpara['state'].copy()
        for key in cstate.keys():
            cstate[key] *= cmask
        cpy = CanopyGrid(cpara, cstate, outputs=False)                             # initialize above ground vegetation hydrology model
        
        for key in org_para.keys():                                                 
            org_para[key] *= cmask
        moss = MossLayer(org_para, outputs=True)  
        print ('Canopy and moss layer hydrology initialized')
    
        #******** Soil and strip parameterization *************************
        stp = StripHydrology(spara)                                                # initialize soil hydrology model
        
        
        ets = np.zeros((length, n))                                                # Evapotranspiration, mm/day
        
        #********initialize result arrays***************************
        scen = spara['scenario name']                                              # scenario name for outputs    
        rounds = len(spara['ditch depth east'])                                    # number of ditch depth scenarios (used in comparison of management)
        
        stpout = stp.create_outarrays(rounds, length, n)    
        intercs, evaps, ETs, transpis, efloors, swes = cpy.create_outarrays(rounds, length, n)
        
        #***********Scenario loop ********************************************************
        
        for r, dr in enumerate(zip(spara['ditch depth west'], spara['ditch depth 20y west'], \
                                   spara['ditch depth east'], spara['ditch depth 20y east'])):   
    
            dwt=spara['initial h']*np.ones(spara['n'])                             # set the initial WT for the scenario
            hdr_west, hdr20y_west,hdr_east, hdr20y_east = dr                       # drain depth [m] in the beginning and after 20 yrs
            h0ts_west = drain_depth_development(length, hdr_west, hdr20y_west)     # compute daily values for drain bottom boundary condition
            h0ts_east = drain_depth_development(length, hdr_east, hdr20y_east)     # compute daily values for drain bottom boundary condition
    
            
            # ---- Initialize integrative output arrays (outputs in nodewise sums) -------------------------------
    
            print ('***********************************')        
            print ('Computing canopy and soil hydrology ', length, ' days', 'scenario:', scen[r])
            
            
            stp.reset_domain()   
            
            d = 0                                                                  # day index
            start = 0                                                              # day counter in annual loop
            # *************** Annual loop *****************************************************************
            for yr in range(start_yr, end_yr + 1):                                 # year loop 
                days = (datetime.datetime(yr,12, 31) - datetime.datetime(yr,1, 1)).days + 1
              
                #**********  Daily loop ************************************************************
                for dd in range(days):                                             # day loop   
                    #-------Canopy hydrology--------------------------            
                    reww = rew_drylimit(dwt)                                       # for each column: moisture limitation from ground water level (Feddes-function)            
                    doy = forc.iloc[d, 14]                                         # day of the year
                    ta =  forc.iloc[d, 4]                                          # air temperature deg C
                    vpd = forc.iloc[d, 13]                                         # vapor pressure deficit
                    rg = forc.iloc[d, 8]                                           # solar radiation
                    par = forc.iloc[d, 10]                                         # photosynthetically active radiation
                    prec=forc.iloc[d, 7]/86400.                                    # precipitation
        
                    potinf, trfall, interc, evap, ET, transpi, efloor, MBE, SWE = cpy.run_timestep(doy, dtc, ta, prec, rg, par, vpd, 
                                                                    hc=hc*np.ones(n), LAIconif=LAI*np.ones(n), Rew=reww, beta=moss.Ree)                       # canopy hydrology computation
                    
                    intercs, evaps, ETs, transpis, efloors, SWEs = cpy.update_outarrays(r, d, interc, evap, ET, transpi, efloor, SWE)
                    
                    potinf, efloor, MBE2 = moss.interception(potinf, efloor)       # ground vegetation and moss hydrology 
                    stpout['deltas'][r, d, :] = potinf - transpi                                   # water flux thru soil surface
                    ets[d] = efloor + transpi + interc                             # evapotranspiration components 
                    
                    
                    #if d%365==0: print ('  - day #', d, ' hdom ', hc, ' m, ',  
                    #                    'LAI ', LAI, ' m2 m-2')
        
                    #--------Soil hydrology-----------------
                    stp.run_timestep(d, h0ts_west[d], h0ts_east[d], stpout['deltas'][r, d,:], moss)  # strip/peat hydrology
                    stpout = stp.update_outarrays(r, d, stpout) 
    
                  
                    swes[r,d] = np.mean(SWE)                                       # snow water equivalent    
                    d += 1
               #******* End of daily loop***************************** 
            
            #----- Hydrology and temperature-related variables to time-indexed dataframes -----------------
                #sday = datetime.datetime(yr, 1, 1)                                 # start day of the year 
                #dfwt = pd.DataFrame(stpout['dwts'][r,start:start+days,:],index=pd.date_range(sday,periods=days))                              # daily water table data frame
                #dfafp = pd.DataFrame(stpout['afps'][r,start:start+days,:],index=pd.date_range(sday,periods=days))                             # air filled porosity
                
        
            self.stpout = stpout
