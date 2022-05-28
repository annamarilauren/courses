# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 14:10:42 2020

@author: alauren
"""
import numpy as np
import pandas as pd
import datetime
from figs import draw_figs
from susi_utils import  read_FMI_weather
from susi_para import get_susi_para
from susi_main import Susi 


#***************** local call for SUSI*****************************************************
folderName=r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_3_py37/outputs/' #'sensitivity/'
susiPath = r'C:/Users/alauren/Documents/Susi_9/'
wpath = r'C:/Users/alauren/Documents/Susi_9/'



wdata='parkano_weather.csv'

start_date = datetime.datetime(2000,1,1)
end_date=datetime.datetime(2000,12,31)
start_yr = start_date.year 
end_yr = end_date.year
yrs = (end_date - start_date).days/365.25
days = (end_date - start_date).days

sarkaSim = 40.                                                                  # strip width, ie distance between ditches, m
n = int(sarkaSim / 2)                                                           # number of computation nodes in the strip
hc = 21.0
LAI = 6.5


sfc =  np.ones(n, dtype=int)*3                                                                        # site fertility class

site = 'develop_scens'

forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
            
wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                          folderName=folderName, hdomSim=hc,  
                                                                          ageSim=None, sarkaSim=sarkaSim, sfc=sfc, 
                                                                          susiPath=susiPath,
                                                                          n=n)


susi = Susi()

susi.run_susi(forc, cpara, org_para, spara, start_yr, end_yr, hc, LAI)

draw_figs(susi.stpout, start_yr, days, start_date)