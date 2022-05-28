# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 14:10:42 2020

@author: alauren
"""
import numpy as np
import pandas as pd
import datetime
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

import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec


def create_profile_line(ax, df, wtmin, cols, ylabel, label, fs, facecolor, colorin, title = '' ,hidey = False, hidex=False):
    wt = df.mean(axis=0)
    sd = df.std(axis=0)
    ax.plot(np.arange(cols)*2, wt, color=colorin, label = label)
    ax.fill_between(np.arange(cols)*2, wt+sd*2, wt-sd, color=colorin, alpha=0.075)
    ax.hlines(y= -0.35, xmin=0, xmax = cols*2, color='red',linestyles='--')
    if hidex: 
        ax.get_xaxis().set_visible(False) 
    else:
        ax.tick_params(axis='x', labelsize=fs)
    
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_ylim([wtmin,0])
    ax.set_ylabel(ylabel, fontsize=fs)
    ax.legend()
    #ax.grid(visible=False)
    if hidey: ax.get_yaxis().set_visible(False)
    
    ax.set_facecolor(facecolor)
    ax.set_title(title)
    return ax


dfwt30 = pd.DataFrame(susi.stpout['dwts'][0,0:days,:],index=pd.date_range(start_date,periods=days))                              # daily water table data frame
mean30 = dfwt30.mean(axis = 0)
gs30 =dfwt30[str(start_yr)+'-05-01': str(start_yr)+'-10-31']    
ls30 =dfwt30[str(start_yr)+'-07-01': str(start_yr)+'-09-30']    

dfwt60 = pd.DataFrame(susi.stpout['dwts'][1,0:days,:],index=pd.date_range(start_date,periods=days))                              # daily water table data frame
mean60 = dfwt60.mean(axis = 0)
gs60 =dfwt60[str(start_yr)+'-05-01': str(start_yr)+'-10-31']      
ls60 =dfwt60[str(start_yr)+'-07-01': str(start_yr)+'-09-30']      


dfwt90 = pd.DataFrame(susi.stpout['dwts'][2,0:days,:],index=pd.date_range(start_date,periods=days))                              # daily water table data frame
mean90 = dfwt90.mean(axis = 0)
gs90 =dfwt90[str(start_yr)+'-05-01': str(start_yr)+'-10-31']      
ls90 =dfwt90[str(start_yr)+'-07-01': str(start_yr)+'-09-30']      


cols =  dfwt30.shape[1]      
facecolor = '#f2f5eb'
fs = 15
fig = plt.figure(num='hydro', figsize=(15,9))   #width, height
gs = gridspec.GridSpec(ncols=12, nrows=12, figure=fig, wspace=0.35, hspace=1.0)

ax = fig.add_subplot(gs[0:4, 0:4])
ax = create_profile_line(ax, dfwt30, -1.0, cols, 'WT m', 'annual', fs, facecolor, 'blue', title= '30 cm', hidex=True)

ax = fig.add_subplot(gs[4:8, 0:4])
ax = create_profile_line(ax, ls30, -1.0, cols, 'WT m', 'late summer', fs, facecolor, 'blue')

ax = fig.add_subplot(gs[0:4, 4:8])
ax = create_profile_line(ax, dfwt60, -1.0, cols, 'WT m', 'annual', fs, facecolor, 'orange', hidey=True , title= '60 cm', hidex=True)

ax = fig.add_subplot(gs[4:8, 4:8])
ax = create_profile_line(ax, ls60, -1.0, cols, 'WT m', 'late summer', fs, facecolor, 'orange', hidey=True)

ax = fig.add_subplot(gs[0:4, 8:12])
ax = create_profile_line(ax, dfwt90, -1.0, cols, 'WT m', 'annual', fs, facecolor, 'green', hidey=True, title= '90 cm', hidex=True)

ax = fig.add_subplot(gs[4:8, 8:12])
ax = create_profile_line(ax, ls90, -1.0, cols, 'WT m', 'late summer', fs, facecolor, 'green', hidey=True)

ax = fig.add_subplot(gs[8:,: ]) 
ax.plot(dfwt30.index, dfwt30.mean(axis=1), color='blue', label='30 cm') 
ax.fill_between(dfwt30.index, dfwt30.mean(axis=1) - dfwt30.std(axis=1), 
                dfwt30.mean(axis=1) + dfwt30.std(axis=1), color='blue', alpha = 0.2) 

ax.plot(dfwt60.index, dfwt60.mean(axis=1), color='orange', label='60 cm') 
ax.fill_between(dfwt60.index, dfwt60.mean(axis=1) - dfwt60.std(axis=1), 
                dfwt60.mean(axis=1) + dfwt60.std(axis=1), color='orange', alpha = 0.2) 


ax.plot(dfwt90.index, dfwt90.mean(axis=1), color='green', label='90 cm') 
ax.fill_between(dfwt90.index, dfwt90.mean(axis=1) - dfwt90.std(axis=1), 
                dfwt90.mean(axis=1) + dfwt90.std(axis=1), color='green', alpha = 0.2) 

ax.set_ylabel('WT m', fontsize=fs)

ax.legend(loc='upper center')
ax.set_facecolor(facecolor) 
ax.tick_params(axis='y', labelsize=fs)
ax.hlines(y= -0.35, xmin=dfwt30.index[0], xmax = dfwt30.index[-1], color='red',linestyles='--')
 