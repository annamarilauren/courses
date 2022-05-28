# -*- coding: utf-8 -*-
"""
Created on Sat May 28 15:20:49 2022

@author: alauren
"""
import numpy as np
import pandas as pd

def draw_figs(stpout, start_yr, days, start_date):
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
    
    
    dfwt30 = pd.DataFrame(stpout['dwts'][0,0:days,:],index=pd.date_range(start_date,periods=days))                              # daily water table data frame
    ls30 =dfwt30[str(start_yr)+'-07-01': str(start_yr)+'-09-30']    
    
    dfwt60 = pd.DataFrame(stpout['dwts'][1,0:days,:],index=pd.date_range(start_date,periods=days))                              # daily water table data frame
    ls60 =dfwt60[str(start_yr)+'-07-01': str(start_yr)+'-09-30']      
    
    
    dfwt90 = pd.DataFrame(stpout['dwts'][2,0:days,:],index=pd.date_range(start_date,periods=days))                              # daily water table data frame
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

#draw_figs(susi.stpout, start_yr)
