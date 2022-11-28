# -*- coding: utf-8 -*-
"""
Created on Sun May 29 20:46:14 2022

@author: alauren
"""
from figs import draw_figs
from susi_utils import  read_FMI_weather
from susi_para import get_susi_para
from susi_main import Susi 
#from google.colab import output
import numpy as np                         # Numerical python, superhero of computing
from datetime import datetime

folderName=r'/content/courses/' 
susiPath = r'/content/courses/'
wpath = r'/content/courses/'


def run_hydro(vuosi = 2000, sarkaleveys = 40.0, turve = 'carex', puusto = 'varttunut kasvatusmetsä'):
  
  wdata='parkano_weather.csv'

  start_date = datetime(int(vuosi),1,1)
  end_date=datetime(int(vuosi),12,31)
  start_yr = start_date.year 
  end_yr = end_date.year
  days = (end_date - start_date).days
  
  sarkaSim = sarkaleveys                                                   # strip width, ie distance between ditches, m
  n = int(sarkaleveys / 2)                                                 # number of computation nodes in the strip
                                                                           # number of computation nodes in the strip
  sfc =  np.ones(n, dtype=int)*3                                           # site fertility class
  site = 'develop_scens'
  hc=10.0

  forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
              
  wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                            folderName=folderName, hdomSim=hc,  
                                                                            ageSim=None, sarkaSim=sarkaSim, sfc=sfc, 
                                                                            susiPath=susiPath,
                                                                            n=n)
  if turve == 'carex':
    spara['peat type'] = ['A','A','A','A','A','A','A','A'] 
    spara['peat type bottom']=['A']
  elif turve =='sphagnum':
    spara['peat type'] = ['S','S','S','S','S','S','S','S'] 
    spara['peat type bottom']=['S']

  if puusto == 'taimikko':
    hc = 3.0; LAI=1.0
  elif puusto == 'nuori kasvatusmetsä':
    hc = 8.0; LAI = 2.5
  elif puusto == 'varttunut kasvatusmetsä':
    hc = 12.0; LAI= 4.0
  elif puusto == 'hakkuukypsä':
    hc = 22; LAI = 6.5

  susi = Susi()
  susi.run_susi(forc, cpara, org_para, spara, start_yr, end_yr, hc, LAI)
  output.clear()
  draw_figs(susi.stpout, start_yr, days, start_date)
  return susi.stpout