# -*- coding: utf-8 -*-
"""
Created on Tue Aug 08 10:38:45 2017

@author: lauren
"""

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline as interS


def peat_hydrol_properties(x, unit='g/cm3', var='bd', ptype='A'):
    """
    Peat water retention and saturated hydraulic conductivity as a function of bulk density
    Päivänen 1973. Hydraulic conductivity and water retention in peat soils. Acta forestalia fennica 129.
    see bulk density: page 48, fig 19; degree of humification: page 51 fig 21
    Hydraulic conductivity (cm/s) as a function of bulk density(g/cm3), page 18, as a function of degree of humification see page 51 
    input:
        - x peat inputvariable in: db, bulk density or dgree of humification (von Post)  as array \n
        - bulk density unit 'g/cm3' or 'kg/m3' \n
        - var 'db' if input variable is as bulk density, 'H' if as degree of humification (von Post) \n
        - ptype peat type: 'A': all, 'S': sphagnum, 'C': Carex, 'L': wood, list with length of x 
    output: (ThetaS and ThetaR in m3 m-3)
        van Genuchten water retention parameters as array [ThetaS, ThetaR, alpha, n] \n
        hydraulic conductivity (m/s)
    """
    #paras is dict variable, parameter estimates are stored in tuples, the model is water content = a0 + a1x + a2x2, where x is
    para={}                                                                     #'bd':bulk density in g/ cm3; 'H': von Post degree of humification
    para['bd'] ={'pF0':(97.95, -79.72, 0.0), 'pF1.5':(20.83, 759.69, -2484.3),
            'pF2': (3.81, 705.13, -2036.2), 'pF3':(9.37, 241.69, -364.6),
            'pF4':(-0.06, 249.8, -519.9), 'pF4.2':(0.0, 174.48, -348.9)}
    para['H'] ={'pF0':(95.17, -1.26, 0.0), 'pF1.5':(46.20, 8.32, -0.54),
            'pF2': (27.03, 8.14, -0.43), 'pF3':(17.59, 3.22, -0.07),
            'pF4':(8.81, 3.03, -0.10), 'pF4.2':(5.8, 2.27, -0.08)}
    
    intp_pF1={}                                                                 # interpolation functions for pF1        
    intp_pF1['bd'] = interp1d([0.04,0.08,0.1,0.2],[63.,84.,86.,80.],fill_value='extrapolate')
    intp_pF1['H'] = interp1d([1.,4.,6.,10.],[75.,84.,86.,80.],fill_value='extrapolate')
    
    #Saturatated hydraulic conductivity parameters
    Kpara ={'bd':{'A':(-2.271, -9.80), 'S':(-2.321, -13.22), 'C':(-1.921, -10.702), 'L':(-1.921, -10.702)}, 
            'H':{'A':(-2.261, -0.205), 'S':(-2.471, -0.253), 'C':(-1.850, -0.278), 'L':(-2.399, -0.124)}}
    
    vg_ini=(0.88,	0.09, 0.03, 1.3)                                              # initial van Genuchten parameters (porosity, residual water content, alfa, n)

    x = np.array(x)
    prs = para[var]; pF1=intp_pF1[var]
    if unit=='kg/m3'and var=='db': x=x/1000.
    if  np.shape(x)[0] >1 and len(ptype)==1:
        ptype=np.repeat(ptype, np.shape(x)[0])        
    vgen = np.zeros((np.size(x),4))
    Ksat = np.zeros((np.size(x)))
    
    #wcont = lambda x, (a0, a1, a2): a0 + a1*x + a2*x**2.
    wcont = lambda x, *a: a[0] + a[1]*x + a[2]*x**2.
    van_g = lambda pot, *p:   p[1] + (p[0] - p[1]) / (1. + (p[2] * pot) **p[3]) **(1. - 1. / p[3])   
    #K = lambda x, (a0, a1): 10.**(a0 + a1*x) / 100.   # to m/s   
    K = lambda x, *a: 10.**(a[0] + a[1]*x) / 100.   # to m/s   
    
    potentials =np.array([0.01, 10.,32., 100.,1000.,10000.,15000. ])
    
    wc = (np.array([wcont(x,*prs['pF0']), pF1(x), wcont(x,*prs['pF1.5']), wcont(x,*prs['pF2']),
               wcont(x,*prs['pF3']), wcont(x,*prs['pF4']),wcont(x,*prs['pF4.2'])]))/100.
        
    for i,s in enumerate(np.transpose(wc)):
        vgen[i],_= curve_fit(van_g,potentials,s, p0=vg_ini)                      # van Genuchten parameters
        
    for i, a, pt in zip(range(len(x)), x, ptype):
        Ksat[i] = K(a, *Kpara[var][pt])                                          # hydraulic conductivity (cm/s -> m/s) 
    
    return vgen, Ksat

def CWTr(nLyrs, z, dz, pF, Ksat, direction='positive'):
    """
    Returns interpolation functions 
        sto=f(gwl)  profile water storage as a function ofground water level
        gwl=f(sto)  ground water level
        tra=f(gwl)  transissivity
    Input:
        nLyrs number of soil layers
        d depth of layer midpoint
        dz layer thickness
        pF van Genuchten water retention parameters: ThetaS, ThetaR, alfa, n
        Ksat saturated hydraulic conductivity in m s-1
        direction: positive or negative downwards
    """    
    #-------Parameters ---------------------
    z = np.array(z)   
    dz =np.array(dz)
    nroot = 8   # 8 good number of layers in rooting zone
    #nroot =10      #for jääli simulations
    nroot2 = 3  #2 10 cm root layers for air-filled porosity

    #--------- Connection between gwl and water storage------------
    d = 6 if direction == 'positive' else -6   
    gwl=np.linspace(0,d,150)
    if direction == 'positive':
        sto = [sum(wrc(pF, x = np.minimum(z-g, 0.0))*dz) for g in gwl]     #equilibrium head m
        storoot = [np.sum(wrc(pF, x = np.minimum(z-g, 0.0))[0:nroot]*dz[0:nroot]) for g in gwl]
        storoot2 = [np.sum(wrc(pF, x = np.minimum(z-g, 0.0))[0:nroot2]*dz[0:nroot2]) for g in gwl]
    else:
        sto = [sum(wrc(pF, x = np.minimum(z+g, 0.0))*dz) for g in gwl]     #equilibrium head m
        storoot = [np.sum(wrc(pF, x = np.minimum(z+g, 0.0))[0:nroot]*dz[0:nroot]) for g in gwl]
        storoot2 = [np.sum(wrc(pF, x = np.minimum(z+g, 0.0))[0:nroot2]*dz[0:nroot2]) for g in gwl]

    gwlToSto = interp1d(np.array(gwl), np.array(sto), fill_value='extrapolate')
    airtot = sto[0]-sto                                                         #m air in the profile
    airroot = storoot[0]-storoot                                                #m air in the rooting zone
    afproot = (storoot2[0]-storoot2)/(sum(dz[:nroot2]))                         #air-filled porosity in root layer
    ratio = airroot[1:]/airtot[1:]                                            #share of air-filled porosity in rooting zone to total air volume
    sto = list(sto); gwl= list(gwl); ratio=list(ratio); afproot = list(afproot)         
    sto.reverse(); gwl.reverse(); ratio.reverse(); afproot.reverse()
    stoToGwl =interp1d(np.array(sto), np.array(gwl), fill_value='extrapolate')
    gwlToRatio = interp1d(np.array(gwl[1:]), np.array(ratio), fill_value='extrapolate' )
    gwlToAfp= interp1d(np.array(gwl), np.array(afproot), fill_value='extrapolate' )
    C = interp1d(np.array(gwl), np.array(np.gradient(gwlToSto(gwl))/np.gradient(gwl)), fill_value='extrapolate')  #storage coefficient function      
    
    del gwl, sto, ratio, afproot
        
    #----------Transmissivity-------------------
    K=np.array(Ksat*86400.)   #from m/s to m/day
    tr =[sum(K[t:]*dz[t:]) for t in range(nLyrs)]        
    if direction=='positive':        
        gwlToTra = interS(z, np.array(tr))            
    else:
        z= list(z);  z.reverse(); tr.reverse()
        gwlToTra = interS(-np.array(z), np.array(tr))                    
    del tr
    return gwlToSto, stoToGwl, gwlToTra, C, gwlToRatio, gwlToAfp

def wrc(pF, x=None, var=None):
    """
    vanGenuchten-Mualem soil water retention curve\n
    IN:
        pF - dict['ThetaS': ,'ThetaR': ,'alpha':, 'n':,] OR
           - list [ThetaS, ThetaR, alpha, n]
        x  - soil water tension [m H2O = 0.1 kPa]
           - volumetric water content [vol/vol]
        var-'Th' is x=vol. wat. cont.
    OUT:
        res - Theta(Psii) or Psii(Theta)
    NOTE:\n
        sole input 'pF' draws water retention curve and returns 'None'. For drawing give only one pF-parameter set. 
        if several pF-curves are given, x can be scalar or len(x)=len(pF). In former case var is pF(x), in latter var[i]=pf[i,x[i]]
               
    Samuli Launiainen, Luke 2/2016
    """
    if type(pF) is dict: #dict input
        #Ts, Tr, alfa, n =pF['ThetaS'], pF['ThetaR'], pF['alpha'], pF['n']
        Ts=np.array(pF['ThetaS'].values()); Tr=np.array( pF['ThetaR'].values()); alfa=np.array( pF['alpha'].values()); n=np.array( pF['n'].values())
        m= 1.0 -np.divide(1.0,n)
    elif type(pF) is list: #list input
        pF=np.array(pF, ndmin=1) #ndmin=1 needed for indexing to work for 0-dim arrays
        Ts=pF[0]; Tr=pF[1]; alfa=pF[2]; n=pF[3] 
        m=1.0 - np.divide(1.0,n)
    elif type(pF) is np.ndarray:
        Ts, Tr, alfa, n = pF.T[0], pF.T[1], pF.T[2], pF.T[3]
        m=1.0 - np.divide(1.0,n)
    else:
        print ('Unknown type in pF')
        
    def theta_psi(x): #'Theta-->Psi'
        x=np.minimum(x,Ts) 
        x=np.maximum(x,Tr) #checks limits
        s= ((Ts - Tr) / (x - Tr))#**(1/m)
        Psi=-1e-2/ alfa*(s**(1/m)-1)**(1/n) # in m
        return Psi
        
    def psi_theta(x): # 'Psi-->Theta'
        x=100*np.minimum(x,0) #cm
        Th = Tr + (Ts-Tr)/(1+abs(alfa*x)**n)**m
        return Th           
 
    if var == 'Th': y=theta_psi(x) #'Theta-->Psi'           
    else: y=psi_theta(x) # 'Psi-->Theta'          
    return y

def hydrCond(pF, x=None, var=None, Ksat=1):
    """
    Hydraulic conductivity following vanGenuchten-Mualem \n
    IN:
        pF - dict or list
        x - Theta [vol/vol] or Psi [m H2O]
        var = 'Th' if x in [vol/vol]\n
        Ksat - saturated hydraulic conductivity [units]\n
    OUT:
        Kh - hydraulic conductivity ( if Ksat ~=1 then in [units], else relative [-]) \n
    """
    import matplotlib.pylab as plt
    if type(pF) is dict: #dict input
        alfa=np.array( pF['alpha']); n=np.array( pF['n'])
        m= 1.0 -np.divide(1.0,n)
        
    else: #list input
        pF=np.array(pF, ndmin=1).T #ndmin=1 needed for indexing of 0-dim arrays
        alfa=pF[2]; n=pF[3]         
        m=1.0 - np.divide(1.0,n)

    def kRel(x):
        nm=(1 - abs(alfa*x)**(n-1) * (1 + abs(alfa*x)**n)**(-m))**2
        dn=(1 + abs(alfa*x)**n)**(m/2.0)
        r=nm/dn
        return r

    if x is None and np.size(alfa)==1:  #draws pf-curve
        xx=-np.logspace(-4,5,100) #cm
        yy=kRel(xx)
        fig=plt.figure()
        fig.suptitle('Hydr. cond. (vanGenuchten-Mualem)', fontsize=16)
        #ttext=str(pF).translate(None,"{}'")
        ttext= r'$K_{sat}=$' +str(Ksat) +r', $\alpha=$'+str(alfa)+ ', n='+str(n)

        plt.title(ttext, fontsize=14)
        plt.semilogx(-xx,yy,'g-')

        plt.ylabel(r'K_{sat}', fontsize=14) 
        plt.xlabel('$\psi$ $(cm)$', fontsize=14)

        del xx, yy
        return None
        
    elif x is None: 
        print ('hydrCond: To draw curve give only one pF -parameter set')
        return None
        
    # this computes and returns    
    x=np.array(x)
    if x is not None and var is 'Th': x=wrc(pF,x=x, var='Th')

    #If psi is positive, set to zero (saturated conductivity)
    if x is not None and var is 'Psii': x[x>0] = 0

    Kh=Ksat*kRel(100.0*x)
    
    return Kh        

def read_FMI_weather(ID, start_date,end_date, sourcefile=None):
    """ 
    reads FMI interpolated daily weather data from file 
    IN: 
        ID - sve catchment ID. set ID=0 if all data wanted
        start_date - 'yyyy-mm-dd'
        end_date - 'yyyy-mm-dd'
    OUT:
        fmi - pd.dataframe with datetimeindex
            fmi columns:['ID','Kunta','aika','lon','lat','T','Tmax','Tmin','Prec','Rg','h2o','dds','Prec_a','Par','RH','esa','VPD','doy']
            units: T, Tmin, Tmax, dds[degC], VPD, h2o,esa[kPa], Prec, Prec_a[mm], Rg,Par[Wm-2],lon,lat[deg]
    """
    
    #OmaTunniste;OmaItä;OmaPohjoinen;Kunta;siteid;vuosi;kk;paiva;longitude;latitude;t_mean;t_max;t_min;
    #rainfall;radiation;hpa;lamposumma_v;rainfall_v;lamposumma;lamposumma_cum
    #-site number
    #-date (yyyy mm dd)
    #-latitude (in KKJ coordinates, metres)
    #-longitude (in KKJ coordinates, metres)
    #-T_mean (degrees celcius)
    #-T_max (degrees celcius)
    #-T_min (degrees celcius)
    #-rainfall (mm)
    #-global radiation (per day in kJ/m2)
    #-H2O partial pressure (hPa)

    ID=0
    print (' + Reading meteorological input file from ')
    print ('    -', sourcefile)
    ID=0
    #import forcing data
    #fmi=pd.read_csv(sourcefile, sep=';', header='infer', usecols=['OmaTunniste','Kunta','aika','longitude','latitude','t_mean','t_max','t_min',\
    #'rainfall','radiation','hpa'],parse_dates='aika')
    #fmi=pd.read_csv(sourcefile, sep=';', header='infer', usecols=['OmaTunniste','Kunta','aika','longitude','latitude','t_mean','t_max','t_min','rainfall','radiation','hpa'])
    #time=pd.to_datetime(fmi['aika'],format='%Y%m%d')
    
    
    fmi=pd.read_csv(sourcefile, sep=';', header='infer', 
                    usecols=['OmaTunniste','Kunta','aika','longitude','latitude','t_mean','t_max','t_min','rainfall',\
                             'radiation','hpa'], encoding= 'ISO-8859-1')
 
    
    #print pd.to_datetime(fmi['aika'][0], format="%Y%m%d")
    #print pd.tseries.tools.to_datetime(fmi['aika'][0], format="%Y%m%d")
    time=pd.to_datetime(fmi['aika'],format='%Y%m%d')
    
    fmi.index=time
    fmi=fmi.rename(columns={'OmaTunniste': 'ID', 'longitude':'lon','latitude':'lat','t_mean':'T','t_max':'Tmax','t_min':'Tmin','rainfall':'Prec',\
        'radiation':'Rg', 'hpa':'h2o'})

    
    fmi['h2o']=1e-1*fmi['h2o'] #hPa-->kPa
    fmi['Rg']=1e3/86400.0*fmi['Rg'] #kJ/m2/d-1 to Wm-2 
    fmi['Par']=0.5*fmi['Rg']

    #saturated vapor pressure    
    esa=0.6112*np.exp((17.67*fmi['T'])/ (fmi['T'] +273.16 -29.66))  #kPa
    vpd=esa - fmi['h2o']; #kPa   
    vpd[vpd<0]=1e-5
    rh=100.0*fmi['h2o']/esa;
    rh[rh<0]=1e-6; rh[rh>100]=100.0
                
    fmi['RH']=rh;
    fmi['esa']=esa;
    fmi['vpd']=vpd
    fmi['doy']=fmi.index.dayofyear
    fmi=fmi.drop(['aika'],axis = 1)



    #replace nan's in prec with 0.0
    fmi['Prec'].fillna(value=0.0)    
    #del dat, fields, n, k, time
    
    #get desired period
    fmi=fmi[(fmi.index >= start_date) & (fmi.index <= end_date)]
    if ID >0:
        fmi=fmi[fmi['ID']==ID]

    return fmi
    
def rew(dwt):
    '''Alla hyväksi havaitut'''
    wt=np.array([-150.0, -1.0, -0.5, -0.3, -0.1, 0.0])   #water table in m
    re=np.array([0.0, 0.2, 0.4 ,1.0, 1.0, 0.7])    #relative water uptake

    frew=interp1d(wt,re,fill_value='extrapolate')    
    return frew(dwt)


def rewFloor(dwt,LAI):          # Leenan testi
    # Restriction for efloor: similar to rew for transpiration but the restriction increases with decreasing LAI (3...1.5)
    LAIarray=np.array([0, 1.5, 3, 5])
    rajKerroin=np.array([1, 1, 1, 1]) 
    
    fLAI = interp1d(LAIarray,rajKerroin,fill_value='extrapolate')  
    
    wt=np.array([-150.0, -1.0, -0.5, -0.3, -0.1, 0.0])   #water table in m
    re=fLAI(LAI) * np.array([1, 1, 1, 1, 1, 1])
    
    frew=interp1d(wt,re,fill_value='extrapolate')   
    
    return frew(dwt)

def rew_drylimit(dwt):
    # Koivusalo et al. 2008 HESS without wet side limit
    wt=np.array([-150.0, -1.2, -0.7, -0.15, 0.0])   #water table in m
    re=np.array([0.0, 0.5, 1.0, 1.0, 1.0])    #relative water uptake
    

    frew=interp1d(wt,re,fill_value='extrapolate')
    return frew(dwt)

    