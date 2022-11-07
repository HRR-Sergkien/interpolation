#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 13:25:17 2022

@author: Paula Izquierdo Sanchez and Sergio H. Ramirez
"""

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import time
import pandas as pd
from os import listdir
import glob
#from PyAstronomy.pyasl import *
from scipy.interpolate import RegularGridInterpolator
import itertools
from functools import lru_cache
from astropy.io import fits

def read_fits(which_params,fits_file):
    hdu=fits.open(fits_file)
    param_vals=[]
    wave_model=hdu[1].data['wavelength']
    flux_model=hdu[1].data['flux']
    model_name=''
    for param in which_params:  #loops over specified parameters in the header of the fits file
        param_val=hdu[0].header[param]
        if param=='teff': #I included this if to be able to apply np.log10(teff)
            param_vals.append(np.log10(param_val))
        else:
            param_vals.append(param_val)
        model_name=model_name+str(param_val)  #I never actually use this name
    hdu.close()
    return param_vals,wave_model,flux_model,model_name  #Returns the values of the specified parameters, and wavelength and flux of each model


#def grid_mods(wwave,fold_name,which_params):
#    return grid_mods_tuple(tuple(wwave), fold_name,tuple(which_params))

#@lru_cache(maxsize=10)
def grid_mods(wwave,fold_name,which_params):

    dic = {}
    master_param=[]
    which_teffs, which_loggs, which_comps, which_siabs = [], [], [], []
    #ficheros=glob
    ext='dk'
    for fichero in glob.glob(fold_name+'*.'+ext):
        tic = time.perf_counter()
        #which_params,wave_model,flux_model,model_name=[],[],[],[]
        print(fichero)
        if ext=='dk':   #This method is still not generalized. It just takes the parameters from the name
            param_vals=[]
            wave_model,flux_model=np.loadtxt(fichero,skiprows=60,unpack=True)
            s=len(fold_name)
            teff_name,logg_name,comps_name,siabs_name=fichero[s:].split('.')[0].split('_')
            name=teff_name+','+logg_name+','+comps_name[1:]+','+siabs_name[2:]
            siabs=float(format(int(siabs_name[2:]),'.2f'))
            teff=int(teff_name)
            comps=float(format(int(comps_name[1:]),'.2f'))
            logg=float(format(int(logg_name),'.2f'))
            param_vals=[teff,logg,comps,siabs]
            int_mod = interpolate.interp1d(wave_model,flux_model,kind='linear')
            dic[tuple(param_vals)] = int_mod(wwave)
            master_param.append(param_vals)

        elif ext=='fits':
            hdu=fits.open(fichero)
            if hdu[0].header['teff']==16100 or hdu[0].header['h']!=-4.1: continue
            else:
                param_vals,wave_model,flux_model,model_name=read_fits(which_params,fichero)
                int_mod = interpolate.interp1d(wave_model,flux_model,kind='linear')
                dic[tuple(param_vals)] = int_mod(wwave)
                master_param.append(param_vals)  #master param is the list of lists of parameters of each file

        else:
            raise Exception(ext+' extension not supported')

        toc = time.perf_counter()
        print('time: '+str(toc-tic)+' sec')

    wh_set,wh_len,wh_index=[],[],[]
    for which in np.array(master_param).transpose():  #This step is required for the itertools.product step
        param_wh=np.array(np.sort(list(set(which))))    #This is the set of all the non-repeated values of each parameter
        wh_set.append(param_wh)
        index=[]
        for k,val in enumerate(param_wh):
            index.append(k)                           #this is to obtain the indeces that supreme matrix will have

        wh_index.append(index)
        wh_len.append(len(index))                       #the len of each list of parameters (also needed to build supreme matrix)

    return(dic,tuple(wh_set),tuple(wh_index),wh_len)

def generator_mods(theta):

    teff, logg, comp,si = theta
    teff_wh=params_wh[0]

    teff_l, teff_u = max(teff_wh[teff_wh < teff]), min(teff_wh[teff_wh > teff])
    #logg_l, logg_u = max(logg_wh[logg_wh < logg]), min(logg_wh[logg_wh > logg])
    #comp_l, comp_u = max(comp_wh[comp_wh < comp]), min(comp_wh[comp_wh > comp])

    x_d = (teff-teff_l)/(teff_u-teff_l)

    c0 = mod_library[(teff_l,logg ,comp,-12.0 )]
    c1 = mod_library[(teff_u,logg ,comp,-12.0 )]

    new_flux= c0*(1-x_d) + c1*x_d

    return(new_flux)


steps_models = 0.15
wave_peo = np.arange(3100,6700,steps_models)
res_uvb, res_vis = 5400, 8900


params=['teff','logg','h']
#numbers = np.array([np.log10(16100),7.5,-4.1])
numbers=np.array([16100,750,410,1200])
#which_fold = 'Modelos_J0827/'
which_fold='InterpolationTest_DBA_2_dk/'
hdus=glob.glob(which_fold+'*.fits')
print('Uploading models')
mod_library, params_wh,index_wh,len_wh= grid_mods(wave_peo, which_fold,params)

print('Models uploaded')
print(params_wh)



perm_index, perm_names,perm = [], [],[]

for v in itertools.product(*index_wh):   #Notice index_wh has been unzipped (it was a tuple object)
    perm_index.append(v)

for v in itertools.product(*params_wh):  #Notice params_wh has been unzipped (it was a tuple object)
    perm_names.append(v)

def f2(j,k):
    value=mod_library[perm_names[j]][k]
    return value
method='linear'
#start = time.time()

if method=='linear':
    len_wh.append(len(wave_peo))
    supreme_matrix = np.zeros(len_wh)
    #supreme_matrix=np.empty((len_wh[0],len_wh[1],len_wh[2],len_wh[3],len(wave_peo)))
    for idx in range(len(perm_index)):
        #if
        print(perm_names[idx])
        supreme_matrix[perm_index[idx]] = mod_library[perm_names[idx]]
    interp=RegularGridInterpolator(params_wh, supreme_matrix, method = method)
    interp_spec=interp(numbers)[0]

elif method=='cubic' or method=='quintic':

    s=len(wave_peo)
    interp_spec=np.zeros(s)
    for k in range(s):
        tic = time.perf_counter()
        supreme_matrix = np.zeros(len_wh)
        for idx in range(len(perm_index)):
            supreme_matrix[perm_index[idx]]=f2(idx,k)
            interp = RegularGridInterpolator(params_wh, supreme_matrix, method = method)
            interp_value=interp(numbers)[0]
        interp_spec[k]=interp_value
        toc = time.perf_counter()
        print("element--"+str(k)+"---"+str(toc-tic)+" sec")
else:
    raise Exception(method+' Method not supported.')

np.savetxt("regulargrid_200.txt",interp_spec)
# PLOT RESULTS
#hdu_model=fits.open(which_fold+'16100_750_H410_Si1200.fits')
#flux_mod=hdu_model[1].data['flux']
#wave_mod=hdu_model[1].data['wavelength']
#int_mod=interpolate.interp1d(wave_mod,flux_mod,kind='linear')
#mod_library[tuple(numbers)]=int_mod(wave_peo)

#paula_interp = generator_mods(numbers)

fig = plt.figure()
gs = fig.add_gridspec(2, hspace=0)
axs = gs.subplots(sharex=True)
#axs[0].plot(wave_peo,mod_4D_teff_logg_comp_si,'r-',label='Paula linear_interp')
axs[0].plot(wave_peo,mod_library[(16200,750,410,1200)], label='16200')
axs[0].plot(wave_peo,mod_library[(16000,750,410,1200)],label='16000')
axs[0].plot(wave_peo,interp_spec,'k-',label='RegularGridInterp linear')
axs[0].plot(wave_peo,mod_library[tuple(numbers)],label='model teff:16100')
#plt.title("LOG10(Teff)")
axs[0].legend()
axs[1].plot(wave_peo,interp_spec-mod_library[tuple(numbers)],'b.',label='difference')
axs[1].legend()
plt.show()

