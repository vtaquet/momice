#############################################################################################
# moman.py - MOMice ANalyser
# reads and analyses the data files from MOMICE and plots the results in postcript files
#############################################################################################
#
from scipy.interpolate import interp1d
from numpy import array, arange, sin
import numpy as np
import math as math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import os
import struct
import read_input
import data_analysis
import plot
from scipy.io import FortranFile
from pathlib import Path
#
print('!! MOMAN: MOMice ANalyser !!')
#
##---------------------------------
##---------------------------------
## input parameters
##---------------------------------
##---------------------------------
#
# # read input file and store input params into dataframe
print('Reading input file...')
Nmodels, Nsp2pl, inp_df = read_input.read_input()
#
##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
## read log files and save physical, surface, and reaction parameters into dataframe
##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
#
print('Reading log files...')
#
if inp_df.typeout[0] == 2:
  Nmodels = read_input.read_logspat(inp_df)
elif inp_df.typeout[0] == 3:
  Nmodels, grid_df = read_input.read_loggrid(inp_df)
elif inp_df.typeout[0] == 4:
  Nmodels = read_input.read_logsens(inp_df)
inp_df['Nmodels'] = Nmodels
#
Nspgas, Nspice, Nsprate, Nrerate, spec_df, reac_df, modparam_df = read_input.read_log(Nmodels, inp_df)
#print(modparam_df.loc[[(i,'nH') for i in range(Nmodels)],'Value'])
#print(modparam_df.loc[(1,'nH'):(3,'nH')])#,'Value'])
#print(modparam_df[modparam_df.index.get_level_values('Param') == 'nH']['Value'])
#print(inp_df)
#exit()
#
# 
#
##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
## read bin files and save temporal evolution of physical and chemical parameters into dataframes
##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
#
print('Reading bin files...')
choice_D, choice_opr, data_df = read_input.read_bin(Nmodels, Nspgas, Nspice, Nsprate, Nrerate, modparam_df, inp_df, spec_df)
inp_df['choice_D'] = choice_D ; inp_df['choice_opr'] = choice_opr
print('Size of all Dataframes:')
print(spec_df.memory_usage(index=True).sum(),\
	  reac_df.memory_usage(index=True).sum(),\
	  data_df.memory_usage(index=True).sum())
#
##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
## perform analysis of temporal evolution of abundances with dataframes
##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
#
print('Performing some data analysis...')
# compute abundance / reference species, deuteration, H2 opr, 
#choice_D, choice_opr, data_df = data_analysis.compute_params(Nmodels, Nspgas, Nspice, inp_df, spec_df, data_df)
#
if inp_df.typeout[0] == 3:
	data_analysis.PCA(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df)
#
##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
## plot results
##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
#
print('Making some plots...')
if inp_df.typeout[0] == 1:
  plot.plot_indmod(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df)
elif inp_df.typeout[0] == 2:
  plot.plot_spat(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df)
elif inp_df.typeout[0] == 3:
  plot.plot_evolparam(Nmodels, inp_df, spec_df, reac_df, data_df, grid_df, modparam_df)
  plot.plot_wholegrid(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df)
elif inp_df.typeout[0] == 4:
  #plot.plot_indmod(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df)
  plot.plot_mean(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df)
  plot.plot_freeparam(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df)
  plot.plot_wholegrid(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df)
#
##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
## for individual simulations, read and plot reaction rates
##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
#
check_rate = Path(inp_df.outdir[0]+'/'+inp_df.outmod[0]+'/reac_numdest.out').is_file()
if check_rate == True and inp_df.plotreac[0] == 'yes' and inp_df.typeout[0] == 1:
  print('Reading and plotting reaction rates... (can take some time due to generation of high number of figures)')
  sp2rate, list_ratform, list_ratdest, list_numform, list_numdest = read_input.read_rates(Nmodels, Nspgas, Nspice, Nsprate, Nrerate, modparam_df, inp_df, spec_df)
  plot.plot_rates(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df, sp2rate, list_ratform, list_ratdest, list_numform, list_numdest)

#
print('Done!')
exit()












