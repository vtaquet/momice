#
##---------------------------------
##---------------------------------
## input parameters
##---------------------------------
##---------------------------------
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
import data_analysis
from scipy.io import FortranFile
from pathlib import Path
#
##################################################################################
##################################################################################
## useful functions
##################################################################################
##################################################################################
#
# check where is the comment starting with #
def checkcomment(inplist):
  #typelist = [type(value) for value in inplist]
  for i, value in enumerate(inplist):
    if value == '#':
      return i
#
# check where is the 1st float in a list
def check1rstfloat(inplist):
  #typelist = [type(value) for value in inplist]
  for i, value in enumerate(inplist):
    try:
      float(value)
      return value
    except:
      pass    
#
# read 1d binary files (time, nH, T, Av)  
def read_bin_1d(fileinp,Nstep,dtype):
    fileid = open(fileinp,'rb')
    init = np.fromfile(fileid, dtype='uint32', count=1)
    data = np.fromfile(fileid, dtype=dtype, count=Nstep)
    fin = np.fromfile(fileid, dtype='uint32', count=1)
    return list(data)
# read 2d binary files (abundances)  
def read_bin_2d(fileinp,Nstep,Nsp,dtype):
    fileid = open(fileinp,'rb')
    init = np.fromfile(fileid, dtype='uint32', count=1)
    data = np.fromfile(fileid, dtype=dtype, count=Nstep*Nsp).reshape(Nstep,Nsp)
    fin = np.fromfile(fileid, dtype='uint32', count=1)
    return (data)
# read 3d binary files (reaction rates)  
def read_bin_3d(fileinp,Nstep,Nreac,Nsprate,dtype):
    fileid = open(fileinp,'rb')
    init = np.fromfile(fileid, dtype='uint32', count=1)
    data = np.fromfile(fileid, dtype=dtype, count=Nstep*Nreac*Nsprate).reshape(Nstep,Nreac,Nsprate)
    fin = np.fromfile(fileid, dtype='uint32', count=1)
    return (data)
#
##################################################################################
##################################################################################
## read input file and store input params into dataframe
##################################################################################
##################################################################################
def read_input():
  #
  inpascii = open("moman.in",'r')
  inpparam = [] ; inpvalue = []
  for line in inpascii:
    s = line.split()
    wherecom = [spl.find('#') for spl in s]
    if s[0] != '#' and s[0] != '*':
      # convert params to float if possible
      try:
        if s[0] != 'outmod':
          s[2:checkcomment(s)] = [float(i) for i in s[2:checkcomment(s)]]
      except:
        pass
      inpparam.append(s[0])
      inpvalue.append(s[2:checkcomment(s)])
    elif s[0] == '*':
      break
  #
  # create input params dataframe
  inp_df = pd.DataFrame({inpparam[i] : [inpvalue[:][i]] for i in range(len(inpparam))},columns=inpparam).loc[0]
  inp_df['species2plot'] = inp_df['species'].copy()
  for isp, sp in enumerate(inp_df['species2plot']):
    for ie in range(15):
      inp_df['species2plot'][isp] = inp_df['species2plot'][isp].replace(str(ie),'$_{'+str(ie)+'}$')
  #
  # useful integers
  Nsp2pl = len(inp_df['species'])
  #
  # invidiual models
  if inp_df.typeout[0] == 1:
    Nmodels = len(inp_df.outmod)
  #
  # spatial evolution
  elif inp_df.typeout[0] == 2:
    Nmodels = 1
  #
  # model grid
  elif inp_df.typeout[0] == 3:
    Nmodels = 1
  #
  # sensitivity analysis
  elif inp_df.typeout[0] == 4:
    Nmodels = 1
  #
  return Nmodels, Nsp2pl, inp_df

  #
#
##################################################################################
##################################################################################
## read spatial log file
##################################################################################
##################################################################################
def read_logspat(inp_df): 
  inp_df.outdir[0] = inp_df.outdir[0]+'/'+inp_df.outmod[0] ; inp_df.outmod = []
  inplogspat = open(inp_df.outdir[0]+'/spatlog.out')
  for il, line in enumerate(inplogspat):
    s = line.split()
    if il == 0:
      Nmodels = int(s[-1])
    elif il >= 1:
      inp_df.outmod.append(s[0])
      
  return Nmodels
#
##################################################################################
##################################################################################
## read grid log file
##################################################################################
##################################################################################
def read_loggrid(inp_df): 
  inp_df.outdir[0] = inp_df.outdir[0]+'/'+inp_df.outmod[0] ; inp_df.outmod = []
  inploggrid = open(inp_df.outdir[0]+'/gridlog.out')
  param = [] ; values = [] ; Nrealparams = 0
  for il, line in enumerate(inploggrid):
    s = line.split()
    if il == 0:
      Nparams, Nmodels = int(s[0]), int(s[1])
    elif il >= 1 and il <= Nparams:
      if len(s[1:]) > 1:
        param.append(s[0]) ; values.append(s[1:])
        Nrealparams += 1
    else:
      inp_df.outmod.append(s[1])
  grid_df = pd.DataFrame(list(map(list, zip(*[param, values]))), columns=['Param', 'Values'])#.set_index(['Model','N','Species'],drop=False)
  
  return Nmodels, grid_df
#
##################################################################################
##################################################################################
## read sensitivity log file
##################################################################################
##################################################################################
def read_logsens(inp_df): 
  inp_df.outdir[0] = inp_df.outdir[0]+'/'+inp_df.outmod[0] ; inp_df.outmod = []
  inplogsens = open(inp_df.outdir[0]+'/sensitivitylog.out') ; freeparam = []
  for il, line in enumerate(inplogsens):
    s = line.split()
    if il == 0:
      Nparams, Nmodels = int(s[0]), int(s[1])
    elif il >= 1 and il <= Nparams:
      freeparam.append(s[0])
      #print(il, Nmodels)
    else:
      inp_df.outmod.append(s[0])
  inp_df['freeparam'] = freeparam
      
  return Nmodels
#
##################################################################################
##################################################################################
## read log files and save physical, surface, and reaction parameters into dataframe
##################################################################################
##################################################################################
def read_log(Nmodels, inp_df): 
  #
  Nspnet = [] ; Nrenet = [] ; s_sp = [] ; s_re = [] ; Nspgas = [] ; Nspice = [] ; Nsprate = [] ; Nrerate = []
  for im in range(Nmodels):
    print(str(im)+'/'+str(Nmodels), end="\r")
    inplog = open(inp_df.outdir[0]+'/'+inp_df.outmod[im]+'/log.out')
    modparam = [] ; modval = [] ; startsp = 0  ; startre = 0 ; checkgr = 0
    ire = 0 ; isp = 0  ; s_sp2 = [] ; s_re2 = [] ; spec_list = []# Nspnet = 0 ; Nrenet = 0 ; 
    for line in inplog:
      s = line.split()
      if s[0][0] != '-' and s[0][0] != '!':
        #
        # physical conditions and grain surface properties
        if s[1] == '!':
          s[2] = s[2].replace(":","")
          modparam.append(s[2]) 
          modval.append(float(s[0]))
        elif line.find('timesteps') > 0:
          modval.append(float(check1rstfloat(s)))
          modparam.append('Nsteps')
        else:
          #
          # Nspecies
          if line.find('SPECIES NETWORK') > 0 and startsp == 0:
            Nspnet.append(int(check1rstfloat(s)))
            Nsprate.append(int(s[-1]))
            startsp = 1
          # species properties
          elif startsp == 1 and isp <= Nspnet[im]:
            if isp == 0:
              spec_head = s #[s2+'_'+str(im) for s2 in s] #['Model'] + 
            else:
              s[0] = int(s[0])
              s[2:] = [float(s2) for s2 in s[2:]]
              s_sp2.append(s)
              spec_list.append(s[1])
              if s[1][0:1] == 'J' and checkgr == 0:
                Nspgas.append(int(s[0])-1)
                Nspice.append(int(Nspnet[im]-Nspgas[im]))
                checkgr = 1
                #exit()
              if s[1][0:1] == 'Q' and checkgr == 1:
                Nspice[im] = int(Nspice[im]/2)
                checkgr += 1
            isp += 1
          #
          # Nreactions
          elif line.find('CHEMICAL NETWORK') > 0 and startre == 0:
            #print(line)
            Nrenet.append(float(check1rstfloat(s)))
            Nrerate.append(int(s[-1]))
            startre = 1
          # reaction properties
          elif startre == 1 and ire <= Nrenet[im]:
            if ire == 0:
              reac_head = s #['Model'] + 
            else:
              s_reac = line[6:51].split() ; s_prod = line[51:113].split() ; s_float = s[-5:]
              s_float = [float(s2) for s2 in s_float]
              if len(s_reac) == 2:
                s_reac = s_reac+[' ']
              elif len(s_reac) == 1:
                s_reac = s_reac+[' ',' ']
              if len(s_prod) == 3:
                s_prod = s_prod+[' ']
              if len(s_prod) == 2:
                s_prod = s_prod+[' ',' ']
              elif len(s_prod) == 1:
                s_prod = s_prod+[' ',' ',' ']
              s_re2.append([int(s[0])]+s_reac+s_prod+s_float)
            ire += 1
    # check if duplicate species 
    sp2rem = -1 ; duplsp = 0
    spec_list = [s_sp2[isp][2] for isp in range(Nspnet[im])]
    if len(spec_list) != len(list(set(spec_list))):
      duplsp = 1
      for isp in range(Nspnet[0]):
        if len(spec_list[0:isp]) != len(list(set(spec_list[0:isp]))):
          s_sp2[isp-1][1] = s_sp2[isp-1][1]+'_2'
          break
    #
    # create dataframe for 1) physical conditions and grain surface properties, 2) species, 3) reactions
    s_par2 = [[im for i in range(len(modval))],modparam,modval]
    modparam_df2 = pd.DataFrame(list(map(list, zip(*s_par2))), columns=['Model','Param','Value']).set_index(['Model','Param'])
    spec_df2 = pd.DataFrame(s_sp2, columns=spec_head).set_index(np.arange(Nspnet[0])+1) #['N'],drop=False)
    reac_df2 = pd.DataFrame(s_re2, columns=reac_head).set_index(['Num'],drop=False)
    #
    if im == 0:
      modparam_df = modparam_df2
      spec_df = spec_df2 ; spec_df = spec_df.rename(columns={'Eb_wat': 'Eb_wat_0', 'Eb_carb': 'Eb_carb_0'})
      reac_df = reac_df2 ; reac_df = reac_df.rename(columns={'C' : 'C_0'})
    else:
      modparam_df = modparam_df.append(modparam_df2) 
      list_inp = ['Eb_carb_'+str(im), 'Eb_wat_'+str(im), 'Xini_'+str(im)]
      spec_df[list_inp] = spec_df2[['Eb_carb', 'Eb_wat', 'Xini']]
      reac_df['C_'+str(im)] = reac_df2['C']
      #spec_df = spec_df.append(spec_df2) 
      #reac_df = reac_df.append(reac_df2)
  #
  return Nspgas, Nspice, Nsprate, Nrerate, spec_df, reac_df, modparam_df

#
##################################################################################
##################################################################################
## read bin files and save temporal evolution of physical and chemical parameters into dataframes
##################################################################################
##################################################################################
def read_bin(Nmodels, Nspgas, Nspice, Nsprate, Nrerate, modparam_df, inp_df, spec_df): 
  #
  # df showing physical properties (nH, Td, Tg, Av, zeta, ...), grain properties (Nlayers, ad, ...) abundances (X_sp, D_sp)
  #
  physfile = ['time','Av','nH','Tdust','Tgas','Nlayers', 'radius']
  chemfile = ['Xgas','Xice','Xsur']
  spgas = [] ; spice = [] ; im2=0
  Nchunks = int(Nmodels/100.)+1 
  num_chunk = [int((ic+1)*Nmodels/Nchunks)-1 for ic in range(Nchunks)] ;num_chunk2 = [int((ic)*Nmodels/Nchunks) for ic in range(Nchunks)]
  #
  # loop over models
  for im in range(Nmodels):
    print(str(im)+'/'+str(Nmodels), end="\r")
    physdata = [] ; chemdata = []
    Nsteps = int(modparam_df.at[(im,'Nstep'),'Value'])
    spgas.append(list(spec_df['Species'])[0:Nspgas[im]])
    spice.append(list(spec_df['Species'])[Nspgas[im]:Nspgas[im]+Nspice[im]])
    spice[im] = [spice2[1:] for spice2 in spice[im]]
    physdata.append([im for i in range(Nsteps)])
    physdata.append([i for i in range(Nsteps)])
    #
    # physical conditions
    for phys in physfile:
      inpfile = inp_df.outdir[0]+'/'+inp_df.outmod[im]+'/'+phys+'.out'
      physdata.append(read_bin_1d(inpfile,Nsteps,'float64'))
    data_df_ind = pd.DataFrame(list(map(list, zip(*physdata))),columns=['Model','Nsteps']+physfile)#.set_index(['Model','time'])
    #
    # chemical properties
    for chem in chemfile:
      inpfile = inp_df.outdir[0]+'/'+inp_df.outmod[im]+'/'+chem+'.out'
      if chem == 'Xgas': 
        Nsp = Nspgas[im]
        sp2plot = spgas[im]
      else:
        Nsp = Nspice[im]
        sp2plot = spice[im]
      #chemdata.append(read_bin_2d(inpfile,Nsteps,Nsp))
      # create dataframe for individual model
      chemdata2 = read_bin_2d(inpfile,Nsteps,Nsp,'float64')
      chem_df_2 = pd.DataFrame(chemdata2,columns=[chem+'_'+sp for sp in sp2plot]) #.set_index(['Model','time'])
      data_df_ind = data_df_ind.join(chem_df_2)
    #print(check_rate == True)
    #data_df_3 = data_df_2#[list_col]
    # save intermediate dataframes in which abundances for input species are saved
    if Nmodels > 1:
      if im in num_chunk and im > 0:
        data_df_interm = data_df_interm.append(data_df_ind)
        choice_D, choice_opr, data_df_interm = data_analysis.compute_params(Nmodels, Nspgas, Nspice, inp_df, spec_df, data_df_interm)
        chemfile2 = chemfile+['Xgasref', 'Xiceref']
        if choice_D >= 1:
          chemfile2 = chemfile2+['Dgas', 'Dice']
        if choice_opr >= 1:
          chemfile2 = chemfile2+['H2opr']
        list_col = ['Model','Nsteps']+physfile+[ord2+'_'+sp2pl for iord, ord2 in enumerate(chemfile2) for isp, sp2pl in enumerate(inp_df['species']) if ord2+'_'+sp2pl in data_df_interm.columns]
        data_df_interm = data_df_interm[list_col]
        if im == num_chunk[0]:
          data_df = data_df_interm
        else:
          data_df = data_df.append(data_df_interm)
        im2=im2+1
      else:
        # create/append full dataframe
        if im in num_chunk2:
          data_df_interm = data_df_ind
        else:
          data_df_interm = data_df_interm.append(data_df_ind)
    else:
      choice_D, choice_opr, data_df_ind = data_analysis.compute_params(Nmodels, Nspgas, Nspice, inp_df, spec_df, data_df_ind)
      data_df = data_df_ind
  #
  # define model number and timesteps as index
  data_df = data_df.set_index(['Model','Nsteps'])
  #print(data_df.head())
  return choice_D, choice_opr, data_df

def read_rates(Nmodels, Nspgas, Nspice, Nsprate, Nrerate, modparam_df, inp_df, spec_df):
  ratefile3d = ['reac_numdest', 'reac_numform', 'reac_ratedest', 'reac_rateform'] #
  ratefile2d = ['reac_contribdest', 'reac_contribform']
  ratefile1d = ['reac_spnum']
  typerate3d = ['uint32', 'uint32', 'float64', 'float64']
  #
  # loop over models
  for im in range(Nmodels):
    Nsteps = int(modparam_df.at[(im,'Nstep'),'Value'])
    # number of species whose formation/destruction have been computed
    for rate in ratefile1d:
      inpfile = inp_df.outdir[0]+'/'+inp_df.outmod[im]+'/'+rate+'.out'
      sp2rate = spec_df.Species[read_bin_1d(inpfile, Nsprate[im],'uint32')].tolist()
    # create two lists of dfs describing formation/destruction rates for each concerned species
    list_rate2 = []
    for irate, rate in enumerate(ratefile3d):
      inpfile = inp_df.outdir[0]+'/'+inp_df.outmod[im]+'/'+rate+'.out'
      list_rate2.append(read_bin_3d(inpfile, Nsteps, Nrerate[im], Nsprate[im],typerate3d[irate])) #Nstep,Nreac,Nsprate
    if im == 0:
      list_ratform = [] ; list_ratdest = [] ; list_numform = [] ; list_numdest = []
    numdest2, numform2 = list_rate2[:2] #pd.DataFrame()
    ratedest2, rateform2 = list_rate2[2:] 
    ratedata = []
    for isp2rate, sp2rate2 in enumerate(sp2rate):
      dfnumform2, dfratform2 = pd.DataFrame(numform2[:,:,isp2rate]), pd.DataFrame(rateform2[:,:,isp2rate])
      dfnumform2['Model'], dfratform2['Model'] = [im for istep in range(Nsteps)], [im for istep in range(Nsteps)]
      dfnumform2['Nsteps'], dfratform2['Nsteps'] = [istep for istep in range(Nsteps)], [istep for istep in range(Nsteps)]
      dfnumdest2, dfratdest2 = pd.DataFrame(numdest2[:,:,isp2rate]), pd.DataFrame(ratedest2[:,:,isp2rate])
      dfnumdest2['Model'], dfratdest2['Model'] = [im for istep in range(Nsteps)], [im for istep in range(Nsteps)]
      dfnumdest2['Nsteps'], dfratdest2['Nsteps'] = [istep for istep in range(Nsteps)], [istep for istep in range(Nsteps)]
      if im == 0:
        list_ratform.append(dfratform2) ; list_numform.append(dfnumform2)
        list_ratdest.append(dfratdest2) ; list_numdest.append(dfnumdest2)
      elif im > 0:
        list_ratform[isp2rate] = list_ratform[isp2rate].append(dfratform2)
        list_ratdest[isp2rate] = list_ratdest[isp2rate].append(dfratform2)
  for isp2rate, sp2rate2 in enumerate(sp2rate):
    list_ratform[isp2rate] =  list_ratform[isp2rate].set_index(['Model','Nsteps'])
    list_ratdest[isp2rate] =  list_ratdest[isp2rate].set_index(['Model','Nsteps'])
  #
  return sp2rate, list_ratform, list_ratdest, list_numform, list_numdest