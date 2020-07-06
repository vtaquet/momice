##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
## plot results from individual models
##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
#
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
from scipy.io import FortranFile
import multiprocessing, subprocess
import time
#
# plot properties
def plot_prop():
	hfont = {'family' : 'DejaVu Sans',
	         'weight' : 'medium',
	         'size'   : 12}
	haxes = {'linewidth' : 1.5}
	hticks = {'major.size' : 6,
	          'major.width' : 1.5,
	          'minor.size' : 3,
	          'minor.width' : 1}
	plt.rc('font', **hfont)
	plt.rc('axes',**haxes)
	plt.rc('xtick',**hticks)
	plt.rc('ytick',**hticks)
	params = {'mathtext.default': 'regular' }       
	colorsc = ['black','salmon','royalblue','seagreen','gold','darkorange','violet','darkturquoise','darkcyan','turquoise']   
	linestyle = ['-','--',':','-.','steps','-','--',':','-.','steps']
	plt.rcParams.update(params)
	return params, colorsc, linestyle
def sp2plot(spname):
	for ie in range(15):
		spname = spname.replace(str(ie),'$_{'+str(ie)+'}$')
	spname = spname.replace('+','$^+$')
	spname = spname.replace('-','$^-$')
	return spname
#
# define properties of axe variables
def axes_var(choice,inp_df, data_df, modparam_df, refspecies2pl):
	if choice == 1:
		# properties of x and y axes for plot  
		colphysplot = ['Param','Title','log','Min','Max']
		physplot_list = [['time','Time [yr]',inp_df['logabs'][0], 1., data_df['time'].max()],\
		                 ['Av','A$_V$ [mag]','log',1.,data_df['Av'].max()],\
		                 ['nH','n$_H$ [cm$^{-3}$]','log',1e3,data_df['nH'].max()],\
		                 ['Tdust','T$_{dust}$ [K]','log',1.,data_df['Tdust'].max()],\
		                 ['Tgas','T$_{gas}$ [K]','log',1.,data_df['Tgas'].max()],\
		                 ['radius','Radius [AU]','log',1.,data_df['radius'].max()],\
		                 ['Nlayers','N$_{layers}$ [MLs]','linear',0.,data_df['Nlayers'].max()]]
		colchemplot = ['Param','Title','log','Min','Max','Choice']
		chemplot_list = [['Xgas','X$_{gas}$(i)','log',1e-12,1e-3,1],\
		                 ['Xice','X$_{ice}$(i)','log',1e-12,1e-3,1],\
		                 ['Xsur','X$_{sur}$(i)','log',1e-6,1e0,1],\
		                 ['Xgasref','X$_{gas}$(i)/X$_{gas}$('+refspecies2pl+')','log',1e-3,1e3,1],\
		                 ['Xiceref','X$_{gas}$(i)/X$_{gas}$('+refspecies2pl+')','log',1e-3,1e3,1],\
		                 ['Dgas','(D/H)$_{gas}$(i) )','log',1e-5,1e1,inp_df.choice_D],\
		                 ['Dice','(D/H)$_{ice}$(i)','log',1e-5,1e1,inp_df.choice_D],\
		                 ['H2opr','o/p(H$_2$)','log',1e-6,1e1,inp_df.choice_opr]]
	elif choice == 2:
		# properties of x and y axes for plot  
		colphysplot = ['Param','Title','log','Min','Max']
		physplot_list = [['time','Time [yr]',inp_df['logabs'][0], 1., data_df['time'].max()],\
		                 ['Av','A$_V$ [mag]','log',1.,data_df['Av'].max()],\
		                 ['nH','n$_H$ [cm$^{-3}$]','log',1e3,data_df['nH'].max()],\
		                 ['Tdust','T$_{dust}$ [K]','log',1.,data_df['Tdust'].max()],\
		                 ['Tgas','T$_{gas}$ [K]','log',1.,data_df['Tgas'].max()],\
		                 ['radius','Radius [AU]','log',1.,data_df['radius'].max()],\
		                 ['Nlayers','N$_{layers}$ [MLs]','linear',0.,data_df['Nlayers'].max()]]
		colchemplot = ['Param','Title','log','Min','Max','Choice']
		chemplot_list = [['Xgas','log( X$_{gas}$(i) )','linear',-12,-3,1],\
		                 ['Xice','log( X$_{ice}$(i) )','linear',-12,-3,1],\
		                 ['Xsur','log( X$_{sur}$(i) )','linear',-6,0,1],\
		                 ['Xgasref','log( X$_{gas}$(i)/X$_{gas}$('+refspecies2pl+') )','linear',-3,3,1],\
		                 ['Xiceref','log( X$_{gas}$(i)/X$_{gas}$('+refspecies2pl+') )','linear',-3,3,1],\
		                 ['Dgas','log( (D/H)$_{gas}$(i) )','linear',-5,1,inp_df.choice_D],\
		                 ['Dice','log( (D/H)$_{ice}$(i) )','linear',-5,1,inp_df.choice_D],\
		                 ['H2opr','log( o/p(H$_2$) )','linear',-6,1,inp_df.choice_opr]]
	elif choice == 3:
		# properties of x and y axes for plot  
		colphysplot = ['Param','Title','log','Min','Max'] 
		physplot_list = [['nH','log(n$_H$ [cm$^{-3}$])','log',
						  np.log10(modparam_df.loc[[(i,'nH') for i in range(inp_df.Nmodels)],'Value'].min()),
						  np.log10(modparam_df.loc[[(i,'nH') for i in range(inp_df.Nmodels)],'Value'].max())],\
		                 ['Tg','T$_{gas}$ [K]','linear',
		                  modparam_df.loc[[(i,'Tg') for i in range(inp_df.Nmodels)],'Value'].min(),
		                  modparam_df.loc[[(i,'Tg') for i in range(inp_df.Nmodels)],'Value'].max()],\
		                 ['Td','T$_{dust}$ [K]','linear',
		                  modparam_df.loc[[(i,'Td') for i in range(inp_df.Nmodels)],'Value'].min(),
		                  modparam_df.loc[[(i,'Td') for i in range(inp_df.Nmodels)],'Value'].max()],\
		                 ['CR','Zeta [s$^{-1}$]','log',
		                  modparam_df.loc[[(i,'zeta') for i in range(inp_df.Nmodels)],'Value'].min(),
		                  modparam_df.loc[[(i,'zeta') for i in range(inp_df.Nmodels)],'Value'].max()],\
		                 ['Av','A$_{V}$ [mag]','linear',0.,modparam_df.loc[[(i,'Av') for i in range(inp_df.Nmodels)],'Value'].max()],\
		                 ['G0','G$_{0}$','linear',0.,modparam_df.loc[[(i,'G0') for i in range(inp_df.Nmodels)],'Value'].max()],\
		                 ['UV','UV Flux [cm$^{-2}$ s$^{-1}$]','log',
		                  modparam_df.loc[[(i,'CRUVf') for i in range(inp_df.Nmodels)],'Value'].min(),
		                  modparam_df.loc[[(i,'CRUVf') for i in range(inp_df.Nmodels)],'Value'].max()]]
		colchemplot = ['Param','Title','log','Min','Max','Choice']
		chemplot_list = [['Xgas','X$_{gas}$(i)','log',1e-12,1e-3,1],\
		                 ['Xice','X$_{ice}$(i)','log',1e-12,1e-3,1],\
		                 ['Xsur','X$_{sur}$(i)','log',1e-6,1e0,1],\
		                 ['Xgasref','X$_{gas}$(i)/X$_{gas}$('+refspecies2pl+')','log',1e-3,1e3,1],\
		                 ['Xiceref','X$_{gas}$(i)/X$_{gas}$('+refspecies2pl+')','log',1e-3,1e3,1],\
		                 ['Dgas','(D/H)$_{gas}$(i) )','log',1e-5,1e1,inp_df.choice_D],\
		                 ['Dice','(D/H)$_{ice}$(i)','log',1e-5,1e1,inp_df.choice_D],\
		                 ['H2opr','o/p(H$_2$)','log',1e-6,1e1,inp_df.choice_opr]]
	#
	physplot_df = pd.DataFrame({colphysplot[i] : list(map(list, zip(*physplot_list)))[i][:] for i in range(len(colphysplot))}).set_index('Param')
	chemplot_df = pd.DataFrame({colchemplot[i] : list(map(list, zip(*chemplot_list)))[i][:] for i in range(len(colchemplot))}).set_index('Param')
	#
	return physplot_df, chemplot_df
#
##################################################################################
##################################################################################
# plot individual models
##################################################################################
##################################################################################
def plot_indmod(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df):
	#
	refspecies2pl = inp_df['refsp'][0]
	for ie in range(15):
	  refspecies2pl = refspecies2pl.replace(str(ie),'$_{'+str(ie)+'}$')
	#
	physplot_df, chemplot_df = axes_var(1,inp_df, data_df, modparam_df, refspecies2pl)
	params, colorsc, linestyle = plot_prop()
	#
	time0 = time.time() ; Nplot = 0
	# create figures folder
	for im in range(Nmodels):
	  #os.system('mkdir '+inp_df.outdir[0]+'/'+inp_df.outmod[im]+'/figures/')
	  subprocess.Popen(['mkdir',inp_df.outdir[0]+'/'+inp_df.outmod[im]+'/figures/'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#
	for iabs, abs2 in enumerate(physplot_df.index.values): 
	  print(abs2)
	  #
	  # physical properties
	  #
	  # plot individual models
	  if Nmodels == 1:
		  for im in range(Nmodels):
		    #
		    plt.figure(1,figsize=[6,6])
		    plt.xlabel(physplot_df.loc[abs2,'Title'])
		    plt.ylabel('n$_H$ [cm$^{-3}$, $T_{dust}$ [K], $T_{gas}$ [K}, A$_V$ [mag]')
		    plt.xscale(physplot_df.loc[abs2,'log'])
		    plt.yscale('log')
		    plt.xlim(data_df.loc[im,abs2].min(),data_df.loc[im,abs2].max())
		    plt.ylim(1e-1,1e2) #data_df.loc[im,'time'].min(),data_df.loc[im,'time'].max())
		    for iph, ph2pl in enumerate(physplot_df.index.values):
		    	if ph2pl != 'nH':
		    		plt.plot(data_df.loc[im,abs2],data_df.loc[im,ph2pl],\
		               		 color=colorsc[iph],lw=1.5,label=physplot_df.ix[iph,'Title'])
		    	elif ph2pl == 'nH':
		    		plt.plot(data_df.loc[im,abs2],data_df.loc[im,ph2pl]/1e5,\
		               		 color=colorsc[iph],lw=1.5,label=physplot_df.ix[iph,'Title'])
		    namefile = inp_df.outdir[0]+'/'+inp_df.outmod[im]+'/figures/physprop_'+abs2
		    plt.legend()
		    plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; Nplot += 1
		    plt.close(1)
	  #
	  # compare individual models 
	  elif Nmodels > 1:
		  plt.figure(1,figsize=[6,6])
		  plt.xlabel(physplot_df.loc[abs2,'Title'])
		  plt.ylabel('n$_H$ [cm$^{-3}$, $T_{dust}$ [K], $T_{gas}$ [K}, A$_V$ [mag]')
		  plt.xscale(physplot_df.loc[abs2,'log'])
		  plt.yscale('log')
		  plt.xlim(data_df.loc[0,abs2].min(),data_df.loc[0,abs2].max())
		  plt.ylim(data_df.loc[0,'time'].min(),data_df.loc[0,'time'].max())
		  for im in range(Nmodels):
		    for iph, ph2pl in enumerate(physplot_df.index.values):
		      plt.plot(data_df.loc[im,abs2],data_df.loc[im,ph2pl],\
		               color=colorsc[iph],lw=1.5,ls='-',label=physplot_df.ix[iph,'Title'])
		  namefile = inp_df.outdir[0]+'/figures/'+inp_df.suffix[0]+'_physprop_'+abs2
		  plt.legend()
		  plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; Nplot += 1
		  plt.close(1)
	  #
	  # abundances
	  for iord, ord2 in enumerate(list(chemplot_df.index)): 
	    if chemplot_df.loc[ord2,'Choice'] >= 1:
	      #
	      # plot individual models
	      if Nmodels == 1:
	      	for im in range(Nmodels):
		        #
		        plt.figure(1,figsize=[6,6])
		        plt.xlabel(physplot_df.loc[abs2,'Title'])
		        plt.ylabel(chemplot_df.loc[ord2,'Title'])
		        plt.xscale(physplot_df.loc[abs2,'log'])
		        plt.yscale(chemplot_df.loc[ord2,'log'])
		        plt.xlim(data_df.loc[im,abs2].min(),data_df.loc[im,abs2].max())
		        plt.ylim(chemplot_df.loc[ord2,'Min'], chemplot_df.loc[ord2,'Max'])
		        for isp, sp2pl in enumerate(inp_df['species']):
		          if ord2+'_'+sp2pl in data_df.columns:
		            plt.plot(data_df.loc[im,abs2],data_df.loc[im,ord2+'_'+sp2pl],\
		                     color=colorsc[isp],lw=1.5,label=inp_df['species2plot'][isp])
		        plt.legend(loc='upper left')
		        namefile = inp_df.outdir[0]+'/'+inp_df.outmod[im]+'/figures/'+ord2+'_'+abs2
		        plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; Nplot += 1
		        plt.close(1)
	      #
	      # compare individual models 
	      if Nmodels > 1:
		      plt.figure(1,figsize=[6,6])
		      plt.xlabel(physplot_df.loc[abs2,'Title'])
		      plt.ylabel(chemplot_df.loc[ord2,'Title'])
		      plt.xscale(physplot_df.loc[abs2,'log'])
		      plt.yscale(chemplot_df.loc[ord2,'log'])
		      plt.xlim(data_df.loc[0,abs2].min(),data_df.loc[0,abs2].max())
		      plt.ylim(chemplot_df.loc[ord2,'Min'], chemplot_df.loc[ord2,'Max'])
		      for im in range(Nmodels):
		        for isp, sp2pl in enumerate(inp_df['species']):
		          if ord2+'_'+sp2pl in data_df.columns:
		          	if im == 0:
		          		plt.plot(data_df.loc[im,abs2],data_df.loc[im,ord2+'_'+sp2pl],\
		                     	 color=colorsc[isp],lw=1.5,ls='-',label=inp_df['species2plot'][isp])
		          	else:
		          		plt.plot(data_df.loc[im,abs2],data_df.loc[im,ord2+'_'+sp2pl],\
		                     	 color=colorsc[isp],lw=1.5,ls='-',label=None)
		      plt.legend([sp2pl for isp, sp2pl in enumerate(inp_df['species2plot'])],loc='upper left')
		      namefile = inp_df.outdir[0]+'/figures/'+inp_df.suffix[0]+'_'+ord2+'_'+abs2
		      plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; Nplot += 1
		      plt.close(1)
	time1 = time.time()
	print(str(iabs)+'/'+str(len(physplot_df.index.values)), end="\r")
#
##################################################################################
##################################################################################
# plot spatial evolution
##################################################################################
##################################################################################
def plot_spat(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df):
	#
	refspecies2pl = inp_df['refsp'][0]
	for ie in range(15):
	  refspecies2pl = refspecies2pl.replace(str(ie),'$_{'+str(ie)+'}$')
	#
	physplot_df, chemplot_df = axes_var(1,inp_df, data_df, modparam_df, refspecies2pl)
	params, colorsc, linestyle = plot_prop()
	#
	time0 = time.time() ; Nplot = 0
	#
	#os.system('mkdir '+inp_df.outdir[0]+'/figures/')
	subprocess.Popen(['mkdir',inp_df.outdir[0]+'/figures/'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#
	# extract the variables for considered times
	for it, time2 in enumerate(inp_df.timeplot):
		#print(data_df[data_df['time'] > time2].index.min()[1])
		data_time2 = data_df[data_df['time'] > time2] 
		data_time2 = data_time2.iloc[data_time2.index.get_level_values('Nsteps') == data_time2.index.get_level_values('Nsteps').min()]
		#print(data_time2.tail())
		#exit()
		if it == 0:
			data_time = data_time2
		else:
			data_time = data_time.append(data_time2)
	#exit()#
	# make the plots
	for iabs, abs2 in enumerate(physplot_df.index.values):
	  print(abs2, end="\r")
	  #
	  # physical properties
	  #
	  plt.figure(1,figsize=[6,6])
	  plt.xlabel(physplot_df.loc[abs2,'Title'])
	  plt.ylabel('n$_H$ [cm$^{-3}$, $T_{dust}$ [K], $T_{gas}$ [K}, A$_V$ [mag]')
	  plt.xscale(physplot_df.loc[abs2,'log'])
	  plt.yscale('log')
	  time2 = data_time.iloc[data_time.index.get_level_values('Nsteps') == it]['time']
	  plt.xlim(physplot_df.loc[abs2,'Min'], physplot_df.loc[abs2,'Max']) #x2.min(),x2.max())
	  plt.ylim(1,1e6) #time2.min(),time2.max())
	  for it, time2 in enumerate(inp_df.timeplot):
	    indextime = data_df[data_df['time'] > time2].index.min()[1]
	    for iph, ph2pl in enumerate(physplot_df.index.values):
	    	x2 = data_df.iloc[data_df.index.get_level_values('Nsteps') == indextime][abs2]
	    	y2 = data_df.iloc[data_df.index.get_level_values('Nsteps') == indextime][ph2pl]
	    	plt.plot(x2,y2,\
	               color=colorsc[iph],lw=1.5,ls=linestyle[it],label=physplot_df.ix[iph,'Title'])
	    if it == 0:
	    	plt.legend()
	  namefile = inp_df.outdir[0]+'/figures/'+inp_df.suffix[0]+'_physprop_'+abs2
	  plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; 
	  Nplot += 1
	  plt.close(1)
	  #
	  # abundances
	  for iord, ord2 in enumerate(list(chemplot_df.index)): 
	    if chemplot_df.loc[ord2,'Choice'] >= 1:
	      #
	      #
	      # compare individual models 
	      plt.figure(1,figsize=[6,6])
	      plt.xlabel(physplot_df.loc[abs2,'Title'])
	      plt.ylabel(chemplot_df.loc[ord2,'Title'])
	      plt.xscale(physplot_df.loc[abs2,'log'])
	      plt.yscale(chemplot_df.loc[ord2,'log'])
	      x2 = data_time.iloc[data_time.index.get_level_values('Nsteps') == it][abs2]
	      plt.xlim(physplot_df.loc[abs2,'Min'], physplot_df.loc[abs2,'Max']) #x2.min(),x2.max())
	      #plt.xlim(x2.min(),x2.max())
	      plt.ylim(chemplot_df.loc[ord2,'Min'], chemplot_df.loc[ord2,'Max'])
	      for it, time2 in enumerate(inp_df.timeplot):
	        indextime = data_df[data_df['time'] > time2].index.min()[1]
	        for isp, sp2pl in enumerate(inp_df['species']):
	          if ord2+'_'+sp2pl in data_df.columns:
	          	x2 = data_df.iloc[data_df.index.get_level_values('Nsteps') == indextime][abs2]
	          	y2 = data_df.iloc[data_df.index.get_level_values('Nsteps') == indextime][ord2+'_'+sp2pl]
	          	plt.plot(x2,y2,color=colorsc[isp],lw=1.5,ls=linestyle[it],label=inp_df['species2plot'][isp])
	        if it == 0:
	        	plt.legend(loc='upper left')
	      namefile = inp_df.outdir[0]+'/figures/'+inp_df.suffix[0]+'_'+ord2+'_'+abs2
	      plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; 
	      Nplot += 1
	      plt.close(1)
	time1 = time.time()
	print(Nplot/(time1-time0))
#
##################################################################################
##################################################################################
# plot whole model grid
##################################################################################
##################################################################################
def plot_wholegrid(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df):
	#
	refspecies2pl = inp_df['refsp'][0]
	for ie in range(15):
	  refspecies2pl = refspecies2pl.replace(str(ie),'$_{'+str(ie)+'}$')
	#
	physplot_df, chemplot_df = axes_var(2,inp_df, data_df, modparam_df, refspecies2pl)
	params, colorsc, linestyle = plot_prop()
	#
	time0 = time.time() ; Nplot = 0
	#
	subprocess.Popen(['mkdir',inp_df.outdir[0]+'/figures/'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#
	# make the plots
	#
	# distribution of all species and all models together
	for iord, ord2 in enumerate(list(chemplot_df.index)): 
	    if chemplot_df.loc[ord2,'Choice'] >= 1:
	      for it, time2 in enumerate(inp_df.timeplot):
	      	plt.figure(1,figsize=[6,6])
	      	plt.ylabel('N$_{models}$')
	      	plt.xlabel(chemplot_df.loc[ord2,'Title'])
	      	plt.xscale(chemplot_df.loc[ord2,'log'])
	      	plt.xlim((chemplot_df.loc[ord2,'Min']), (chemplot_df.loc[ord2,'Max']))
	      	indextime = data_df[data_df['time'] > time2].index.min()[1]
	      	for isp, sp2pl in enumerate(inp_df['species']):
	        	if ord2+'_'+sp2pl in data_df.columns:
	        		y2 = np.log10(data_df.iloc[data_df.index.get_level_values('Nsteps') == indextime][ord2+'_'+sp2pl])
	        		plt.hist(y2, range=[chemplot_df.loc[ord2,'Min'], chemplot_df.loc[ord2,'Max']],
	        				 bins=int(np.sqrt(Nmodels)),lw=1.5,color=colorsc[isp],ls=linestyle[it], alpha=0.3,label=inp_df['species2plot'][isp]) # x2,y2,color=colorsc[isp],lw=1.5,ls=linestyle[it],label=inp_df['species2plot'][isp])
	      	plt.legend(loc='upper left')
	      	namefile = inp_df.outdir[0]+'/figures/'+inp_df.suffix[0]+'_'+ord2+'_'+"{:.2e}".format(time2)+'_hist'
	      	plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; 
	      	Nplot += 1
	      	plt.close(1)
	time1 = time.time()
	print(Nplot/(time1-time0))
#
##################################################################################
##################################################################################
# plot whole model grid
##################################################################################
##################################################################################
def plot_evolparam(Nmodels, inp_df, spec_df, reac_df, data_df, grid_df, modparam_df):
	#
	refspecies2pl = inp_df['refsp'][0]
	for ie in range(15):
	  refspecies2pl = refspecies2pl.replace(str(ie),'$_{'+str(ie)+'}$')
	#
	physplot_df, chemplot_df = axes_var(3,inp_df, data_df, modparam_df, refspecies2pl)
	params, colorsc, linestyle = plot_prop()
	#
	time0 = time.time() ; Nplot = 0
	#
	subprocess.Popen(['mkdir',inp_df.outdir[0]+'/figures/'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#
	# make the plots
	#
	# distribution of all species and all models together
	for iord, ord2 in enumerate(list(chemplot_df.index)): # loop over abundances and related parameters
		if chemplot_df.loc[ord2,'Choice'] >= 1:
			for it, time2 in enumerate(inp_df.timeplot): # loop over fixed times
				indextime = data_df[data_df['time'] > time2].index.min()[1]
				# if only one free param: normal plot 
				if len(grid_df.index) == 1:
					param = grid_df.Param[0]
					#if param == 'nH':
					#modparam2 = modparam_df.iloc[modparam_df.index.get_level_values('Param') == 'nH']
					modparam2 = modparam_df.loc[[(i,param) for i in range(Nmodels)],'Value']
					plt.figure(1,figsize=[6,6])
					plt.ylabel(chemplot_df.loc[ord2,'Title'])
					plt.xlabel(physplot_df.loc[param,'Title'])
					plt.xscale(physplot_df.loc[param,'log'])
					plt.yscale(chemplot_df.loc[ord2,'log'])
					plt.xlim(physplot_df.loc[param,'Min'],physplot_df.loc[param,'Max'])
					plt.ylim(chemplot_df.loc[ord2,'Min'], chemplot_df.loc[ord2,'Max'])
					#indextime = data_df[data_df['time'] > time2].index.min()[1]
					for isp, sp2pl in enumerate(inp_df['species']):
						if ord2+'_'+sp2pl in data_df.columns:
							y2 = (data_df.iloc[data_df.index.get_level_values('Nsteps') == indextime][ord2+'_'+sp2pl])
							plt.plot(modparam2, y2, '-',
									 lw=1.5,color=colorsc[isp],label=inp_df['species2plot'][isp]) # x2,y2,color=colorsc[isp],lw=1.5,ls=linestyle[it],label=inp_df['species2plot'][isp])
					plt.legend(loc='upper left')
					namefile = inp_df.outdir[0]+'/figures/'+inp_df.suffix[0]+'_'+ord2+'_'+"{:.2e}".format(time2)+'_'+param
					plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; 
					Nplot += 1
					plt.close(1)
				# if two free parameters: 2D plot
				elif len(grid_df.index) == 2:
					#print(len(grid_df.Values[0]))
					#exit()
					param1 = grid_df.Param[0] ; param2 = grid_df.Param[1]
					#indextime = data_df[data_df['time'] > time2].index.min()[1]
					for isp, sp2pl in enumerate(inp_df['species']):
						if ord2+'_'+sp2pl in data_df.columns:
							y2 = (data_df.iloc[data_df.index.get_level_values('Nsteps') == indextime][ord2+'_'+sp2pl]).tolist()
							chimaplog = np.reshape(np.log10(y2), (len(grid_df.Values[0]), len(grid_df.Values[1])))
							plt.figure(1,figsize=[6,6])
							plt.set_cmap('rainbow')
							plt.ylabel(physplot_df.loc[param1,'Title'])
							plt.xlabel(physplot_df.loc[param2,'Title'])
							plt.ylim(physplot_df.loc[param1,'Min'],physplot_df.loc[param1,'Max'])
							plt.xlim(physplot_df.loc[param2,'Min'],physplot_df.loc[param2,'Max'])
							chimap = plt.imshow(chimaplog,extent=(physplot_df.loc[param2,'Min'],physplot_df.loc[param2,'Max'],\
                        						physplot_df.loc[param1,'Max'],physplot_df.loc[param1,'Min']), vmin=-7, vmax=-4,\
                        						aspect='auto')#, interpolation='spline36')#,origin='lower', cmap='seismic',interpolation='spline36', aspect='auto')
							plt.colorbar(chimap,label=chemplot_df.loc[ord2,'Title']) 
							namefile = inp_df.outdir[0]+'/figures/'+inp_df.suffix[0]+'_'+ord2+'_'+"{:.2e}".format(time2)+'_'+sp2pl+'_'+param1+'_'+param2
							plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; 
							plt.close(1)
							#CS = plt.contour(np.log10((cdmol13d[0:ncd+1])), temp3d[0:ntemp+1], chimaplog, \
							#                 levels=levels,colors='k',lw=2)

							#plt.plot([np.log10(cdmol13d[0]),np.log10(cdmol13d[ncd])],[10,10],'w--')
							#plt.plot([np.log10(cdmol13d[0]),np.log10(cdmol13d[ncd])],[12,12],'w--')
	time1 = time.time()
	print(Nplot/(time1-time0))
#
##################################################################################
##################################################################################
# plot free param histogram
##################################################################################
##################################################################################
def plot_freeparam(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df):
	#
	refspecies2pl = inp_df['refsp'][0]
	for ie in range(15):
	  refspecies2pl = refspecies2pl.replace(str(ie),'$_{'+str(ie)+'}$')
	#
	physplot_df, chemplot_df = axes_var(1,inp_df, data_df, modparam_df, refspecies2pl)
	params, colorsc, linestyle = plot_prop()
	#
	time0 = time.time() ; Nplot = 0
	#
	subprocess.Popen(['mkdir',inp_df.outdir[0]+'/figures/'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#
	# make the plots
	#
	# distribution of all species and all models together
	for ip, param in enumerate(inp_df['freeparam']):
		if param == 'Ed':
			modparam2 = np.array(modparam_df.iloc[modparam_df.index.get_level_values('Param') == 'EdHs'])
		elif param == 'Eb':
			modparam2 = np.array(spec_df[spec_df['Species'] == 'JH'][['Eb_wat_'+str(im) for im in range(Nmodels)]].T.values.tolist() )
		elif param == 'Ea':
			firstreac = reac_df[(reac_df['Type'] == 14) & (reac_df['C_0'] > 0)].loc[:,'Num'].tolist()[0]
			modparam2 = np.array(reac_df[reac_df['Num'] == firstreac][['C_'+str(im) for im in range(Nmodels)]].T.values.tolist()) 
		#
		# histogram of free parameter
		plt.figure(1,figsize=[6,6])
		plt.xscale('linear')
		plt.yscale('linear')
		plt.xlabel(param)
		plt.ylabel('N$_{models}$')
		plt.xlim(np.array(modparam2).min(),np.array(modparam2).max())
		plt.hist(np.array(modparam2),range=[np.array(modparam2).min(),np.array(modparam2).max()],
				 bins=int(np.sqrt(Nmodels)),lw=1.5, alpha=0.3) #,color=colorsc[0]
		plt.legend(loc='upper left')
		namefile = inp_df.outdir[0]+'/figures/distrib_'+param
		plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; 
		Nplot += 1
		plt.close(1)
		#
		# plot abundances vs free param values as scatter plot
		for iord, ord2 in enumerate(list(chemplot_df.index)): 
		    if chemplot_df.loc[ord2,'Choice'] >= 1:
		      for it, time2 in enumerate(inp_df.timeplot):
		      	plt.figure(1,figsize=[6,6])
		      	plt.ylabel(chemplot_df.loc[ord2,'Title'])
		      	plt.xlabel(param)
		      	plt.xscale('linear')
		      	plt.yscale(chemplot_df.loc[ord2,'log'])
		      	plt.xlim(np.array(modparam2).min(),np.array(modparam2).max())
		      	plt.ylim(1e-7,2e-4) #chemplot_df.loc[ord2,'Min'], chemplot_df.loc[ord2,'Max'])
		      	indextime = data_df[data_df['time'] > time2].index.min()[1]
		      	for isp, sp2pl in enumerate(inp_df['species']):
		        	if ord2+'_'+sp2pl in data_df.columns:
		        		y2 = (data_df.iloc[data_df.index.get_level_values('Nsteps') == indextime][ord2+'_'+sp2pl])
		        		plt.plot(modparam2, y2, '+',
		        				 lw=1.5,color=colorsc[isp],label=inp_df['species2plot'][isp]) # x2,y2,color=colorsc[isp],lw=1.5,ls=linestyle[it],label=inp_df['species2plot'][isp])
		      	plt.legend(loc='upper left')
		      	namefile = inp_df.outdir[0]+'/figures/'+inp_df.suffix[0]+'_'+ord2+'_t='+"{:.2e}".format(time2)+'yr_'+param
		      	plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; 
		      	Nplot += 1
		      	plt.close(1)
		#
		# plot abundances vs free param values with mean and std
		for iord, ord2 in enumerate(list(chemplot_df.index)): 
		    if chemplot_df.loc[ord2,'Choice'] >= 1:
		      for it, time2 in enumerate(inp_df.timeplot):
		      	bins = np.arange(np.float(modparam2.min()),np.float(modparam2.max()),(np.float(modparam2.max())-np.float(modparam2.min()))/10., dtype='float')
		      	plt.figure(1,figsize=[6,6])
		      	plt.ylabel(chemplot_df.loc[ord2,'Title'])
		      	plt.xlabel(param)
		      	plt.xscale('linear')
		      	plt.yscale(chemplot_df.loc[ord2,'log'])
		      	plt.xlim(np.array(modparam2).min(),np.array(modparam2).max())
		      	plt.ylim(1e-7,2e-4) #chemplot_df.loc[ord2,'Min'], chemplot_df.loc[ord2,'Max'])
		      	indextime = data_df[data_df['time'] > time2].index.min()[1]
		      	for isp, sp2pl in enumerate(inp_df['species']):
		        	if ord2+'_'+sp2pl in data_df.columns:
		        		y2 = pd.DataFrame(data_df.iloc[data_df.index.get_level_values('Nsteps') == indextime][ord2+'_'+sp2pl])
		        		y2['modparam'] = modparam2 #list(modparam2.Value)
		        		y2_bins = y2.groupby(np.digitize(y2.modparam, bins))
		        		plt.fill_between(y2_bins.mean()['modparam'],y2_bins.mean()[ord2+'_'+sp2pl]-y2_bins.std()[ord2+'_'+sp2pl],\
	      								 y2_bins.mean()[ord2+'_'+sp2pl]+y2_bins.std()[ord2+'_'+sp2pl],\
	      								 color=colorsc[isp],alpha=0.3)
		        		plt.plot(y2_bins.mean()['modparam'], y2_bins.mean()[ord2+'_'+sp2pl], '-',
		        				 lw=1.5,color=colorsc[isp],label=inp_df['species2plot'][isp]) # x2,y2,color=colorsc[isp],lw=1.5,ls=linestyle[it],label=inp_df['species2plot'][isp])
		      	plt.legend(loc='upper left')
		      	namefile = inp_df.outdir[0]+'/figures/'+inp_df.suffix[0]+'_'+ord2+'_t='+"{:.2e}".format(time2)+'yr_'+param+'_mean'
		      	plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; 
		      	Nplot += 1
		      	plt.close(1)
	#
	# plot abundances vs free param values
	for iord, ord2 in enumerate(list(chemplot_df.index)): 
	    if chemplot_df.loc[ord2,'Choice'] >= 1:
	      for it, time2 in enumerate(inp_df.timeplot):
	      	plt.figure(1,figsize=[6,6])
	      	plt.ylabel(chemplot_df.loc[ord2,'Title'])
	      	plt.xlabel(chemplot_df.loc[ord2,'Title']+'('+inp_df.refsp[0]+')')
	      	plt.xscale(chemplot_df.loc[ord2,'log'])
	      	plt.yscale(chemplot_df.loc[ord2,'log'])
	      	plt.xlim(1e-7,2e-4)#chemplot_df.loc[ord2,'Min'], chemplot_df.loc[ord2,'Max'])
	      	plt.ylim(1e-7,2e-4)#chemplot_df.loc[ord2,'Min'], chemplot_df.loc[ord2,'Max'])
	      	indextime = data_df[data_df['time'] > time2].index.min()[1]
	      	for isp, sp2pl in enumerate(inp_df['species']):
	        	if ord2+'_'+sp2pl in data_df.columns:
	        		x2 = (data_df.iloc[data_df.index.get_level_values('Nsteps') == indextime][ord2+'_'+inp_df.refsp[0]])
	        		y2 = (data_df.iloc[data_df.index.get_level_values('Nsteps') == indextime][ord2+'_'+sp2pl])
	        		plt.plot(x2, y2, '+',
	        				 lw=1.5,color=colorsc[isp],label=inp_df['species2plot'][isp]) # x2,y2,color=colorsc[isp],lw=1.5,ls=linestyle[it],label=inp_df['species2plot'][isp])
	      	plt.legend(loc='upper left')
	      	namefile = inp_df.outdir[0]+'/figures/'+inp_df.suffix[0]+'_'+ord2+'_t='+"{:.2e}".format(time2)+'yr_'+inp_df.refsp[0]
	      	plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; 
	      	Nplot += 1
	      	plt.close(1)
	    #
	    # plot abundances vs free param values
	    for iord, ord2 in enumerate(list(chemplot_df.index)): 
		    if chemplot_df.loc[ord2,'Choice'] >= 1:
		      for it, time2 in enumerate(inp_df.timeplot):
		      	bins = np.arange(np.float(modparam2.min()),np.float(modparam2.max()),(np.float(modparam2.max())-np.float(modparam2.min()))/10., dtype='float')
		      	plt.figure(1,figsize=[6,6])
		      	plt.ylabel(chemplot_df.loc[ord2,'Title'])
		      	plt.xlabel(param)
		      	plt.xscale('linear')
		      	plt.yscale(chemplot_df.loc[ord2,'log'])
		      	plt.xlim(np.array(modparam2).min(),np.array(modparam2).max())
		      	plt.ylim(1e-7,2e-4) #chemplot_df.loc[ord2,'Min'], chemplot_df.loc[ord2,'Max'])
		      	indextime = data_df[data_df['time'] > time2].index.min()[1]
		      	for isp, sp2pl in enumerate(inp_df['species']):
		        	if ord2+'_'+sp2pl in data_df.columns:
		        		y2 = pd.DataFrame(data_df.iloc[data_df.index.get_level_values('Nsteps') == indextime][ord2+'_'+sp2pl])
		        		y2['modparam'] = modparam2 #list(modparam2.Value)
		        		y2_bins = y2.groupby(np.digitize(y2.modparam, bins))
		        		plt.fill_between(y2_bins.mean()['modparam'],y2_bins.mean()[ord2+'_'+sp2pl]-y2_bins.std()[ord2+'_'+sp2pl],\
	      								 y2_bins.mean()[ord2+'_'+sp2pl]+y2_bins.std()[ord2+'_'+sp2pl],\
	      								 color=colorsc[isp],alpha=0.3)
		        		plt.plot(y2_bins.mean()['modparam'], y2_bins.mean()[ord2+'_'+sp2pl], '-',
		        				 lw=1.5,color=colorsc[isp],label=inp_df['species2plot'][isp]) # x2,y2,color=colorsc[isp],lw=1.5,ls=linestyle[it],label=inp_df['species2plot'][isp])
		      	plt.legend(loc='upper left')
		      	namefile = inp_df.outdir[0]+'/figures/'+inp_df.suffix[0]+'_'+ord2+'_t='+"{:.2e}".format(time2)+'yr_'+param+'_mean'
		      	plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; 
		      	Nplot += 1
		      	plt.close(1)
	time1 = time.time()
	print(Nplot/(time1-time0))
#
##################################################################################
##################################################################################
# plot mean and std evolution
##################################################################################
##################################################################################
def plot_mean(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df):
	#
	refspecies2pl = inp_df['refsp'][0]
	for ie in range(15):
	  refspecies2pl = refspecies2pl.replace(str(ie),'$_{'+str(ie)+'}$')
	#
	physplot_df, chemplot_df = axes_var(1,inp_df, data_df, modparam_df, refspecies2pl)
	params, colorsc, linestyle = plot_prop()
	#
	time0 = time.time() ; Nplot = 0
	#
	subprocess.Popen(['mkdir',inp_df.outdir[0]+'/figures/'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#
	# make the plots
	#
	mean_df = (data_df).groupby(level=1).mean().fillna(0)
	std_df = (data_df).groupby(level=1).std().fillna(0)
	#
	for iabs, abs2 in enumerate(physplot_df.index.values):
	  print(abs2+'   ', end="\r")
	  #
	  # abundances
	  for iord, ord2 in enumerate(list(chemplot_df.index)): 
	    if chemplot_df.loc[ord2,'Choice'] >= 1:
	      # 
	      #
	      # compare individual models 
	      plt.figure(1,figsize=[6,6])
	      plt.xlabel(physplot_df.loc[abs2,'Title'])
	      plt.ylabel(chemplot_df.loc[ord2,'Title'])
	      plt.xscale(physplot_df.loc[abs2,'log'])
	      plt.yscale(chemplot_df.loc[ord2,'log'])
	      plt.xlim(physplot_df.loc[abs2,'Min'], physplot_df.loc[abs2,'Max']) #x2.min(),x2.max())
	      plt.ylim(chemplot_df.loc[ord2,'Min'], chemplot_df.loc[ord2,'Max'])
	      for isp, sp2pl in enumerate(inp_df['species']):
	      	if ord2+'_'+sp2pl in data_df.columns:
	      		plt.fill_between(mean_df[abs2],(mean_df[ord2+'_'+sp2pl])-(std_df[ord2+'_'+sp2pl]),\
	      			(mean_df[ord2+'_'+sp2pl])+(std_df[ord2+'_'+sp2pl]),\
	      			color=colorsc[isp],alpha=0.3)
	      		plt.plot(mean_df[abs2],(mean_df[ord2+'_'+sp2pl]),\
	                     color=colorsc[isp],lw=1.5,label=inp_df['species2plot'][isp])
	      plt.legend(loc='upper left')
	      namefile = inp_df.outdir[0]+'/figures/'+inp_df.suffix[0]+'_'+ord2+'_'+abs2+'_mean'
	      plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; 
	      Nplot += 1
	      plt.close(1)
	time1 = time.time()
	print(Nplot/(time1-time0))
#
##################################################################################
##################################################################################
# plot the reaction rates
##################################################################################
##################################################################################
def plot_rates(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df, sp2rate, list_ratform, list_ratdest, list_numform, list_numdest):
	list_reacs = ['React1','React2','Prod1','Prod2']
	#
	refspecies2pl = inp_df['refsp'][0]
	for ie in range(15):
	  refspecies2pl = refspecies2pl.replace(str(ie),'$_{'+str(ie)+'}$')
	#
	physplot_df, chemplot_df = axes_var(1,inp_df, data_df, modparam_df, refspecies2pl)
	params, colorsc, linestyle = plot_prop()
	#
	time0 = time.time() ; Nplot = 0
	# create figures folder
	for im in range(Nmodels):
	  #os.system('mkdir '+inp_df.outdir[0]+'/'+inp_df.outmod[im]+'/figures/')
	  subprocess.Popen(['mkdir',inp_df.outdir[0]+'/'+inp_df.outmod[im]+'/figures_rates/'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#
	for iabs, abs2 in enumerate(physplot_df.index.values): 
	  #
	  # abundances
	  for isp, sp in enumerate(sp2rate):
	  	print(abs2+' '+sp, end="\r")
	  	for im in range(Nmodels):
	  		#
	  		plt.figure(1,figsize=[6,6])
	  		plt.xlabel(physplot_df.loc[abs2,'Title'])
	  		plt.ylabel('Rate [cm$^{-3}$ s$^{-1}$]')
	  		plt.xscale(physplot_df.loc[abs2,'log'])
	  		plt.yscale('log')
	  		plt.xlim(data_df.loc[im,abs2].min(),data_df.loc[im,abs2].max())
	  		plt.ylim(1e-30,1e-5)
	  		for ir, rate2pl in enumerate(list_ratform[isp]):
	  			ireac = int(list_numform[isp].loc[im][ir])
	  			if ireac > 0:
		  			list_reacprod = reac_df.loc[ireac,list_reacs].tolist()
		  			for imol,mol in enumerate(list_reacprod):
		  				list_reacprod[imol] = sp2plot(list_reacprod[imol])
		  			#str_reacprod = ('%s+%s $\rightarrow$ %s+%s' % (list_reacprod[0], list_reacprod[1], list_reacprod[2], list_reacprod[3] ))
		  			str_reacprod = (list_reacprod[0]+'+'+list_reacprod[1]+'->'+list_reacprod[2]+'+'+list_reacprod[3])
		  			plt.plot(data_df.loc[im,abs2],list_ratform[isp].loc[im][ir],lw=1.5,label=str_reacprod)
		  		plt.legend(loc=[1,0], fontsize=7)#'upper left')
	  		namefile = inp_df.outdir[0]+'/'+inp_df.outmod[im]+'/figures_rates/'+sp+'_'+abs2+'_form'
	  		plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; Nplot += 1
	  		plt.close(1)
	  		#
	  		plt.figure(1,figsize=[6,6])
	  		plt.xlabel(physplot_df.loc[abs2,'Title'])
	  		plt.ylabel('Rate [cm$^{-3}$ s$^{-1}$]')
	  		plt.xscale(physplot_df.loc[abs2,'log'])
	  		plt.yscale('log')
	  		plt.xlim(data_df.loc[im,abs2].min(),data_df.loc[im,abs2].max())
	  		plt.ylim(1e-30,1e-5)
	  		for ir, rate2pl in enumerate(list_ratdest[isp]):
	  			ireac = int(list_numdest[isp].loc[im][ir])
	  			if ireac > 0:
		  			list_reacprod = reac_df.loc[ireac,list_reacs].tolist()
		  			for imol,mol in enumerate(list_reacprod):
		  				list_reacprod[imol] = sp2plot(list_reacprod[imol])
		  			#str_reacprod = ('%s+%s $\rightarrow$ %s+%s' % (list_reacprod[0], list_reacprod[1], list_reacprod[2], list_reacprod[3] ))
		  			str_reacprod = (list_reacprod[0]+'+'+list_reacprod[1]+'->'+list_reacprod[2]+'+'+list_reacprod[3] )
		  			plt.plot(data_df.loc[im,abs2],list_ratdest[isp].loc[im][ir],lw=1.5,label=str_reacprod)
		  		plt.legend(loc=[1,0], fontsize=7)#'upper left')
	  		namefile = inp_df.outdir[0]+'/'+inp_df.outmod[im]+'/figures_rates/'+sp+'_'+abs2+'_dest'
	  		plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True) ; Nplot += 1
	  		plt.close(1)