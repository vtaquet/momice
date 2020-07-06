##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
## perform analysis of the data
##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
#
from scipy.interpolate import interp1d
from numpy import array, arange, sin
import numpy as np
import math as math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import multiprocessing, subprocess
import os
import struct
from scipy.io import FortranFile
#
##################################################################################
##################################################################################
# compute abundance / reference species, deuteration, H2 opr, 
##################################################################################
##################################################################################
#
def compute_params(Nmodels, Nspgas, Nspice, inp_df, spec_df, data_df):
	#
	spgas = [] ; spice = []
	for im in range(Nmodels):
		spgas.append(list(spec_df['Species'])[0:Nspgas[im]])
		spice.append(list(spec_df['Species'])[Nspgas[im]:Nspgas[im]+Nspice[im]])
		spice[im] = [spice2[1:] for spice2 in spice[im]]
    #
	#print(data_df.head())
	# abundance / reference species
	#chemfile = chemfile+['Xgasref','Xiceref']
	Xgasref_list = ['Xgasref_'+spgas2 for spgas2 in spgas[0]]
	Xgas_list = ['Xgas_'+spgas2 for spgas2 in spgas[0]]
	for isp, spgas2 in enumerate(spgas[0]):
		data_df['Xgasref_'+spgas2] = data_df['Xgas_'+spgas2]/data_df['Xgas_'+inp_df['refsp'][0]]
	for isp, spice2 in enumerate(spice[0]):
	  data_df['Xiceref_'+spice2] = data_df['Xice_'+spice2]/data_df['Xice_'+inp_df['refsp'][0]]
	#
	# deuteration
	choice_D = spec_df['Species'].str.contains("D").sum() #== 'D'
	if choice_D >= 1:
	  for isp, spgas2 in enumerate(spgas[0]):
	    if spgas2.find('D') >= 0:
	      for isp2 in range(isp):
	        if spgas[0][isp-isp2].find('D') < 0:
	          data_df['Dgas_'+spgas[0][isp]] = data_df['Xgas_'+spgas2]/data_df['Xgas_'+spgas[0][isp-isp2]]
	          break
	  for isp, spice2 in enumerate(spice[0]):
	    if spice2.find('D') >= 0:
	      for isp2 in range(isp):
	        if spice[0][isp-isp2].find('D') < 0:
	          data_df['Dice_'+spice[0][isp]] = data_df['Xice_'+spice2]/data_df['Xice_'+spice[0][isp-isp2]]
	          break
	#
	# H2 opr
	choice_opr = spec_df['Species'].str.contains("oH2").sum()
	if choice_opr >= 1:
	  data_df['H2opr'] = data_df['Xgas_oH2']/data_df['Xgas_pH2']
	#
	return choice_D, choice_opr, data_df


##################################################################################
# perform PCA analysis
##################################################################################
#
def PCA(Nmodels, inp_df, spec_df, reac_df, data_df, modparam_df):
	#
	print('performing pca analysis...')
	#
	# import packages
	from sklearn.preprocessing import StandardScaler
	from sklearn.pipeline import make_pipeline
	from sklearn.decomposition import PCA
	#
	list_chem = ['Xice']
	namecol_sel = [ord2+'_'+sp2pl for iord, ord2 in enumerate(list_chem) for isp, sp2pl in enumerate(inp_df['species'])] 
	#print(namecol_sel)
	#exit()
	#
	# scale the data to get same mean and standard deviation for each feature
	Nsteps = modparam_df.loc[(0,'Nstep'),'Value']
	data_pca = data_df.loc[[(i,Nsteps-1) for i in range(inp_df.Nmodels)],namecol_sel].fillna(0) #df_obsdata2[df_obsdata2 == np.nan] = 0
	#print(data_pca[namecol_sel])
	scaler = StandardScaler()
	data_scaled = scaler.fit_transform(data_pca[namecol_sel]) #df_obsdata2[namecol2[0*Nspecies:1*Nspecies]] #s
	#print(data_scaled)
	#exit()
	#
	#
	pca = PCA(n_components=len(namecol_sel))
	pca.fit(data_scaled)
	features = range(pca.n_components_)
	variance = pca.explained_variance_ratio_
	components = pca.components_
	loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
	singular_values = pca.singular_values_
	transformed = pca.transform(data_scaled)
	#
	for ipca in range(1, pca.n_components_):
		data_pca['PC'+str(ipca)] = [float(dat) for dat in np.nditer(transformed[:,ipca-1])]  #list_val 
		#df_obsdata2 = df_obsdata[df_obsdata[nameTmb] > 3*row['rms']]
		#df_obsdata2 = df_obsdata[df_obsdata['H2cd'] > 3*1.9e21]
		#corr = df_obsdata2['PC'+str(ipca)].corr(df_obsdata2['H2cd'],method='spearman')
		#corr_list2.append(corr)
		#print(ipca,corr)
	#print(data_pca)

	subprocess.Popen(['mkdir',inp_df.outdir[0]+'/figures/'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	# plot the variance of each PC
	fig = plt.figure(0,figsize=[5,5])
	plt.bar(features, variance, edgecolor='k')
	plt.xticks(features)
	plt.xlabel('PC feature')
	plt.ylabel('Variance')
	namefile = inp_df.outdir[0]+'/figures/pca_variance'
	plt.savefig(namefile+'.pdf',transparent=True,bbox_inches='tight')
	plt.close(0)

	# plot the contribution of each species to PCs 
	for ipca in range(1, pca.n_components_):
	  fig = plt.figure(ipca,figsize=[5,2])
	  plt.bar(range(len(inp_df['species'])), pca.components_[ipca-1,:], edgecolor='k')#, alpha=0.3)
	  #plt.xticks(inpparam_df['species2plot'])
	  plt.xticks(range(len(inp_df['species'])), inp_df['species'], rotation=90)
	  plt.title('PC '+str(ipca))
	  plt.ylim(-1,1)
	  plt.xlabel(' ')
	  plt.ylabel(' ')
	  namefile = inp_df.outdir[0]+'/figures/pca_contribution_'+str(ipca)
	  plt.savefig(namefile+'.pdf',transparent=True,bbox_inches='tight')
	  #plt.savefig('../../figures/pca/contribution_'+str(ipca)+'.pdf',bbox_inches='tight',transparent=True)
	  plt.close(ipca)
	  #exit()

	  fig = plt.figure(ipca,figsize=[4,4])
	  for isp, sp2pl in enumerate(inp_df['species']):
	    plt.plot([0,pca.components_[ipca-1,isp]], [0,pca.components_[ipca,isp]], '-k')
	    plt.text(pca.components_[ipca-1,isp],pca.components_[ipca,isp],sp2pl)
	  #plt.xticks(features)
	  plt.xlabel('PC '+str(ipca))
	  plt.ylabel('PC '+str(ipca+1))
	  plt.xlim([-1,1])
	  plt.ylim([-1,1])
	  #plt.show()
	  namefile = inp_df.outdir[0]+'/figures/pca_correlationwheel_pc'+str(ipca)
	  plt.savefig(namefile+'.pdf',bbox_inches='tight',transparent=True)
	  #plt.savefig('../../figures/pca/correlationwheel_pc'+'_'+str(ipca)+'.pdf',bbox_inches='tight',transparent=True)
	  plt.close(ipca)

	#exit()

