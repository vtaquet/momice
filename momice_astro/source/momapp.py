###############################################################################################
# momapp.py: MOMice APPlier
# runs MoMICE for multiple simulations (spatial evolution, model grid, or statistical analysis)
###############################################################################################
#
import numpy as np
import pandas as pd
import os
import multiprocessing, subprocess
import time
#
print('MOMAPP: MOMICE APPLIER')
#
# check where is the comment starting with #
def checkcomment(inplist):
  #typelist = [type(value) for value in inplist]
  for i, value in enumerate(inplist):
    if value == '#' or value == '!':
      return i
#
# function that checks number of species and reactions of input files
def check_network():
	inpascii = open('../input/'+inp_df['inpdir'][0]+'/input_parameters.in','r')
	line2 = []
	for il, line in enumerate(inpascii):
	    s = line.split()
	    if il == 10:
	    	filesp = s[0]
	    elif il == 11:
	    	filere = s[0]
	inpsp = open('../input/'+inp_df['inpdir'][0]+'/'+filesp,'r')
	for il, line in enumerate(inpsp):		
		s = line.split()
		if s[0] == 'Nspecies':
			Nspecies = int(s[-1])
	inpre = open('../input/'+inp_df['inpdir'][0]+'/'+filere,'r')
	for il, line in enumerate(inpre):		
		s = line.split()
		if s[0] == 'Nreactions':
			Nreacs = int(s[-1])
		elif s[0] == '*':
			break
	return Nspecies, Nreacs
#
#
# worker function that runs MOMICE
def mp_worker(args):
	#
    global Ndone, Nmodels
    Ndone += 1 
    irun = args[0]
    #
    print("Process %s run..." % (irun+1))
    #
    cmd = ['./momice.gfort', str(inp_df.typesimu[0]), inp_df.outdir[0]] + args[1:] #inp_df.outdir[0], namemod, namefile]
    print(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (output, err) = p.communicate()  
    p_status = p.wait()
    output_list = str(output).split("\\n") ; err_list = str(err).split("\\n")
    #
    try:
    	outlogterm = open(inp_df.outdir[0]+'/'+args[1]+'/'+'logterm.out','w')
    	for line in output_list:
    		outlogterm.write(line+' \n')
    	for line in err_list:
    		outlogterm.write(line+' \n')
    	outlogterm.close()
    except:
    	pass
    print("Process %s done!" % (irun+1))
    print(str(Ndone)+'/'+str(Nmodels)+' done!', end="\r")
#
# handler function for multiprocessing computation
def mp_handler(Ncpus, data):
    p = multiprocessing.Pool(Ncpus)
    p.map(mp_worker, data)
#
# main program
if __name__=='__main__':
	#
	global Ndone, Nmodels
	Ndone = 0
	#
	# read grapp input file
	inpascii = open("momapp.in",'r')
	inpparam = [] ; inpvalue = []
	for line in inpascii:
	    s = line.split()
	    wherecom = [spl.find('#') for spl in s]
	    if s[0] != '#' and s[0] != '*':
	      # convert params to float if possible
	      try:
	        s[2:checkcomment(s)] = [int(i) for i in s[2:checkcomment(s)]]
	      except:
	        pass
	      inpparam.append(s[0])
	      inpvalue.append(s[2:checkcomment(s)])
	    elif s[0] == '*':
	      break
	#
	# create input params dataframe
	inp_df = pd.DataFrame({inpparam[i] : [inpvalue[:][i]] for i in range(len(inpparam))},columns=inpparam).loc[0]  
	subprocess.Popen(['mkdir',inp_df.outdir[0]], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#print(inp_df.head())
	#
	# read number of species and reactions
	Nspecies, Nreacs = check_network()
	#
	# modify name input directory
	with open('../input/input_file.in', "r+") as f:
		data = f.read().split('\n')
		data[1] = inp_df.inpdir[0]
		f.seek(0)
		for item in data:
			f.write("%s\n" % item)
		f.truncate()
	#
	# spatial evolution
	if inp_df.typesimu[0] == 2:
		#
		## read spatial log file
		logfile = []
		inplogspat = open('../input/'+inp_df.inpdir[0]+'/'+inp_df.nameinp2[0], 'r')
		for il, line in enumerate(inplogspat):
			s = line.split()
			if s[0] != '#' and s[0] != '*':
				if il == 0:
					spatdir = s[0]
				elif il >= 1:
					logfile.append([il, s[0], spatdir+'/'+s[0]+'.dat']) 
			elif s[0] == '*':
				break
		Nmodels = len(logfile)
		#
		# write spatial log file in output directory
		subprocess.Popen(['mkdir',inp_df.outdir[0]], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		outlogspat = open(inp_df.outdir[0]+'/spatlog.out','w')
		outlogspat.write(str(Nmodels)+'\n')
		for line2 in logfile:
			outlogspat.write(line2[1]+'\n')
	#
	# model grid
	elif inp_df.typesimu[0] == 3:
		#
		## read model grid file
		gridfile = [] ; gridparam = [] ; gridvalue = []
		inploggrid = open('../input/'+inp_df.inpdir[0]+'/'+inp_df.nameinp3[0], 'r')
		for il, line in enumerate(inploggrid):
			s = line.split()
			s = [s[0]] + [float(s2) for s2 in s[1:]]
			if s[0] != '#' and s[0] != '*':
				gridparam.append(s[0])
				gridvalue.append(s[:])
			elif s[0] == '*':
				break
		#
		# create input params dataframe
		grid_df = pd.DataFrame(gridvalue, columns=['param', 'Nsteps', 'min', 'max'])#.set_index('Param')
		Nparams = len(grid_df)
		#Nmodels = len(spatfile)
		# 
		# make grid of parameters
		paramvalue = []
		for irow, row in grid_df.iterrows():
			if row['param'] in ['nH', 'CR', 'UV']:
				paramvalue2 = np.logspace(np.log10(row['min']), np.log10(row['max']), int(row['Nsteps']))
			else:
				paramvalue2 = np.linspace(row['min'], row['max'], int(row['Nsteps']))
			paramvalue.append(paramvalue2)
		#print(len(paramvalue))
		if len(paramvalue) == 1:
			list_grid = [[p1] for p1 in paramvalue[0]]
		elif len(paramvalue) == 2:
			list_grid = [[p1,p2] for p1 in paramvalue[0] for p2 in paramvalue[1]]
		elif len(paramvalue) == 3:
			list_grid = [[p1,p2,p3] for p1 in paramvalue[0] for p2 in paramvalue[1] for p3 in paramvalue[2]]
		elif len(paramvalue) == 4:
			list_grid = [[p1,p2,p3,p4] for p1 in paramvalue[0] for p2 in paramvalue[1] for p3 in paramvalue[2] \
						 for p4 in paramvalue[3]]
		elif len(paramvalue) == 5:
			list_grid = [[p1,p2,p3,p4,p5] for p1 in paramvalue[0] for p2 in paramvalue[1] for p3 in paramvalue[2] \
						 for p4 in paramvalue[3] for p5 in paramvalue[4]]
		elif len(paramvalue) == 6:
			list_grid = [[p1,p2,p3,p4,p5,p6] for p1 in paramvalue[0] for p2 in paramvalue[1] for p3 in paramvalue[2] \
						 for p4 in paramvalue[3] for p5 in paramvalue[4] for p6 in paramvalue[5]]
		elif len(paramvalue) == 7:
			list_grid = [[p1,p2,p3,p4,p5,p6,p7] for p1 in paramvalue[0] for p2 in paramvalue[1] for p3 in paramvalue[2] \
						 for p4 in paramvalue[3] for p5 in paramvalue[4] for p6 in paramvalue[5] for p7 in paramvalue[6]]
		elif len(paramvalue) == 8:
			list_grid = [[p1,p2,p3,p4,p5,p6,p7,p8] for p1 in paramvalue[0] for p2 in paramvalue[1] for p3 in paramvalue[2] \
						 for p4 in paramvalue[3] for p5 in paramvalue[4] for p6 in paramvalue[5] for p7 in paramvalue[6] \
						 for p8 in paramvalue[7]]
		Nmodels = len(list_grid)
		#
		# write grid log file in output directory
		outloggrid = open(inp_df.outdir[0]+'/gridlog.out','w')
		outloggrid.write(" ".join([str(Nparams), str(Nmodels)])+'\n') ; logfile = []
		for irow, row in grid_df.iterrows():
			outloggrid.write(" ".join([row['param']]+[str(param) for param in paramvalue[irow]]+['\n']))
		for il, line in enumerate(list_grid):
			namemod = "_".join(["{:.2e}".format(line2) for line2 in line])[:]
			logfile.append([il, namemod] + ["{:.2e}".format(line2) for line2 in line])
			outloggrid.write(" ".join([str(il), namemod, '\n']))
	#
	# sensitivity analysis
	elif inp_df.typesimu[0] == 4:
		#
		## read model grid file
		sensfile = [] ; sensparam = [] ; sensvalue = [] ; Nmodels = 0
		inplogsens = open('../input/'+inp_df.inpdir[0]+'/'+inp_df.nameinp4[0], 'r')
		for il, line in enumerate(inplogsens):
			s = line.split()
			if s[0] != '#' and s[0] != '*':
				if il == 0:
					Nmodels = int(s[0])
				elif il == 1:
					distrib = s[0]
				else:
					s = [s[0], s[1]] + [float(s2) for s2 in s[2:]]
					sensparam.append(s[0])
					sensvalue.append(s[:])
			elif s[0] == '*':
				break
		#
		# create input params dataframe
		sens_df = pd.DataFrame(sensvalue, columns=['param', 'selection', 'mean', 'std'])#.set_index('Param')
		Nparams = len(sens_df) 
		#
		# create distribution of input parameters
		Nvalues = {'Ed' : 1 , 'Eb' : Nspecies , 'Ea' : Nreacs}
		for irow, row in sens_df.iterrows():
			outinpparam = open('../input/'+inp_df.inpdir[0]+'/'+row['param']+'.in','w')
			distribvalues = []
			for iv in range(Nvalues[row['param']]):
				if distrib == 'normal':
					distribvalues.append(list(np.random.normal(row['mean'], row['std'],Nmodels)))
					#print(len(list(np.random.normal(row['mean'], row['std'],Nmodels))))
				elif distrib == 'uniform':
					distribvalues.append(list(np.random.uniform(row['mean']-row['std'],row['mean']+row['std'],Nmodels)))
			distribvalues_rev = np.array(distribvalues).T.tolist()
			for im in range(Nmodels):
				outinpparam.write(" ".join(["{:.5e}".format(distribvalues2) for distribvalues2 in distribvalues_rev[im]]+[' \n']))
		#		
		# write sensitivity log file in output directory
		outlogsens = open(inp_df.outdir[0]+'/sensitivitylog.out','w')
		outlogsens.write(" ".join([str(Nparams), str(Nmodels)])+'\n') ; logfile = [] ; randinp = [] ; logfile2 = []
		for irow, row in sens_df.iterrows():
			outlogsens.write(" ".join([row['param'],row['selection'],str(row['mean']),str(row['std']),'\n']))
			for row2 in [row['param'],row['selection']]:
				logfile2.append(row2) 
		for im in range(Nmodels):
			logfile.append([(im),str(im),distrib]+logfile2)
			outlogsens.write(str(im)+'\n')
	#
	## run MOMICE on several CPUs
	#print(logfile)
	mp_handler(inp_df.Ncpus[0], logfile)
	#
	exit()
