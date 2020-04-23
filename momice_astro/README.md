This README describes how to compile and execute the code, the input files needed by MOMICE_ASTRO, and the routine that analyzes and visualises the results.

## 1. Compilation and Execution

The code is written in Fortran 90, it is therefore necessary to compile its source files located in the `source` directory before running the model. 

To compile the source code: 
- go to the `source` directory. 
- type "make gfort" or "make ifort" depending on the compiler (gfortran or ifortran) you wish to use. 
- a `momice.gfort` or `momice.ifort` executable file is created.

The fortran executable is called by the MOMAPP python routine (`momapp.py`) also located in the `source` directory. MOMAPP reads a first input file `momapp.in` in which the user defines the type of simulation to be run and the location of input files read by MOMICE. See the next section for a description of `momapp.in` and all the input files read by MOMICE.

Once the input files are configured, the simulations can be run by executing the `momapp.py` python routine. 


## 2. Input files

### `momapp.in`

MOMAPP first reads `momapp.in` to define the type of simulation. A few parameters need to be chosen:

- type of simulation: the user can choose to run 1) individual simulations with constant physical conditions, 2) simulations for evolving physical conditions as function of distance of a central source, 3) model grids in which (constant) physical conditions are explored, 4) sensitivity analyses in which the distributions of key input parameters impacts the uncertainties of the model predictions. 

- Input directory: name of the directory within the `input` folder that contains all the input files, described below.

- Output directory path with respect to the source folder.

- For options 2 to 4, one needs to set up a file located in the input directory and that contains the parameters of the simulations (depending on the type of simulation specified above). 

- Number of CPUs to be used for running the simulations (-1 means that all available CPUs are used).


### `input_parameters.in`

It is the main and first input file to check. It allows you to specify the physical conditions, the grain and ice properties, the names of the chemical networks, and to choose whether you wish to run an individual model or a model grid. 

#### DATA OPTIONS:

- Saving results: Whether or not you want to save the results. 
- Formation/destruction rates: Whether or not you want to save the formation and destruction rates of selected species. If so, please define the list of selected species in `reacrates.in`. 

#### INPUT/OUTPUT FILES:

- Location of the output directory with respect to the source folder. You already defined the folder in `momapp.in`. 
- Name of directory: By default, the name of the directory is the date and the time at which the model has been run (ex: 20130211_0956) if "date" is chosen. If you want to specify another name, please edit this line.
- The file specifying the list of species and their initial abundances (see below for more details).
- The file specifying the list of reactions and their properties (see below for more details).
- The file specifying the list of binding energies of surface species (see below for more details).
- The file specifying the list of files that gives the temporal evolution of physical conditions as function of the distance from the central source (see below for more details).

#### PHYSICAL CONDITIONS:

- Constant physical conditions. Whether or not you want to apply evolving physical conditions. If physical conditions are constant (ch_ph=1), then the following physical conditions are specified here: 
	- total time of integrations (in years)
	- total density of H nuclei n<sub>H</sub> (cm<sup>-3</sup>)
	- gas temperature (K) and grain/ice temperature (K)
	- cosmic ray ionization rate zeta (s-1)
	- visual extinction Av (mag)
	- scaling factor for external radiation field G0
	- "secondary" UV flux generated from the H<sub>2</sub> ionization by cosmic rays (cm<sup>-2</sup> s<sup>-1</sup>)


#### GRAIN AND ICE PROPERTIES

This section allows you to fix the main properties of interstellar grains and ices that affect the gas-grain process and surface chemistry.

- Grain properties: 
	- dust-to-gas mass ratio R<sub>dg</sub>
	- grain size a<sub>d</sub>
	- volumic mass of grains \rho<sub>d</sub>
From these parameters, MOMICE determines the grain abundance and the grain cross section. In this model, only one constant grain size can be specify.

- Ice properties: 
	- site size d<sub>s</sub>
	- sticking coefficient for species heavier than H and H<sub>2</sub>
	- diffusion-to-binding energy ratios by distinguishing the surface and the bulk of the ice, and "light" (H, H<sub>2</sub>) and "heavy" species
	- number of active surface monolayers and approximate number of timesteps needed to fill one monolayer
The values of these ice parameters are extensively discussed in my thesis and in the papers. Please refer to them for the best choice.

- Porosity parameters. A 3D porosity treatment is yet to be included. has been recently included. Please disregard the porosity parameters at the moment. %From these three input parameters (size of each square pore, fraction of the sphere occupied by pore, vacuum in the grain), MOMICE shall compute the area of the grain surface and the exchange rate between the pores and the non-porous surface. If the grain is porous, you need to add the porous species and the exchange processes in your chemical network (see next section). 

#### MODEL SWITCHES

This section allows you to switch between different chemistry options:
- Ice formation: old 2-phase "bulk" or "new" 3-phase multi-phase approach. In this version, the 3-phase approach follows the method introduced by Hasegawa & Herbst (1993) that allows a faster computation with respect to the approach by Taquet et al. (2012). 
- Modified rates: One can include the modified rate approach proposed by Garrod (2008). At this moment, it's still in beta, please do not use it yet. 
- Grain size evolution due to the growth of ice, assuming that the thickness of each monolayer is equal to the size of a site. 
- Binding energy evolution at the surface and in the bulk as a function of the coverage of bare grain, H2, and water. 
- Atomic hydrogen in gas phase: One can choose to keep nH(H) or X(H) constants if, for example, you're using a small chemical network which does not allow you to self consistently compute the abundance of H. 
- H$_22$ ortho/para ratio: With the last version of the chemical network, the H2 ortho/para ratio evolves with time due to gas phase reactions. But if you use a simple chemical network, you can also choose to keep it constant in order to study its influence on the chemistry. 
- Reaction probability: Choose "1" if the reaction rate is directly computed from the transmission probability of the reaction (Hasegawa et al. 1992 approach for rate equations). Choose "2" if the reaction rate is computed from the competition process between the transmission probability of the reaction and the diffusion rate of the reactants (Chang et al. 2007 approach for microscopic Monte-Carlo calculations).  
- Transmission probability: Choose "1" to compute the transmission probability by assuming a rectangular barrier for all reactions (in that case, the activation energies need to be specified in the reaction network and the width of the barrier is specified just below). Choose "2" to compute the transmission probability by using the Eckart model which fits an approximate Potential Energy Surface for the reactions available in the data file "eckart\_data.in". 
- Methanol network: Choose options "2" or "3" to use the rates deduced by Watanabe et al. (2008) from their experiments or by Caselli et al. (2002) from theoretical calculations. Read my thesis for more details. 

#### NUMERICAL PARAMETERS



### 1.3 species_file.in

 This file contains the species network. Each line refers to a specific species and includes: i) its name (15 caracters), its charge, its number of elements, and its initial abundance relative to H nuclei. If you add a species, please do not forget to specify the right number of elements and the charge. GRAINOBLE uses this data to check the conservation of mass and charge. Please do not forget to specify the number of species at the beginning of the file.
The species network must respect several rules to be read correctly by the fortran code:
i) the name of non-porous surface species start with a J and the name of porous surface or bulk species start with a Q. At this moment, one cannot combine the 3-phase and the porosity treatment simultaneously. 
ii) gas phase species are listed first, then the external surface species, and then the porous or bulk species. 
iii) The list of gaseous neutral species, surface species, and porous/bulk species need to be identical. 
iv) Gas phase species that are not considered in surface (e.g. ions) need to be included after the neutral gas phase species.  
Example of a species network included 4 neutral species and 2 ions: 
H 
D 
O 
CO 
H+ 
D+ 
JH 
JD  
JO 
JCO 
QH 
QD 
QO 
QCO
 
 
### 1.4 reaction_file.in

 This file contains the reaction network. Each line refers to a specific reaction and includes i) up to 3 reactants (15 caracters each), ii) up to 4 products (15 caracters each); iii) 3 parameters A, B, C; iv) the reaction type, v) the range of temperatures where the reaction is considered, vi) the type of the formula. Please don't forget to specify the number of reactions at the beginning of the file.

For gas phase reactions (from 1 to 8), the reaction rates are computed following the methodology adopted in KIDA. A description of the types and formulas can be found \href{hhttp://kida.obs.u-bordeaux1.fr/help}{here}. In addition to the standard gas-phase reactions, GRAINOBLE takes the following processes into account: \\
0: Gas-grain interaction and electron-grain recombination ($Rate(s^{-1}) = A \times \zeta$) (from Flower and Pineau des Forets 2003). 
20: Accretion from gas phase to grain surfaces. A=pre-factor constant. 
14: Langmuir-Hinshelwood grain-surface chemical reactions. 
21: Thermal evaporation. A=pre-factor constant. 
22: Cosmic Ray induced general desorption following Hasegawa \& Herbst 93. A=branching ratio. 
23: Photo-desorption by induced-CR photons + background photons coming from wavelength-dependent experimental studies (see Taquet et al. 2013). A=branching ratio. 
26: Photodissociation/photodesorption processes coming from molecular dynamics simulations (see Taquet et al. 2012c). A=branching ratio. 
30: Non-porous surace->Pores rate exchange. A=branching ratio. 
31: Pores->Non-porous rate following the method introduced by Taquet et al. (2012). A=branching ratio. 
32: Surface<->mantle exchange rate following the method introduced by Hasegawa \& Herbst (1993).


### 1.5 energy_file.in

 This file contains the list of binding energies of surface species relative to a bare grain surface (amorphous carbon and/or silicate), amorphous water ice, H2 ice, and as a pure substrate. For most species, it also contains the references where the values are taken from. This file should be considered as a data file since it contains all possible species which are not necessarily included in the species network. However, every surface species included in the species network must be located somewhere in this file. 


### 1.6 model_grid.in

 If a model grid is desired in input_parameters.in, one needs to specify this file. 

- Initial abundances: Instead of using the same initial abundances for all the model grid given in the species network, the results of another model grid can be used as initial abundances for the new model grid. To use the results of another model grid, choose option "2" and specify the location and the name of input grid (relative to the root directory). Of course, the number of free parameters and the values between the input model grid and the new one need to identical.

- Input parameters: 14 input parameters can be chosen as free parameters. For each parameter, one needs to specify the number of parameter steps as well as its extremal values. GRAINOBLE will compute all the values following linear (or logarithmic) steps. For example, if you choose 4 steps for $n_H$ between $10^3$ and $10^6$ cm$^{-3}$, GRAINOBLE will specify $n_H=10^3$, $10^4$, $10^5$, $10^6$ cm$^{-3}$. If you don't want to use a specific parameter as free, choose "0" for its number of steps and specify its value in input_parameters.in.


### 1.7 reacrates.in

This file contains the list of species whose formation and destruction rates will be computed and saved as binary files. For each species, the rates of the 30 most important reactions of formation and destruction are saved. The amount of data tends to quickly increase with the number of species, one should limit the number of species to 50.


### 1.8 spat_evol.in

If a spatial evolution has been chosen in input\_parameters.in, the folder and the list of files giving the one-dimensional evolution of the physical conditions have to be specified in this file. The format of the input files follows the format given in the examples located in the "phys" folders. 



## 3. OUTPUT ANALYSIS AND VISUALISATION

The output files are saved in binary files located in the selected output directory. An IDL code, called GRAN (for GRainoble ANalyzer) is provided in the output directory. It reads and analyzes the binary files, and plots the results in postcript files. To run GRAN, go to the output directory, run IDL and type "GRAN".

### 3.1 Individual model

The input_mod.in file, located in the output directory, allows you to analyze and plot the results of individual models. 
First, you need to specify the location and the name of the models you want to read. When several models are read in one call, their results are also compared in the same plots located in the figures folder within the directory where the models are located. You can specify a prefix to the postcript files which show the comparison between the different models. 
Then, you need to specify different options for plotting the results (plotting the formation and destruction rates, log or linear axes, plotting observational data, or the name of the reference species). Finally, the list of species that will be showed has to be given in one row. each species has to be followed by a blank character. 

GRAN creates the postcript files in the model directory (yes). They show the evolution of various variables (physical conditions, gas and ice abundances, gas and ice deuterium fractionation, CO depletion, fractional composition within the ices, ice thickness, H2 o/p ratio, formation rates) as function of time, CO depletion, and ice thickness, and physical parameters. If chosen, GRAN also presents the contribution of the major reactions to the formation and the destruction of specified species.


### 3.2 Model grid

The input_grid.in file allows you to analyse a model grid. First, you need to specify the location of the grid directory and its name. For one grid, several gridlog.out files can be edited, depending or not if you want to analyze the entire grid of just a part of it. If you want to read a modified gridlog file, please edit the suffix of the gridlog''...''.out file. 
Before creating the postcript files, GRAN saves the analyzed data. If you want to restore some data already saved in a previous run of GRAN, type "yes" to "Restore the results data already saved by GRAN?". For more details about the analysis of the model grid, please take a look at my thesis. GRAN can plot the results of the model grid via distributions of abundances/fractionations at specified times, or time-dependent evolutions of abundances/fractionations. You can decide what you want to show by typing "yes" or "no" to the two concerned questions. Last, you can specify the species that are compared in the same plot. 

GRAN creates the postcript files in the figures directory located in model grid directory. The directory evoltime contains the plots showing the abundances as function of time while distribs contains the distribution plots. GRAN compares the abundances and distributions for each value of every free parameter by varying the values of other parameters. For instance, a model grid is built by choosing four values for the density, the temperature, and the grain size. The file Dgas_time_CH2DOH_Tg=Td_sd.ps shows the mean and the standard deviation of the deuterium fractionation of CH2DOH as function of time for the different values of the temperature and induced by the variation of the values of other parameters.


### 3.3 Spatial evolution

The input_spat.in file allows you to analyse a set of models in order to follow the one-dimensional evolution of the chemistry. As usual, you need to specify the location of the set of output folders where the "spatlog.out" file, generated by GRAINOBLE, is located. In addition to the usual different options and the list of species to plot, a list of time steps has to be specified.

GRAN generates a set of postcript and ascii files in different folders showing the evolution of the chemical abundances as function of the radius for the different times chosen in input_spat.in. 
