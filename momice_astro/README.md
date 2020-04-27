This README describes how to compile and execute the code, the input files needed by MOMICE_ASTRO, and the routine that analyzes and visualises the results.

## 1. Compilation and Execution

The code is written in Fortran 90, it is therefore necessary to compile its source files located in the `source` directory before running the model. 

To compile the source code: 
- go to the `source` directory. 
- type "make gfort" or "make ifort" depending on the compiler (gfortran or ifortran) you wish to use. 
- a `momice.gfort` or `momice.ifort` executable file is created.

The fortran executable is called by the MOMAPP python routine (`momapp.py`) also located in the `source` directory. MOMAPP reads a first input file `momapp.in` in which the user defines the type of simulation to be run and the location of input files read by MOMICE. See the next section for a description of `momapp.in` and all the input files read by MOMICE.

Once the input files are configured, the simulations can be run by executing the `momapp.py` python routine. 

***

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
	- dust-to-gas mass ratio R<sub>dg</sub>, the default value is 1%. 
	- grain size a<sub>d</sub>, the default value is 5x10<sup>-2</sup> um.
	- volumic mass of grains rho<sub>d</sub>, the default value is 3 g cm<sup>-3</sup>
From these parameters, MOMICE determines the grain abundance and the grain cross section. Here, only one constant grain size can be specify.

- Ice properties: 
	- site size d<sub>s</sub>
	- sticking coefficient for species heavier than H and H<sub>2</sub>
	- diffusion-to-binding energy ratios by distinguishing the surface and the bulk of the ice, and "light" (H, H<sub>2</sub>) and "heavy" species
	- number of active surface monolayers and approximate number of timesteps needed to fill one monolayer
	- barrier width (in angstroms) for surface reactions and for diffusion through tunnelling if they are activated (see model switches below)
The values of these ice parameters are extensively discussed in my thesis and in the papers listed in the README of the main MOMICE directory. Please refer to them for the best choice.

- Porosity parameters: a 3D porosity treatment is yet to be included. Please disregard the porosity parameters at the moment. 

#### MODEL SWITCHES

This section allows you to switch between different chemistry options:
- Phase: 2-phase or 3-phase approach. In the "old" 2-phase approach, the ice is considered as homogeneous whereas the 3-phase approach distinguishes the ice surface (with a specified number of monolayers) and the ice bulk mantle. The 3-phase approach follows the method introduced by Hasegawa & Herbst (1993) that allows a faster computation with respect to the approach by Taquet et al. (2012). 
- Rmod: Modified rates. One can include the modified rate approach proposed by Garrod (2008). At this moment, it's still in beta, please do not use it yet. 
- GSevo: Grain size evolution due to the growth of ice, assuming that the thickness of each monolayer is equal to the site size of a site. 
- Ebevo: Binding energy evolution at the surface and in the bulk as a function of the coverage of bare grain, H2, and water. 
- surfc: Turn on/off chemistry at the ice surface.
- bulkc: Turn on/off chemistry in the ice bulk.
- DiffH: Choose whether the diffusion of light particles (H) is governed by thermal hopping, quantum tunnelling, or both.
- DiffO: Choose whether the diffusion of heavy atoms (C, N, O) is governed by thermal hopping, quantum tunnelling, or both.
- XgasH: One can choose to keep the density or the abundance of atomic hydrogen constant if, for example, you’re using a small chemical network which does not allow you to self consistently compute the abundance of H.
- H2opr: Whether or not keep the H<sub>2</sub> ortho/para ratio constant (if ortho/param states are included in the network).
- Pcomp: Choose "1" if the reaction rate is directly computed from the transmission probability of the reaction (Hasegawa et al. 1992 approach for rate equations). Choose "2" if the reaction rate is computed from the competition process between the transmission probability of the reaction and the diffusion rate of the reactants (Chang et al. 2007 approach for microscopic Monte-Carlo calculations).  
- Ptran: Choose "1" to compute the transmission probability by assuming a rectangular barrier for all reactions (in that case, the activation energies need to be specified in the reaction network and the width of the barrier is specified just below). Choose "2" to compute the transmission probability by using the Eckart model which fits an approximate Potential Energy Surface for the reactions available in the data file `eckart_data.in`. 

#### ODE SOLVER PARAMETERS

- ATOL0: Initial absolute tolerance set to the DLSODES solver, default value is 1e-20
- RTOL0: Initial relative tolerance set to the DLSODES solver, default value is 1e-4
- ITOL: default value is 2
- ITASK: default value is 1


### `species_file.in`

The species file defined in the `input_parameters.in` file lists the species included in the network, together with their number of elements, their mass, and their charge, and their initial abundance relative to the total number of H nuclei. 
MOMICE uses this data to check the conservation of mass and charge. Please do not forget to specify the number of species at the beginning of the file.

The species network must respect several rules to be read correctly by the fortran code:
- Neutral gas phase species are listed first, followed by ions, ice surface, and ice bulk species. Neutral gas phase species need to have a icy counterpart (at the surface and in the bulk if 3-phase approach is used). In other words, the list of gaseous neutral species, surface species, and porous/bulk species need to be identical. 
- Names of surface species start with a J character and the name of bulk species start with a Q.

We advise the user to use the species files available here as templates to build their own species network.

 
### `reaction_file.in`

 This file contains the network of chemical reactions. Each line refers to a specific reaction and includes:
 - 1 to 3 reactants (with 15 caracters each)
 - 1 to 4 products (with 15 caracters each)
 - 3 parameters A, B, C
 - the type of reaction
 - the range of temperatures where the reaction is considered
 - the type of the formula used to compute the rate. 

For gas phase reactions (from 1 to 8), the reaction rates are computed following the methodology adopted in KIDA. A description of the types and formulas can be found [here](http://kida.astrophy.u-bordeaux.fr/help.html). In addition to the standard gas-phase reactions, GRAINOBLE takes the following processes into account: 
- 0: Gas-grain interaction and electron-grain recombination (Rate(s<sup>-1</sup>) = A x zeta) (from Flower and Pineau des Forets 2003). 
- 20: Accretion from gas phase to grain surfaces. A is the pre-factor constant. 
- 14: Langmuir-Hinshelwood grain-surface chemical reactions. A is the pre-factor constant, B is the efficiency induced by the chemical desorption, C is the activation barrier for computing the transmission probability through a rectangular barrier. 
- 21: Thermal evaporation. A is the pre-factor constant. 
- 22: Cosmic Ray induced general desorption following Hasegawa & Herbst 93. A=branching ratio. 
- 23: Photo-desorption by secondary CR-induced photons + background photons coming from wavelength-dependent experimental studies (see Taquet et al. 2013). A=branching ratio. 
- 26: Photodissociation/photodesorption processes coming from molecular dynamics simulations (see Taquet et al. 2013). A is the absorption cross section, B is the branching ratio between different products. 
- 27: Ice sputtering through cosmic rays, scaled on water ice efficiency measured experimentally by Dartois et al. (2015).
- 30: Non-porous surace->Pores rate exchange. A=branching ratio. 
- 31: Pores->Non-porous rate following the method introduced by Taquet et al. (2012). A=branching ratio. 
- 32: Surface<->mantle exchange rate following the method introduced by Hasegawa & Herbst (1993).


### `energy_file.in`

 This file contains the list of binding energies of surface species relative to a bare grain surface (amorphous carbon and/or silicate), amorphous water ice, H2 ice, and as a pure substrate. For most species, it also contains the references where the values are taken from. This file should be considered as a data file since it contains all possible species which are not necessarily included in the species network. However, every surface species included in the species network must be located somewhere in this file. 

### `reacrates.in`

This file contains the list of species whose formation and destruction rates will be computed and saved as binary files. For each species, the rates of the 30 most important reactions of formation and destruction are saved. The amount of data tends to quickly increase with the number of species, one should limit the number of species to 50.


### `spat_evol.in`

If a spatial evolution is chosen in the `momapp.in` input file (as type of simulation 2), the folder and the list of files giving the one-dimensional evolution of the physical conditions have to be specified in this file. The format of the input files follows the format given in the examples located in the "phys" folders. 


### `model_grid.in`

 If a model grid is chosen in the `momapp.in` input file (as type of simulation 3), one needs to specify a file with the name defined in `momapp.in`. This file lists the physical parameters whose ranges of values are explored, with four columns: 
 - name of parameter among the followinng physical parameters: nH (total density in cm<sup>-3</sup>), Tg (gas temperature in K), Td (dust/ice temperature in K), Av (visual extinction in mag), UV (flux of secondary photons in cm<sup>-2</sup>s<sup>-1</sup>), CR (cosmic ray ionization rate in s<sup>-1</sup>)
- number of values 
- minimal value
- maximal value

### `sensitivity.in`

If a sensitivity analysis in the `momapp.in` input file (as type of simulation 4), one needs to specify a file with the name defined in `momapp.in`. This file first sets the number of models to be run, and the distribution followed by the parameters (either `uniform` or `normal` distribution). The user then defines the parameters whose values are explored with four columns:
- name of parameter among the following chemical and surface parameters: Ed (diffusion-to-binding energy ratio), Eb (binding energy factor with respect to those defined in the energy data file), Ea (activation barrier factor of surface reactions with respect to those defined in the reaction file)
- the mean of the distribution
- the standard deviation (for a gaussian distribution) or the range (for a uniform distribution)



## 3. Output analysis and Visualisation

The output files are saved in binary files located in the selected output directory. A Python routine, called `moman.py` (MOMAN for MOMice ANalyzer, in Python 3) is provided in the output directory. It reads and analyzes the binary files generated by MOMICE, and plots the results in pdf files. Before running MOMAN, the user needs to define the options in the input `moman.in` file:

- typeout: Type of output: 1) individual models, 2) spatial evolution, 3) model grid, 4) sensitivity analysis 
- outdir: Location of the output directory where model folders are saved relative to current one
- outmod: Name of the folder(s) to read (one or several folders if typeout=1 chosen, only one otherwise)
- prefix: Name of the prefix of figures that compare different models if typeout=1:  
- plotreac: Whether or not plot the rates of reactions ('yes' or 'no')
- refsp: Reference species for abundance ratios
- species: Name of the species shown in figures
- invtime: Whether or not inverse time on x-axis ('yes' or 'no')
- timein: Initial time if different from 0 (in yr)
- timeplot: Timesteps for figures showing distributions, spatial or parameter evolutions (for typeout 2 to 4)
- logabs: Logarithmic (log) or linear (lin) evolution of time on x-axis


### Individual model(s)

If typeout=1, MOMAN creates the pdf files in the directory of each selected model. They show the evolution of various parameters (physical conditions, gas and ice abundances, gas and ice deuterium fractionation, CO depletion, fractional composition within the ices, ice thickness, H2 o/p ratio) as function of time, CO depletion, and ice thickness, and physical parameters. If chosen, MOMAN also presents the contribution of the major reactions to the formation and the destruction of specified species.

### Spatial evolution

If typeout=2, MOMAN allows you to follow the chemical evolution in a one-dimensional dimension. outdir, specified in `moman.in`, determines the location of the set of output folders and where the `spatlog.out` file is located. 
MOMICE generates a set of pdf and ascii files in different folders showing the evolution of the chemical abundances as function of the distance from the central source for the different times chosen in `moman.in`. 


### Model grid

If typeout=3, MOMAN allows you to analyse a model grid in which the values of physical conditions are explored. This allows you to assess the importance of physical conditions on the chemical abundances. First, you need to specify the location of the grid directory and its name, through the outdir parameter in `moman.in`, that contains all the models of the grid and the `gridlog.out` file. 

MOMAN creates the pdf files in the `figures` folder located in model grid directory. MOMAN generates distributions of chemical abundances for various timesteps (defined as timeplot in `moman.in`). If values of one physical parameter are explored, MOMAN shows the evolution of abundances as function of the parameter values. If the values of two physical parameters are explored, MOMAN shows the evolution of abundances as function of parameter values as 2D images, with one figure per species. 

### Sensitivity

If typeout=4, MOMAN allows you to visualise the result of a sensitivity analysis in which the values of surface and chemical parameters are explored. This allows you to explore how these parameters affect the "uncertainty" on the chemical predictions. Several types of figures are generated: 1) like for individual models, figures showing the evolution of chemical abundances as function of time of all models, but they get messy when the number of models is high, 2) figures the evolution of averaged abundances as function of and their standard deviation induced by the variation of input parameters, 3) distribution of chemical abundances for different timesteps (defined as timeplot in `moman.in`), 4) figures showing the chemical abundances for different timesteps as function of the input parameter values. 

In addition, figures showing the distributions of input parameter values are shown in order to check that they correctly follow the desired distribution (normal or uniform). 


## 4. Data files

In addition to the input files described above and located in a specific folder within the `input` directory, MOMICE reads several files located in the `data` folder for various purposes. 

### `eckart_data.in`

This file lists reactions with parameters that are used to compute the transmission probability with the Eckart model. The three parameters are: the activation barrier of the forward reaction (Vf), the activation barrier of the reverse reaction (Vr), and the frequency of the imaginary mode of the transition state (wf). Here, most reactions are involed in the methanol and water networks but feel free to update the list if you are more data. 


### `self_shielding.in`

`self_shielding.in` and `self_shielding_PDRmeudon.in` list the self-shielding function as function of the species column density for several key species: H<sub>2<\sub>, HD, and CO. 
