*******************************************************************
HISTORY OF GRAINOBLE VERSIONS 
*******************************************************************
*
*******************************************************************
GRAINOBLE_2013.06
*******************************************************************
This version is essentially the version developed during my PhD and 
should be considered as the IPAG version. 
A full description of the processes, equations, and chemical 
networks is provided here: 
http://adsabs.harvard.edu/abs/2012PhDT........23T. 
*
INPUT FILES AND CHEMICAL NETWORKS:
Chemical networks provided here are those from Taquet et al. 
(2012a, 2012b, 2013a). An "old" ASCII format is used for the 
species, and reaction files (species with 10 caracters, no 
temperature range, no choice of formula).
The 2012a one considers grain surface chemistry only (+gas-grain 
interactions) to follow a simple ice formation. 
The 2013a is based on KIDA and updated to follow the evolution of 
atomic D in the gas phase. Surface network follows the formation 
and deuteration of main icy species.
Photodissociation of icy species is followed following the results 
of MD simulations.
*
CODE AND PROCESSES:
The "official" multilayer method, tested against several networks 
and physical cases, is the one described in Taquet et al. (2012).
A porosity treatment is included in this version but has been 
abandoned in the subsequent versions. 
*
PHYSICAL EVOLUTION:
This version only considers only pseudo-time dependent evolution of
the chemistry: the physical conditions remain constant at 0D. 
*
*******************************************************************
GRAINOBLE_2014.05 (Taquet, Charnley, Sipila 2014, ApJ)
*******************************************************************
This version has been developed in 2013 during the first half of my 
post-doc at NASA Goddard with S. Charnley. The main novelties 
concern the multilayer treatment that follow Hasegawa et al. (1993) 
to consider the evaporation of ices, and the incorporation of a 1D 
physical collapse model. 
*
INPUT FILES AND CHEMICAL NETWORKS:
The new chemical network used in Taquet et al. (2014) is called: 
"re_Dnetwork_short_new.in"). 
The gas phase chemical network based on the 2013 KIDA network has 
been "deuterated" for hydrogenated species with four or less atoms 
and includes the spin states of H2, H2+, and H3+ with routines from 
O. Sipila (Sipila et al. 2013).
The surface network includes the formation of COMs following the 
networks from Hasegawa et al. (1992, 1993) and Garrod et al. (2006, 
2008).
*
CODE AND PROCESSES:
From this version, the "official" multilayer method is the one 
described in Hasegawa et al. (1993). Although it adds a third set 
of chemical reactions, it allows us to follow the bulk chemistry 
and the ice evaporation. 
The parameters for the porosity treatment are still shown, but
one should take this option with caution since it hasn't been 
tested on this version. 
*
PHYSICAL EVOLUTION:
This version has been applied against a 1D physical collapse model 
to follow the ice formation during the static formation of a 
prestellar core followed by the ice warm-up and evaporation during 
the free-fall collapse of the protostellar envelope. The dynamical 
evolution follow the equations from Whitworth & Ward-Thompson 
(2000), while the dust temperature profiles are computed with the 
DUSTY radiative transfer code (Izevic & Ellitzur 1997).
The "phys_time" folder contains all the models considered in 
Taquet et al. (2014).
*
*******************************************************************
GRAINOBLE_2016.02 (Taquet, Wirstrom, Charnley 2016, ApJ)
*******************************************************************
This version of the code is the same as GRAINOBLE_2014.05. Special 
attention has been paid for ion-neutral reactions to follow hot core
chemistry during protostellar outbursts. 
*
INPUT FILES AND CHEMICAL NETWORKS:
All input parameters are now located within the same directory 
specified in "input_files.in". 
The chemical network called "re_HCchem_*.in" used here is based on 
the Hot-Core network developed by S. Charnley in the 90s and updated 
by E. Wirstrom. It only contains gas phase reactions and gas-grain 
processes.
* 
CODE AND PROCESSES:
No major modifications have been implemented in this version.
*
PHYSICAL EVOLUTION:
This version was used to follow the Hot-Core chemistry during 
protostellar outbursts. A set of files describes the physical 
evolution during typical outbursts as presented in Taquet et al. 
(2016a). 
*
*******************************************************************
GRAINOBLE_2016.08 (Taquet, Furuya, Walsh, van Dishoeck 2016, MNRAS)
*******************************************************************
This version has been developed mostly in 2015 during the first 
half of my post-doc in Leiden. It contains several modifications 
with respect fo 2016.02, mostly regarding the efficiency of the 
code and the input parameters. 
*
INPUT FILES AND CHEMICAL NETWORKS:
The "input_parameters.in" layout has been rearranged with clear 
subsections. Ice parameters between the surface and the bulk have 
been distinguished. 
The chemical network used in this version is a combination of the 
non-deuterated network from Taquet et al. (2014) and the Hot-Core 
network from Taquet et al. (2016a) supplemented by reactions from 
Balucani et al. (2015) and Charnley (2009). 
* 
CODE AND PROCESSES:
The routines describing the ODEs have been rearranged. In particular, 
individual routines have been added for the computation of the 
jacobian matrices for each type of reaction in order to optimize 
the computation time for ODE solving by the DLSODE solver. To this 
aim, reactions are now listed type by type in the reaction file and 
can no longer be mixed up. 
Bulk chemistry has been added and tested. Bulk species are now 
explicitely mentioned in the species and reaction files starting 
with "Q". Surface-bulk exchange processes are also explicitely 
written in the reaction file with reaction type 32. 
Chemical desorption probability is now explicitely defined in 
the reaction file with the B parameter. 
*
PHYSICAL EVOLUTION:
Similar to 2016.02, several files describing protostellar outbursts 
but for disk conditions are located in "phys_time".
*
*******************************************************************
GRAINOBLE_2017.09 (van't Hoff et al. 2018)
*******************************************************************
GRAINOBLE was used to predict the HCO+ and DCO+ abundance profiles 
for the IRAS4A physical profile. 
*
INPUT FILES AND CHEMICAL NETWORKS:
Same type of input parameters and chemical networks as 
GRAINOBLE_2016.08 were used in this version.
* 
CODE AND PROCESSES:
The code is essentially the same as GRAINOBLE_2016.08. The 
"non-diffusion" reaction type was implemented but not used in the 
chemical network.
*
PHYSICAL EVOLUTION:
The chemical evolution was followed for the IRAS4A physical profile, 
by assuming conditions physical conditions.
*
*******************************************************************
GRAINOBLE_2018.08 
*******************************************************************
*
INPUT FILES AND CHEMICAL NETWORKS:
- Re-arrangement of reactions so types are in good order and 
double-checking of "similar" reactions between different types
- Add reactions from Skouteris et al. (2017)
- Non-diffusion reaction type (18) added
- TO FINISH: "Thermal" reactions studied by Theule et al. (2013) and
others to be added (NH4+OCN-)
- TO FINISH: Abstraction reactions with H and OH from Garrod (2013) 
and Lamberts et al.
* 
CODE AND PROCESSES:
- Implementation of "non-diffusive" reaction rate and test it for 
abundant reactants studied in experiments.
- Implementation of formula by Minissale et al. (2016) and Vasyunin 
et al. (2017) estimating the chemical desorption probability from 
the surface composition.
