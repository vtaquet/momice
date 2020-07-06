MODULE VARIABLES

  IMPLICIT NONE
  SAVE

  !----------------------------------------
  ! Types of structures used in all modules
  !----------------------------------------

  ! kind of variables
  INTEGER, PARAMETER :: DP=selected_real_kind(P=15) ! real kind
  INTEGER, PARAMETER :: LONG=KIND(1)                ! integer kind

  ! parameters for size of vectors
  INTEGER(KIND=LONG), PARAMETER :: Nelements=14
  INTEGER(KIND=LONG), PARAMETER :: Neckartreacs = 77
  INTEGER(KIND=LONG), PARAMETER :: Nenthalpy = 51
  INTEGER(KIND=LONG), PARAMETER :: NCOnetreacs = 25
  INTEGER(KIND=LONG), PARAMETER :: Ntype = 99

  ! constants
  REAL(KIND=DP), PARAMETER :: Pi = 3.14159265
  REAL(KIND=DP), PARAMETER :: G = 6.67384d-11
  REAL(KIND=DP), PARAMETER :: Mpr = 1.672623D-27
  REAL(KIND=DP), PARAMETER :: Me = 9.10938188D-30
  REAL(KIND=DP), PARAMETER :: kb = 1.3806504D-23
  REAL(KIND=DP), PARAMETER :: Rys = 31536000
  REAL(KIND=DP), PARAMETER :: hb = 6.62606896D-34/(2*Pi)
  REAL(KIND=DP), PARAMETER :: cmtos = 2.99792e10
  REAL(KIND=DP), PARAMETER :: mdyneatonm=1e-18
  REAL(KIND=DP), PARAMETER :: conv=1.69733088E+6

  ! loops, logicals, ...
  INTEGER(KIND=LONG) :: ok
  LOGICAL :: precision
  CHARACTER :: toto

  ! output directories
  TYPE :: outputdirectory
     CHARACTER(LEN=8) :: date
     CHARACTER(LEN=8) :: date_ini, date_fin
     CHARACTER(LEN=10) :: realtime
     CHARACTER(LEN=4) :: realtime2
     CHARACTER(LEN=10) ::  time_ini, time_fin, test
     CHARACTER(len=300) :: direct_mod
     CHARACTER(len=300) :: direct_grid
     CHARACTER(len=300) :: direct_gridlog
  END TYPE outputdirectory

  ! parameters of input_parameters.in
  TYPE :: inputparam
     ! data options
     CHARACTER(len=80) :: comments, folder
     INTEGER(KIND=LONG) :: chdata, chgrid, chreacrates
     CHARACTER(len=100) :: name_output, location_output
     ! physical conditions
     REAL(KIND=DP) :: tmax, tmaxini, nHini, nHfin, Tgini, Tdini, Tini, Tfin
     REAL(KIND=DP) :: Avini, zetaini, UVfluxCRini, G0ext, tff2
     CHARACTER(len=3) :: chphys, chinphys
     CHARACTER(len=200) :: inputphys
     REAL(KIND=DP) :: dens_evol, temp_evol
     ! grain parameters
     REAL(KIND=DP) :: Rdg, acst, rhod
     ! ice properties
     INTEGER(KIND=LONG) :: chlayer, chmulti, chmodifrates, chgrainevol, Nsteps_lay
     INTEGER(KIND=LONG) :: chevolen, chHtunnel, chCtunnel, chsurfchem, chbulkchem
     REAL(KIND=DP) :: stick_coef, REdH, REdoth, REdbulkH, REdbulkoth, ds, laywid, Nsurf
     CHARACTER(len=3) :: chlayer2
     ! grain porosity
     REAL(KIND=DP) :: Slat, Fporsur, Fvacuum, Fedgesur, Fnpsur
     ! chemical parameters
     CHARACTER(len=65) :: filesp, filereac, fileen
     INTEGER(KIND=LONG) :: chH, chH2op
     INTEGER(KIND=LONG) :: chreacproba, checkart, chDnetwork
     REAL(KIND=DP) :: rectwid, H2opr, diffwid, Fnd
     INTEGER(KIND=LONG) :: Nspecies, Nreactions, Nspgr, Nspgr2, Nspgas, Nspnd
     INTEGER(KIND=LONG) :: Nsprate, Nspeciestot, Neq
     ! ODE parameters
     REAL(KIND=DP) :: ATOL0, RTOL0
     INTEGER(KIND=LONG) :: ITASK0, ITOL0
     ! random parameters
     CHARACTER(len=10) :: distrib
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: randomEb, randomEa
  END TYPE inputparam

  TYPE :: inputphysparam
     INTEGER(KIND=LONG) :: Nphysteps
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: time
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: rad
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: veloc
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: nH
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Tg
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Td
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Av
  END TYPE inputphysparam


  ! parameters of model_grid.in
  TYPE :: gridini
     INTEGER(KIND=LONG) :: NT, NnH, NUV, Nzeta, NAv
     INTEGER(KIND=LONG) :: Nacst, Ntreat, Nenratio, Nporo, Nds, NEbH
     INTEGER(KIND=LONG) :: NEa, NXO, Nop, Nmodels
     INTEGER(KIND=LONG) :: NT2, NnH2, NUV2, Nzeta2, NAv2
     INTEGER(KIND=LONG) :: Nacst2, Ntreat2, Nenratio2, Nporo2, Nds2
     INTEGER(KIND=LONG) :: NEbH2, NEa2, NXO2, Nop2
     INTEGER(KIND=LONG) :: choicetreat, Nparam, chabuini
     REAL(KIND=DP) :: Tlow, Tup, nHlow, nHup, EbHlow, EbHup
     REAL(KIND=DP) :: acstlow, acstup
     REAL(KIND=DP) :: enratiolow, enratioup, porolow, poroup, dslow
     REAL(KIND=DP) :: dsup, Ealow, Eaup, XOlow, XOup
     REAL(KIND=DP) :: UVlow, UVup, oplow, opup, zetalow, zetaup
     REAL(KIND=DP) :: Avlow, Avup
     CHARACTER(len=10) :: char_nH, char_T, char_dist, char_treat, char_EbH
     CHARACTER(len=10) :: char_UV, char_en, char_poro, char_ds
     CHARACTER(len=10) :: char_Ea, char_XO, char_op, char_Av, char_zeta
     CHARACTER(len=100) :: input_abuini
  END TYPE gridini

  TYPE :: spatialevol
     INTEGER(KIND=LONG) :: Nspat
     CHARACTER(len=100), DIMENSION(:), ALLOCATABLE :: file
     CHARACTER(len=100) :: folder
  END TYPE spatialevol

  
  ! physical properties 
  TYPE :: physparam
     REAL(KIND=DP) :: nH, nH_prev, dnH, dnH_prev, Tg, Tg_prev, Td, Td_prev
     REAL(KIND=DP) :: Av, Av_prev, zeta, zeta_prev, UVfluxCR, UVfluxCR_prev, v, v_prev
     REAL(KIND=DP) :: dtime, time, time_prev, DTOGN, Nsndtot, signdtot, rad, rad_prev
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Vth
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: stick
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Kreac
  END TYPE physparam


  ! species properties 
  TYPE :: species
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: num
     CHARACTER(LEN=15), DIMENSION(:), ALLOCATABLE :: name
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: charge
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: mass
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: element
     CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE :: name_elem
     CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE :: rate_name
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: rate_num
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Eb
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Ed
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Eb_carb
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Ed_carb
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Eb_wat
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Ed_wat
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Eb_hh93
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Eb_H2
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Ed_H2
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Eb_pure
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Ed_pure
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Eb_lay
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Eb_ini
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Ed_lay
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Ed_ini
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: cr_coef
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xini
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: dHf
     INTEGER(KIND=LONG), DIMENSION(:,:), ALLOCATABLE :: grproc
     INTEGER(KIND=LONG), DIMENSION(:,:), ALLOCATABLE :: formreac
     INTEGER(KIND=LONG), DIMENSION(:,:), ALLOCATABLE :: destreac
     INTEGER(KIND=LONG), DIMENSION(:,:), ALLOCATABLE :: maxreac
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: rateform
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ratedest
     INTEGER(KIND=LONG), DIMENSION(:,:), ALLOCATABLE :: numform
     INTEGER(KIND=LONG), DIMENSION(:,:), ALLOCATABLE :: numdest
     INTEGER(KIND=LONG), DIMENSION(:,:), ALLOCATABLE :: degform
     INTEGER(KIND=LONG), DIMENSION(:,:), ALLOCATABLE :: degdest
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: Nreacform
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: Nreacdest
     INTEGER(KIND=LONG) :: maxNreacform, maxNreacdest, maxNreac, maxNreac2
     CHARACTER(len=173) :: fmt 
  END TYPE species
  TYPE :: specnum
     INTEGER(KIND=LONG) :: H, D, He, JOH, JOD
     INTEGER(KIND=LONG) :: H2, oH2, pH2, D2, oD2, pD2
     INTEGER(KIND=LONG) :: HD, CO, O, O2, H2O, CO2, C, N, CH3OH
     INTEGER(KIND=LONG) :: H2DPLUS, E, H3PL 
     INTEGER(KIND=LONG) :: JH, JD, JHD, JCO, JH2O, JCN
     INTEGER(KIND=LONG) :: JCO2, JO, JH2, JoH2, JpH2, JC, JN, JCH3OH
     INTEGER(KIND=LONG) :: JD2, JoD2, JpD2
     INTEGER(KIND=LONG) :: G0, GM, GP, grinput
  END TYPE specnum

  ! reaction properties
  TYPE :: reactions
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: num_react1
     CHARACTER(LEN=15), DIMENSION(:), ALLOCATABLE :: react1
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: num_react2
     CHARACTER(LEN=15), DIMENSION(:), ALLOCATABLE :: react2
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: num_react3
     CHARACTER(LEN=15), DIMENSION(:), ALLOCATABLE :: react3
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: num_prod1
     CHARACTER(LEN=15), DIMENSION(:), ALLOCATABLE :: prod1
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: num_prod2
     CHARACTER(LEN=15), DIMENSION(:), ALLOCATABLE :: prod2
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: num_prod3
     CHARACTER(LEN=15), DIMENSION(:), ALLOCATABLE :: prod3
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: num_prod4
     CHARACTER(LEN=15), DIMENSION(:), ALLOCATABLE :: prod4
     CHARACTER(LEN=15), DIMENSION(:,:), ALLOCATABLE :: ndreact1
     INTEGER(KIND=LONG), DIMENSION(:,:), ALLOCATABLE :: num_ndreact1
     CHARACTER(LEN=15), DIMENSION(:), ALLOCATABLE :: ndreact2
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: num_ndreact2
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: A
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: B
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: C
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: type
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: num
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: formula
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: Tmin, Tmax
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: Nrate
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: sim
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: act_en
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: proba
     CHARACTER(len=69) :: fmt
  END TYPE reactions
  TYPE :: reactionnumber
     INTEGER(KIND=LONG) :: acc_ini, acc_fin, reac_ini, reac_fin
     INTEGER(KIND=LONG) :: evth_ini, evth_fin, nondiffreac_ini, nondiffreac_fin
     INTEGER(KIND=LONG) :: evCR_ini, evCR_fin, ic_accH2, evcol_ini, evcol_fin
     INTEGER(KIND=LONG) :: evUVCR_ini, evUVCR_fin
     INTEGER(KIND=LONG) :: disUVCR_ini, disUVCR_fin
     INTEGER(KIND=LONG) :: COH, ref_photo, reacCO2, accH2, reacref
     INTEGER(KIND=LONG), DIMENSION(0:Ntype) :: istart, iend
  END TYPE reactionnumber
 
  ! grain properties
  TYPE :: grainparam
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Nssur
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Nstot
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Npore
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Sdep
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Spore
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Fnptot
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Fportot
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Fedgetot
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Vd
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Vvacuum
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: a
     REAL(KIND=DP) :: X, ad
!!$     REAL(KIND=DP) :: aini, Xini, Nssurini, Nstotini
!!$     REAL(KIND=DP) :: Nporeini, Sdepini, Sporeini
!!$     REAL(KIND=DP) :: Fnptotini, Fportotini, Fedgetotini
!!$     REAL(KIND=DP) :: Vdini, Vvacuumini
!!$     REAL(KIND=DP) :: afin, Xfin, Ns_surfin, Nstotfin
!!$     REAL(KIND=DP) :: Nporefin, Sdepfin, Sporefin
!!$     REAL(KIND=DP) :: Fnptotfin, Fportotfin, Fedgetotfin
!!$     REAL(KIND=DP) :: Vdfin, Vvacuumfin
!!$     REAL(KIND=DP) :: Xgrain
  END TYPE grainparam

  ! abundances in gas phase and populations on grains
  TYPE :: abundances
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xgas_ini
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xgas
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xgas_prev
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xice
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xice_prev
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xsur
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xsur_prev
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xbulk
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xbulk_prev
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xtot
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xtot_prev
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xelem
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xelem_prev
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xelem_ini
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xelem_lay
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xcharge
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xcharge_prev
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Psur
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Psur_prev
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Pbulk
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Pbulk_prev
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xsurml
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xsurml_prev
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xbulkml
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Xbulkml_prev
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Ptotml
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Ptotml_prev
     !REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Pfracml
     !REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Pfracml_prev
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Pice
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Pice_prev
     !REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pop
     !REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pop_prev
     !REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pop_lay
     !REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pop_lay_prev
     !REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pop_tot
     !REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pop_tot_prev
     !REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: pop_spec
     REAL(KIND=DP) :: pop_int, abuH2, abuH2bulk, abuH, abuCO, Nlay, Nlay_prev
     REAL(KIND=DP) :: k_evgr_tot, max_kgr, maxlayers, Psurtot, Pbulktot, Ptot
  END TYPE abundances

  ! output vectors
  TYPE :: outputparam
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: v
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: rad
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: time
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: timefix
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: nH
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Av
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Tg
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Td
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Nlay
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Xgas
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Xice
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Xsur
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Xsurml
     !REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: pop_spec
     !REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: pop_tot
     INTEGER(KIND=LONG), DIMENSION(:,:,:), ALLOCATABLE :: numform
     INTEGER(KIND=LONG), DIMENSION(:,:,:), ALLOCATABLE :: numdest
     INTEGER(KIND=LONG), DIMENSION(:,:,:), ALLOCATABLE :: numform2
     INTEGER(KIND=LONG), DIMENSION(:,:,:), ALLOCATABLE :: numdest2
     REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: rateform
     REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ratedest
     REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: rateform2
     REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ratedest2
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: contribform
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: contribdest
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: rateformtot
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: ratedesttot
     INTEGER(KIND=LONG) :: Ninisteps, Nplotsteps, Nrealsteps, Ninilayer
  END TYPE outputparam

  ! model grid properties
  TYPE :: gridparam
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: nH
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: T
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: UV
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: zeta
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Av
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: acst
     CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: treatment
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: treatment2
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: enratio
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Fpor
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ds
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Ea
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XO
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: EbH
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: op
     CHARACTER(len=200), DIMENSION(:), ALLOCATABLE :: direct
     CHARACTER(len=200), DIMENSION(:), ALLOCATABLE :: Xini_file
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Xgasini
     REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Xiceini
  END TYPE gridparam
  
  ! self shielding parameters
  TYPE :: selfshielding
     REAL(KIND=DP), DIMENSION(52) :: NCO, T_CO, NH2_1, T_H2_1, AV2, T_AV
     REAL(KIND=DP), DIMENSION(105) :: NH2_2, T_H2_2
     REAL(KIND=DP), DIMENSION(203) :: NH2_3, T_H2_3, T_HD
  END TYPE selfshielding
  
  ! eckart model parameters
  TYPE :: eckartmod
     CHARACTER(len=10), DIMENSION(Neckartreacs) :: react1
     CHARACTER(len=10), DIMENSION(Neckartreacs) :: react2
     CHARACTER(len=10), DIMENSION(Neckartreacs) :: react3
     CHARACTER(len=10), DIMENSION(Neckartreacs) :: prod1
     CHARACTER(len=10), DIMENSION(Neckartreacs) :: prod2
     REAL(KIND=DP), DIMENSION(Neckartreacs) :: Vf, Vr, wf
  END TYPE eckartmod
  
  ! enthalpy parameters
  TYPE :: enthalpy
     CHARACTER(len=10), DIMENSION(Nenthalpy) :: species
     REAL(KIND=DP), DIMENSION(Neckartreacs) :: dHf
  END TYPE enthalpy

  ! CO network parameters
  TYPE :: COnetwork
     CHARACTER(len=10), DIMENSION(NCOnetreacs) :: react1, react2, react3
     CHARACTER(len=10), DIMENSION(NCOnetreacs) :: prod1, prod2
     REAL(KIND=DP), DIMENSION(NCOnetreacs) :: A, B
  END TYPE COnetwork

  ! ODE parameters
  TYPE :: varode
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Y
     INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK, IPAR, ISTATS, IOUT
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RWORK, RPAR, RSTATS, RTOL
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ATOL, ATOLini, ROUT
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Ysave, IWK, Yini
     INTEGER :: ITASK, ISTATE, nonzerocount
     INTEGER :: ITOL, MF, LRW, LIW, IOPT
     CHARACTER(20) :: JAC
     INTEGER(KIND=LONG), DIMENSION(:,:),ALLOCATABLE :: nonzero
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: IA
     INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: JA
  END TYPE varode

  ! date and time to create the folder



END MODULE VARIABLES

MODULE SHARED_VARIABLES

  !----------------------------------------
  ! Variables shared between main and Fode
  !----------------------------------------

  USE VARIABLES 

  IMPLICIT NONE 

  TYPE(inputparam) :: inp
  TYPE(physparam) :: ph
  TYPE(species) :: sp
  TYPE(specnum) :: ns
  TYPE(reactions) :: re
  TYPE(reactionnumber) :: nr
  TYPE(grainparam) :: gr
  TYPE(selfshielding) :: ss
  TYPE(abundances) :: ab
  TYPE(inputphysparam) :: iph
  TYPE(eckartmod) :: ec
  TYPE(enthalpy) :: ent
  TYPE(COnetwork) :: con
  TYPE(outputparam) :: out
  TYPE(gridini) :: gi
  TYPE(gridparam) :: gp
  TYPE(spatialevol) :: spat
  TYPE(outputdirectory) :: di
  TYPE(varode) :: ode
  TYPE(varode) :: ode2
  INTEGER(KIND=LONG) :: il, itest, itest2, nonzerocount
  CHARACTER(len=80) :: ERR_MESS
  REAL(KIND=DP) :: alphacc, alphades, F2tot, nstot, nmtot, dXsurf

END MODULE SHARED_VARIABLES


MODULE SHARED_VARIABLES2

  !----------------------------------------
  ! Variables shared between main and Fode
  !----------------------------------------

  USE VARIABLES 

  IMPLICIT NONE 
!!$
  INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: inzI
  INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: inzJ
  INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: inzIA
  INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: inzJA
  INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: inzIord
  INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: inzIord2
  INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: inzIordini
  INTEGER(KIND=LONG), DIMENSION(:), ALLOCATABLE :: inzIord2ini
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: inzNZord
  INTEGER(KIND=LONG) :: inzNZ
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Pbulkact
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: Pbulkact0

END MODULE SHARED_VARIABLES2
