MODULE OUTPUT

  !-----------------------------------------------------
  ! This module contains several subroutines that 
  ! create directories
  ! write the main results in the terminal
  ! write the main results in a ASCII log file
  ! write the results in binary files
  !-----------------------------------------------------

  USE VARIABLES
  USE PHYSPROP
  USE CHEMPROP


CONTAINS


  !--------------------!
  !--------------------!
  ! CREATE DIRECTORIES !
  !--------------------!
  !--------------------!


  SUBROUTINE OUTPUT_DIRECTORY

    !--------------------------------------------------------
    ! create the directory and the folder of the output data
    !--------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, di

    CHARACTER(len=80) :: loc_output, inputphys

    ! select the location of the output folder depending on the input choice
    IF (trim(inp%location_output) .eq. 'output') THEN
       loc_output = '../output/results'
    ELSE IF (trim(inp%location_output) .eq. 'forcoms') THEN
       loc_output = '/data/forcoms/results'
    ELSE IF (trim(inp%location_output) .ne. 'output' .and. &
         trim(inp%location_output) .ne. 'forcoms' ) THEN
       loc_output = trim(inp%location_output)
    END IF

    ! create the location if not already created
    CALL system("mkdir "//trim(loc_output)//" 2> /dev/null")

    ! specify the name of the output folder
    IF (trim(inp%name_output) .eq. 'date') THEN
       di%direct_mod = trim(loc_output)//'/'//di%date//"_"//di%realtime2
    ELSE
       di%direct_mod = trim(loc_output)//'/'//trim(inp%name_output)
    END IF

    ! create the output folder
    CALL system("mkdir "//trim(di%direct_mod)//" 2> /dev/null")

  END SUBROUTINE OUTPUT_DIRECTORY


  !--------------------------!
  !--------------------------!
  ! OPEN AND WRITE LOG FILES !
  !--------------------------!
  !--------------------------!


  SUBROUTINE GRIDLOG_PARAM

    !---------------------------------------------------------------------
    ! open the gridlog file and write the grid of parameters
    !---------------------------------------------------------------------

    USE SHARED_VARIABLES, only: di, gi

    CHARACTER(len=15), DIMENSION(14) :: name_param
    INTEGER(KIND=LONG) :: file_gridlog

    ! specify the name of free parameters
    name_param(1) = '            nH'; name_param(2) = '         Tg=Td'
    name_param(3) = '      UVfluxCR'; name_param(4) = '          zeta'
    name_param(5) = '            Av'; name_param(6) = '            ad'
    name_param(7) = ' ice_behaviour'; name_param(8) = '         Ed_Eb'
    name_param(9) = '         Eb(H)'; name_param(10) = '           Fin'
    name_param(11) = '            ds'; name_param(12) = '            Ea'
    name_param(13) = '       X(O)ini'; name_param(14) = '        H2_opr'

    ! open gridlog.out file
    di%direct_gridlog=trim(di%direct_mod)//"/gridlog.out"
    file_gridlog = 495
    OPEN(file_gridlog,file=TRIM(di%direct_gridlog), status='REPLACE',&
         access='SEQUENTIAL',form='FORMATTED',action='WRITE')
       
    ! write the extremal values of the free parameters
    WRITE(file_gridlog,*) '--------------------------------------------',&
         '----------------------'
    WRITE(file_gridlog,'(I2,1x,A11)') gi%Nparam,' parameters'
    WRITE(file_gridlog,'(A53)') &
         '     Parameter    Number   Lower value   Higher value'
    IF (gi%NnH .ne. 0) THEN
       WRITE(file_gridlog,'(A15,7x,I2,7x,ES7.1,8x,ES7.1)') &
            name_param(1),gi%NnH,gi%nHlow,gi%nHup
    END IF
    IF (gi%NT .ne. 0) THEN
       WRITE(file_gridlog,'(A15,7x,I2,7x,ES7.1,8x,ES7.1)') &
            name_param(2),gi%NT,gi%Tlow,gi%Tup
    END IF
    IF (gi%NUV .ne. 0) THEN
       WRITE(file_gridlog,'(A15,7x,I2,7x,ES7.1,8x,ES7.1)') &
            name_param(3),gi%NUV,gi%UVlow,gi%UVup
    END IF
    IF (gi%Nzeta .ne. 0) THEN
       WRITE(file_gridlog,'(A15,7x,I2,7x,ES7.1,8x,ES7.1)') &
            name_param(4),gi%Nzeta,gi%zetalow,gi%zetaup
    END IF
    IF (gi%NAv .ne. 0) THEN
       WRITE(file_gridlog,'(A15,7x,I2,7x,ES7.1,8x,ES7.1)') &
            name_param(5),gi%NAv,gi%Avlow,gi%Avup
    END IF
    IF (gi%Nacst .ne. 0) THEN
       WRITE(file_gridlog,'(A15,7x,I2,7x,ES7.1,8x,ES7.1)') &
            name_param(6),gi%Nacst,gi%acstlow,gi%acstup
    END IF
    IF (gi%choicetreat .ne. 0) THEN
       WRITE(file_gridlog,'(A15,7x,I2,13x,I1)') &
            name_param(7),gi%Ntreat,gi%choicetreat
    END IF
    IF (gi%Nenratio .ne. 0) THEN
       WRITE(file_gridlog,'(A15,7x,I2,7x,ES7.1,8x,ES7.1)') &
            name_param(8),gi%Nenratio,gi%enratiolow,gi%enratioup
    END IF
    IF (gi%NEbH .ne. 0) THEN
       WRITE(file_gridlog,'(A15,7x,I2,7x,ES7.1,8x,ES7.1)') &
            name_param(9),gi%NEbH,gi%EbHlow,gi%EbHup
    END IF
    IF (gi%Nporo .ne. 0) THEN
       WRITE(file_gridlog,'(A15,7x,I2,7x,ES7.1,8x,ES7.1)') &
            name_param(10),gi%Nporo,gi%porolow,gi%poroup
    END IF
    IF (gi%Nds .ne. 0) THEN
       WRITE(file_gridlog,'(A15,7x,I2,7x,ES7.1,8x,ES7.1)') &
            name_param(11),gi%Nds,gi%dslow,gi%dsup
    END IF
    IF (gi%NEa .ne. 0) THEN
       WRITE(file_gridlog,'(A15,7x,I2,7x,ES7.1,8x,ES7.1)') &
            name_param(12),gi%NEa,gi%Ealow,gi%Eaup
    END IF
    IF (gi%NXO .ne. 0) THEN
       WRITE(file_gridlog,'(A15,7x,I2,7x,ES7.1,8x,ES7.1)') &
            name_param(13),gi%NXO,gi%XOlow,gi%XOup
    END IF
    IF (gi%Nop .ne. 0) THEN
       WRITE(file_gridlog,'(A15,7x,I2,7x,ES7.1,8x,ES7.1)') &
            name_param(14),gi%Nop, gi%oplow, gi%opup
    END IF
    
    ! write the number of models
    WRITE(file_gridlog,'(I6,A7)') &
         gi%Nmodels, ' models'
    WRITE(file_gridlog,*) '--------------------------------------------',&
         '----------------------'
    WRITE(file_gridlog,*) 'Number    Folder name'


  END SUBROUTINE GRIDLOG_PARAM


  SUBROUTINE SPATLOG_PARAM

    !---------------------------------------------------------------------
    ! open the spatlog file and write the set of streamline files
    !---------------------------------------------------------------------

    USE SHARED_VARIABLES, only: di, spat

    INTEGER(KIND=LONG) :: file_spatlog, i
    CHARACTER(len=200) :: direct_spatlog

    ! open gridlog.out file
    direct_spatlog=trim(di%direct_mod)//"/spatlog.out"
    file_spatlog = 495
    OPEN(file_spatlog,file=TRIM(direct_spatlog), status='REPLACE',&
         access='SEQUENTIAL',form='FORMATTED',action='WRITE')
       
    ! write the extremal values of the free parameters
    WRITE(file_spatlog,*) '--------------------------------------------',&
         '----------------------'
    WRITE(file_spatlog,'(A23,I4)') 'Number of streamlines: ',spat%Nspat
    WRITE(file_spatlog,*) '--------------------------------------------',&
         '----------------------'
    DO i=1,spat%Nspat
       WRITE(file_spatlog,*) spat%file(i)
    END DO

    CLOSE(file_spatlog)

  END SUBROUTINE SPATLOG_PARAM


  SUBROUTINE GRIDLOG_DIRECT(im)

    !---------------------------------------------------------------------
    ! write the name of each model of the grid
    !---------------------------------------------------------------------

    USE SHARED_VARIABLES, only: di, gi, gp

    INTEGER(KIND=LONG), INTENT(in) :: im
    INTEGER(KIND=LONG) :: len_mod, file_gridlog
    CHARACTER(len=200) :: form_mod, char_mod
    CHARACTER(len=200) :: direct_mod2

    gp%direct(im) = ' '

    IF (gi%NnH .ne. 0) gp%direct(im) = &
         trim(gp%direct(im))//trim(gi%char_nH)//'_'
    IF (gi%NT .ne. 0) gp%direct(im) = &
         trim(gp%direct(im))//trim(gi%char_T)//'_'
    IF (gi%NUV .ne. 0) gp%direct(im) = &
         trim(gp%direct(im))//trim(gi%char_UV)//'_'
    IF (gi%Nzeta .ne. 0) gp%direct(im) = &
         trim(gp%direct(im))//trim(gi%char_zeta)//'_'
    IF (gi%NAv .ne. 0) gp%direct(im) = &
         trim(gp%direct(im))//trim(gi%char_Av)//'_'
    IF (gi%Nacst .ne. 0) gp%direct(im) = &
         trim(gp%direct(im))//trim(gi%char_dist)//'_'
    IF (gi%choicetreat .ne. 0) gp%direct(im) = &
         trim(gp%direct(im))//trim(gi%char_treat)//'_'
    IF (gi%Nenratio .ne. 0) gp%direct(im) = &
         trim(gp%direct(im))//trim(gi%char_en)//'_'
    IF (gi%NEbH .ne. 0) gp%direct(im) = &
         trim(gp%direct(im))//trim(gi%char_EbH)//'_'
    IF (gi%Nporo .ne. 0) gp%direct(im) = &
         trim(gp%direct(im))//trim(gi%char_poro)//'_'
    IF (gi%Nds .ne. 0) gp%direct(im) = &
         trim(gp%direct(im))//trim(gi%char_ds)//'_'
    IF (gi%NEa .ne. 0) gp%direct(im) = &
         trim(gp%direct(im))//trim(gi%char_Ea)//'_'
    IF (gi%NXO .ne. 0) gp%direct(im) = &
         trim(gp%direct(im))//trim(gi%char_XO)//'_'
    IF (gi%Nop .ne. 0) gp%direct(im) = &
         trim(gp%direct(im))//trim(gi%char_op)//'_'

    direct_mod2=gp%direct(im)
    len_mod=len(trim(direct_mod2))
    gp%direct(im) = direct_mod2(1:len_mod-1)

    CALL system("mkdir "//trim(di%direct_mod)//'/'//trim(gp%direct(im))//" 2> /dev/null")

    ! open gridlog.out file
    !di%direct_gridlog=trim(di%direct_mod)//"/gridlog.out"
    file_gridlog = 495
    OPEN(file_gridlog,file=TRIM(di%direct_gridlog), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='WRITE')
       
    WRITE(char_mod,'(I3)') len_mod
    form_mod = '(I6,3x,A'//char_mod
    WRITE(file_gridlog,form_mod//')') im, gp%direct(im)
    
    IF(im .eq. gi%Nmodels) CLOSE(file_gridlog)

  END SUBROUTINE GRIDLOG_DIRECT

  SUBROUTINE OUTPUT_RATES(it,il,im)

    !---------------------------------------------------------------------
    ! write the contribution of formation and destruction of each species
    !---------------------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, sp, re, gr, ph, ab, di, gp, out, spat
    
    CHARACTER(len=52) :: fmt_reac1
    CHARACTER(len=35) :: fmt_reac2
    CHARACTER(len=71) :: fmt_reac_output
    CHARACTER(len=300) :: direct_rates
    INTEGER(KIND=LONG) :: log, ic, ie, is, ic2
    INTEGER(KIND=LONG), INTENT(in) :: it, il, im

    fmt_reac1 = "(I5,1x,A10,A3,A10,A3,A10,A4,A10,A3,A10,A3,A10,A3,A10"
    fmt_reac2 = ",1x,ES12.3,ES12.3)"
    fmt_reac_output = fmt_reac1//fmt_reac2

    ! open rates.out file depending on the folder location
    IF (inp%chgrid .eq. 1) THEN
       direct_rates=trim(di%direct_mod)//"/log_rates.out"
    ELSE IF (inp%chgrid .eq. 2) THEN
       direct_rates=trim(di%direct_mod)//'/'//trim(gp%direct(im))//"/log_rates.out"
    ELSE IF (inp%chgrid .eq. 3) THEN
       direct_rates=trim(di%direct_mod)//'/'//trim(spat%file(im))//"/log_rates.out"
    END IF

    log = 357
    OPEN(log,file=TRIM(direct_rates), status='NEW',&
         access='SEQUENTIAL',form='FORMATTED',action='WRITE')

    WRITE(log,*) "----------------------------------------",&
         "------------------------------------------------"
    WRITE(log,*) "List of reactions which contribute to the formation",&
         " and the destruction of each species"
    WRITE(log,*) "----------------------------------------",&
         "------------------------------------------------"

    DO is=1,inp%Nsprate
       WRITE(log,*) 
       WRITE(log,*) sp%rate_name(is)
       WRITE(log,*) 'FORMATION:'
       WRITE(log,*) " Num React1     + React2     + React3     -> Prod1      ",&
            "+ Prod2      + Prod3      + Prod4       Contrib. (%) Total rate" 
       DO ic=1,sp%Nreacform(is) 
          ic2 = out%numform2(is,ic,1)
          WRITE(log,fmt_reac_output) &
               ic2, re%react1(ic2),' + ', re%react2(ic2),' + ',&
               re%react3(ic2),' -> ', re%prod1(ic2),' + ',re%prod2(ic2),' + ',&
               re%prod3(ic2),' + ',re%prod4(ic2),&
               out%contribform(is,ic)*1d2,out%rateformtot(is,ic)
       END DO
       WRITE(log,*) 'DESTRUCTION:'
       WRITE(log,*) " Num React1     + React2     + React3     -> Prod1      ",&
            "+ Prod2      + Prod3      + Prod4       Contrib. (%) Total rate" 
       DO ic=1,sp%Nreacdest(is) 
          ic2 = out%numdest2(is,ic,1)
          WRITE(log,fmt_reac_output) &
               ic2, re%react1(ic2),' + ', re%react2(ic2),' + ',&
               re%react3(ic2),' -> ', re%prod1(ic2),' + ',re%prod2(ic2),' + ',&
               re%prod3(ic2),' + ',re%prod4(ic2),&
               out%contribdest(is,ic)*1d2,out%ratedesttot(is,ic)
       END DO

       WRITE(log,*) "---------------------------------------------------------",&
            "-------------------------------------------------------------"
    END DO


  END SUBROUTINE OUTPUT_RATES

  SUBROUTINE OUTPUT_LOG(it,il,im)

    !----------------------------------------------------------
    ! write the physical parameters of the model in a log file
    !----------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, sp, re, gr, ph, ab, di, gp, out, spat

    CHARACTER(len=34) :: fmt_reac1
    CHARACTER(len=43) :: fmt_reac2
    CHARACTER(len=78) :: fmt_reac_output
    CHARACTER(len=300) :: direct_log
    INTEGER(KIND=LONG) :: log, ic, ie, is
    INTEGER(KIND=LONG), INTENT(in) :: it, il, im

    fmt_reac1 = "(I5,1x,A15,A15,A15,A15,A15,A15,A15"
    fmt_reac2 = ",1x,ES9.2,1x,ES9.2,1x,ES9.2,1x,ES9.2,1x,I2)"
    fmt_reac_output = fmt_reac1//fmt_reac2

    ! open log.out file depending on the folder location
    IF (inp%chgrid .eq. 1) THEN
       direct_log=trim(di%direct_mod)//"/log.out"
    ELSE IF (inp%chgrid .eq. 2) THEN
       direct_log=trim(di%direct_mod)//'/'//trim(gp%direct(im))//"/log.out"
    ELSE IF (inp%chgrid .eq. 3) THEN
       direct_log=trim(di%direct_mod)//'/'//trim(spat%file(im))//"/log.out"
    END IF

    log = 258
    OPEN(log,file=TRIM(direct_log), status='NEW',&
         access='SEQUENTIAL',form='FORMATTED',action='WRITE')

    ! write input parameters
    WRITE(log,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(log,*) "NUMERICAL PARAMETERS:"
    WRITE(log,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(log,'(ES9.3,5x,A41)') inp%tmax, "! time : Time of Prestellar phase (years)"
    WRITE(log,'(I6,8x,A29)') out%Nrealsteps, "! Nstep: Number of time steps"
    WRITE(log,'(I6,8x,A33)') il, "! Nlay : Maximal number of layers"

    WRITE(log,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(log,*) "PHYSICAL PARAMETERS:"
    WRITE(log,*) "----------------------------------------",&
         "---------------------------------------"
    IF (inp%chphys .eq. 'no') THEN 
       WRITE(log,*) 'Constant physical conditions: '
       WRITE(log,'(ES9.3,5x,A35)') inp%nHini, "! nH   : Constant density nH (cm-3)"
       WRITE(log,'(ES9.3,5x,A36)') inp%Tgini,   "! Tg   : Constant temperature Tg (K)"
       WRITE(log,'(ES9.3,5x,A36)') inp%Tdini, "! Td   : Constant temperature Td (K)"
       WRITE(log,'(ES9.3,5x,A35)') inp%Avini, "! Av   : Visual extinction Av (mag)"
    ELSE IF (inp%chphys .eq. 'yes') THEN
       WRITE(log,*) 'Evolving physical conditions: '
       IF (inp%chinphys .ne. 'yes') THEN
          WRITE(log,'(ES9.3,5x,A7)') inp%dens_evol, "! Alpha"
       ELSE 
          WRITE(log,'(A37)') 'Density evolution given by input file'
       END IF
       WRITE(log,'(ES9.3,5x,A38)') ph%Td,   "! Tini : Initial temperature Tg=Td (K)"
       IF (inp%chinphys .ne. 'yes') THEN
          WRITE(log,'(ES9.3,5x,A6)') inp%temp_evol, "! Beta"
       ELSE 
          WRITE(log,'(A37)') 'Temperature evolution given by input file'
       END IF
       WRITE(log,'(ES9.3,5x,A35)') ph%Av, &
            "! Avini: Initial visual extinction Av (mag)"
    END IF
    WRITE(log,'(ES9.3,5x,A41)') inp%zetaini, &
         "! zeta : Cosmic ray ionization rate (s-1)"
    WRITE(log,'(ES9.3,5x,A36)') inp%G0ext, "! G0   : External radiation field G0"
    WRITE(log,'(ES9.3,5x,A38)') inp%UVfluxCRini, &
         "! CRUVf: CR-induced UV flux (cm-2.s-1)"

    WRITE(log,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(log,*) "GRAIN AND ICE PROPERTIES:"
    WRITE(log,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(log,'(ES9.3,5x,A28)') inp%Rdg,"! Rdg  : Dust/Gas mass ratio"
    WRITE(log,'(ES9.3,5x,A24)') inp%acst, "! ad   : Grain size (um)"
    WRITE(log,'(ES9.3,5x,A37)') gr%X,"! Xd   : Grain abundance wrt H nuclei"
    WRITE(log,'(ES9.3,5x,A33)') inp%ds, "! ds   : Size of a grain site (A)"
    WRITE(log,'(ES9.3,5x,A32)') gr%Nstot(1), "! Nsite: Initial number of sites"
    WRITE(log,'(ES9.3,5x,A40)') inp%rhod,'! rhod : Volumic mass of grains (g.cm-3)'
    WRITE(log,'(ES9.3,5x,A34)') inp%REdH, '! EdHs : Ed/Eb(H) ratio on surface'
    WRITE(log,'(ES9.3,5x,A39)') inp%REdoth, '! EdOts: Ed/Eb(others) ratio on surface'
    WRITE(log,'(ES9.3,5x,A31)') inp%REdbulkH, '! EdHb : Ed/Eb(H) ratio on bulk'
    WRITE(log,'(ES9.3,5x,A36)') inp%REdbulkoth, '! EdOtb: Ed/Eb(others) ratio on bulk'
    WRITE(log,'(ES9.3,5x,A26)') inp%Fporsur, "! Fpor : Porosity Sin/Stot"
    WRITE(log,'(ES9.3,5x,A38)') inp%Slat, "! Slat : Size of a square pore (pores)"
    WRITE(log,'(ES9.3,5x,A62)') inp%rectwid, &
         "! Abarr: Barrier width (in A) if rectangular barrier is chosen"
           
    WRITE(log,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(log,*) "MODEL SWITCHES:"
    WRITE(log,*) "----------------------------------------",&
         "---------------------------------------"

    WRITE(log,'(I1,13x,A48)') inp%chlayer, &
         '! Phase: Formation of ices: 1=2-phase, 2=3-phase'
    WRITE(log,'(I1,13x,A69)') inp%chreacproba, &
         "! Proba: Reaction proba given by 1=the barrier, 2=competition process"
    WRITE(log,'(I1,13x,A72)') inp%checkart, &
         "! Trans: Transmission proba computed with 1=rect barrier, 2=Eckart model"
    WRITE(log,'(I1,13x,A71)') inp%chDnetwork, &
         "! Meth : Methanol network from: 1=input file, 2=Watanabe08, 3=Caselli02"


    WRITE(log,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(log,*) "CHEMICAL PARAMETERS:"
    WRITE(log,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(log,'(A18,I4,A22,I4)') " SPECIES NETWORK: ",inp%Nspecies,&
         " species     Nsprate: ", inp%Nsprate
    WRITE(log,*) "  N Species   Mass  Eb_carb   Eb_wat    Xini" 
    DO is=1,inp%Nspecies
       WRITE(log,'(I4,A1,A15,F4.0,A1,ES8.2,A1,ES8.2,A1,ES9.3)') &
            sp%num(is), ' ', sp%name(is), &
            sp%mass(is), ' ', sp%Eb_ini(is), ' ',sp%Eb_wat(is), &
            ' ', sp%Xini(is)
    END DO

    WRITE(log,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(log,'(A19,I5,A26,I4)') " CHEMICAL NETWORK: ",inp%Nreactions, &
         " reactions     Nreacrate: ", sp%maxNreac2
    WRITE(log,'(A50,A67,A67,A72)')" Num React1         React2         React3         ",&
                          "Prod1          Prod2          Prod3          Prod4",&
                          "            A         B         C         Ea       Type" 
    DO ic=1,inp%Nreactions
       IF (re%type(ic) .le. 10) re%proba(ic) = 0.
       WRITE(log,fmt_reac_output) &
            ic, re%react1(ic), re%react2(ic), re%react3(ic), &
            re%prod1(ic), re%prod2(ic), re%prod3(ic), re%prod4(ic), &
            re%A(ic), re%B(ic), re%C(ic), re%proba(ic), re%type(ic)
    END DO
    

    ! Results at t = tmax_prest
    WRITE(log,*) "--------------------------------------------"
    WRITE(log,*) 'ABUNDANCES AT THE END OF THE SIMULATION:'
    WRITE(log,'(A7,ES9.3,A6)') ' Time: ', ph%time, ' years'
    WRITE(log,'(A22,I6)') ' Number of timesteps: ', out%Nrealsteps
    WRITE(log,*) "--------------------------------------------"
    DO is=1,inp%Nspgas
       WRITE(log,'(A6,A10,A4,ES12.5)') " Xgas(",TRIM(sp%name(is)), ") = ",&
            ab%Xgas(is)
    END DO
    WRITE(log,*) "--------------------------------------------"
    DO is=1,inp%Nspgr
       WRITE(log,'(A8,A10,A4,ES12.5)') " Xice(",TRIM(sp%name(inp%Nspgas+is)),&
            ") = ",ab%Xice(is)
    END DO
    WRITE(log,*) "--------------------------------------------"
    DO is=1,inp%Nspgas
       WRITE(log,'(A6,A10,A4,ES12.5)') " Xtot(",TRIM(sp%name(is)), ") = ",&
            ab%Xtot(is)
    END DO
    WRITE(log,*) "--------------------------------------------"
    DO ie=1,Nelements-1
       WRITE(log,'(A7,A10,A4,ES12.5)') " Xelem(",TRIM(sp%name_elem(ie)), &
            ") = ", ab%Xelem(ie)

    END DO
    WRITE(log,'(A7,A10,A4,ES11.5,2x,A10,ES12.4)') " Xelem(",&
         TRIM(sp%name_elem(Nelements)), ") = ", ab%Xelem(Nelements), &
         "Xcharge = ",ab%Xcharge(1)
    WRITE(log,*) "--------------------------------------------"
    WRITE(log,'(A33,I4)') ' Number of monolayers on grains: ', il
    DO is=1,inp%Nspgr
       WRITE(log,'(A8,A10,A4,ES12.5)') " Ngrain(",&
            TRIM(sp%name(inp%Nspgas+is)), ") = ",&
            ab%Pice(is)
    END DO
    WRITE(log,'(A8,ES11.5)') " Ntot = ", ab%pop_int
    WRITE(log,*) "--------------------------------------------",&
         "-----------------------------------"
    CLOSE(log)

  END SUBROUTINE OUTPUT_LOG



  !-----------------------------!
  !-----------------------------!
  ! OPEN AND WRITE BINARY FILES !
  !-----------------------------!
  !-----------------------------!

  SUBROUTINE WRITE_BINARY(it,il,im)
    
    !---------------------------------------------------------------------
    ! write the physical properties during the prestellar phase
    !---------------------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, sp, re, gr, ph, ab, di, gp, out, spat
 
    INTEGER(KIND=LONG), INTENT(in) :: it, il, im
    INTEGER(KIND=LONG) :: j
    CHARACTER(LEN=4) :: realtime2
    CHARACTER(LEN=200) :: direct_model
    CHARACTER(len=300) :: direct_rad
    CHARACTER(len=300) :: direct_time
    CHARACTER(len=300) :: direct_Tgas
    CHARACTER(len=300) :: direct_Tdust
    CHARACTER(len=300) :: direct_nH
    CHARACTER(len=300) :: direct_Av
    CHARACTER(len=300) :: direct_Xgas
    CHARACTER(len=300) :: direct_Xice
    CHARACTER(len=300) :: direct_Xsur
    CHARACTER(len=300) :: direct_Nlay
    CHARACTER(len=300) :: direct_poplay
    CHARACTER(len=300) :: direct_grain
    CHARACTER(len=300) :: direct_input
    CHARACTER(len=300) :: direct_reacs1
    CHARACTER(len=300) :: direct_reacs2
    CHARACTER(len=300) :: direct_reacs3
    CHARACTER(len=300) :: direct_reacs4
    CHARACTER(len=300) :: direct_reacs5
    CHARACTER(len=300) :: direct_reacs6
    CHARACTER(len=300) :: direct_reacs7
    INTEGER(KIND=LONG) :: file_rad
    INTEGER(KIND=LONG) :: file_time
    INTEGER(KIND=LONG) :: file_Tgas
    INTEGER(KIND=LONG) :: file_Tdust
    INTEGER(KIND=LONG) :: file_nH
    INTEGER(KIND=LONG) :: file_Av
    INTEGER(KIND=LONG) :: file_Xgas
    INTEGER(KIND=LONG) :: file_Xice
    INTEGER(KIND=LONG) :: file_Xsur
    INTEGER(KIND=LONG) :: file_Nlay
    INTEGER(KIND=LONG) :: file_poplay
    INTEGER(KIND=LONG) :: file_grain
    INTEGER(KIND=LONG) :: file_input
    INTEGER(KIND=LONG) :: file_reacs1
    INTEGER(KIND=LONG) :: file_reacs2
    INTEGER(KIND=LONG) :: file_reacs3
    INTEGER(KIND=LONG) :: file_reacs4
    INTEGER(KIND=LONG) :: file_reacs5
    INTEGER(KIND=LONG) :: file_reacs6
    INTEGER(KIND=LONG) :: file_reacs7

    direct_Xgas = ' '
    direct_Xice = ' '
    direct_Xsur = ' '
    direct_Nlay = ' '
    direct_poplay = ' '
    direct_grain = ' '
    direct_rad = ' '
    direct_time = ' '
    direct_Tgas = ' '
    direct_Tdust = ' '
    direct_nH = ' '
    direct_Av = ' '
    direct_input = ' '
    direct_reacs1 = ' '
    direct_reacs2 = ' '
    direct_reacs3 = ' '
    direct_reacs4 = ' '
    file_Xgas = 0
    file_Xice = 0
    file_Xsur = 0
    file_Nlay = 0
    file_poplay = 0
    file_grain = 0
    file_rad = 0
    file_time = 0
    file_Tgas = 0
    file_Tdust = 0
    file_nH = 0
    file_Av = 0
    file_input = 0
    file_reacs1 = 0
    file_reacs2 = 0
    file_reacs3 = 0
    file_reacs4 = 0
    file_reacs5 = 0
    file_reacs6 = 0
    file_reacs7 = 0

    realtime2 = di%realtime

    ! select the right folder
    IF (inp%chgrid .eq. 1) THEN
       direct_model = trim(di%direct_mod)
    ELSE IF (inp%chgrid .eq. 2) THEN
       direct_model =  trim(di%direct_mod)//'/'//trim(gp%direct(im))
    ELSE IF (inp%chgrid .eq. 3) THEN
       direct_model =  trim(di%direct_mod)//'/'//trim(spat%file(im))
    END IF

    ! specify the name of binary file
    !IF (inp%chgrid .eq. 3) 
    direct_rad=trim(direct_model)//"/radius.out"
    direct_time=trim(direct_model)//"/time.out"
    direct_Tgas=trim(direct_model)//"/Tgas.out"
    direct_Tdust=trim(direct_model)//"/Tdust.out"
    direct_nH=trim(direct_model)//"/nH.out"
    direct_Av=trim(direct_model)//"/Av.out"
    direct_Xgas=trim(direct_model)//"/Xgas.out"
    direct_Xice=trim(direct_model)//"/Xice.out"
    direct_Xsur=trim(direct_model)//"/Xsur.out"
    direct_Nlay=trim(direct_model)//"/Nlayers.out"
    direct_poplay=trim(direct_model)//"/populationlay.out"
    direct_grain=trim(direct_model)//"/grain.out"
    direct_input=trim(direct_model)//"/input.out"
    direct_reacs1=trim(direct_model)//"/reac_numform.out"
    direct_reacs2=trim(direct_model)//"/reac_numdest.out"
    direct_reacs3=trim(direct_model)//"/reac_rateform.out"
    direct_reacs4=trim(direct_model)//"/reac_ratedest.out"
    direct_reacs5=trim(direct_model)//"/reac_contribform.out"
    direct_reacs6=trim(direct_model)//"/reac_contribdest.out"
    direct_reacs7=trim(direct_model)//"/reac_spnum.out"

    file_Xgas = 401
    file_Xice = 402
    file_Xsur = 403
    file_Nlay = 501
    file_poplay = 601
    file_grain = 701
    file_rad = 751
    file_time = 801
    file_Tgas = 802
    file_Tdust = 803
    file_nH = 804
    file_Av = 805
    file_input = 901
    file_reacs1 = 1001
    file_reacs2 = 1002
    file_reacs3 = 1003
    file_reacs4 = 1004
    file_reacs5 = 1005
    file_reacs6 = 1006
    file_reacs7 = 1007

    ! create binary files
    OPEN(file_Xgas,file=TRIM(direct_Xgas),form='UNFORMATTED')
    OPEN(file_Xice,file=TRIM(direct_Xice),form='UNFORMATTED')
    OPEN(file_Xsur,file=TRIM(direct_Xsur),form='UNFORMATTED')
    OPEN(file_Nlay,file=TRIM(direct_Nlay),form='UNFORMATTED')
    !OPEN(file_poplay,file=TRIM(direct_poplay),form='UNFORMATTED')
    OPEN(file_grain,file=TRIM(direct_grain),form='UNFORMATTED')
    !IF (inp%chgrid .eq. 3 .or. (inp%chphys .eq. 'yes' .and. inp%chinphys .eq. 'yes')) &
    OPEN(file_rad,file=TRIM(direct_rad),form='UNFORMATTED')
    OPEN(file_time,file=TRIM(direct_time),form='UNFORMATTED')
    OPEN(file_Tgas,file=TRIM(direct_Tgas),form='UNFORMATTED')
    OPEN(file_Tdust,file=TRIM(direct_Tdust),form='UNFORMATTED')
    OPEN(file_nH,file=TRIM(direct_nH),form='UNFORMATTED')
    OPEN(file_Av,file=TRIM(direct_Av),form='UNFORMATTED')
    OPEN(file_input,file=TRIM(direct_input),form='UNFORMATTED')
    IF (inp%chreacrates .eq. 1) THEN
       OPEN(file_reacs1,file=TRIM(direct_reacs1),form='UNFORMATTED')
       OPEN(file_reacs2,file=TRIM(direct_reacs2),form='UNFORMATTED')
       OPEN(file_reacs3,file=TRIM(direct_reacs3),form='UNFORMATTED')
       OPEN(file_reacs4,file=TRIM(direct_reacs4),form='UNFORMATTED')
       OPEN(file_reacs5,file=TRIM(direct_reacs5),form='UNFORMATTED')
       OPEN(file_reacs6,file=TRIM(direct_reacs6),form='UNFORMATTED')
       OPEN(file_reacs7,file=TRIM(direct_reacs7),form='UNFORMATTED')
    END IF

    ! write in the binary files
    !IF (inp%chgrid .eq. 3 .or. (inp%chphys .eq. 'yes' .and. inp%chinphys .eq. 'yes')) &
    WRITE(file_rad) out%rad(1:out%Nrealsteps), out%v(1:out%Nrealsteps)
    !write(6,*) out%rad(1:out%Nrealsteps)
    WRITE(file_time) out%time(1:out%Nrealsteps), out%Nrealsteps
    WRITE(file_Tgas) out%Tg(1:out%Nrealsteps)
    WRITE(file_Tdust) out%Td(1:out%Nrealsteps)
    WRITE(file_nH) out%nH(1:out%Nrealsteps)
    WRITE(file_Av) out%Av(1:out%Nrealsteps)
    WRITE(file_Xgas) out%Xgas(1:inp%Nspgas,1:out%Nrealsteps)
    WRITE(file_Xice) out%Xice(1:inp%Nspgr,1:out%Nrealsteps)
    WRITE(file_Xsur) out%Xsurml(1:inp%Nspgr,1:out%Nrealsteps)
    WRITE(file_Nlay) out%Nlay(1:out%Nrealsteps)
    WRITE(file_grain) gr%a(1),gr%X,gr%Nstot(1),gr%a(il),gr%X,gr%Nstot(il)
    WRITE(file_input) out%Nrealsteps, sp%name(1:inp%Nspecies)

    IF (inp%chreacrates .eq. 1) THEN
       !write(6,*) inp%Nsprate, sp%maxNreac2, out%Nrealsteps, sp%maxNreac
       WRITE(file_reacs1) out%numform2(1:inp%Nsprate,1:sp%maxNreac2,&
            1:out%Nrealsteps)
       WRITE(file_reacs2) out%numdest2(1:inp%Nsprate,1:sp%maxNreac2,&
            1:out%Nrealsteps)
       WRITE(file_reacs3) out%rateform2(1:inp%Nsprate,1:sp%maxNreac2,&
            1:out%Nrealsteps)
       WRITE(file_reacs4) out%ratedest2(1:inp%Nsprate,1:sp%maxNreac2,&
            1:out%Nrealsteps)
       WRITE(file_reacs5) out%contribform(1:inp%Nsprate,1:sp%maxNreac2)
       WRITE(file_reacs6) out%contribdest(1:inp%Nsprate,1:sp%maxNreac2)
       WRITE(file_reacs7) sp%rate_num(1:inp%Nsprate)
    END IF

    CLOSE(file_Xgas)
    CLOSE(file_Xice)
    CLOSE(file_Xsur)
    CLOSE(file_Nlay)
    CLOSE(file_grain)
    CLOSE(file_rad)
    CLOSE(file_time)
    CLOSE(file_Tgas)
    CLOSE(file_Tdust)
    CLOSE(file_nH)
    CLOSE(file_Av)
    CLOSE(file_input)
    IF (inp%chreacrates .eq. 1) THEN
       CLOSE(file_reacs1)
       CLOSE(file_reacs2)
       CLOSE(file_reacs3)
       CLOSE(file_reacs4)
       CLOSE(file_reacs5)
       CLOSE(file_reacs6)
       CLOSE(file_reacs7)
    END IF

  END SUBROUTINE WRITE_BINARY



  !----------------------------------------------
  !----------------------------------------------
  ! WRITE DATA IN THE TERMINAL
  !----------------------------------------------
  !----------------------------------------------


  SUBROUTINE TERM_BEGINNING

    WRITE(6,*) "-----------------------------------------------",&
         "--------------------------------"
    WRITE(6,*) "                       GRAINOBLE astrochemical model"
    WRITE(6,*) "  gas-grain, multilayer, ",&
         "time-dependent based on the rate equations"
    WRITE(6,*) "-----------------------------------------------",&
         "--------------------------------"
    

  END SUBROUTINE TERM_BEGINNING


  SUBROUTINE TERM_INPUT(im) 

    !-----------------------------------------------
    ! write the input parameters in the terminal
    !-----------------------------------------------
    
    USE SHARED_VARIABLES, only: inp, out, ph, gr, di, gi

    INTEGER(KIND=LONG), INTENT(in) :: im

    IF (inp%chgrid .eq. 2) THEN
       WRITE(6,'(I7,A11,I7)') im, 'th model on ', gi%Nmodels
   ELSE IF (inp%chgrid .eq. 3) THEN
       WRITE(6,'(I7,A11,I7)') im, 'th model on ', gi%Nmodels
    END IF

    ! write input parameters
    WRITE(6,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(6,*) "NUMERICAL PARAMETERS:"
    WRITE(6,'(ES9.3,5x,A34)') inp%tmax, "! Time of Prestellar phase (years)"
    WRITE(6,'(I6,8x,A22)') out%Ninisteps, "! Number of time steps"
    WRITE(6,'(I5,10x,A26)') out%Ninilayer, "! Maximal number of layers"

    WRITE(6,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(6,*) "PHYSICAL PARAMETERS:"
    IF (inp%chphys .eq. 'no') THEN 
       WRITE(6,*) 'Constant physical conditions: '
       WRITE(6,'(ES9.3,5x,A28)') inp%nHini, "! Constant density nH (cm-3)"
       WRITE(6,'(ES9.3,5x,A29)') inp%Tgini,   "! Constant temperature Tg (K)"
       WRITE(6,'(ES9.3,5x,A29)') inp%Tdini, "! Constant temperature Td (K)"
       WRITE(6,'(ES9.3,5x,A28)') inp%Avini, "! Visual extinction Av (mag)"
    ELSE IF (inp%chphys .eq. 'yes') THEN
       WRITE(6,*) 'Evolving physical conditions: '
       WRITE(6,'(ES9.3,5x,A27)') ph%nH, "! Initial density nH (cm-3)"
       IF (inp%chinphys .ne. 'yes') THEN
          WRITE(6,'(ES9.3,5x,A7)') inp%dens_evol, "! Alpha"
       ELSE 
          WRITE(6,'(A37)') 'Density evolution given by input file'
       END IF
       WRITE(6,'(ES9.3,5x,A31)') ph%Td,   "! Initial temperature Tg=Td (K)"
       IF (inp%chinphys .ne. 'yes') THEN
          WRITE(6,'(ES9.3,5x,A6)') inp%temp_evol, "! Beta"
       ELSE 
          WRITE(6,'(A37)') 'Temperature evolution given by input file'
       END IF
       WRITE(6,'(ES9.3,5x,A28)') ph%Av, &
            "! Initial visual extinction Av (mag)"
    END IF
    WRITE(6,'(ES9.3,5x,A34)') inp%zetaini, "! Cosmic ray ionization rate (s-1)"
    WRITE(6,'(ES9.3,5x,A29)') inp%G0ext, "! External radiation field G0"
    WRITE(6,'(ES9.3,5x,A31)') inp%UVfluxCRini,"! CR-induced UV flux (cm-2.s-1)"

    WRITE(6,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(6,*) "GRAIN PROPERTIES:"
    WRITE(6,'(ES9.3,5x,A17)') inp%acst, "! Grain size (um)"
    WRITE(6,'(ES9.3,5x,A30)') gr%X,"! Grain abundance wrt H nuclei"
    WRITE(6,'(ES9.3,5x,A21)') inp%Rdg,"! Dust/Gas mass ratio"
    WRITE(6,'(ES9.3,5x,A33)') inp%rhod,'! Volumic mass of grains (g.cm-3)'
           
    WRITE(6,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(6,*) "ICE PROPERTIES:"
    IF (inp%chmulti .eq. 1) THEN
       WRITE(6,'(A3,11x,A54)') inp%chlayer2, &
            '! Multilayer ice formation (Taquet et al. 2012 method)'
    ELSE IF (inp%chmulti .eq. 2) THEN
       WRITE(6,'(A3,11x,A56)') inp%chlayer2, &
            '! Multilayer ice formation (Hasegawa et al. 1993 method)'
    END IF
    WRITE(6,'(ES9.3,5x,A7)') inp%REdH, '! Ed/Eb'
    WRITE(6,'(ES9.3,5x,A26)') inp%ds, "! Size of a grain site (A)"
    WRITE(6,'(ES9.3,5x,A25)') gr%Nstot(1), "! Initial number of sites"

    WRITE(6,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(6,*) "GRAIN POROSITY:"
    WRITE(6,'(ES9.3,5x,A19)') inp%Fporsur, "! Porosity Sin/Stot"
    WRITE(6,'(ES9.3,5x,A31)') inp%Slat, "! Size of a square pore (pores)"
    WRITE(6,'(ES9.3,5x,A19)') inp%Fporsur, "! Porosity Sin/Stot"

    WRITE(6,*) "----------------------------------------",&
         "---------------------------------------"
    WRITE(6,*) "CHEMICAL PARAMETERS:"
    WRITE(6,'(I1,13x,A62)') inp%chreacproba, &
         "! Reaction proba given by 1=the barrier, 2=competition process"
    WRITE(6,'(I1,13x,A65)') inp%checkart, &
         "! Transmission proba computed with 1=rect barrier, 2=Eckart model"
    WRITE(6,'(ES9.3,5x,A55)') inp%rectwid, &
         "! Barrier width (in A) if rectangular barrier is chosen"
    WRITE(6,'(I1,13x,A64)') inp%chDnetwork, &
         "! Methanol network from: 1=input file, 2=Watanabe08, 3=Caselli02"

    WRITE(6,*) "----------------------------------------",&
         "---------------------------------------"

  END SUBROUTINE TERM_INPUT



  SUBROUTINE TERM_LOG_LAYER(checklay,it,il) 

    !-----------------------------------------------
    ! write the layer properties in the terminal
    !-----------------------------------------------

    USE SHARED_VARIABLES, only: ph, sp, ns, ab

    INTEGER(KIND=LONG), INTENT(in) :: it, il, checklay

    CHARACTER(LEN=8) :: date_lay
    CHARACTER(LEN=10) ::  time_lay

    CALL DATE_AND_TIME(date_lay,time_lay) 

    ! 1st layer
    IF (it == 1) THEN 
       WRITE(6,*) "-----------------------------------",&
            "--------------------------------------------"
       WRITE(6,'(A33,A6,A6,ES9.3,A6)') &
            " Beginning of the calculation at ", time_lay, ", t = ", ph%time, &
            " years"
       WRITE(6,'(A9,ES9.3,A12,ES9.3,A14,ES9.3,A9,ES9.3,A4)') &
            " dtime = ", ph%dtime/Rys, " years, n = ", ph%nH, &
            " cm-3, Tgas = ",ph%Tg, " K, Av = ", ph%Av, " mag"
       IF (ns%H .ne. 0 .and. ns%O .ne. 0 .and. ns%CO .ne. 0) THEN
          WRITE(6,'(A8,ES9.3,A9,ES9.3,A10,ES9.3,A12,ES9.3,A4)') &
               " X(H) = ", sp%Xini(ns%H), &
               ", X(O) = ", sp%Xini(ns%O), &
               ", X(CO) = ", sp%Xini(ns%CO)
       END IF
    ! other layers
    ELSE IF (it .gt. 1) THEN
       WRITE(6,*) "-----------------------------------",&
            "--------------------------------------------"
       IF (checklay .eq. 1) THEN 
          WRITE(6,'(I4,A21,A6,A6,ES9.3,A26,I6)') &
               il, "th monolayer full at ", time_lay, ", t = ", ph%time_prev, &
               " years, timestep number = ", it
       ELSE IF (checklay .eq. 2) THEN 
          WRITE(6,'(I4,A22,A6,A6,ES9.3,A26,I6)') &
               il, "th monolayer empty at ", time_lay, ", t = ", ph%time_prev, &
               " years, timestep number = ", it
       END IF
       WRITE(6,'(A9,ES9.3,A12,ES9.3,A11,ES9.3,A9,ES9.3,A4)') &
            " dtime = ", ph%dtime/Rys*100, " years, n = ", ph%nH, &
            " cm-3, T = ",ph%Tg, " K, Av = ", ph%Av, " mag"
       WRITE(6,'(A12,ES9.3,A13,ES9.3,A13,ES9.3)') &
            " Xelem(H) = ", ab%Xelem(1), ", Xelem(O) = ", ab%Xelem(6), &
            ", Xelem(C) = ", ab%Xelem(4)
    END IF

  END SUBROUTINE TERM_LOG_LAYER


  SUBROUTINE FIN_TERM_LOG(it,il)
    
    !-----------------------------------------------
    ! write the final abundances in the terminal
    !-----------------------------------------------

    USE SHARED_VARIABLES, only: inp, ph, sp, ab, di
    
    INTEGER(KIND=LONG), INTENT(in) :: it, il
  
!!$    WRITE(6,*) 'PARAMETERS AT THE END OF CALCULATION:'
!!$    WRITE(6,'(A7,ES9.3,A6)') ' Time: ', ph%time, ' years'
!!$    WRITE(6,*) 'Number of timesteps: ', it
!!$    WRITE(6,*) "--------------------------------------------"
!!$    DO is=1,inp%Nspgas
!!$       WRITE(6,'(A6,A10,A4,ES12.5)') " Xgas(",TRIM(sp%name(is)), ") = ",&
!!$            ab%Xgas(is)
!!$    END DO
!!$    WRITE(6,*) "--------------------------------------------"
!!$    DO is=1,inp%Nspgr
!!$       WRITE(6,'(A9,A10,A4,ES12.5,2x,ES12.5)') &
!!$            " Xice(",TRIM(sp%name(inp%Nspgas+is)),") = ",ab%Xice(is)
!!$    END DO
!!$    WRITE(6,*) "--------------------------------------------"
!!$    DO is=1,inp%Nspgas
!!$       WRITE(6,'(A6,A10,A4,ES12.5)') &
!!$            " Xtot(",TRIM(sp%name(is)), ") = ",ab%Xtot(is)
!!$    END DO
!!$    WRITE(6,*) "--------------------------------------------"
!!$    DO ie=1,Nelements
!!$       WRITE(6,'(A7,A10,A4,ES12.5)') " Xelem(",TRIM(sp%name_elem(ie)), &
!!$            ") = ", ab%Xelem(ie)
!!$    END DO
!!$    WRITE(6,'(A11,ES11.4)') " Xcharge = '",ab%Xcharge(1)
!!$    WRITE(6,*) "--------------------------------------------"
!!$    WRITE(6,*) 'Number of monolayers on grains: ', il
!!$    DO is=1,inp%Nspgr
!!$       WRITE(6,'(A8,A10,A4,ES12.5,2x,ES12.5)') &
!!$            " Ngrain(",TRIM(sp%name(is+inp%Nspgas)),") = ",ab%Pice(is)
!!$    END DO
!!$    WRITE(6,'(A8,ES11.5)') " Ntot = ", ab%pop_int
!!$    WRITE(6,*) "--------------------------------------------",&
!!$         "-----------------------------------"
!!$    
    CALL DATE_AND_TIME(di%date_fin,di%time_fin) 

    WRITE(6,'(A31,A8,A2,A6)') " Beginning of calculations at: ", di%date_ini, &
         '  ', di%time_ini
    WRITE(6,'(A25,A8,A2,A6)') " End of calculations at: ", di%date_fin, &
         '  ', di%time_fin 
    WRITE(6,*) "--------------------------------------------",&
         "-----------------------------------"
    WRITE(6,*) "--------------------------------------------",&
         "-----------------------------------"

  END SUBROUTINE FIN_TERM_LOG


END MODULE OUTPUT
