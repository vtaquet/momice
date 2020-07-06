MODULE INPUT

  USE VARIABLES
  USE ALLOCATION

  IMPLICIT NONE

  !----------------------------------------------------------------------------
  ! This module contains several subroutines that: 
  ! read the various input files
  !----------------------------------------------------------------------------

CONTAINS 

  SUBROUTINE INPUT_FILE
    
    USE SHARED_VARIABLES, only: inp 

    INTEGER(KIND=LONG) :: inputfile
    CHARACTER(len=50) :: direct_inputfile

    ! open file
    direct_inputfile="../input/input_file.in"
    inputfile = 101
    OPEN(inputfile,file=TRIM(direct_inputfile), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    REWIND(inputfile)

    READ(inputfile,*) 
    READ(inputfile,'(A80)') inp%folder
    inp%folder = TRIM(inp%folder)

    CLOSE(inputfile)
    

  END SUBROUTINE INPUT_FILE


  SUBROUTINE INPUT_PARAMETERS

    !-----------------------------
    ! read input_parameters.in
    !-----------------------------
    
    USE SHARED_VARIABLES, only: inp 

    INTEGER(KIND=LONG) :: inputparameters, i, ir, check_J, j
    CHARACTER(len=84) :: direct_inputparameters
    CHARACTER(len=5) :: Nspeciestotchar, Nspecieschar, Nreactionstotchar
    INTEGER(KIND=LONG) :: file_species, file_reactions
    INTEGER(KIND=LONG) :: file_energies, file_rate
    CHARACTER(len=100) :: direct_species, direct_reactions
    CHARACTER(len=100) :: direct_energies, direct_rate
    CHARACTER(len=7) :: fmt_species, char
    CHARACTER(len=1) :: spname
    REAL(KIND=DP) :: RrhonH, Mpr2, G2, RcmAU, rhoflat, Ryrsec

    ! open file
    direct_inputparameters="../input/"//TRIM(inp%folder)//"/input_parameters.in"
    inputparameters = 101
    OPEN(inputparameters,file=TRIM(direct_inputparameters), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    REWIND(inputparameters)

    !--- DATA OPTIONS ---
    READ(inputparameters,*)
    READ(inputparameters,*)
    READ(inputparameters,*)
    READ(inputparameters,'(8x,I1)') inp%chdata
    READ(inputparameters,'(8x,I1)') inp%chreacrates
    
    !--- INPUT AND OUTPUT FILES ---
    READ(inputparameters,*)
    READ(inputparameters,*)
    READ(inputparameters,*)
    READ(inputparameters,'(A60)') inp%location_output
    READ(inputparameters,'(A60)') inp%name_output
    READ(inputparameters,'(A60)') inp%filesp
    READ(inputparameters,'(A60)') inp%filereac
    READ(inputparameters,'(A60)') inp%fileen
    READ(inputparameters,'(A60)') inp%inputphys

    !--- CONSTANT PHYSICAL CONDITIONS ---
    READ(inputparameters,*)
    READ(inputparameters,*)
    READ(inputparameters,*)
    READ(inputparameters,'(8x,A1)') inp%chphys
    READ(inputparameters,'(8x,ES13.7)') inp%tmax
    inp%tmaxini = inp%tmax
    READ(inputparameters,'(8x,ES13.7)') inp%nHini
    READ(inputparameters,'(8x,ES13.7)') inp%Tgini
    READ(inputparameters,'(8x,ES13.7)') inp%Tdini
    READ(inputparameters,'(8x,ES13.7)') inp%zetaini
    READ(inputparameters,'(8x,ES13.7)') inp%Avini
    READ(inputparameters,'(8x,ES13.7)') inp%G0ext
    READ(inputparameters,'(8x,ES13.7)') inp%UVfluxCRini
    !READ(inputparameters,'(A20)') inp%inputphys
    IF (inp%chphys .eq. '1') THEN 
      inp%chphys = 'no'
    ELSE IF (inp%chphys .eq. '2') THEN
      inp%chphys = 'yes'
    END IF

    !--- EVOLVING PHYSICAL CONDITIONS ---
    !READ(inputparameters,*)
    !READ(inputparameters,*)
    !READ(inputparameters,*)
    !IF (trim(inp%chphys) .eq. 'yes') THEN
    !   READ(inputparameters,'(A3)') inp%chinphys
    !   READ(inputparameters,'(ES13.7)') inp%dens_evol
    !   READ(inputparameters,'(ES13.7)') inp%tff2
    !   READ(inputparameters,'(ES13.7)') inp%nHini
    !   READ(inputparameters,'(ES13.7)') inp%nHfin
    !   READ(inputparameters,'(ES13.7)') inp%temp_evol
    !   READ(inputparameters,'(ES13.7)') inp%Tini
    !   READ(inputparameters,'(ES13.7)') inp%Tfin
    !ELSE 
    !   DO i=1,8
    !      READ(inputparameters,*)
    !   END DO       
    !END IF

    !--- GRAIN AND ICE PROPERTIES ---
    READ(inputparameters,*)
    READ(inputparameters,*)
    READ(inputparameters,*)
    READ(inputparameters,'(8x,ES13.7)') inp%Rdg
    READ(inputparameters,'(8x,ES13.7)') inp%acst
    READ(inputparameters,'(8x,ES13.7)') inp%rhod
    READ(inputparameters,'(8x,ES13.7)') inp%ds
    READ(inputparameters,'(8x,ES13.7)') inp%stick_coef
    READ(inputparameters,'(8x,ES13.7)') inp%REdH
    READ(inputparameters,'(8x,ES13.7)') inp%REdbulkH
    READ(inputparameters,'(8x,ES13.7)') inp%REdoth
    READ(inputparameters,'(8x,ES13.7)') inp%REdbulkoth
    READ(inputparameters,'(8x,ES13.7)') inp%Nsurf
    READ(inputparameters,'(8x,A3)') char
    READ(char,*) inp%Nsteps_lay
    READ(inputparameters,'(8x,ES13.7)') inp%rectwid
    READ(inputparameters,'(8x,ES13.7)') inp%diffwid
    READ(inputparameters,'(8x,ES13.7)') inp%Slat
    READ(inputparameters,'(8x,ES13.7)') inp%Fporsur
    READ(inputparameters,'(8x,ES13.7)') inp%Fvacuum

    !--- MODEL SWITCHES ---
    READ(inputparameters,*)
    READ(inputparameters,*)
    READ(inputparameters,*)
    READ(inputparameters,'(8x,I1)') inp%chlayer
    READ(inputparameters,'(8x,I1)') inp%chmodifrates
    READ(inputparameters,'(8x,I1)') inp%chgrainevol
    READ(inputparameters,'(8x,I1)') inp%chevolen
    READ(inputparameters,'(8x,I1)') inp%chsurfchem
    READ(inputparameters,'(8x,I1)') inp%chbulkchem
    READ(inputparameters,'(8x,I1)') inp%chHtunnel
    READ(inputparameters,'(8x,I1)') inp%chCtunnel
    READ(inputparameters,'(8x,I1)') inp%chH
    READ(inputparameters,'(8x,I1)') inp%chH2op
    READ(inputparameters,'(8x,I1)') inp%chreacproba
    READ(inputparameters,'(8x,I1)') inp%checkart

    !--- ODE SOLVER PARAMETERS ---
    READ(inputparameters,*)
    READ(inputparameters,*) 
    READ(inputparameters,*) 
    READ(inputparameters,'(8x,ES13.7)') inp%ATOL0
    READ(inputparameters,'(8x,ES13.7)') inp%RTOL0
    READ(inputparameters,'(8x,I1)') inp%ITOL0
    READ(inputparameters,'(8x,I1)') inp%ITASK0

    inp%Fnd = 1
    inp%chmulti = 2
    
    ! Porosity parameters
    inp%Fnpsur = 1 - inp%Fporsur
    IF (inp%Slat .gt. 2) THEN
       inp%Fedgesur = 4*inp%Fporsur*(inp%Slat-1)/inp%Slat**2
    ELSE IF (inp%Slat .le. 2 .and. inp%Slat .gt. 0) THEN
       inp%Fedgesur = inp%Fporsur
    ELSE IF (inp%Slat .lt. 0) THEN
       inp%Fedgesur = 0
    END IF

    ! layer thickness = site size
    inp%laywid = inp%ds

    IF (inp%chlayer .eq. 2) inp%chlayer2 = 'yes' 
    IF (inp%chlayer .eq. 1) inp%chlayer2 = 'no '

    CLOSE(inputparameters)

    ! read Nspecies
    direct_species="../input/"//TRIM(inp%folder)//"/"//TRIM(inp%filesp)
    file_species = 107
    OPEN(file_species,file=TRIM(direct_species), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    REWIND(file_species)
    READ(file_species,*)
    READ(file_species,*)
    READ(file_species,*)
    READ(file_species,'(A10,A5)') toto, Nspecieschar
    READ(Nspecieschar, '(I5)' )  inp%Nspecies
    READ(file_species,*)
    !write(6,*) inp%Nspecies
    !stop
    
    fmt_species="(5x,A1)"
    check_J = 1
    DO j=1,inp%Nspecies
       READ(file_species,fmt_species) spname
       ! count the number of gaseous species 
       IF (spname .eq. 'J' .and. check_J .eq. 1) THEN
          inp%Nspgas = j-1
          inp%Nspgr = (inp%Nspecies - inp%Nspgas)
          inp%Nspgr2 = inp%Nspgr
          check_J = 2
       END IF
       ! recompute Nspgr if 3-phase model
       IF (spname .eq. 'Q' .and. check_J .eq. 2) THEN
          inp%Nspgr = inp%Nspgr/2.
          check_J = 3
       END IF
    END DO
    CLOSE(file_species)
    
    ! read Nenergies
    direct_energies="../input/"//TRIM(inp%folder)//"/"//TRIM(inp%fileen)
    file_energies = 165
    OPEN(file_energies,file=TRIM(direct_energies), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    REWIND(file_energies)
    READ(file_energies,*)
    READ(file_energies,*)
    READ(file_energies,*)
    READ(file_energies,'(A10,A5)') toto, Nspeciestotchar
    READ(Nspeciestotchar, '(I5)' )  inp%Nspeciestot
    CLOSE(file_energies)

    ! read Nreactions
     direct_reactions="../input/"//TRIM(inp%folder)//"/"//TRIM(inp%filereac)
    file_reactions = 108
    OPEN(file_reactions,file=TRIM(direct_reactions), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    REWIND(file_reactions)
    READ(file_reactions,*)
    READ(file_reactions,*)
    READ(file_reactions,*)
    READ(file_reactions,'(A12,A6)') toto, Nreactionstotchar
    READ(Nreactionstotchar, '(I5)' )  inp%Nreactions
    CLOSE(file_reactions)

    ! read Nsprate
    IF (inp%chreacrates .eq. 1) THEN
       direct_rate = "../input/"//TRIM(inp%folder)//"/reacrates.in"
       file_rate = 109
       OPEN(file_rate,file=TRIM(direct_rate), status='OLD',&
            access='SEQUENTIAL',form='FORMATTED',action='READ')
       REWIND(file_rate)
       READ(file_rate,*)
       READ(file_rate,'(I2)') inp%Nsprate
       CLOSE(file_rate)
    END IF

    ! change tmax to tff from nH
    rhoflat = inp%nHini * RrhonH*Mpr2
    inp%tmax = min(inp%tmaxini,sqrt(3*Pi/(32*G2*rhoflat))/Ryrsec)
    
    inp%Neq = inp%Nspecies

  END SUBROUTINE INPUT_PARAMETERS
  

  SUBROUTINE READ_TERMINAL

    !-----------------------------
    ! read terminal
    !-----------------------------
    
    USE SHARED_VARIABLES, only: inp, sp, ns

    CHARACTER(len=100) :: arg, direct_param
    CHARACTER(len=10) :: param_sens(100)
    CHARACTER(len=300000) :: toto, format_toto
    INTEGER(KIND=LONG) :: i, i2, j, intfile, io, num_mod, im
    !INTEGER(KIND=LONG), PARAMETER :: n = 
    REAL(KIND=DP) :: array(inp%Nspecies), mean, std

    !CHARACTER(len=84) :: direct_inputparameters
    !CHARACTER(len=5) :: Nspeciestotchar, Nspecieschar, Nreactionstotchar
    !INTEGER(KIND=LONG) :: file_species, file_reactions
    !INTEGER(KIND=LONG) :: file_energies, file_rate
    !CHARACTER(len=100) :: direct_species, direct_reactions
    !CHARACTER(len=100) :: direct_energies, direct_rate
    !CHARACTER(len=7) :: fmt_species, char
    !CHARACTER(len=1) :: spname
    !REAL(KIND=DP) :: RrhonH, Mpr2, G2, RcmAU, rhoflat, Ryrsec

    i = 0 ; i2 = 1

    ALLOCATE(inp%randomEb(inp%Nspecies), inp%randomEa(inp%Nreactions))
    inp%randomEb(1:inp%Nspecies) = 1. ; inp%randomEa(1:inp%Nreactions) = 1.

    DO
      CALL get_command_argument(i, arg)
      IF (LEN_TRIM(arg) == 0) EXIT
      !WRITE (*,*) TRIM(arg)

      ! common parameters
      ! 1st parameter: type of simulation
      IF (i .eq. 1) THEN
        READ(arg, *) inp%chgrid
      ! 2nd parameter: output directory location
      ELSE IF (i .eq. 2 .and. inp%chgrid .gt. 1) THEN
        inp%location_output = TRIM(arg)
      ! 3rd parameter: output directory
      ELSE IF (i .eq. 3 .and. inp%chgrid .gt. 1) THEN
        inp%name_output = TRIM(arg)
        IF (inp%chgrid .eq. 4) READ(arg,*) num_mod 
        num_mod = num_mod+1
      END IF
      ! specific parameters to each simulation type
      ! spatial evolution
      IF (inp%chgrid .eq. 2) THEN
        IF (i .eq. 4) THEN
          inp%inputphys = TRIM(arg)
        END IF
      ! model grid
      ELSE IF (inp%chgrid .eq. 3) THEN
        IF (i .eq. 4) THEN
          READ(arg,*) inp%nHini
        ELSE IF (i .eq. 5) THEN 
          READ(arg,*) inp%Tgini
        ELSE IF (i .eq. 6) THEN 
          READ(arg,*) inp%Tdini
          inp%Tgini = inp%Tdini
        ELSE IF (i .eq. 7) THEN
          READ(arg,*) inp%zetaini
        ELSE IF (i .eq. 8) THEN
          READ(arg,*) inp%Avini
        ELSE IF (i .eq. 9) THEN
          READ(arg,*) inp%G0ext
        ELSE IF (i .eq. 10) THEN
          READ(arg,*) inp%UVfluxCRini
        END IF
      ! sensitivity analysis
      ELSE IF (inp%chgrid .eq. 4) THEN
        IF (i .ge. 5) THEN
          READ(arg,*) param_sens(i2)
          !WRITE(6,*) param_sens(i2)
          i2 = i2+1
        ELSE IF (i .eq. 4) THEN
          READ(arg,*) inp%distrib
        END IF
      END IF

      i = i+1
    END DO

  ! compute random number array
  DO i=1,i2-1
    IF (TRIM(param_sens(i)) .eq. 'Ed') THEN
      IF (TRIM(param_sens(i+1)) .eq. 'all') THEN
        !READ(param_sens(i+2),*) mean ; READ(param_sens(i+3),*) std
        !CALL RANDOM(inp%Nspecies, (inp%distrib), mean, std, array(1:inp%Nspecies))
        direct_param="../input/"//TRIM(inp%folder)//"/Ed.in"
        intfile = 101
        OPEN(intfile,file=TRIM(direct_param), status='OLD',&
            access='SEQUENTIAL',form='FORMATTED',action='READ')
        DO im=1,num_mod
          READ(intfile,*,iostat=io) toto
          IF (io/=0) EXIT
        END DO
        READ(toto,*) inp%REdH ; inp%REdoth = inp%REdH
        !inp%REdH = array(num_mod) ; inp%REdoth = array(num_mod)
        CLOSE(intfile)
      END IF
    ELSE IF (TRIM(param_sens(i)) .eq. 'Eb') THEN
      IF (TRIM(param_sens(i+1)) .eq. 'all') THEN
        !READ(param_sens(i+2),*) mean ; READ(param_sens(i+3),*) std
        !CALL RANDOM(inp%Nspecies, (inp%distrib), 1d0, std, inp%randomEb)
        direct_param="../input/"//TRIM(inp%folder)//"/Eb.in"
        intfile = 102
        OPEN(intfile,file=TRIM(direct_param), status='OLD',&
            access='SEQUENTIAL',form='FORMATTED',action='READ')
        format_toto = '('
        DO i2=1,inp%Nspecies-1
          format_toto = trim(format_toto)//'ES11.5,1x,'
        END DO
        format_toto = trim(format_toto)//'ES11.5)'
        DO im=1,num_mod
          READ(intfile,format_toto,iostat=io) inp%randomEb(1:inp%Nspecies)
          IF (io/=0) EXIT
        END DO
        inp%randomEb(1:inp%Nspgas) = 0.
        CLOSE(intfile)
      END IF
    ELSE IF (TRIM(param_sens(i)) .eq. 'Ea') THEN
      IF (TRIM(param_sens(i+1)) .eq. 'all') THEN
        !READ(param_sens(i+2),*) mean ; READ(param_sens(i+3),*) std
        !CALL RANDOM(inp%Nreactions, (inp%distrib), 1d0, std, inp%randomEa)
        direct_param="../input/"//TRIM(inp%folder)//"/Ea.in"
        intfile = 103
        OPEN(intfile,file=TRIM(direct_param), status='OLD',&
            access='SEQUENTIAL',form='FORMATTED',action='READ')
        format_toto = '('
        DO i2=1,inp%Nreactions-1
          format_toto = trim(format_toto)//'ES11.5,1x,'
        END DO
        format_toto = trim(format_toto)//'ES11.5)'
        DO im=1,num_mod
          READ(intfile,format_toto,iostat=io) inp%randomEa(1:inp%Nreactions)
          IF (io/=0) EXIT
        END DO
        CLOSE(intfile)
      END IF

    END IF

  END DO

  !write(6,*) inp%randomEb(1:inp%Nspecies)
  !write(6,*) inp%location_output
  !write(6,*) inp%name_output
  !write(6,*) inp%chgrid, inp%nHini, inp%Tdini, inp%zetaini, inp%Avini, inp%G0ext, inp%UVfluxCRini
  !write(6,*) inp%REdH, inp%REdoth
  !stop
  
  IF (inp%chgrid .eq. 2) THEN
    inp%chphys = 'yes'
    inp%chinphys = 'yes'
  ELSE IF (inp%chgrid .ne. 2) THEN
    inp%chphys = 'no'
    inp%chinphys = 'no'
  END IF
  inp%chgrid = 1
  !inp%inputphys = TRIM(inp%inputphys)

  !write(6,*) inp%chgrid
  !write(6,*) TRIM(inp%location_output)
  !write(6,*) TRIM(inp%name_output)
  !write(6,*) TRIM(inp%inputphys)
  !stop

  END SUBROUTINE READ_TERMINAL


  SUBROUTINE  RANDOM(n, distrib, meaninp, stdinp, array)

    INTEGER(KIND=LONG), INTENT(in) :: n
    CHARACTER(len=10), INTENT(in) :: distrib ! = 1000: 
    INTEGER(KIND=LONG) :: i, size, values(1:8)
    REAL(KIND=DP), INTENT(inout) :: array(n)
    INTEGER(KIND=LONG), allocatable :: seed(:)
    REAL(KIND=DP), INTENT(in) :: meaninp, stdinp
    REAL(KIND=DP) :: temp, mean, pi = 1.0, sd = 0.5
 
    pi = 4.0*ATAN(1.0)
    ! generate an array of random values
    CALL DATE_AND_TIME(values=values)
    CALL RANDOM_SEED(size=size)
    ALLOCATE(seed(size))
    seed(:) = values(8)
    CALL RANDOM_SEED(put=seed)
    CALL RANDOM_NUMBER(array) ! Uniform distribution
    DEALLOCATE(seed)
 
    ! Now convert to normal distribution
    IF (TRIM(distrib) .eq. 'normal') THEN
      DO i = 1, n-1, 2
        temp = stdinp * SQRT(-2.0*LOG(array(i))) * COS(2*pi*array(i+1)) + meaninp
        array(i+1) = stdinp * SQRT(-2.0*LOG(array(i))) * SIN(2*pi*array(i+1)) + meaninp
        array(i) = temp
      END DO
    ELSE IF (TRIM(distrib) .eq. 'uniform') THEN
      array = (array-0.5)*2*stdinp+meaninp
    END IF

    ! Check mean and standard deviation
    mean = SUM(array)/n
    sd = SQRT(SUM((array - mean)**2)/n)
 
    WRITE(*, "(A,F8.6)") "Mean = ", mean
    WRITE(*, "(A,F8.6)") "Standard Deviation = ", sd
    !write(6,*) array!(1), meaninp, stdinp!, MIN(array), MAX(array)!:10)!, n2
    
  END SUBROUTINE RANDOM


  SUBROUTINE INPUT_SPECIES

    !---------------------
    ! read species.in
    !---------------------

    USE SHARED_VARIABLES, only: inp, sp, ns
    
    INTEGER(KIND=LONG) :: file_species, check_J, check_Q
    CHARACTER(len=100) :: direct_species
    INTEGER(KIND=LONG) :: is
    CHARACTER(len=45) :: fmt_species1
    CHARACTER(len=61) :: fmt_species2
    CHARACTER(len=58) :: fmt_species3
    CHARACTER(len=48+61+58) :: fmt_species
    CHARACTER(len=5) :: Nspecieschar
    CHARACTER(len=10) :: toto
    CHARACTER(len=1) :: letter_Q, letter_J
    REAL(KIND=DP) :: Xcat, Xan

    WRITE(6,*) "Reading the species file..."


    check_J = 1
    check_Q = 1
    ns%grinput = 0
    Xcat = 0d0; Xan = 0d0

    ! open file
    sp%name(1:inp%Nspecies) = ' '
    direct_species="../input/"//TRIM(inp%folder)//"/"//TRIM(inp%filesp)
    file_species = 107
    OPEN(file_species,file=TRIM(direct_species), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    ! format
    fmt_species1="(I4,1x,A15,F2.0,2x,F2.0,1x,F2.0,1x,F2.0,1x,"
    fmt_species2="F2.0,1x,F2.0,1x,F2.0,1x,F2.0,1x,F2.0,1x,F2.0,1x,F2.0,1x,F2.0,"
    fmt_species3="1x,F2.0,1x,F2.0,1x,F2.0,2x,ES9.3)"
    fmt_species=fmt_species1//fmt_species2//fmt_species3
    
    REWIND(file_species)
    
    READ(file_species,*)
    READ(file_species,*)
    READ(file_species,*)
    READ(file_species,*)

    ! read the name of elements
    READ(file_species,'(A25,A3,A3,A3,A3,A3,A3,A3,A3,A3,A3,A3,A3,A3,A3)') &
         toto, &
         sp%name_elem(1), sp%name_elem(2), sp%name_elem(3), &
         sp%name_elem(4), sp%name_elem(5), sp%name_elem(6), &
         sp%name_elem(7), sp%name_elem(8), sp%name_elem(9), &
         sp%name_elem(10), sp%name_elem(11), sp%name_elem(12), &
         sp%name_elem(13), sp%name_elem(14)
    
    sp%element(1:Nelements,inp%Nspecies+1) = 0d0
    sp%mass(inp%Nspecies+1) = 0d0
    sp%charge(inp%Nspecies+1) = 0d0

    ! read the number of elements and initial abundance of each species
    DO is=1,inp%Nspecies
       
       READ(file_species,fmt_species) sp%num(is), sp%name(is), &
            sp%charge(is), sp%element(1,is), &
            sp%element(2,is), sp%element(3,is), sp%element(4,is), &
            sp%element(5,is), sp%element(6,is), sp%element(7,is), &
            sp%element(8,is), sp%element(9,is), sp%element(10,is), &
            sp%element(11,is), sp%element(12,is), sp%element(13,is),&
            sp%element(14,is), sp%Xini(is)

       sp%mass(is) = sp%element(1,is)*1d0 + &
            sp%element(2,is)*2d0 + &
            sp%element(3,is)*4d0 + &
            sp%element(4,is)*12d0 + &
            sp%element(5,is)*14d0 + &
            sp%element(6,is)*16d0 + &
            sp%element(7,is)*28d0 + &
            sp%element(8,is)*32d0 + &
            sp%element(9,is)*55.6d0 + &
            sp%element(10,is)*23d0 + &
            sp%element(11,is)*24.3d0 + &
            sp%element(12,is)*35.4d0 + &
            sp%element(13,is)*31d0 + &
            sp%element(14,is)*19d0
       
       sp%num(is) = is


       IF (sp%charge(is) .gt. 0) THEN
          Xcat = Xcat + sp%charge(is)*sp%Xini(is)
       ELSE IF (sp%charge(is) .lt. 0) THEN
          Xan = Xan + sp%charge(is)*sp%Xini(is)
       END IF
       
       WRITE(letter_J,'(A1)') sp%name(is)
       WRITE(letter_Q,'(A1)') sp%name(is)
              
       ! save the number of specific species 
       IF (TRIM(sp%name(is)) .eq. 'H') ns%H = is
       IF (TRIM(sp%name(is)) .eq. 'He') ns%He = is
       IF (TRIM(sp%name(is)) .eq. 'D') ns%D = is
       IF (TRIM(sp%name(is)) .eq. 'H2') ns%H2 = is
       IF (TRIM(sp%name(is)) .eq. 'H2') ns%pH2 = is
       IF (TRIM(sp%name(is)) .eq. 'oH2') ns%oH2 = is
       IF (TRIM(sp%name(is)) .eq. 'pH2') ns%H2 = is
       IF (TRIM(sp%name(is)) .eq. 'pH2') ns%pH2 = is
       IF (TRIM(sp%name(is)) .eq. 'CO') ns%CO = is
       IF (TRIM(sp%name(is)) .eq. 'O') ns%O = is
       IF (TRIM(sp%name(is)) .eq. 'C') ns%C = is
       IF (TRIM(sp%name(is)) .eq. 'N') ns%N = is
       IF (TRIM(sp%name(is)) .eq. 'O2') ns%O2 = is
       IF (TRIM(sp%name(is)) .eq. 'CO2') ns%CO2 = is
       IF (TRIM(sp%name(is)) .eq. 'H2O') ns%H2O = is
       IF (TRIM(sp%name(is)) .eq. 'CH3OH') ns%CH3OH = is
       IF (TRIM(sp%name(is)) .eq. 'JO') ns%JO = is
       IF (TRIM(sp%name(is)) .eq. 'JC') ns%JC = is
       IF (TRIM(sp%name(is)) .eq. 'JN') ns%JN = is
       IF (TRIM(sp%name(is)) .eq. 'JCO') ns%JCO = is
       IF (TRIM(sp%name(is)) .eq. 'JCO2') ns%JCO2 = is
       IF (TRIM(sp%name(is)) .eq. 'JCH3OH') ns%JCH3OH = is
       IF (TRIM(sp%name(is)) .eq. 'JCH4O') ns%JCH3OH = is
       IF (TRIM(sp%name(is)) .eq. 'JOH') ns%JOH = is
       IF (TRIM(sp%name(is)) .eq. 'JOD') ns%JOD = is
       IF (TRIM(sp%name(is)) .eq. 'H2D+') ns%H2DPLUS = is
       IF (TRIM(sp%name(is)) .eq. 'HD') ns%HD = is
       IF (TRIM(sp%name(is)) .eq. 'D2') ns%D2 = is
       IF (TRIM(sp%name(is)) .eq. 'oD2') ns%oD2 = is
       IF (TRIM(sp%name(is)) .eq. 'pD2') ns%D2 = is
       IF (TRIM(sp%name(is)) .eq. 'pD2') ns%pD2 = is
       IF (TRIM(sp%name(is)) .eq. 'E') ns%E = is
       IF (TRIM(sp%name(is)) .eq. 'H3+') ns%H3PL = is
       IF (TRIM(sp%name(is)) .eq. 'CO') ns%CO = is
       IF (TRIM(sp%name(is)) .eq. 'G0') ns%G0 = is
       IF (TRIM(sp%name(is)) .eq. 'G0') ns%grinput = 1
       IF (TRIM(sp%name(is)) .eq. 'G+') ns%GP = is
       IF (TRIM(sp%name(is)) .eq. 'G-') ns%GM = is
       IF (TRIM(sp%name(is)) .eq. 'JH') ns%JH = is
       IF (TRIM(sp%name(is)) .eq. 'JD') ns%JD = is
       IF (TRIM(sp%name(is)) .eq. 'JH2') ns%JH2 = is
       IF (TRIM(sp%name(is)) .eq. 'JH2') ns%JpH2 = is
       IF (TRIM(sp%name(is)) .eq. 'JoH2') ns%JoH2 = is
       IF (TRIM(sp%name(is)) .eq. 'JpH2') ns%JpH2 = is
       IF (TRIM(sp%name(is)) .eq. 'JpH2') ns%JH2 = is
       IF (TRIM(sp%name(is)) .eq. 'JHD') ns%JHD = is
       IF (TRIM(sp%name(is)) .eq. 'JD2') ns%JD2 = is
       IF (TRIM(sp%name(is)) .eq. 'JoD2') ns%JoD2 = is
       IF (TRIM(sp%name(is)) .eq. 'JpD2') ns%JpD2 = is
       IF (TRIM(sp%name(is)) .eq. 'JpD2') ns%JD2 = is
       IF (TRIM(sp%name(is)) .eq. 'JH2O') ns%JH2O = is
       IF (TRIM(sp%name(is)) .eq. 'JH2O') ns%JH2O = is
       IF (TRIM(sp%name(is)) .eq. 'JCN') ns%JCN = is
       
    END DO
    
    IF (inp%chH2op .eq. 1 .and. ns%oH2 .ne. 0) THEN
       inp%H2opr = sp%Xini(ns%oH2)/sp%Xini(ns%pH2)
    END IF

    ! specific case in gas phase chemistry only
    IF (inp%Nspgr .eq. 0) THEN
       inp%Nspgas = inp%Nspecies
       inp%Nspgr = 0
    END IF

    IF (Xcat .ne. Xan) THEN 
       IF ((Xcat+Xan) .gt. 0) sp%Xini(ns%E) = Xcat+Xan
    END IF

    !WRITE(6,*) 'Nspecies=',inp%Nspecies,', Nspgas=',inp%Nspgas,', Nspgr=',inp%Nspgr

    
  END SUBROUTINE INPUT_SPECIES


  SUBROUTINE INPUT_ENERGIES

    !---------------------
    ! read energies.in
    !---------------------

    !IMPLICIT NONE

    USE SHARED_VARIABLES, only: inp, sp, ns
    
    INTEGER(KIND=LONG) :: file_energies, check_J
    CHARACTER(len=100) :: direct_energies
    INTEGER(KIND=LONG) :: is, j
    CHARACTER(len=173) :: fmt_energies 
    CHARACTER(len=43) :: fmt_energies1
    CHARACTER(len=61) :: fmt_energies2
    CHARACTER(len=58) :: fmt_energies3
    CHARACTER(len=1) :: letter_J
    CHARACTER(len=5) :: Nspeciestotchar
    CHARACTER(len=10) :: toto
    REAL(KIND=DP) :: spec_Eb_carb2, spec_Eb_wat2, spec_Eb_pure2
    REAL(KIND=DP) :: spec_cr_coef2, spec_Eb_HH932, spec_Eb_H22
    INTEGER(KIND=LONG) :: spec_num2, check_spec
    CHARACTER(len=23) :: spec_name2
            
    check_J = 1
    letter_J = ' '

    ! open file
    direct_energies="../input/"//TRIM(inp%folder)//"/"//TRIM(inp%fileen)
    file_energies = 165
    OPEN(file_energies,file=TRIM(direct_energies), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    ! format
    fmt_energies1="(I4,1x,A14,1x,F5.0,1x,F5.0,1x,F5.0,1x,F5.0)"
    fmt_energies=fmt_energies1
    
    REWIND(file_energies)
    
    READ(file_energies,*)
    READ(file_energies,*)
    READ(file_energies,*)
    READ(file_energies,*)
    READ(file_energies,*)

!!$    DO j=inp%Nspgas+1,inp%Nspecies
!!$       spec_Eb_wat2 = 2d3
!!$       spec_Eb_carb2 = 2d3
!!$       spec_Eb_H22 = 2d3
!!$       spec_Eb_pure2 = 2d3
!!$       sp%Eb_wat(j) = spec_Eb_wat2
!!$       sp%Eb_carb(j) = spec_Eb_carb2
!!$       sp%Eb_HH93(j) = spec_Eb_HH932
!!$       sp%Eb_H2(j) = spec_Eb_H22
!!$       sp%Eb_H2(j) = sp%Eb_wat(j)/sp%Eb_wat(ns%H)*sp%Eb_H2(ns%H)
!!$       sp%Eb_pure(j) = spec_Eb_pure2
!!$       sp%Eb_pure(j) = sp%Eb_wat(j)/1.15d3*8.5d2
!!$       sp%Eb(j) = sp%Eb_wat(j) ! initially, Eb wrt water is assumed
!!$       sp%Eb_ini(j) = sp%Eb(j) ! initially, Eb wrt water is assumed
!!$       sp%cr_coef(j) = spec_cr_coef2
!!$       sp%Ed_carb(j) = inp%EdEbsurf*sp%Eb_carb(j)
!!$       sp%Ed_wat(j) = inp%EdEbsurf*sp%Eb_wat(j)
!!$       sp%Ed_H2(j) = inp%EdEbsurf*sp%Eb_H2(j)
!!$       sp%Ed_pure(j) = inp%EdEbsurf*sp%Eb_pure(j)
!!$       sp%Ed(j) = inp%EdEbsurf*sp%Eb(j)
!!$    END DO

    sp%num = 0
    ! read and attribute energies for each species
    DO is=1,inp%Nspeciestot

       READ(file_energies,fmt_energies) spec_num2, spec_name2, &
            spec_Eb_carb2, spec_Eb_wat2, spec_Eb_H22, spec_Eb_pure2

       check_spec = 0
       DO j=inp%Nspgas+1,inp%Nspgas+inp%Nspgr
                 
          IF (trim(sp%name(j)) .eq. trim(spec_name2)) THEN

          sp%num(j) = j

             check_spec = check_spec + 1
             
             sp%Eb_wat(j) = spec_Eb_wat2*inp%randomEb(j)
             sp%Eb_carb(j) = spec_Eb_carb2*inp%randomEb(j)
             sp%Eb_HH93(j) = spec_Eb_HH932*inp%randomEb(j)
             sp%Eb_H2(j) = spec_Eb_H22*inp%randomEb(j)
             sp%Eb_H2(j) = sp%Eb_wat(j)/sp%Eb_wat(ns%JH)*&
                  sp%Eb_H2(ns%JH)
             sp%Eb_pure(j) = spec_Eb_pure2*inp%randomEb(j)
             !sp%Eb_pure(j) = sp%Eb_wat(j)/1.15d3*8.5d2
             sp%Eb(j) = sp%Eb_carb(j)
             sp%Eb_ini(j) = sp%Eb(j)
             sp%cr_coef(j) = spec_cr_coef2
             !write(6,*) sp%name(j), sp%Eb(j)
             
             IF (j .eq. ns%JH) THEN 
                sp%Ed(j) = inp%REdH*sp%Eb(j)
             ELSE
                sp%Ed(j) = inp%REdoth*sp%Eb(j)
             END IF

             IF (inp%chlayer .eq. 2) THEN

                sp%Eb_wat(j+inp%Nspgr) = sp%Eb_wat(j)
                sp%Eb_carb(j+inp%Nspgr) = sp%Eb_carb(j)
                sp%Eb_HH93(j+inp%Nspgr) = sp%Eb_HH93(j)
                sp%Eb_H2(j+inp%Nspgr) = sp%Eb_H2(j)
                sp%Eb_pure(j+inp%Nspgr) = sp%Eb_pure(j)
                sp%Eb(j+inp%Nspgr) = sp%Eb(j) 
                sp%Eb_ini(j+inp%Nspgr) = sp%Eb(j)
                sp%cr_coef(j+inp%Nspgr) = spec_cr_coef2
                
                IF (j .eq. ns%JH) THEN 
                   sp%Ed(j+inp%Nspgr) = inp%REdbulkH*sp%Eb(j)
                ELSE
                   sp%Ed(j+inp%Nspgr) = inp%REdbulkoth*sp%Eb(j)
                END IF
                

             END IF
             
          END IF
             
       END DO
       DO j=1,inp%Nspecies
          IF (trim(sp%name(j)) .eq. trim(spec_name2)) THEN
             sp%Eb_H2(j) = spec_Eb_H22*inp%randomEb(j)
             sp%Eb_H2(j) = sp%Eb_wat(j)/sp%Eb_wat(ns%JH)*sp%Eb_H2(ns%JH)
             IF (j .eq. ns%JH) THEN 
                sp%Ed_H2(j) = inp%REdH*sp%Eb_H2(j)
             ELSE
                sp%Ed_H2(j) = inp%REdoth*sp%Eb_H2(j)
             END IF
          END IF
       END DO
    END DO
    
    ! attribute energies to species in pores
    IF (inp%Nspgr .eq. 2*inp%Nspgr2) THEN
       DO j=1,inp%Nspgr2
          sp%Eb_wat(j+inp%Nspgas+inp%Nspgr2) = sp%Eb_wat(j+inp%Nspgas) 
          sp%Eb_carb(j+inp%Nspgas+inp%Nspgr2) = sp%Eb_carb(j+inp%Nspgas) 
          sp%Eb_HH93(j+inp%Nspgas+inp%Nspgr2) = sp%Eb_HH93(j+inp%Nspgas) 
          sp%Eb_H2(j+inp%Nspgas+inp%Nspgr2) = sp%Eb_H2(j+inp%Nspgas) 
          sp%Eb(j+inp%Nspgas+inp%Nspgr2) = sp%Eb_wat(j+inp%Nspgas) 
          sp%Eb_ini(j+inp%Nspgas+inp%Nspgr2) = sp%Eb(j+inp%Nspgas)  
          sp%cr_coef(j+inp%Nspgas+inp%Nspgr2) = sp%cr_coef(j+inp%Nspgas) 
          IF (j .eq. ns%JH) THEN 
             sp%Ed(j+inp%Nspgas+inp%Nspgr2) = inp%REdH*sp%Eb(j+inp%Nspgas+inp%Nspgr2)
          ELSE
             sp%Ed(j+inp%Nspgas+inp%Nspgr2) = inp%REdoth*sp%Eb(j+inp%Nspgas+inp%Nspgr2)
          END IF
       END DO
    END IF

    ! attribute binding energies to species without data
    DO j=inp%Nspgas+1,inp%Nspgas+inp%Nspgr
       IF (sp%num(j) .ne. j) THEN
          
          sp%Eb_wat(j) = sp%Eb_wat(j-1)
          sp%Eb_carb(j) = sp%Eb_carb(j-1)
          sp%Eb_HH93(j) = sp%Eb_HH93(j-1)
          sp%Eb_H2(j) = sp%Eb_H2(j-1)
          sp%Eb_H2(j) = sp%Eb_wat(j)/sp%Eb_wat(ns%JH)*sp%Eb_H2(ns%JH)
          sp%Eb_pure(j) = sp%Eb_pure(j-1)
          sp%Eb(j) = sp%Eb_wat(j) 
          sp%Eb_ini(j) = sp%Eb(j)
          sp%cr_coef(j) = sp%cr_coef(j-1)

          IF (j .eq. ns%JH) THEN 
             sp%Ed(j) = inp%REdH*sp%Eb(j)
          ELSE
             sp%Ed(j) = inp%REdoth*sp%Eb(j)
          END IF

       END IF
    END DO


    ! set Eb=0 for gas phase species
    DO j=1,inp%Nspgas
          sp%Eb_wat(j) = 0. ; sp%Eb_carb(j) = 0. ; sp%Eb_HH93(j) = 0.
          sp%Eb_H2(j) = 0. ; sp%Eb_H2(j) = 0. ; sp%Eb_pure(j) = 0.
          sp%Eb(j) = 0. ; sp%Eb_ini(j) = 0. ;  sp%Ed(j) = 0.
    END DO
    
   ! DO j=inp%Nspgas+1,inp%Nspgas+inp%Nspgr
   !   IF (sp%num(j) .ne. j)write(6,*) sp%name(j), sp%Eb(j), sp%Eb(j+inp%Nspgr)
   ! END DO


  END SUBROUTINE INPUT_ENERGIES


  SUBROUTINE INPUT_REACTIONS

    !----------------------
    ! reading reactions.in
    !----------------------

    USE SHARED_VARIABLES, only: inp, sp, ns, re, nr

    INTEGER(KIND=LONG) :: file_reactions, ir1, ir2, ir3, ip1, ip2, ip3, ip4, indsp
    CHARACTER(len=100) :: direct_reactions
    INTEGER(KIND=LONG) :: ic, is, is2, log_reac, log_nondiffreac, inz, icheck, itype
    INTEGER(KIND=LONG) :: log_acc, log_evth, log_evCR, log_evUVCR, log_disUVCR, log_evcol
    CHARACTER(len=29) :: fmt_reac1
    CHARACTER(len=36) :: fmt_reac2
    CHARACTER(len=24) :: fmt_reac3
    CHARACTER(len=92) :: fmt_reac
    CHARACTER(len=5) :: Nreactionschar
    CHARACTER(len=10) :: toto
    REAL(KIND=DP), DIMENSION(14,inp%Nreactions) :: elemreac, elemprod
    REAL(KIND=DP), DIMENSION(inp%Nreactions) :: chargereac, chargeprod
    
    WRITE(6,*) "Reading the reactions file..."

    fmt_reac1="(A15,A15,A15,A15,A15,A15,A15,"
    fmt_reac2="ES9.2,1x,ES9.2,1x,ES9.2,1x,I2,1x,I5,"
    fmt_reac3="1x,I5,1x,I2,1x,I5,1x,I2)"
    fmt_reac=fmt_reac1//fmt_reac2//fmt_reac3

    re%react1(1:inp%Nreactions) = ' '
    re%react2(1:inp%Nreactions) = ' '
    re%react3(1:inp%Nreactions) = ' '
    re%prod1(1:inp%Nreactions) = ' '
    re%prod2(1:inp%Nreactions) = ' '
    re%prod3(1:inp%Nreactions) = ' '
    re%prod4(1:inp%Nreactions) = ' '
    re%act_en(1:inp%Nreactions) = 0d0
    log_acc = 1
    log_evth = 1
    log_evCR = 1
    log_evUVCR = 1
    log_evcol = 1
    sp%maxreac(1:2,1:inp%Nspecies) = 0d0
    indsp = 1

    ! open file
    direct_reactions="../input/"//TRIM(inp%folder)//"/"//TRIM(inp%filereac)
    file_reactions = 108
    OPEN(file_reactions,file=TRIM(direct_reactions), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    REWIND(file_reactions)

    READ(file_reactions,*)
    READ(file_reactions,*)
    READ(file_reactions,*)
    READ(file_reactions,*)
    READ(file_reactions,*)

    re%num_prod2 = 0d0

    re%num_react1 = inp%Neq+1 ; re%num_react2 = inp%Neq+1 ; re%num_react3 = inp%Neq+1
    re%num_prod1 = inp%Neq+1 ; re%num_prod2 = inp%Neq+1 ; re%num_prod3 = inp%Neq+1 ; re%num_prod4 = inp%Neq+1

    DO ic=1,inp%Nreactions

       ! read reaction properties
       READ(file_reactions,fmt_reac) re%react1(ic), re%react2(ic),&
            re%react3(ic), re%prod1(ic), re%prod2(ic), &
            re%prod3(ic), re%prod4(ic), re%A(ic), &
            re%B(ic), re%C(ic), re%type(ic), re%Tmin(ic), re%Tmax(ic), &
            re%formula(ic), re%num(ic), re%Nrate(ic)
       re%num(ic) = ic
       IF (re%type(ic) .ge. 14) re%C(ic) = re%C(ic)*inp%randomEa(ic)
       !IF (re%type(ic) .ne. 13) re%C(ic) = 0.

!!$       IF (re%C(ic) .lt. 0 .and. re%type(ic) .lt. 14) THEN 
!!$          write(6,*) re%react1(ic), re%react2(ic),&
!!$               re%react3(ic), re%prod1(ic), re%prod2(ic), &
!!$               re%prod3(ic), re%prod4(ic), re%A(ic), &
!!$               re%B(ic), re%C(ic), re%type(ic), re%Tmin(ic), re%Tmax(ic), &
!!$               re%formula(ic), re%num(ic), re%Nrate(ic)
!!$       END IF
       
       IF (re%Tmin(ic) .eq. 10 .and. re%Tmax(ic) .gt. 10) re%Tmin(ic) = 0d0
       IF (ic .gt. 1) THEN
          IF (re%Tmax(ic) .le. 300 .and. re%Tmax(ic) .ge. 200 .and. &
               (re%Tmax(ic) .ne. re%Tmin(ic-1) .or. re%Tmax(ic) .ne. re%Tmin(ic+1))) re%Tmax(ic) = 9999.
          IF (re%Tmax(ic) .ge. 200 .and. &
               (re%Tmax(ic) .ne. re%Tmin(ic-1) .or. re%Tmax(ic) .ne. re%Tmin(ic+1))) re%Tmax(ic) = 9999.
       END IF

       
       DO is=1,inp%Nspecies
          IF (re%react1(ic) == sp%name(is)) re%num_react1(ic) = is
          IF (re%react2(ic) == sp%name(is)) re%num_react2(ic) = is
          IF (re%react3(ic) == sp%name(is)) re%num_react3(ic) = is
          IF (re%prod1(ic) == sp%name(is)) re%num_prod1(ic) = is
          IF (re%prod2(ic) == sp%name(is)) re%num_prod2(ic) = is
          IF (re%prod3(ic) == sp%name(is)) re%num_prod3(ic) = is
          IF (re%prod4(ic) == sp%name(is)) re%num_prod4(ic) = is
       END DO
       
       ! attribute right number for non-species
       IF (re%num_react1(ic) .le. 0 .or. re%num_react1(ic) .ge. inp%Neq+1) &
            re%num_react1(ic) = inp%Neq+1
       IF (re%num_react2(ic) .le. 0 .or. re%num_react2(ic) .ge. inp%Neq+1) &
            re%num_react2(ic) = inp%Neq+1
       IF (re%num_react3(ic) .le. 0 .or. re%num_react3(ic) .ge. inp%Neq+1) &
            re%num_react3(ic) = inp%Neq+1
       IF (re%num_prod1(ic) .le. 0 .or. re%num_prod1(ic) .ge. inp%Neq+1) &
            re%num_prod1(ic) = inp%Neq+1
       IF (re%num_prod2(ic) .le. 0 .or. re%num_prod2(ic) .ge. inp%Neq+1) &
            re%num_prod2(ic) = inp%Neq+1
       IF (re%num_prod3(ic) .le. 0 .or. re%num_prod3(ic) .ge. inp%Neq+1) &
            re%num_prod3(ic) = inp%Neq+1
       IF (re%num_prod4(ic) .le. 0 .or. re%num_prod4(ic) .ge. inp%Neq+1) &
            re%num_prod4(ic) = inp%Neq+1

       ! save the formation reactions for each species
       IF (re%num_prod1(ic) .le. inp%Nspecies) THEN
          sp%maxreac(1,re%num_prod1(ic)) = sp%maxreac(1,re%num_prod1(ic))+1
          !sp%formreac(re%num_prod1(ic),sp%maxreac(1,re%num_prod1(ic))) = ic
       END IF
       IF (re%num_prod2(ic) .le. inp%Nspecies) THEN
          !sp%maxreac(1,re%num_prod2(ic)) = sp%maxreac(1,re%num_prod2(ic))+1
          !write(6,*) re%num_prod2(ic),sp%maxreac(1,re%num_prod2(ic)), inp%Nspecies
          !sp%formreac(re%num_prod2(ic),sp%maxreac(1,re%num_prod2(ic))) = ic
       END IF

       ! save the destruction reactions for each species
       IF (re%num_react1(ic) .le. inp%Nspecies) THEN
          sp%maxreac(2,re%num_react1(ic)) = sp%maxreac(2,re%num_react1(ic))+1
          !sp%destreac(re%num_react1(ic),sp%maxreac(2,re%num_react1(ic))) = ic
       END IF
       IF (re%num_react2(ic) .le. inp%Nspecies) THEN
          sp%maxreac(2,re%num_react2(ic)) = sp%maxreac(2,re%num_react2(ic))+1
          !sp%destreac(re%num_react2(ic),sp%maxreac(2,re%num_react2(ic))) = ic
       END IF
       IF (re%num_react3(ic) .le. inp%Nspecies) THEN
          sp%maxreac(2,re%num_react3(ic)) = sp%maxreac(2,re%num_react3(ic))+1
          !sp%destreac(re%num_react3(ic),sp%maxreac(2,re%num_react3(ic))) = ic
       END IF

       ! count the number of grain surface reactions
       IF (re%type(ic) .eq. 14) THEN
          IF (log_reac .eq. 1) THEN
             nr%reac_ini = re%num(ic)
             log_reac = 2
          ELSE IF (log_reac .eq. 2) THEN
             nr%reac_fin = re%num(ic)
          END IF
          IF (re%num_prod1(ic).le.inp%Nspgas .and. re%B(ic).ne.1.) THEN
             !re%B(ic) = re%B(ic)/512.
             !if (re%B(ic).ge.1) re%B(ic) = 1.
             !re%B(ic) = 0.012
             !re%B(ic-1) = 1.-re%B(ic)
             !write(6,*) ic,re%prod1(ic),re%B(ic),re%B(ic-1)
          END IF
       END IF
       ! non-diffusion reactions
       IF (re%type(ic) .eq. 18) THEN
          IF (log_nondiffreac .eq. 1) THEN
             nr%nondiffreac_ini = re%num(ic)
             log_nondiffreac = 2
             re%ndreact1(1,1) = re%react1(ic)
             re%num_ndreact1(1,1) = re%num_react1(ic)
             re%ndreact1(1,2) = sp%name(re%num_react1(ic)-1)
             re%num_ndreact1(1,2) = re%num_react1(ic)-1
             re%ndreact2(1) = re%react2(ic)
             re%num_ndreact2(1) = re%num_react2(ic)
             indsp = 2
          ELSE IF (log_nondiffreac .eq. 2) THEN
             nr%nondiffreac_fin = re%num(ic)
             IF (re%react1(ic).ne.re%react1(ic-1).and.&
                  re%react1(ic).ne.re%react1(ic-2).and.&
                  re%react1(ic).ne.re%react1(ic-3)) THEN
                re%ndreact1(indsp,1) = re%react1(ic)
                re%num_ndreact1(indsp,1) = re%num_react1(ic)
                re%ndreact1(indsp,2) = sp%name(re%num_react1(ic)-1)
                re%num_ndreact1(indsp,2) = re%num_react1(ic)-1
                re%ndreact2(indsp) = re%react2(ic)
                re%num_ndreact2(indsp) = re%num_react2(ic)
                indsp = indsp+1
             END IF
             !write(6,*) ic, indsp-1, re%num_ndreact1(indsp-1,1), re%num_ndreact1(indsp-1,2)
          END IF
       END IF
       ! count the number of accretion reactions
       IF (re%type(ic) .eq. 20) THEN
          IF (log_acc .eq. 1) THEN
             nr%acc_ini = re%num(ic)
             log_acc = 2
          ELSE IF (log_acc .eq. 2) THEN
             nr%acc_fin = re%num(ic)
          END IF
          IF (TRIM(re%react1(ic)) .eq. 'H2' .or. &
               TRIM(re%react1(ic)) .eq. 'oH2' ) THEN
             nr%accH2 = re%num(ic)
          END IF
          !write(6,*) re%num_react1(ic)
          sp%grproc(re%num_react1(ic),1) = ic
       END IF
       ! count the number of desorption reactions
       IF (re%type(ic) .eq. 21) THEN
          IF (log_evth .eq. 1) THEN
             nr%evth_ini = re%num(ic)
             log_evth = 2
          ELSE IF (log_evth .eq. 2) THEN
             nr%evth_fin = re%num(ic)
          END IF
          !write(6,*) re%react1(ic)
          sp%grproc(re%num_react1(ic),2) = ic
       END IF          
       IF (re%type(ic) .eq. 22) THEN
          IF (log_evCR .eq. 1) THEN
             nr%evCR_ini = re%num(ic)
             log_evCR = 2
          ELSE IF (log_evCR .eq. 2) THEN
             nr%evCR_fin = re%num(ic)
          END IF
          sp%grproc(re%num_react1(ic),3) = ic
       END IF
       IF (re%type(ic) .eq. 23) THEN
          IF (log_evUVCR .eq. 1) THEN
             nr%evUVCR_ini = re%num(ic)
             log_evUVCR = 2
          ELSE IF (log_evUVCR .eq. 2) THEN
             nr%evUVCR_fin = re%num(ic)
          END IF
          sp%grproc(re%num_react1(ic),4) = ic
       END IF
       IF (re%type(ic) .eq. 26) THEN
          IF (log_disUVCR .eq. 1) THEN
             nr%disUVCR_ini = re%num(ic)
             log_disUVCR = 2
          ELSE IF (log_disUVCR .eq. 2) THEN
             nr%disUVCR_fin = re%num(ic)
          END IF
          sp%grproc(re%num_react1(ic),5) = ic
       END IF
       IF (re%type(ic) .eq. 28) THEN
          IF (log_evcol .eq. 1) THEN
             nr%evcol_ini = re%num(ic)
             log_disUVCR = 2
          ELSE IF (log_evcol .eq. 2) THEN
             nr%evcol_fin = re%num(ic)
          END IF
          sp%grproc(re%num_react1(ic),5) = ic
       END IF

       ! count the number of several reference reactions
       IF (trim(re%react1(ic)) .eq. 'JCO' .and. &
            trim(re%react2(ic)) .eq. 'JH' .and. &
            trim(re%prod1(ic)) .eq. 'JHCO') THEN
          nr%COH = ic
       END IF
       IF (trim(re%react1(ic)) .eq. 'JH2O' .and. &
            trim(re%prod1(ic)) .eq. 'H2O' .and. &
            (re%type(ic) .eq. 26 .or. re%type(ic) .eq. 27)) THEN
          nr%ref_photo = ic
       END IF
       IF (trim(re%react1(ic)) .eq. 'JCO' .and. &
            trim(re%react2(ic)) .eq. 'JH' .and. & 
            trim(re%prod1(ic)) .eq. 'JHCO') THEN
          nr%reacref = ic
       END IF    
       IF (trim(re%react1(ic)) .eq. 'JCO' .and. &
            trim(re%react2(ic)) .eq. 'JO' .and. & 
            trim(re%react3(ic)) .eq. 'JH') THEN
          nr%reacCO2 = ic
       END IF

       ! check if previous reaction is a counterpart
       IF (ic .gt. 1) THEN
          IF (re%react1(ic).eq.re%react1(ic-1).and.&
               re%react2(ic).eq.re%react2(ic-1).and.&
               re%prod1(ic).eq.re%prod1(ic-1).and.&
               re%prod2(ic).eq.re%prod2(ic-1).and.&
               re%prod3(ic).eq.re%prod3(ic-1).and.&
               re%prod4(ic).eq.re%prod4(ic-1).and.&
               re%Tmin(ic).eq.re%Tmax(ic-1)) THEN
             re%sim(ic-1)=ic
             re%sim(ic)=ic-1
          END IF
       END IF
       
    END DO

    ! recompute number of species for non diffusion reactions
    inp%Nspnd = indsp-1   
    
    ! dummy reaction
    re%num_react1(inp%Nreactions+1) = 1; re%num_react2(inp%Nreactions+1) = 1
    re%num_react3(inp%Nreactions+1) = 1
    re%num_prod1(inp%Nreactions+1) = 1; re%num_prod2(inp%Nreactions+1) = 1
    re%num_prod3(inp%Nreactions+1) = 1; re%num_prod4(inp%Nreactions+1) = 1
    re%A(inp%Nreactions+1) = 1.; re%B(inp%Nreactions+1) = 1.
    re%C(inp%Nreactions+1) = 1.; re%type(inp%Nreactions+1) = 1.
    re%Tmin(inp%Nreactions+1) = 1.; re%Tmax(inp%Nreactions+1) = 1.
    re%formula(inp%Nreactions+1) = 1.; re%num(inp%Nreactions+1) = 1.
    re%Nrate(inp%Nreactions+1) = 1.
    
    ! compute the conservation of elements for all reactions
    DO ic=1,inp%Nreactions

       DO is=1,Nelements
          elemreac(is,ic) = sp%element(is,re%num_react1(ic))+&
               sp%element(is,re%num_react2(ic))+&
               sp%element(is,re%num_react3(ic))
          
          elemprod(is,ic) = sp%element(is,re%num_prod1(ic))+&
               sp%element(is,re%num_prod2(ic))+&
               sp%element(is,re%num_prod3(ic))+&
               sp%element(is,re%num_prod4(ic))
          
          IF ((elemreac(is,ic) .ne. elemprod(is,ic) .and. &
               elemreac(is,ic) .gt. 1d-40)) THEN
             WRITE(6,*) re%react1(ic), re%react2(ic),&
                  re%react3(ic), re%prod1(ic), re%prod2(ic), &
                  re%prod3(ic), re%prod4(ic)
             WRITE(6,*) ic, sp%name_elem(is), &
                  re%num_react1(ic), &
                  re%num_react2(ic), &
                  re%num_react3(ic), &
                  re%num_prod1(ic), &
                  re%num_prod2(ic), &
                  re%num_prod3(ic), &
                  re%num_prod4(ic)
             WRITE(6,*) ic, sp%name_elem(is), &
                  sp%element(is,re%num_react1(ic)), &
                  sp%element(is,re%num_react2(ic)), &
                  sp%element(is,re%num_react3(ic)), &
                  sp%element(is,re%num_prod1(ic)), &
                  sp%element(is,re%num_prod2(ic)), &
                  sp%element(is,re%num_prod3(ic)), &
                  sp%element(is,re%num_prod4(ic)), &
                  elemreac(is,ic), elemprod(is,ic)
             !write(6,*) re%num_react1(ic), re%num_react2(ic), re%num_prod1(ic)
             WRITE(6,*) ic, ' Elements are not conservative !!!!!!'
          END IF
       END DO


       chargereac(ic) = sp%charge(re%num_react1(ic))+&
            sp%charge(re%num_react2(ic))+&
            sp%charge(re%num_react3(ic))
       chargeprod(ic) = sp%charge(re%num_prod1(ic))+&
            sp%charge(re%num_prod2(ic))+&
            sp%charge(re%num_prod3(ic))+&
            sp%charge(re%num_prod4(ic))

       IF (chargereac(ic) .ne. chargeprod(ic)) THEN
          WRITE(6,*) re%react1(ic), re%react2(ic),&
               re%react3(ic), re%prod1(ic), re%prod2(ic), &
               re%prod3(ic), re%prod4(ic)
          WRITE(6,*) ic, ' Charge is not conservative !!!!!!'
       END IF
       
    END DO

    ! locate the first and last reactions of each type
    DO itype=0,Ntype
       nr%istart(itype) = inp%Nreactions+1
       nr%iend(itype) = inp%Nreactions+1
       DO ic=1,inp%Nreactions
          IF (itype .eq. re%type(ic) .and. nr%istart(itype) .eq. inp%Nreactions+1) &
               nr%istart(itype) = ic
          IF (itype .eq. re%type(ic)) nr%iend(itype) = ic
       END DO
    END DO
    
!!$    DO ic=1,inp%Nreactions
!!$       IF (re%sim(ic).eq.0) THEN
!!$          re%Tmin(ic) = -9999.
!!$          re%Tmax(ic) = 9999.
!!$       END IF
!!$    END DO
    
    WRITE(6,*) 'Nreactions=',inp%Nreactions
    
  END SUBROUTINE INPUT_REACTIONS
  

  SUBROUTINE INPUT_RATES!(inp,sp,ns,re)

    !---------------------
    ! read formrates.in
    !---------------------

    !IMPLICIT NONE

    USE SHARED_VARIABLES, only: inp, sp, ns, re
    
    INTEGER(KIND=LONG) :: ir, is, is2, is3, ic
    INTEGER(KIND=LONG) :: file_rate
    CHARACTER(len=100) :: direct_rate
    INTEGER(KIND=LONG), DIMENSION(inp%Nspecies) :: ispecform, ispecdest
    
    ispecform = 0d0
    ispecdest = 0d0

    ! open file
    direct_rate = "../input/"//TRIM(inp%folder)//"/reacrates.in"
    file_rate = 102
    OPEN(file_rate,file=TRIM(direct_rate), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    REWIND(file_rate)
    
    ! read species names to plot
    READ(file_rate,*)
    READ(file_rate,*)
    READ(file_rate,*)
    DO ir=1,inp%Nsprate
       READ(file_rate,'(A10)') sp%rate_name(ir)
    END DO

    ! find the numbers of each species
    is3 = 1
    DO is2=1,inp%Nsprate
       DO is=1,inp%Nspecies
          IF (trim(sp%name(is)) .eq. trim(sp%rate_name(is2))) THEN
             sp%rate_num(is2) = is
          END IF
       END DO
    END DO

    ! select the reactions involved in the form/dest of each species
    DO is=1,inp%Nsprate
       DO ic=1,inp%Nreactions
       
          ! select the reactions involved in the form/dest of each species
          IF (re%num_react1(ic) .eq. sp%rate_num(is) .or. &
               re%num_react2(ic) .eq. sp%rate_num(is) .or. &
               re%num_react3(ic) .eq. sp%rate_num(is)) THEN
             
             ispecdest(is) = ispecdest(is) + 1
             sp%numdest(is,ispecdest(is)) = ic
             sp%degdest(is,ispecdest(is)) = 1
             
             IF ((re%num_react1(ic) .eq. sp%rate_num(is) .and. &
                  re%num_react2(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_react1(ic) .eq. sp%rate_num(is) .and. &
                  re%num_react3(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_react2(ic) .eq. sp%rate_num(is) .and. &
                  re%num_react3(ic) .eq. sp%rate_num(is))) THEN
                sp%degdest(is,ispecdest(is)) = 2
             END IF
             IF (re%num_react1(ic) .eq. sp%rate_num(is) .and. &
                  re%num_react2(ic) .eq. sp%rate_num(is) .and. &
                  re%num_react3(ic) .eq. sp%rate_num(is)) THEN
                sp%degdest(is,ispecdest(is)) = 3
             END IF
             
          ELSE IF (re%num_prod1(ic) .eq. sp%rate_num(is) .or. &
               re%num_prod2(ic) .eq. sp%rate_num(is) .or. &
               re%num_prod3(ic) .eq. sp%rate_num(is) .or. &
               re%num_prod4(ic) .eq. sp%rate_num(is)) THEN
             
             ispecform(is) = ispecform(is) + 1
             sp%numform(is,ispecform(is)) = ic
             sp%degform(is,ispecform(is)) = 1
             
             IF ((re%num_prod1(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod2(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_prod1(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod3(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_prod1(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod4(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_prod2(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod3(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_prod2(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod4(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_prod2(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod1(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_prod3(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod1(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_prod3(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod2(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_prod3(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod4(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_prod4(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod1(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_prod4(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod2(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_prod4(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod3(ic) .eq. sp%rate_num(is))) THEN
                sp%degform(is,ispecform(is)) = 2
             END IF
             IF ((re%num_prod1(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod2(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod3(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_prod1(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod2(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod4(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_prod1(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod3(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod4(ic) .eq. sp%rate_num(is)) .or. &
                  (re%num_prod2(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod3(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod4(ic) .eq. sp%rate_num(is))) THEN
                sp%degform(is,ispecform(is)) = 3
             END IF
             IF ((re%num_prod1(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod2(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod3(ic) .eq. sp%rate_num(is) .and. &
                  re%num_prod4(ic) .eq. sp%rate_num(is))) THEN
                sp%degform(is,ispecform(is)) = 4
             END IF
             
          END IF
       
       END DO
    END DO

    ! compute the number of reactions involved in the form/dest of each species
    DO is=1,inp%Nsprate
       sp%Nreacform(is) = ispecform(is)
       sp%Nreacdest(is) = ispecdest(is)
    END DO
    sp%maxNreacform = maxval(sp%Nreacform(1:inp%Nsprate),inp%Nsprate)
    sp%maxNreacdest = maxval(sp%Nreacdest(1:inp%Nsprate),inp%Nsprate)
    sp%maxNreac = max(sp%maxNreacform,sp%maxNreacdest)
    sp%maxNreac2 = 30 !min(sp%maxNreac,inp%Nreactions)

  END SUBROUTINE INPUT_RATES



  SUBROUTINE INPUT_GRID
    
    !----------------------
    ! read model_grid.in
    !----------------------

    USE SHARED_VARIABLES, only: inp, gi
    
    INTEGER(KIND=LONG) :: file_grid, i
    CHARACTER(len=41) :: direct_grid

    ! open file
    direct_grid="../input/"//TRIM(inp%folder)//"/model_grid.in"
    file_grid = 111
    OPEN(file_grid,file=TRIM(direct_grid), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    REWIND(file_grid)

    ! read file
    READ(file_grid,*) 
    READ(file_grid,'(I1)') gi%chabuini
    READ(file_grid,*) 
    READ(file_grid,'(A100)') gi%input_abuini
    READ(file_grid,*) 
    READ(file_grid,'(I2)') gi%NnH
    READ(file_grid,'(ES13.7)') gi%nHlow
    READ(file_grid,'(ES13.7)') gi%nHup
    READ(file_grid,*)
    READ(file_grid,'(I2)') gi%NT
    READ(file_grid,'(ES13.7)') gi%Tlow
    READ(file_grid,'(ES13.7)') gi%Tup
    READ(file_grid,*) 
    READ(file_grid,'(I2)') gi%NUV
    READ(file_grid,'(ES13.7)') gi%UVlow
    READ(file_grid,'(ES13.7)') gi%UVup
    READ(file_grid,*)
    READ(file_grid,'(I2)') gi%Nzeta
    READ(file_grid,'(ES13.7)') gi%zetalow
    READ(file_grid,'(ES13.7)') gi%zetaup
    READ(file_grid,*)
    READ(file_grid,'(I2)') gi%NAv
    READ(file_grid,'(ES13.7)') gi%Avlow
    READ(file_grid,'(ES13.7)') gi%Avup
    READ(file_grid,*) 
    READ(file_grid,'(I2)') gi%Nacst
    READ(file_grid,'(ES13.7)') gi%acstlow
    READ(file_grid,'(ES13.7)') gi%acstup
    READ(file_grid,*) 
    READ(file_grid,'(I1)') gi%choicetreat
    READ(file_grid,*) 
    READ(file_grid,'(I2)') gi%Nenratio
    READ(file_grid,'(ES13.7)') gi%enratiolow
    READ(file_grid,'(ES13.7)') gi%enratioup
    READ(file_grid,*) 
    READ(file_grid,'(I2)') gi%NEbH
    READ(file_grid,'(ES13.7)') gi%EbHlow
    READ(file_grid,'(ES13.7)') gi%EbHup
    READ(file_grid,*) 
    READ(file_grid,'(I2)') gi%Nporo
    READ(file_grid,'(ES13.7)') gi%porolow
    READ(file_grid,'(ES13.7)') gi%poroup
    READ(file_grid,*) 
    READ(file_grid,'(I2)') gi%Nds
    READ(file_grid,'(ES13.7)') gi%dslow
    READ(file_grid,'(ES13.7)') gi%dsup
    READ(file_grid,*) 
    READ(file_grid,'(I2)') gi%NEa
    READ(file_grid,'(ES13.7)') gi%Ealow
    READ(file_grid,'(ES13.7)') gi%Eaup
    READ(file_grid,*) 
    READ(file_grid,'(I2)') gi%NXO
    READ(file_grid,'(ES13.7)') gi%XOlow
    READ(file_grid,'(ES13.7)') gi%XOup
    READ(file_grid,*) 
    READ(file_grid,'(I2)') gi%Nop
    READ(file_grid,'(ES13.7)') gi%oplow
    READ(file_grid,'(ES13.7)') gi%opup

    IF (gi%choicetreat .eq. 0) gi%Ntreat = 1
    IF (gi%choicetreat .eq. 1) gi%Ntreat = 1
    IF (gi%choicetreat .eq. 2) gi%Ntreat = 1
    IF (gi%choicetreat .eq. 3) gi%Ntreat = 2
   
    
    ! compute the number of free parameters
    gi%NnH2 = gi%NnH
    gi%NT2 = gi%NT
    gi%NUV2 = gi%NUV
    gi%Nzeta2 = gi%Nzeta
    gi%NAv2 = gi%NAv
    gi%Nacst2 = gi%Nacst
    gi%Ntreat2 = gi%Ntreat
    gi%Nenratio2 = gi%Nenratio
    gi%NEbH2 = gi%NEbH
    gi%Nporo2 = gi%Nporo
    gi%Nds2 = gi%Nds
    gi%NEa2 = gi%NEa
    gi%NXO2 = gi%NXO
    gi%Nop2 = gi%Nop

    IF (gi%NnH .eq. 0) gi%NnH2 = 1
    IF (gi%NT .eq. 0) gi%NT2 = 1
    IF (gi%NUV .eq. 0) gi%NUV2 = 1
    IF (gi%Nzeta .eq. 0) gi%Nzeta2 = 1
    IF (gi%NAv .eq. 0) gi%NAv2 = 1
    IF (gi%Nacst .eq. 0) gi%Nacst2 = 1
    IF (gi%Nenratio .eq. 0) gi%Nenratio2 = 1
    IF (gi%NEbH .eq. 0) gi%NEbH2 = 1
    IF (gi%Nporo .eq. 0) gi%Nporo2 = 1
    IF (gi%Nds .eq. 0) gi%Nds2 = 1
    IF (gi%NEa .eq. 0) gi%NEa2 = 1
    IF (gi%NXO .eq. 0) gi%NXO2 = 1
    IF (gi%Nop .eq. 0) gi%Nop2 = 1

    gi%Nparam=0
    IF (gi%NnH .ne. 0) gi%Nparam=gi%Nparam+1
    IF (gi%NT .ne. 0) gi%Nparam=gi%Nparam+1
    IF (gi%NUV .ne. 0) gi%Nparam=gi%Nparam+1
    IF (gi%Nzeta .ne. 0) gi%Nparam=gi%Nparam+1
    IF (gi%NAv .ne. 0) gi%Nparam=gi%Nparam+1
    IF (gi%Nacst .ne. 0) gi%Nparam=gi%Nparam+1
    IF (gi%Ntreat .ne. 0) gi%Nparam=gi%Nparam+1
    IF (gi%Nenratio .ne. 0) gi%Nparam=gi%Nparam+1
    IF (gi%NEbH .ne. 0) gi%Nparam=gi%Nparam+1
    IF (gi%Nporo .ne. 0) gi%Nparam=gi%Nparam+1
    IF (gi%Nds .ne. 0) gi%Nparam=gi%Nparam+1
    IF (gi%NEa .ne. 0) gi%Nparam=gi%Nparam+1
    IF (gi%NXO .ne. 0) gi%Nparam=gi%Nparam+1
    IF (gi%Nop .ne. 0) gi%Nparam=gi%Nparam+1

    ! compute the number of runs
    gi%Nmodels = gi%NnH2*gi%NT2*gi%NUV2*gi%Nzeta2*gi%NAv2*&
         gi%Nacst2*gi%Ntreat2*gi%Nenratio2*gi%NEbH2*gi%Nporo2*&
         gi%Nds2*gi%NEa2*gi%NXO2*gi%Nop2

    CLOSE(file_grid)

  END SUBROUTINE INPUT_GRID


  SUBROUTINE INPUT_SPAT

    USE SHARED_VARIABLES, only: inp, spat
    
    INTEGER(KIND=LONG) :: i
    CHARACTER(len=100) :: direct_spat
    INTEGER(KIND=LONG) :: file_spat
    CHARACTER(len=1) :: checkchar

    direct_spat = "../input/"//TRIM(inp%folder)//"/spat_evol.in"

    file_spat = 1
    OPEN(file_spat,file=TRIM(direct_spat), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    REWIND(file_spat) 
    checkchar = 'A'
    READ(file_spat,'(A100)') spat%folder
    DO WHILE (checkchar .ne. '*')
      READ(file_spat,'(A1)') checkchar
      spat%Nspat = spat%Nspat+1
    END DO
    spat%Nspat = spat%Nspat-1

    ALLOCATE (spat%file(spat%Nspat), stat=ok)
    REWIND(file_spat) 
    READ(file_spat,'(A100)') spat%folder

    DO i=1,spat%Nspat
       READ(file_spat,'(A80)') spat%file(i)
    END DO

    CLOSE(file_spat)

  END SUBROUTINE INPUT_SPAT


  SUBROUTINE INPUT_PHYSPARAM

    !---------------------------------------------------------------------
    ! read inphysparam.in:
    ! data table coming from physical model giving the evolution of 
    ! radius, nH, T, Av along one streamline
    ! ---------------------------------------------------------------------

    !IMPLICIT NONE

    USE SHARED_VARIABLES, only: inp, iph
    
    CHARACTER(len=200) :: direct_inphys
    CHARACTER(len=100) :: formatphys, toto
    INTEGER(KIND=LONG) :: file_inphys, ipt, checkneg
    REAL(KIND=DP) :: timeneg

    direct_inphys = TRIM(inp%inputphys) !"../input/"//TRIM(inp%folder)//"/"//
    formatphys = '(2x,ES16.9,2x,ES16.9,2x,ES16.9,2x,ES16.9,2x,ES16.9,2x,ES16.9,2x,ES16.9)'
    checkneg = 1
    timeneg = 0d0

    file_inphys = 1
    OPEN(file_inphys,file=TRIM(direct_inphys), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    
    REWIND(file_inphys)
    READ(file_inphys,*) 
    READ(file_inphys,'(A23,I7)') toto, iph%Nphysteps
    READ(file_inphys,*)
    
    ALLOCATE(iph%time(iph%Nphysteps), iph%rad(iph%Nphysteps), stat=ok)
    ALLOCATE(iph%veloc(iph%Nphysteps), iph%nH(iph%Nphysteps), iph%Tg(iph%Nphysteps), stat=ok)
    ALLOCATE(iph%Td(iph%Nphysteps), iph%Av(iph%Nphysteps), stat=ok)

    DO ipt=1,iph%Nphysteps
       READ(file_inphys,formatphys) iph%time(ipt), iph%rad(ipt), iph%veloc(ipt), &
            iph%nH(ipt), iph%Td(ipt), iph%Tg(ipt), iph%Av(ipt)
       
       IF (ipt .eq. 1 .and. iph%time(1) .lt. 0) THEN
          timeneg = iph%time(ipt)
          checkneg = 2
       END IF
       IF (checkneg .eq. 2) iph%time(ipt) = iph%time(ipt) - timeneg
       !iph%Av(ipt) = iph%Av(ipt)/2.
       !iph%Td(ipt) = iph%Td(ipt) + 3.
       !iph%Tg(ipt) = iph%Tg(ipt) + 3.
       !iph%nH(ipt) = iph%nH(ipt)*3.33
       !write(6,*) ipt,timeneg,iph%rad(ipt), iph%time(ipt), &
       !     iph%nH(ipt), iph%Tg(ipt), iph%Td(ipt), iph%Av(ipt)
    END DO
    
    inp%tmax = iph%time(iph%Nphysteps)

    CLOSE(file_inphys)


  END SUBROUTINE INPUT_PHYSPARAM

  
  SUBROUTINE INPUT_SELFSHIELDING

    !------------------------------------------------------
    ! read self_shielding.in:
    ! data files for the self-shielding of H2, CO, and HD
    !------------------------------------------------------

    USE SHARED_VARIABLES, only: ss
    
    INTEGER(KIND=LONG) :: file_shield1, file_shield2, i
    CHARACTER(len=31) :: direct_shield1
    CHARACTER(len=41) :: direct_shield2
    CHARACTER(len=52) :: format1
    CHARACTER(len=16) :: format2 
    CHARACTER(len=25) :: format3

    direct_shield1="../input/data/self_shielding.in"
    direct_shield2="../input/data/self_shielding_PDRmeudon.in"

    ! read the CO and H2 self-shielding parameters from Lee et al. (1996) 
    file_shield1 = 113
    OPEN(file_shield1,file=TRIM(direct_shield1), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    ! format
    format1='(ES9.3,3x,ES9.3,3x,ES9.3,3x,ES9.3,3x,ES9.3,3x,ES9.3)'
    format2='(ES9.3,3x,ES9.3)'
    format3='(ES9.3,2x,ES9.3,2x,ES9.3)'

    REWIND(file_shield1)
                   
    READ(file_shield1,*)                                                   
    DO i=1,43                                                       
       READ(file_shield1,format1) &
            ss%NCO(i),ss%T_CO(i),ss%NH2_1(i),&
            ss%T_H2_1(i),ss%AV2(i),ss%T_AV(i)
    END DO
    DO i=44,52                                                      
       READ(file_shield1,format2) ss%NCO(i),ss%T_CO(i) 
    END DO            
    READ(file_shield1,*)                                                
    DO i=1,105                                                      
       READ(file_shield1,format2) ss%NH2_2(i),ss%T_H2_2(i) 
    ENDDO
    CLOSE(file_shield1)

    ! read the H2 and HD self-shielding parameters from the Meudon PDR code
    file_shield2 = 114
    OPEN(file_shield2,file=TRIM(direct_shield2), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    REWIND(file_shield2)
                 
    READ(file_shield2,*)                                           
    DO i=1,203                                                      
       READ(file_shield2,format3) ss%NH2_3(i),ss%T_H2_3(i),ss%T_HD(i) 
    ENDDO
    CLOSE(file_shield2)

  END SUBROUTINE INPUT_SELFSHIELDING


  SUBROUTINE INPUT_ECKART

    !--------------------------------------------------------------------------
    ! read eckart_data.in: 
    ! data file for computing transmission probabilities with the Eckart model
    !--------------------------------------------------------------------------

    USE SHARED_VARIABLES, only: ec
    
    INTEGER(KIND=LONG) :: file_eckart, i
    CHARACTER(len=35) :: direct_eckart
    CHARACTER(len=45) :: format1
    CHARACTER(len=16) :: format2 
    
    ! format
    format1='(A10,A10,A10,A10,A10,ES9.3,1x,ES9.3,1x,ES9.3)'

    ! open file
    direct_eckart="../input/data/eckart_data.in"
    !direct_eckart="../input/data/eckart_data_rimola.in"
    file_eckart = 116
    OPEN(file_eckart,file=TRIM(direct_eckart), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    REWIND(file_eckart)
                 
    READ(file_eckart,*)
    DO i=1,Neckartreacs
       READ(file_eckart,trim(format1)) ec%react1(i), ec%react2(i), &
            ec%react3(i), ec%prod1(i), ec%prod2(i), &
            ec%Vf(i), ec%Vr(i), ec%wf(i)
      ! write(6,*) i,ec%react1(i),ec%react2(i),ec%prod1(i),ec%prod2(i)
    END DO
    !stop
    CLOSE(file_eckart)
    
  END SUBROUTINE INPUT_ECKART

  
  SUBROUTINE INPUT_ENTHALPY

    !--------------------------------------------------------------------------
    ! read enthalpy_data.in: 
    ! data file for computing chemical desorption probabilities
    !--------------------------------------------------------------------------

    USE SHARED_VARIABLES, only: ent, sp, inp
    
    INTEGER(KIND=LONG) :: file_enthalpy, i, j
    CHARACTER(len=35) :: direct_enthalpy
    CHARACTER(len=45) :: format1
    CHARACTER(len=16) :: format2 
    
    ! format
    format1='(A10,F6.1)'

    ! open file
    direct_enthalpy="../input/data/enthalpy_data.in"
    file_enthalpy = 116
    OPEN(file_enthalpy,file=TRIM(direct_enthalpy), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    sp%dHf(1:inp%Neq) = -1e4
    
    REWIND(file_enthalpy)
                 
    READ(file_enthalpy,*)
    DO i=1,Nenthalpy
       READ(file_enthalpy,trim(format1)) ent%species(i), ent%dHf(i)

       DO j=inp%Nspgas+1,inp%Nspgas+inp%Nspgr
          
          sp%num(j) = j
       
          IF (trim(sp%name(j)) .eq. 'J'//trim(ent%species(i))) THEN
             sp%dHf(j) = ent%dHf(i)
          END IF
       END DO
    END DO

    CLOSE(file_enthalpy)
    
  END SUBROUTINE INPUT_ENTHALPY


  SUBROUTINE INPUT_CONETWORK

    !------------------------------------------------------------
    ! read COnetwork_data.in: 
    ! data to compute the reaction rates of the methanol network
    !------------------------------------------------------------

    USE SHARED_VARIABLES, only: con

    INTEGER(KIND=LONG) :: file_network, i
    CHARACTER(len=31) :: direct_COnetwork
    CHARACTER(len=45) :: format1
    CHARACTER(len=16) :: format2 
    
    format1='(A10,A10,A10,A10,A10,ES9.3,1x,ES9.3)'

    direct_COnetwork="../input/data/COnetwork_data.in"

    file_network = 116
    OPEN(file_network,file=TRIM(direct_COnetwork), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    REWIND(file_network)
                 
    READ(file_network,*)
    READ(file_network,*)
    READ(file_network,*)
    READ(file_network,*)
    DO i=1,NCOnetreacs
       READ(file_network,trim(format1)) con%react1(i), con%react2(i), &
            con%react3(i), con%prod1(i), con%prod2(i), &
            con%A(i), con%B(i)
    END DO

    CLOSE(file_network)


  END SUBROUTINE INPUT_CONETWORK

END MODULE
    
