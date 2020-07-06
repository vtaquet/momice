PROGRAM MAIN

  USE VARIABLES
  USE OUTPUT
  USE NUMPROP
  USE PHYSPROP
  USE CHEMPROP
  USE INPUT
  USE ODES
  USE ALLOCATION

  !!--------------------------------------------------------------------------
  !!--------------------------------------------------------------------------
  !! MAIN PROGRAM OF THE MOMICE CODE:
  !! Please read the readme file for a description of the physical and chemical 
  !! processes, the parameters, and the options
  !!
  !! Steps of the code:
  !! - call input files and create directories
  !! - allocate memory for species, reaction, and grid vectors
  !! - read input files: species, reactions, data
  !! - compute the initial conditions and initialize the vectors
  !! - TEMPORAL LOOP
  !!   - LAYER LOOP
  !!     - compute physical parameters and rates
  !!     - compute the timestep based on rate evaluation
  !!     - call ODE integrator
  !!     - save the parameters at the previous time
  !!     - if calculation does not converge, need to modify ATOL
  !!     - conversion and save of parameters
  !!     - save output parameters at specific timesteps
  !!     - check if the system is conservative by looking 
  !!       the elemental abundance and exit the layer loop if not
  !!     - leave the loop when t = tmax
  !!   - layer is filled
  !!   - system not conservative
  !! - calculation did not converge
  !! - layer is filled
  !! - write results in the terminal/files at the end of calculations
  !! - deallocate vectors
  !!
  !!--------------------------------------------------------------------------
  !!--------------------------------------------------------------------------


  !--------------------------------
  ! declaration of local variables
  !--------------------------------
  
  ! several structures are shared with Fode subroutine
  USE SHARED_VARIABLES
  USE SHARED_VARIABLES2

  IMPLICIT NONE

  INTEGER(KIND=LONG) :: error, im, it, ir, is, ic, ie, k
  INTEGER(KIND=LONG) :: i, j, iplot, iplot2, it2
  REAL(KIND=DP) :: timein, timeout, timein2, time_in, time_out, dtime_out, layloop
  REAL(KIND=DP), DIMENSION(2,10000) :: time_lay
  REAL(KIND=DP), DIMENSION(1) :: dturb
  REAL(KIND=DP) :: lay_full, time_phase, XHwtH2, XHwtH2_prev, grab
  REAL(KIND=DP) :: Nsites_int, Nlay_max, Nlayer, Nsteps_prest, Nsteps_plot
  INTEGER(KIND=LONG) :: checkcons, iplot_lay, iplot_lay2, checklay
  REAL(KIND=DP) :: A, B, C, D, E, time_min, Xsurf, Xmant
  REAL(KIND=DP) :: dNgr, dNacc, dNev, div_dNgr, try_div
  CHARACTER(LEN=8) :: actualdate
  CHARACTER(LEN=10) :: actualtime
  CHARACTER(len=300) :: direct_logtime, test_read
  CHARACTER(len=100) :: arg
  INTEGER(KIND=LONG) :: logtime

!  TYPE (VODE_OPTS) :: OPTIONS

!!$  TYPE(outputparam) :: out
!!$  TYPE(gridini) :: gi
!!$  TYPE(gridparam) :: gp
!!$  TYPE(spatialevol) :: spat
!!$  TYPE(eckartmod) :: ec
!!$  TYPE(COnetwork) :: con
!!$  TYPE(outputdirectory) :: di
!!$  TYPE(varode) :: ode
!!$  TYPE(varode) :: ode2

!!!------------------------------------------------
!!!------------------------------------------------
!!! call input files and create directories
!!!------------------------------------------------
!!!------------------------------------------------


  ! read input files
  CALL INPUT_FILE ! in INPUT module
  CALL INPUT_PARAMETERS ! in INPUT module

  ! read input arguments from terminal
  CALL READ_TERMINAL

  gi%Nmodels = 1

  ! read date and time and create output directory
  CALL DATE_AND_TIME(di%date,di%realtime)
  di%date_ini = di%date; di%time_ini = di%realtime
  di%realtime2 = di%realtime

  IF (inp%chdata .eq. 1) CALL OUTPUT_DIRECTORY ! in OUTPUT module

!!!---------------------------------------------------------------
!!!---------------------------------------------------------------
!!! allocate memory for species, reaction, and grid vectors
!!!---------------------------------------------------------------
!!!---------------------------------------------------------------
  
  CALL ALLO_SPECIES ! in ALLOCATION module
  CALL ALLO_REACTIONS ! in ALLOCATION module

 
!!!-----------------------
!!!-----------------------
!!! MODEL GRID LOOP
!!!-----------------------
!!!-----------------------

  DO im=1,gi%Nmodels 

     !!---------------------------------------------------
     !!---------------------------------------------------
     !! read input files: species, reactions, data
     !!---------------------------------------------------
     !!---------------------------------------------------

     !------------------
     ! read input files
     !------------------

     CALL INPUT_SPECIES ! in INPUT module
     CALL INPUT_REACTIONS ! in INPUT module
     IF (inp%chreacrates .eq. 1) THEN
        CALL INPUT_RATES ! in INPUT module
     END IF
     IF (TRIM(inp%chinphys) .eq. 'yes') THEN
        IF (inp%chgrid .eq. 3) &
             inp%inputphys = TRIM(spat%folder)//'/'//TRIM(spat%file(im))//'.dat'
        CALL INPUT_PHYSPARAM ! in INPUT module
     END IF

     CALL INPUT_SELFSHIELDING ! in INPUT module
     CALL INPUT_ECKART ! in INPUT module
     CALL INPUT_CONETWORK ! in INPUT module
     
     IF (inp%Nspgr .gt. 0) THEN
        CALL INPUT_ENERGIES ! in INPUT module
        CALL INPUT_ENTHALPY
     END IF



     !!----------------------------------------------------------------
     !!----------------------------------------------------------------
     !! compute the initial conditions and initialize the vectors
     !!----------------------------------------------------------------
     !!----------------------------------------------------------------
     
     !---------------------------------------------------------------
     ! create directories for each model and write the grid log file
     !---------------------------------------------------------------


     inp%chgrid = 1
     IF (inp%chgrid .eq. 1) THEN
        direct_logtime=trim(di%direct_mod)//"/logtime.out"
     !   direct_logtime=trim(di%direct_mod)//'/'//trim(spat%file(im))//"/logtime.out"
     END IF
     
     IF (inp%chdata .eq. 1) THEN
        logtime = 258
        OPEN(logtime,file=TRIM(direct_logtime), status='REPLACE',&
             access='SEQUENTIAL',form='FORMATTED',action='READWRITE')
     END IF

     !-------------------------------------------------------------------------
     ! compute the physical/numerical parameters at the beginning of the phase
     !-------------------------------------------------------------------------

     CALL ALLO_PHYS ! in ALLOCATION module
     out%Ninilayer = 0.; ph%nH = inp%nHini
     
     CALL ALLO_GRAIN(1) ! in ALLOCATION module
     CALL GRAIN_PROP(ns%grinput,grab) ! in PHYSPROP module
     CALL PHYS_PROP(1,ph%time) ! in PHYSPROP module
     CALL DEALLO_GRAIN ! in ALLOCATION module


     !----------------------------------------------
     ! compute time for outputs
     !----------------------------------------------

     CALL COMPUTE_NSTEPS ! in NUMPROP module
     CALL ALLO_PLOT ! in ALLOCATION module
     CALL COMPUTE_TIME ! in NUMPROP module

     
     !---------------------------
     ! allocate vectors
     !---------------------------

     CALL ALLO_GRAIN(out%Ninilayer) ! in ALLOCATION module
     CALL ALLO_ABU ! in ALLOCATION module
     CALL GRAIN_PROP(ns%grinput,grab) ! in PHYSPROP module
     IF (ns%G0 .gt. 0) sp%Xini(ns%G0) = grab


     !---------------------------------------
     ! Set initial options for DVODE routine
     !---------------------------------------

     CALL INI_ODEPARAMS ! in NUMPROP module
     

     !----------------------------------------------------------------
     ! compute the reaction probabilities
     !----------------------------------------------------------------

     CALL CHEM_DESORPTION
     CALL REAC_PROBAS
     
     !----------------------------------------------------------------------
     ! initialize the vectors at the beginning of each phase
     !----------------------------------------------------------------------

     ! integers for loops
     iplot = 1; iplot2 = 1; iplot_lay = 1; iplot_lay2 = 1
     it = 1; il = 1
     div_dNgr=1d0; try_div = 0d0; ph%time = 0d0
     timein = 0d0; timeout = 0d0; time_lay = 0d0
     layloop = 0d0; ERR_MESS = ' '; checkcons = 1

     CALL INI_PHYSCHEMPARAMS(iplot,il) ! in NUMPROP module

     CALL PHYS_PROP(il,ph%time) ! in PHYSPROP module
     CALL GRAIN_EVOLUTION(il)!,ab%Ptotml(il)/gr%Nstot(il)) ! in PHYSPROP module

     CALL REAC_RATES(il,ode%Y) ! in CHEMPROP module


     !------------------------------------------------------------------
     ! write the input parameters in the terminal
     !------------------------------------------------------------------

     ! CALL TERM_INPUT(im) ! in OUTPUT module

     ! write the progress of the computation into a log file
     IF (inp%chdata .eq. 1) THEN
        CALL DATE_AND_TIME(actualdate,actualtime)
        WRITE(logtime,'(A8,1x,A4,A5,I6,A7,ES9.3,A7,ES9.3,A10,ES9.3,A2)') &
             actualdate,actualtime,', it=',it,'time=',ph%time,'yr, nH=',ph%nH,&
             ' cm-3, Tg=',ph%Tg,' K' 
        CLOSE(logtime)
     END IF

!!!----------------------
!!!----------------------
!!! TEMPORAL LOOP
!!!----------------------
!!!----------------------

     chemiloop : DO WHILE (ph%time .le. inp%tmax)

        !------------------------------------------------------------------
        ! condition for stopping calculation, depending on 3-phase/2-phase
        !------------------------------------------------------------------

        IF(inp%chlayer .eq. 1) lay_full = 1d40
        IF(inp%chlayer .eq. 2) lay_full = gr%Nstot(il)
        lay_full = 1d40
        
!!!--------------------
!!!--------------------
!!! LAYER LOOP
!!!--------------------
!!!--------------------

        layer : DO WHILE (layloop .le. lay_full .and. layloop .ge. 0d0)


           !--------------------------------------------------
           ! compute physical parameters and rates
           !--------------------------------------------------

           !CALL PHYS_PROP(il,ph%time,inp,iph,sp,ns,gr,ph)
           !CALL GRAIN_EVOLUTION(il,inp,gr,ph,ab%Ptotml(il)/gr%Nstot(il))
           !CALL REAC_RATES(il,inp,ph,sp,ns,re,nr,gr,ss,ab,ode%Y)


           !---------------------------------------------------------
           ! compute the timestep based on rate evaluation
           !---------------------------------------------------------

           CALL COMPUTE_DTIME(it,iplot,il,timeout,timein)
           
           !IF (it == 1) CALL TERM_LOG_LAYER(1,it,il) ! in OUTPUT module 
           
           !-------------------------------
           ! call ODE integrator
           !-------------------------------
           
           timein2 = timein
           timein = 0d0
           timeout = timeout-timein2

           itest = 0
           
           ! loop to re-run calculation if problem in the calculation
           DO WHILE (timein .lt. timeout .and. itest .lt. 1000) 
              
              ! input parameters for DLSODES
              CALL COMPUTE_IAJA(it) ! in NUMPROP module
              ode%IWORK(30+1:30+inp%Neq+1)=ode%IA(1:inp%Neq+1)
              ode%IWORK(31+inp%Neq+1:31+inp%Neq+ode%nonzerocount) = &
                   ode%JA(1:ode%nonzerocount)
              DO is=1,inp%Neq
                 !ode%ATOL(is)=max(ode%ATOLini(1),inp%ATOL0/1d4*ode%Y(is))
                 ode%ATOL(is)=ode%ATOLini(1)
              END DO
              ode%ISTATE = 1!; itest = 0

              ! decrease timestep if solver failed to converge
              IF (itest .gt. 0) THEN
                 timein = 0.
                 timeout = timeout/10.
                 ph%dtime = ph%dtime/10.
                 ode%Y(1:inp%Neq) = ode%Ysave(1:inp%Neq)
              ELSE IF (itest .eq. 0) THEN
                 ode%Ysave(1:inp%Neq) = ode%Y(1:inp%Neq) 
              END IF
              IF (inp%chH .eq. 1 .and. ns%H .ne. 0) ode%Y(ns%H) = ode%Yini(ns%H)
              
              write(6,'(A3,I7,A8,I7,A7,ES15.9,A11,ES9.3,A8,ES9.3,A10,ES9.3,'//&
                   'A7,ES9.3,A12,ES9.3,A12,ES9.3)') &
                   'it=',it,', iplot=',iplot,', time=',ph%time,' yr, dtime=',&
                   ph%dtime/Rys,' yr, nH=',ph%nH,' cm-3, Tg=',ph%Tg,' K, Av=',&
                   ph%Av,' mag, N(ML)=',ab%Nlay,', Xice(H2O)=', ab%Xice(ns%H2O)
              !write(6,*) sp%Eb(ns%JCO2), sp%Eb(ns%JH2O) !ph%rad, out%rad!(iplot2)
              ! DLSODES
              CALL DLSODES(Fode,inp%Neq,ode%Y,timein,timeout,ode%ITOL,&
                   ode%RTOL,ode%ATOL,ode%ITASK,ode%ISTATE,ode%IOPT,&
                   ode%RWORK,ode%LRW,ode%IWORK,ode%LIW,JACode,ode%MF)

              itest = itest + 1

           END DO

           IF (itest .ge. 1000) THEN
              WRITE(6,*) 'Calculation did not converge after 1000 attempts, program stopped'
              STOP
           END IF
           
           timein = timein2 + timein
           timeout = timein2 + timeout
           
           ! write the progress of the computation into a log file
           IF (inp%chdata .eq. 1) THEN
              CALL DATE_AND_TIME(actualdate,actualtime)
              OPEN(logtime,file=TRIM(direct_logtime), status='OLD',&
                   access='SEQUENTIAL',form='FORMATTED',action='READWRITE',position='APPEND')
              WRITE(logtime,'(A8,1x,A4,A5,I6,A8,ES9.3,A8,ES9.3,A10,ES9.3,A14,ES9.3)') &
                   actualdate,actualtime,', it=',it,'time=',ph%time,' yr, nH=',ph%nH,&
                   ' cm-3, Tg=',ph%Tg,' K, Xgas(H2O)=', ab%Xgas(ns%H2O)
              CLOSE(logtime)
           END IF

           !----------------------------------------------------
           ! save the parameters at the previous time
           !----------------------------------------------------

           CALL SAVE_PREVPARAMS(1,il) ! in NUMPROP module
           it = it + 1
           
           !-----------------------------------------------------------------
           ! if calculation does not converge, need to modify ATOL
           !-----------------------------------------------------------------

           IF (ERR_MESS .ne. ' ' .or. ode%ISTATE .eq. -4) THEN
              WRITE(6,*) 'ERROR MESSAGE: ', ERR_MESS
              write(6,*) ode%ISTATE
              EXIT chemiloop
           END IF

           !---------------------------------------------
           ! conversion and save of parameters
           !---------------------------------------------

           ! compute integration time
           ph%time = timeout/Rys

           
           ! compute physical conditions for new t
           CALL GRAIN_EVOLUTION(il)
           CALL PHYS_PROP(il,ph%time) ! in PHYSPROP module
           CALL REAC_RATES(il,ode%Y) ! in CHEMPROP module

           ! convert chemical parameters
           CALL CONVERT_CHEMPARAMS(il,layloop) ! in NUMPROP module


           !-----------------------------------------------------------------
           ! save output parameters at specific timesteps
           !-----------------------------------------------------------------

           IF (abs(ph%time-out%timefix(iplot)) .lt. 1/Rys .or. &
                int(ab%Nlay) .ne. int(ab%Nlay_prev)) THEN
              CALL SAVE_OUTPARAMS(1,il,iplot,iplot2) ! in NUMPROP module
           END IF


           !---------------------------------------------------------------
           ! check if the system is conservative by looking 
           ! the elemental abundance and exit the layer loop if not
           !---------------------------------------------------------------

           !CALL CHECK_ELEMS(it,checkcons) ! in CHEMPROP module
           !IF (checkcons .eq. 2) STOP
           checkcons = 1

           IF (ab%Ptotml(il) .le. 1d0 .and. ab%Ptotml_prev(il) .gt. 1d0) THEN
              EXIT layer
           END IF

           !---------------------------------------------------------------
           ! leave the loop when t = tmax
           !---------------------------------------------------------------

           !IF (iplot .eq. 731) EXIT chemiloop
           IF (ph%time .ge. inp%tmax*0.9999999) THEN
           !IF (ph%time .ge. 1e6) THEN
              WRITE(6,'(A9,ES16.10)') 't(fin) = ', ph%time
              EXIT chemiloop
           END IF

        END DO layer



     END DO chemiloop


     !!----------------------------------------------------------------
     !!----------------------------------------------------------------
     !! save the number of species in the last ML on the surface 
     !!----------------------------------------------------------------
     !!----------------------------------------------------------------

     ! save the population at the end of the calculation
     out%Nrealsteps = iplot2-1

     ! write parameters into output vectors at tmax
     CALL SAVE_OUTPARAMS(2,il,iplot,iplot2) ! in NUMPROP module


     !!-----------------------------------------------------------------
     !!-----------------------------------------------------------------
     !! order the reaction rates from its integrated contribution
     !!-----------------------------------------------------------------
     !!-----------------------------------------------------------------

     IF (inp%chreacrates .eq. 1) CALL ORDER_RATES(ode%Y) ! in CHEMPROP module


     !!------------------------------------------------------------------------
     !!------------------------------------------------------------------------
     !! write results in the terminal/files at the end of calculations
     !!------------------------------------------------------------------------
     !!------------------------------------------------------------------------

     ! in terminal
     CALL FIN_TERM_LOG(it,il) ! in OUTPUT module

     ! write data and log in output files
     IF (inp%chdata .eq. 1) THEN
        CALL OUTPUT_LOG(it,il,im) ! in OUTPUT module
        CALL WRITE_BINARY(it,il,im) ! in OUTPUT module
        IF (inp%chreacrates .eq. 1) &
             CALL OUTPUT_RATES(it,il,im) ! in OUTPUT module
     END IF
    

     !----------------------------
     !----------------------------
     ! deallocate vectors
     !----------------------------
     !----------------------------

     ! deallocate local vectors
     CALL DEALLO_VARED ! in ALLOCATION module

     ! deallocate physical/chemical parameters
     CALL DEALLO_ABU ! in ALLOCATION module
     CALL DEALLO_PHYS ! in ALLOCATION module
     IF (trim(inp%chinphys) .eq. 'yes') CALL DEALLO_INPHYS ! in ALLOCATION module
     CALL DEALLO_GRAIN ! in ALLOCATION module
     CALL DEALLO_PLOT ! in ALLOCATION module

    DEALLOCATE(inp%randomEb, inp%randomEa)
     
  END DO
  

  ! deallocate grid, species, and reactions variables
  CALL DEALLO_SPECIES ! in ALLOCATION module
  CALL DEALLO_REACTIONS ! in ALLOCATION module
  


END PROGRAM
