MODULE GRIDS

  USE VARIABLES

  IMPLICIT NONE

  !----------------------------------------------------------------------------
  ! This module contains several subroutines that: 
  ! compute and assign the values of the free parameters 
  ! read the final abundances of another model grid as initial abundances
  !----------------------------------------------------------------------------

CONTAINS 

  SUBROUTINE COMPUTE_GRID

    !------------------------------------------------------------
    ! compute the input parameters of each run of the model grid
    !------------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, gi, gp

    INTEGER(KIND=LONG) :: idens, itemp, idistrib, itreat, iratio, iporo, im
    INTEGER(KIND=LONG) :: isize, iea, imodel, ioxab, iEbH, iUV, iop, iAv, izeta
    REAL(KIND=DP) :: lognH, logoxab

    DO idens=1,gi%NnH2
    DO itemp=1,gi%NT2
    DO iUV=1,gi%NUV2
    DO izeta=1,gi%Nzeta2
    DO iAv=1,gi%NAv2
    DO idistrib=1,gi%Nacst2
    DO itreat=1,gi%Ntreat2
    DO iratio=1,gi%Nenratio2
    DO iEbH=1,gi%NEbH2
    DO iporo=1,gi%Nporo2
    DO isize=1,gi%Nds2
    DO iea=1,gi%NEa2
    DO ioxab=1,gi%NXO2
    DO iop=1,gi%Nop2

       im= iop+&
            (ioxab-1)*gi%Nop2+&
            (iea-1)*gi%NXO2*gi%Nop2+&
            (isize-1)*gi%NEa2*gi%NXO2*gi%Nop2+&
            (iporo-1)*gi%Nds2*gi%NEa2*gi%NXO2*gi%Nop2+&
            (iEbH-1)*gi%Nporo2*gi%Nds2*gi%NEa2*gi%NXO2*gi%Nop2+&
            (iratio-1)*gi%NEbH2*gi%Nporo2*gi%Nds2*gi%NEa2*gi%NXO2*gi%Nop2+&
            (itreat-1)*gi%Nenratio2*gi%NEbH2*gi%Nporo2*gi%Nds2*gi%NEa2*&
            gi%NXO2*gi%Nop2+&
            (idistrib-1)*gi%Ntreat2*gi%Nenratio2*gi%NEbH2*gi%Nporo2*gi%Nds2*&
            gi%NEa2*gi%NXO2*gi%Nop2+&
            (iAv-1)*gi%Nacst2*gi%Ntreat2*gi%Nenratio2*gi%NEbH2*gi%Nporo2*&
            gi%Nds2*gi%NEa2*gi%NXO2*gi%Nop2+&
            (izeta-1)*gi%NAv2*gi%Nacst2*gi%Ntreat2*gi%Nenratio2*gi%NEbH2*&
            gi%Nporo2*gi%Nds2*gi%NEa2*gi%NXO2*gi%Nop2+&
            (iUV-1)*gi%Nzeta2*gi%NAv2*gi%Nacst2*gi%Ntreat2*gi%Nenratio2*&
            gi%NEbH2*gi%Nporo2*gi%Nds2*gi%NEa2*gi%NXO2*gi%Nop2+&
            (itemp-1)*gi%NUV2*gi%Nzeta2*gi%NAv2*gi%Nacst2*gi%Ntreat2*&
            gi%Nenratio2*gi%NEbH2*gi%Nporo2*gi%Nds2*gi%NEa2*gi%NXO2*gi%Nop2+&
            (idens-1)*gi%NT2*gi%NUV2*gi%Nzeta2*gi%NAv2*gi%Nacst2*gi%Ntreat2*&
            gi%Nenratio2*gi%NEbH2*gi%Nporo2*gi%Nds2*gi%NEa2*gi%NXO2*gi%Nop2
     
       ! density
       IF (gi%NnH .eq. 0) THEN
          gp%nH(im) = inp%nHini
       ELSE IF (gi%NnH .eq. 1) THEN
          gp%nH(im) = gi%nHlow
       ELSE IF (gi%NnH .gt. 1) THEN
          lognH = log10(gi%nHlow)+(idens-1)*(log10(gi%nHup)-log10(gi%nHlow))&
               /(gi%NnH-1)
          gp%nH(im) = 10**lognH
       END IF
       ! temperature
       IF (gi%NT .eq. 0) THEN 
          gp%T(im) = inp%Tgini
       ELSE IF (gi%NT .eq. 1) THEN 
          gp%T(im) = gi%Tlow
       ELSE IF (gi%NT .gt. 1) THEN 
          gp%T(im) = gi%Tlow+(itemp-1)*(gi%Tup-gi%Tlow)/(gi%NT-1)
       END IF
       ! UV flux induced by CR ionization of H2
       IF (gi%NUV .eq. 0) THEN
          gp%UV(im) = inp%UVfluxCRini
       ELSE IF (gi%NUV .eq. 1) THEN
          gp%UV(im) = gi%UVlow
       ELSE IF (gi%NUV .gt. 1) THEN
          gp%UV(im) = 10**(log10(gi%UVlow)+(iUV-1)*&
               (log10(gi%UVup)-log10(gi%UVlow))/(gi%NUV-1))
       END IF
       ! CR ionization rate
       IF (gi%Nzeta .eq. 0) THEN
          gp%zeta(im) = inp%zetaini
       ELSE IF (gi%Nzeta .eq. 1) THEN
          gp%zeta(im) = gi%zetalow
       ELSE IF (gi%Nzeta .gt. 1) THEN
          gp%zeta(im) = 10**(log10(gi%zetalow)+(izeta-1)*&
               (log10(gi%zetaup)-log10(gi%zetalow))/(gi%Nzeta-1))
       END IF
       ! Visual extinction
       IF (gi%NAv .eq. 0) THEN
          gp%Av(im) = inp%Avini
       ELSE IF (gi%NAv .eq. 1) THEN
          gp%Av(im) = gi%Avlow
       ELSE IF (gi%NAv .gt. 1) THEN
          gp%Av(im) = gi%Avlow+(iAv-1)*(gi%Avup-gi%Avlow)/(gi%NAv-1)
       END IF
       ! grain size
       IF (gi%Nacst .eq. 0) THEN
          gp%acst(im) = inp%acst
       ELSE IF (gi%Nacst .eq. 1) THEN
          gp%acst(im) = gi%acstlow
       ELSE IF (gi%Nacst .gt. 1) THEN
          gp%acst(im) = gi%acstlow+(idistrib-1)*(gi%acstup-gi%acstlow)&
               /(gi%Nacst-1)
       END IF
       ! mantle behaviour
       IF (gi%choicetreat .eq. 0) THEN 
          gp%treatment2(im) =  inp%chlayer
          IF (inp%chlayer .eq. 1) THEN
             gp%treatment(im) = 'no'
          ELSE IF (inp%chlayer .eq. 2) THEN
             gp%treatment(im) = 'yes'
          END IF
       ELSE IF (gi%choicetreat .eq. 1) THEN 
          gp%treatment(im) = 'no'
          gp%treatment2(im) =  1d0
       ELSE IF (gi%choicetreat .eq. 2) THEN
          gp%treatment(im) = 'yes '
          gp%treatment2(im) =  2d0
       ELSE IF (gi%choicetreat .eq. 3) THEN
          IF (itreat .eq. 1) THEN
             gp%treatment(im) = 'yes'
             gp%treatment2(im) = 2d0
          ELSE IF (itreat .eq. 2) THEN
             gp%treatment(im) = 'no '
             gp%treatment2(im) = 1d0
          END IF
       END IF
       ! Ed/Eb ratio
       IF (gi%Nenratio .eq. 0) THEN
          gp%enratio(im) = inp%REdH
       ELSE IF (gi%Nenratio .eq. 1) THEN
          gp%enratio(im) = gi%enratiolow
       ELSE IF (gi%Nenratio .gt. 1) THEN
          gp%enratio(im) = gi%enratiolow+(iratio-1)*&
               (gi%enratioup-gi%enratiolow)/(gi%Nenratio-1)
       END IF
       ! binding energy of volatile species
       IF (gi%NEbH .eq. 0) THEN
          gp%EbH(im) = 0d0
       ELSE IF (gi%NEbH .eq. 1) THEN
          gp%EbH(im) = gi%EbHlow
       ELSE IF (gi%NEbH .gt. 1) THEN
          gp%EbH(im) = gi%EbHlow+(iEbH-1)*(gi%EbHup-gi%EbHlow)/(gi%NEbH-1)
       END IF
       ! porosity factor
       IF (gi%Nporo .eq. 0) THEN
          gp%Fpor(im) = inp%Fporsur
       ELSE IF (gi%Nporo .eq. 1) THEN
          gp%Fpor(im) = gi%porolow
       ELSE IF (gi%Nporo .gt. 1) THEN
          gp%Fpor(im) = gi%porolow+(iporo-1)*(gi%poroup-gi%porolow)/&
               (gi%Nporo-1)
       END IF
       ! site size
       IF (gi%Nds .eq. 0) THEN
          gp%ds(im) = inp%ds
       ELSE IF (gi%Nds .eq. 1) THEN
          gp%ds(im) = gi%dslow
       ELSE IF (gi%Nds .gt. 1) THEN
          gp%ds(im) = gi%dslow+(isize-1)*(gi%dsup-gi%dslow)&
               /(gi%Nds-1)
       END IF
       ! activation energy of the H+CO reaction
       IF (gi%NEa .eq. 0) THEN
          gp%Ea(im) = 0d0! reac_C(num_reac_COH)
       ELSE IF (gi%NEa .eq. 1) THEN
          gp%Ea(im) = gi%Ealow
       ELSE
          gp%Ea(im) = gi%Ealow+(iea-1)*(gi%Eaup-gi%Ealow)/(gi%NEa-1)
       END IF
       ! initial oxygen abundance
       IF (gi%NXO .eq. 0) THEN
          gp%XO(im) = 0d0!spec_abu_ini(num_O)
       ELSE IF (gi%NXO .eq. 1) THEN
          gp%XO(im) = gi%XOlow
       ELSE IF (gi%NXO .gt. 1) THEN
          logoxab = log10(gi%XOlow)+(ioxab-1)*(log10(gi%XOup)-log10(gi%XOlow))&
               /(gi%NXO-1)
          gp%XO(im) = 10**logoxab
       END IF
       ! ortho/para ratio of H2
       IF (gi%Nop .eq. 0) THEN
          gp%op(im) = 1
       ELSE IF (gi%Nop .eq. 1) THEN
          gp%op(im) = gi%oplow
       ELSE IF (gi%Nop .gt. 1) THEN
          gp%op(im) = 10**(log10(gi%oplow)+(iop-1)*&
               (log10(gi%opup)-log10(gi%oplow))/(gi%Nop-1))
       END IF

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE COMPUTE_GRID


  SUBROUTINE INPUT_GRID_INPUTGRIDLOG

    !------------------------------------------------------------------------
    ! read the log file of an other model grid to get the initial abundances
    !------------------------------------------------------------------------

    USE SHARED_VARIABLES, only: gi, gp 

    CHARACTER(len=75) :: toto_abuini
    CHARACTER(len=6) :: toto_abuini2
    INTEGER(KIND=LONG) :: i, log_inputgrid
    CHARACTER(len=130) :: direct_inputgrid

    ! open log file
    direct_inputgrid = trim(gi%input_abuini)//'/gridlog.out'

    log_inputgrid = 414
    OPEN(log_inputgrid,file=TRIM(direct_inputgrid), status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    ! read log file
    DO i=1,6+gi%Nparam
       READ(log_inputgrid,*) 
    END DO

    DO i=1,gi%Nmodels
       READ(log_inputgrid,'(A75)') toto_abuini
       gp%Xini_file(i) = toto_abuini(10:75)
    END DO

    CLOSE(log_inputgrid)


  END SUBROUTINE INPUT_GRID_INPUTGRIDLOG


  SUBROUTINE INPUT_GRID_ABUINI(num_model)

    !----------------------------------------------------------------------
    ! read the final abundances of an other model grid to get the 
    ! initial abundances
    !----------------------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, gi, gp

    INTEGER(KIND=LONG) :: i, Nrealsteps
    CHARACTER(len=200) :: direct_Xgas
    CHARACTER(len=200) :: direct_Xice
    CHARACTER(len=200) :: direct_inp
    INTEGER(KIND=LONG) :: file_Xgas2
    INTEGER(KIND=LONG) :: file_Xice2
    INTEGER(KIND=LONG) :: file_inp
    REAL(KIND=DP) :: abuOi, abuO2i
    INTEGER(KIND=LONG), INTENT(in) :: num_model
    REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: Xgas2
    REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: Xice2
    
    direct_inp = trim(gi%input_abuini)//'/'//trim(gp%Xini_file(num_model))&
         //'/input.out'
    file_inp = 414
    OPEN(file_inp,file=TRIM(direct_inp), status='OLD',&
         form='UNFORMATTED',action='READ')
    READ(file_inp) Nrealsteps
    CLOSE(file_inp)

    ALLOCATE(Xgas2(inp%Nspgas,Nrealsteps,1), stat=ok)
    ALLOCATE(Xice2(inp%Nspgr,Nrealsteps,1), stat=ok)

    ! open binary files
    direct_Xgas = trim(gi%input_abuini)//'/'//trim(gp%Xini_file(num_model))&
         //'/Xgas.out'
    direct_Xice = trim(gi%input_abuini)//'/'//&
         trim(gp%Xini_file(num_model))//'/Xice.out'
    file_Xgas2 = 415
    file_Xice2 = 416
    OPEN(file_Xgas2,file=TRIM(direct_Xgas), status='OLD',&
         form='UNFORMATTED',action='READ')
    OPEN(file_Xice2,file=TRIM(direct_Xice), status='OLD',&
         form='UNFORMATTED',action='READ')

    ! read binary files and use them for gas and ice initial abundances
    READ(file_Xgas2) Xgas2
    gp%Xgasini(num_model,1:inp%Nspgas) = Xgas2(1:inp%Nspgas,Nrealsteps,1)
    READ(file_Xice2) Xice2
    gp%Xiceini(num_model,1:inp%Nspgr) = Xice2(1:inp%Nspgr,Nrealsteps,1)
    
    CLOSE(file_Xgas2)
    CLOSE(file_Xice2)
    
    DEALLOCATE(Xgas2,Xice2)

  END SUBROUTINE INPUT_GRID_ABUINI


  SUBROUTINE ASSIGN_PARAMS(im)

    USE SHARED_VARIABLES, only: gp, gi, inp, sp, ns, re, nr 

    INTEGER(KIND=LONG), INTENT(in) :: im
    INTEGER(KIND=LONG) :: ic
    INTEGER(KIND=LONG) :: log_inputgrid
    REAL(KIND=DP) :: RrhonH, Mpr2, G2, RcmAU, rhoflat, Ryrsec

    G2 = 6.67d-11*1d6*1d-3
    Mpr2 = 1.67d-24                
    RcmAU = 149598000d5          
    Ryrsec = 31536000d0           
    RrhonH = 1.36          

    !----------------------------------------------------------------------
    ! assign the input parameters chosen in model_grid.in 
    !----------------------------------------------------------------------

    ! density
    IF (gi%NnH .ne. 0) THEN
       inp%nHini = gp%nH(im)
       WRITE(gi%char_nH,'(ES7.1)') gp%nH(im)

       ! change tmax to tff from nH
       rhoflat = inp%nHini * RrhonH*Mpr2
       inp%tmax = min(inp%tmaxini,sqrt(3*Pi/(32*G2*rhoflat))/Ryrsec)
    END IF
    ! temperature
    IF (gi%NT .ne. 0) THEN
       inp%Tini = gp%T(im)
       inp%Tgini = inp%Tini
       inp%Tdini = inp%Tini
       WRITE(gi%char_T,'(ES7.1)') gp%T(im)
    END IF
    ! UV flux induced by CRs
    IF (gi%NUV .ne. 0) THEN
       inp%UVfluxCRini = gp%UV(im)
       WRITE(gi%char_UV,'(ES7.1)') gp%UV(im)
    END IF
    ! CR ionization rate
    IF (gi%Nzeta .ne. 0) THEN
       inp%zetaini = gp%zeta(im)
       WRITE(gi%char_zeta,'(ES7.1)') gp%zeta(im)
    END IF
    ! Visual extinction
    IF (gi%NAv .ne. 0) THEN
       inp%Avini = gp%Av(im)
       WRITE(gi%char_Av,'(ES7.1)') gp%Av(im)
    END IF
    ! grain size
    IF (gi%Nacst .ne. 0) THEN
       inp%acst = gp%acst(im)
       WRITE(gi%char_dist,'(ES7.1)') gp%acst(im)
    END IF
    ! mantle behaviour
    IF (gi%choicetreat .ne. 0) THEN
       inp%chlayer = gp%treatment2(im)
       inp%chlayer2 = gp%treatment(im)
       inp%Neq = inp%Nspecies+1
       IF (inp%chlayer .eq. 2 .and. inp%chmulti .eq. 2) THEN
          inp%Neq = inp%Nspgas+2*inp%Nspgr+1
       END IF
       WRITE(gi%char_treat,'(ES7.1)') gp%treatment2(im)
    END IF
    ! energy ratio
    IF (gi%Nenratio .ne. 0) THEN
       inp%REdH = gp%enratio(im)
       !inp%REdoth = gp%enratio(im)
       inp%REdbulkH = gp%enratio(im)*2.
       !inp%REdbulkoth = gp%enratio(im)*2.
       WRITE(gi%char_en,'(ES7.1)') gp%enratio(im)
    END IF
    ! binding energy
    IF (gi%NEbH .gt. 0) THEN
       ! binding energy of H relative to ASW
       sp%Eb(ns%JH) = gp%EbH(im)
       sp%Eb_ini(ns%JH) = gp%EbH(im)
       sp%Ed(ns%JH) = inp%REdH*sp%Eb(ns%JH)
       IF (inp%Nspgr2 .eq. inp%Nspgr/2d0) THEN
          sp%Eb(ns%JH+inp%Nspgr) = gp%EbH(im)
          sp%Eb_ini(ns%JH+inp%Nspgr) = gp%EbH(im)
          sp%Ed(ns%JH+inp%Nspgr) = inp%REdH*sp%Eb(ns%JH)
       END IF
       ! binding energy of D relative to ASW
       sp%Eb(ns%JD) = gp%EbH(im)!+20
       sp%Eb_ini(ns%JD) = sp%Eb(ns%JD)
       sp%Ed(ns%JD) = inp%REdH*sp%Eb(ns%JD)
       IF (inp%Nspgr2 .eq. inp%Nspgr/2d0) THEN
          sp%Eb(ns%JD+inp%Nspgr) = sp%Eb(ns%JD)
          sp%Eb_ini(ns%JD+inp%Nspgr) = sp%Eb(ns%JD)
          sp%Ed(ns%JD+inp%Nspgr) = inp%REdH*sp%Eb(ns%JD)
       END IF
       ! binding energy of H2 relative to ASW
       sp%Eb(ns%JH2) = gp%EbH(im)
       sp%Eb_ini(ns%JH2) = sp%Eb(ns%JH2)
       sp%Ed(ns%JH2) = inp%REdH*sp%Eb(ns%JH2)
       sp%Eb(ns%JoH2) = gp%EbH(im)
       sp%Eb_ini(ns%JoH2) = sp%Eb(ns%JH2)
       sp%Ed(ns%JoH2) = inp%REdH*sp%Eb(ns%JH2)
       IF (inp%Nspgr2 .eq. inp%Nspgr/2d0) THEN
          sp%Eb(ns%JH2+inp%Nspgr) = sp%Eb(ns%JH2)
          sp%Eb_ini(ns%JH2+inp%Nspgr) = sp%Eb(ns%JH2)
          sp%Ed(ns%JH2+inp%Nspgr) = inp%REdH*sp%Eb(ns%JH2)
          sp%Eb(ns%JoH2+inp%Nspgr) = sp%Eb(ns%JH2)
          sp%Eb_ini(ns%JoH2+inp%Nspgr) = sp%Eb(ns%JH2)
          sp%Ed(ns%JoH2+inp%Nspgr) = inp%REdH*sp%Eb(ns%JH2)
       END IF
       ! binding energy of HD relative to ASW
       sp%Eb(ns%JHD) = gp%EbH(im)!+15
       sp%Eb_ini(ns%JHD) = sp%Eb(ns%JHD)
       sp%Ed(ns%JHD) = inp%REdH*sp%Eb(ns%JHD)
       IF (inp%Nspgr2 .eq. inp%Nspgr/2d0) THEN
          sp%Eb(ns%JHD+inp%Nspgr) = sp%Eb(ns%JHD)
          sp%Eb_ini(ns%JHD+inp%Nspgr) = sp%Eb(ns%JHD)
          sp%Ed(ns%JHD+inp%Nspgr) = inp%REdH*sp%Eb(ns%JHD)
       END IF
       ! binding energy of D2 relative to ASW
       sp%Eb(ns%JD2) = gp%EbH(im)!+30
       sp%Eb_ini(ns%JD2) = sp%Eb(ns%JD2)
       sp%Ed(ns%JD2) = inp%REdH*sp%Eb(ns%JD2)
       sp%Eb(ns%JoD2) = gp%EbH(im)+30
       sp%Eb_ini(ns%JoD2) = sp%Eb(ns%JD2)
       sp%Ed(ns%JoD2) = inp%REdH*sp%Eb(ns%JD2)
       IF (inp%Nspgr2 .eq. inp%Nspgr/2d0) THEN
          sp%Eb(ns%JD2+inp%Nspgr) = sp%Eb(ns%JD2)
          sp%Eb_ini(ns%JD2+inp%Nspgr) = sp%Eb(ns%JD2)
          sp%Ed(ns%JD2+inp%Nspgr) = inp%REdH*sp%Eb(ns%JD2)
          sp%Eb(ns%JoD2+inp%Nspgr) = sp%Eb(ns%JD2)
          sp%Eb_ini(ns%JoD2+inp%Nspgr) = sp%Eb(ns%JD2)
          sp%Ed(ns%JoD2+inp%Nspgr) = inp%REdH*sp%Eb(ns%JD2)
       END IF
       WRITE(gi%char_EbH,'(ES7.1)') gp%EbH(im)
    ELSE IF (gi%NEbH .eq. 0) THEN
       gp%EbH(im) = sp%Eb(ns%JH)
    END IF
    ! porosity factor
    IF (gi%Nporo .ne. 0) THEN
       gp%Fpor = gp%Fpor(im)
       WRITE(gi%char_poro,'(ES7.1)') gp%Fpor(im)
    END IF
    ! site size
    IF (gi%Nds .ne. 0) THEN
       inp%ds = gp%ds(im)
       inp%laywid = inp%ds
       inp%Fnpsur = 1 - inp%Fporsur
       IF (inp%Slat .gt. 2) THEN
          inp%Fedgesur = 4*inp%Fporsur*(inp%Slat-1)/inp%Slat**2
       ELSE IF (inp%Slat .le. 2 .and. inp%Slat .gt. 0) THEN
          inp%Fedgesur = inp%Fporsur
       ELSE IF (inp%Slat .le. 0) THEN
          inp%Fedgesur = 0
       END IF
       WRITE(gi%char_ds,'(ES7.1)') gp%ds(im)
    END IF
    ! CO and H2CO hydrogenation activation energies
    IF (gi%NEa .ne. 0) THEN
       DO ic=1,inp%Nreactions
          IF (trim(re%react1(ic)) .eq. 'JCO' .and. &
               trim(re%react2(ic)) .eq. 'JH') THEN
             re%C(ic) = gp%Ea(im)
          END IF
          !reac_proba(ic) = proba_reac(ic)
       END DO
       WRITE(gi%char_Ea,'(ES7.1)') gp%Ea(im)
    END IF
    ! initial abundance of atomic oxygen
    IF (gi%NXO .ne. 0) THEN
       sp%Xini(ns%O) = gp%XO(im)
       WRITE(gi%char_XO,'(ES7.1)') gp%XO(im)
       !sp%Xini(ns%O2) = (2.56d-4 - sp%Xini(ns%CO) - &
       !     sp%Xini(ns%O))/2d0
    END IF
    ! ortho/para ratio of H2
    IF (gi%Nop .ne. 0) THEN
       sp%Xini(ns%pH2) = 0.5/(1+gp%op(im))
       sp%Xini(ns%oH2) = 0.5 - sp%Xini(ns%pH2)
       inp%H2opr = gp%op(im)
       WRITE(gi%char_op,'(ES7.1)') gp%op(im)
    END IF

  END SUBROUTINE ASSIGN_PARAMS

END MODULE
    
