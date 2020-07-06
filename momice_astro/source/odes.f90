MODULE ODES

  USE VARIABLES
  USE INPUT
  USE PHYSPROP
  USE CHEMPROP

  IMPLICIT NONE

  !---------------------------------------------------------------------------
  ! This module contains the subroutine which computes differential equations
  !---------------------------------------------------------------------------

CONTAINS

!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------
!! DIFFERENTIAL EQUATIONS
!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------

  SUBROUTINE Fode(Neq,t,Y,F)

    !---------------------------------------------------------------
    ! main routine to compute the ODEs
    !---------------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, ns, nr, il, ph, itest
    USE SHARED_VARIABLES2

    INTEGER(KIND=LONG), INTENT(in) :: Neq
    REAL(KIND=DP), DIMENSION(Neq), INTENT(in) :: Y
    REAL(KIND=DP), DIMENSION(Neq), INTENT(out) :: F
    REAL(KIND=DP), INTENT(in) :: t
    REAL(KIND=DP), DIMENSION(inp%Nreactions) :: Rreactot
    REAL(KIND=DP), DIMENSION(Neq+1) :: Rform, Rdest, F2
    !REAL(KIND=DP) :: nstot,nmtot,alphacc,alphades
    INTEGER(KIND=LONG) :: j, ic, ic2, ir1, ir2, ir3, ip1, ip2, ip3, ip4

    ! compute reaction rates
    CALL REAC_RATES(il,Y)

    ! compute derivatives due to "normal" (type<26) reactions
    CALL Freacs(Neq,Y,F)

    ! compute derivatives due to mantle/surface exchange rates
    IF (inp%chlayer .eq. 2 .and. nr%istart(32) .le. inp%Nreactions) &
         CALL Fsurfmant(Neq,Y,F)

    ! fix the H2 opr as constant if needed
    IF (inp%chH2op .eq. 1 .and. ns%oH2 .ne. 0) THEN
       F(ns%H2) = (F(ns%H2) + F(ns%oH2))/(1+inp%H2opr)
       F(ns%oH2) = inp%H2opr*F(ns%pH2)
    END IF

    IF (inp%chH .eq. 1 .and. ns%H .ne. 0) F(ns%H) = 0.
    !write(6,*) t/Rys, ph%dtime/Rys, F(ns%JCO2), Y(ns%JCO2)

    
  END SUBROUTINE Fode


  
  SUBROUTINE Fode0(Neq,t,Y,F)

    !------------------------------------------------------------------------------
    ! compute the ODEs given only by chemical processes (no bulk-surface exchange)
    !------------------------------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, ns, nr, il
    USE SHARED_VARIABLES2

    INTEGER(KIND=LONG), INTENT(in) :: Neq
    REAL(KIND=DP), DIMENSION(Neq), INTENT(in) :: Y
    REAL(KIND=DP), DIMENSION(Neq), INTENT(out) :: F
    REAL(KIND=DP), INTENT(in) :: t
    REAL(KIND=DP), DIMENSION(inp%Nreactions) :: Rreactot
    REAL(KIND=DP), DIMENSION(Neq+1) :: Rform, Rdest, F2
    !REAL(KIND=DP) :: nstot,nmtot,alphacc,alphades
    INTEGER(KIND=LONG) :: j, ic, ic2, ir1, ir2, ir3, ip1, ip2, ip3, ip4

    
    ! compute reaction rates
    CALL REAC_RATES(il,Y)

    ! compute derivatives due to "normal" (type<26) reactions
    CALL Freacs(Neq,Y,F)

    ! fix the H2 opr as constant if needed
    IF (inp%chH2op .eq. 1 .and. ns%oH2 .ne. 0) THEN
       F(ns%H2) = (F(ns%H2) + F(ns%oH2))/(1+inp%H2opr)
       F(ns%oH2) = inp%H2opr*F(ns%pH2)
    END IF

    IF (inp%chH .eq. 1 .and. ns%H .ne. 0) THEN
       F(ns%H) = 0.
    END IF

  END SUBROUTINE Fode0


  
  SUBROUTINE Freacs(Neq,Y,F)

    !---------------------------------------------------------------
    ! subroutine to compute the ODEs from "standard" reactions
    !---------------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, ph, ns, re, nr, ab, sp
    USE SHARED_VARIABLES2

    INTEGER(KIND=LONG), INTENT(in) :: Neq
    REAL(KIND=DP), DIMENSION(Neq), INTENT(in) :: Y
    REAL(KIND=DP), DIMENSION(Neq), INTENT(out) :: F
    REAL(KIND=DP), DIMENSION(inp%Nreactions) :: Rreactot
    REAL(KIND=DP), DIMENSION(Neq+1) :: Rform, Rdest, F2
    REAL(KIND=DP) :: Ytotreact, Rreactot2, COfact
    INTEGER(KIND=LONG) :: j, ic, ic2, ir1, ir2, ir3, ip1, ip2, ip3, ip4, checkCO2reac

    Rreactot = 0d0; Rform = 0d0; Rdest = 0d0; COfact = 1d0; F = 0d0; F2 = 0d0

    ! compute the rate of each reaction of "normal" types (type<28)
    DO ic=1,nr%iend(28)

       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       Ytotreact = 0d0; Rreactot2 = 0d0;  checkCO2reac = 0
   
       ! assign numbers to reactants and products
       ir1 = re%num_react1(ic)
       ir2 = re%num_react2(ic)
       ir3 = re%num_react3(ic)
       ip1 = re%num_prod1(ic)
       ip2 = re%num_prod2(ic)
       ip3 = re%num_prod3(ic)
       ip4 = re%num_prod4(ic)

       ! compute the COfact factor to compute the CO2 formation from Garrod11
       IF ((re%num_react1(ic) .eq. ns%JO .and. &
            re%num_react2(ic) .eq. ns%JH) .or. &
            (re%num_react1(ic) .eq. ns%JO .and. &
            re%num_react2(ic) .eq. ns%JD) .and. &
            nr%reacCO2 .ne. 0) THEN
          COfact = 1-ab%abuCO
          ir3 = inp%Nspecies+1
       ELSE IF ((re%num_react1(ic) .eq. ns%JCO .and. &
            re%num_react2(ic) .eq. ns%JO .and. &
            re%num_react3(ic) .eq. ns%JH) .or. &
            (re%num_react1(ic) .eq. ns%JCO .and. &
            re%num_react2(ic) .eq. ns%JO .and. &
            re%num_react3(ic) .eq. ns%JD) .and. &
            nr%reacCO2 .ne. 0) THEN
          COfact = ab%abuCO
          !ir1 = re%num_react2(ic)
          !ir2 = re%num_react3(ic)
          !ir3 = inp%Nspecies+1
          checkCO2reac = 1
          !write(6,*) COfact,sp%name(ir1),sp%name(ir2)
       END IF

       ! compute formation and destruction rates
       IF (checkCO2reac .eq. 0) THEN
          IF (ir3 .ne. inp%Neq+1) THEN
             Ytotreact = Y(ir1)*Y(ir2)*Y(ir3)*ph%nH*ph%nH
          ELSE IF (ir3 .eq. inp%Neq+1) THEN
             IF (ir2 .eq. inp%Neq+1) THEN 
                !write(6,*) ic, ir1, re%react1(ic)
                Ytotreact = Y(ir1)
             ELSE IF (ir2 .ne. inp%Nspecies+1) THEN
                Ytotreact = Y(ir1)*Y(ir2)*ph%nH
             END IF
          END IF
          
       ELSE IF (checkCO2reac .eq. 1) THEN
          
          Ytotreact = Y(ir2)*Y(ir3)*ph%nH
          
       END IF

       Rreactot(ic) = COfact*ph%Kreac(ic)*Ytotreact
!!$
!!$       ! adjust rate for non-diffusion reactions
!!$       IF (re%type(ic).eq.14) THEN 
!!$          DO isnd=1,inp%Nspnd
!!$             IF (ip1.eq.re%num_ndreact1(isnd)-1) THEN
!!$                Rreactot(ic) = Rreactot(ic)*ab%F(re%num_ndreact2(isnd))
!!$                F2(ip1+1) = F2(ip1+1) + Rreactot2*(1-abu(isnd))
!!$             ELSE IF (ip2.eq.re%num_ndreact1(1)) THEN
!!$                F2(ip2) = F2(ip2) + Rreactot2*(1-abu(isnd))
!!$                F2(ip2) = F2(ip2) + Rreactot2*(abu(isnd))
!!$             END IF
!!$          END DO
!!$       END IF

       Rreactot2 = Rreactot(ic)
       
       F2(ir1) = F2(ir1) - Rreactot2
       F2(ir2) = F2(ir2) - Rreactot2
       F2(ir3) = F2(ir3) - Rreactot2
       F2(ip1) = F2(ip1) + Rreactot2
       F2(ip2) = F2(ip2) + Rreactot2
       F2(ip3) = F2(ip3) + Rreactot2
       F2(ip4) = F2(ip4) + Rreactot2

       
       !IF (re%type(ic).ge.14 .and. re%type(ic).le.18) THEN !ph%time.ge.7.852356346E+03 .and. 
       !   WRITE(6,*) ic, re%type(ic), Rreactot(ic)
       !END IF
       
    END DO

    F(1:Neq) = F2(1:Neq)
    
  END SUBROUTINE Freacs


  SUBROUTINE Fsurfmant(Neq,Y,F)

    !---------------------------------------------------------------
    ! routine to compute the ODEs from surface/mantle exchange
    !---------------------------------------------------------------

    USE SHARED_VARIABLES, only: &
         inp, re, nr, ab, gr, F2tot, nstot, nmtot, alphacc, alphades, il, ph, ns
    USE SHARED_VARIABLES2

    INTEGER(KIND=LONG), INTENT(in) :: Neq
    REAL(KIND=DP), DIMENSION(Neq), INTENT(in) :: Y
    REAL(KIND=DP), DIMENSION(Neq), INTENT(out) :: F
    REAL(KIND=DP), DIMENSION(inp%Nreactions) :: Rsurfmant
    REAL(KIND=DP), DIMENSION(Neq+1) :: Rform, Rdest, F2
    REAL(KIND=DP) :: Ytotreact, dXsurf
    INTEGER(KIND=LONG) :: j, ic, ic2, ir1, ir2, ir3, ip1, ip2, ip3, ip4

    F2 = 0d0
    
    ! intermediate derivatives
    F2tot = 0d0; alphacc = 0d0; alphades = 0d0
    nstot = 0d0; nmtot = 0d0; dXsurf = 0d0
    DO j=inp%Nspgas+1,inp%Nspgas+inp%Nspgr
       F2tot = F2tot + F(j)
       alphacc = alphacc + Y(j)/(gr%Nstot(il)*gr%X*inp%Nsurf)
       alphades = alphades + Y(j+inp%Nspgr)/(gr%Nstot(il)*gr%X*inp%Nsurf)
       nstot = nstot + Y(j)
       nmtot = nmtot + Y(j+inp%Nspgr)
    END DO
    IF (alphacc .ge. 1) alphacc = 1d0
    IF (alphades .ge. 1) alphades = 1d0
    
    Pbulkact = 0d0
    
    ! exchange between the surface layer and the mantle
    DO ic=nr%istart(32),nr%iend(32)
       
       ! assign numbers to reactants and products
       ir1 = re%num_react1(ic)
       ip1 = re%num_prod1(ic)

       ! net accretion: surface->bulk exchange
       IF (F2tot .ge. 0d0 .and. nstot .gt. 0d0) THEN
          Rsurfmant(ic) = alphacc*F2tot*Y(ir1)/nstot
          F(ir1) = F(ir1) -  Rsurfmant(ic)
          F(ip1) = F(ip1) + Rsurfmant(ic)
          dXsurf = dXsurf + Rsurfmant(ic)
          ! net desorption: bulk->surface exchange
       ELSE IF (F2tot .lt. 0d0 .and. nmtot .gt. 0d0) THEN
          Rsurfmant(ic) = alphades*F2tot*Y(ip1)/nmtot
          F(ir1) = F(ir1) - Rsurfmant(ic)
          F(ip1) = F(ip1) + Rsurfmant(ic)        
       END IF
       
    END DO

    
    !IF (ph%time.ge.1.7e4) & ! .and. re%num_react1(ic).eq.ns%JN
    !     write(6,*) Rsurfmant(10276)/Rsurfmant(10275), Y(re%num_react1(10276))/Y(re%num_react1(10275))

    
  END SUBROUTINE Fsurfmant

  
  SUBROUTINE INI_ODEPARAMS

    !---------------------------------------
    ! Set initial options for DVODE routine
    !---------------------------------------
    
    USE SHARED_VARIABLES, only: inp, ph, re, ode 
    USE SHARED_VARIABLES2


    TYPE(varode) :: ode2
    INTEGER(KIND=LONG) :: is

     ode%ITOL = 2
     ode%ITASK = 1
     ode%ISTATE = 1
     ode%IOPT = 1
     ode%MF = 21
     ode%LRW = 10*(20+9*inp%Neq+2000*inp%Neq) 
     ode%LIW = 10*(31+inp%Neq+2000*inp%Neq) 
     CALL ALLO_VARED(inp%Neq)
     ode%IWORK(5) = 5
     !ode%IWORK(6) = 5000. !100000
     ode%IWORK(7) = 2
     ode%RWORK(6) = 3.154E+15
     precision = .false.

     ode%ATOLini(1) = inp%ATOL0
     ode%ATOL(1:inp%Nspecies) = ode%ATOLini(1)
     ode%RTOL = inp%RTOL0

   END SUBROUTINE INI_ODEPARAMS

   
   SUBROUTINE COMPUTE_IAJA(it)

     USE SHARED_VARIABLES, only: inp, ph, re, ode 
     USE SHARED_VARIABLES2
     
     IMPLICIT NONE
     
     INTEGER(KIND=LONG), INTENT(in) :: it
     REAL(KIND=DP), DIMENSION(inp%Neq) :: F
     INTEGER(KIND=LONG) :: ir1, ir2, ir3, ip1, ip2, ip3, ip4, ic, is, is2, inz, inz2
     INTEGER(KIND=LONG), DIMENSION(inp%Neq+1) :: IAN
     INTEGER(KIND=LONG), DIMENSION(inp%Neq*inp%Neq) :: JAN
     INTEGER(KIND=LONG), DIMENSION(inp%Neq+1) :: IA2
     INTEGER(KIND=LONG), DIMENSION(inp%Neq*inp%Neq) :: JA2
     REAL(KIND=DP), DIMENSION(inp%Neq) :: PDJ, countNZ, countJ, NZord
     INTEGER(KIND=LONG), DIMENSION(inp%Neq) :: Iord
     INTEGER(KIND=LONG), DIMENSION(inp%Neq,inp%Neq) :: nonzero2
     REAL(KIND=DP) :: eps
     

     ode%nonzero = 0d0
     ode%IA = 0d0; ode%JA = 0d0
     PDJ = 0d0
     NZord(1:inp%Neq) = 0d0
     inzIord = inzIordini
     inzIord2 = inzIord2ini
     inzNZ = 0d0
          
     ! count the non-zero elements from the jacobian
     ode%nonzero(1:inp%Neq,1:inp%Neq) = 0d0
     CALL Fode(inp%Neq,ph%time*Rys,ode%Y,F)
     DO is2=1,inp%Neq
        CALL JACode(inp%Neq,ph%time*Rys,ode%Y,is2,IAN,JAN,PDJ)
        DO is=1,inp%Neq
           IF (abs(PDJ(is)) .gt. 1d-40) ode%nonzero(is,is2) = 1
        END DO
        ode%nonzero(is2,is2) = 1
     END DO

     ode%nonzerocount=0; inz = 1; inz2 = 1
     
     DO is2=1,inp%Neq
        ode%IA(is2) = inz
        DO is=1,inp%Neq
           IF (ode%nonzero(is,is2) .gt. 1d-20) THEN
              ode%JA(inz) = is
              inz=inz+1
              ode%nonzerocount = ode%nonzerocount + 1
              countJ(is2) = countJ(is2)+1
           END IF
        END DO
     END DO

     ode%IA(inp%Neq+1) = inz
     inzNZ = inz-1
     
     !inzIA(1:inp%Neq+1) = ode%IA(1:inp%Neq+1) 
     !inzJA(1:inzNZ) = ode%JA(1:inzNZ)
     
   END SUBROUTINE COMPUTE_IAJA



!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------
!! JACOBIAN
!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------

  SUBROUTINE JACode(Neq,t,Y,JD,IA,JA,PDJ)

    !---------------------------------------------------------------
    ! compute the differential equations 
    !---------------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, ns, nr, il, ph, &
         F2tot, nstot, nmtot, alphacc, alphades
    USE SHARED_VARIABLES2

    INTEGER, INTENT(in) :: Neq, JD
    REAL(KIND=DP), INTENT(in) :: t
    INTEGER(KIND=LONG), DIMENSION(Neq+1), INTENT(inout) :: IA
    INTEGER(KIND=LONG), DIMENSION(inzNZ), INTENT(inout) :: JA
    REAL(KIND=DP), DIMENSION(Neq), INTENT(in) :: Y
    REAL(KIND=DP), DIMENSION(Neq), INTENT(out) :: PDJ
    REAL(KIND=DP), DIMENSION(Neq+1) :: PDJ2
    REAL(KIND=DP), DIMENSION(Neq) :: F
    REAL(KIND=DP) :: Ytotreact, COfact, Rreactot2
    !REAL(KIND=DP) :: nstot,nmtot,alphacc,alphades
    INTEGER(KIND=LONG) :: i, j, j2, k, ir1, ir2, ir3, ip1, ip2, ip3, ip4, ic, inz, choice

    PDJ = 0d0; PDJ2 = 0d0
    COfact = 1d0
    !itest=itest+1
    choice = 1

    ! 0 Gas-grain interaction, electron-grain recombination (Flower-PdF 03)
    DO ic=nr%istart(0),nr%iend(0)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       !WRITE(6,*) ic, JD, inp%Nreactions
       CALL JACtwobody1(Neq,ic,JD,Y,PDJ2)
    END DO

    ! 1 Direct cosmic-ray processes
    DO ic=nr%istart(1),nr%iend(1)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JAConebody(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 2 Photo-processes induced by cosmic-rays (secondary photons)
    DO ic=nr%istart(2),nr%iend(2)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JAConebody(Neq,ic,JD,Y,PDJ2)
    END DO

    ! 3 Photo-processes
    DO ic=nr%istart(3),nr%iend(3)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JAConebody(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 4 Bimolecular reactions
    DO ic=nr%istart(4),nr%iend(4)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JACtwobody1(Neq,ic,JD,Y,PDJ2)
    END DO

    ! 5 Charge exchange reactions
    DO ic=nr%istart(5),nr%iend(5)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JACtwobody1(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 6 Radiative associations
    DO ic=nr%istart(6),nr%iend(6)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JACtwobody1(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 7 Associative detachment
    DO ic=nr%istart(7),nr%iend(7)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JACtwobody1(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 8 Electronic recombination and attachment
    DO ic=nr%istart(8),nr%iend(8)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JACtwobody1(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 9 Third-body assisted association
    DO ic=nr%istart(9),nr%iend(9)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JACthreebody(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 10 Electron detachment from grains
    DO ic=nr%istart(10),nr%iend(10)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JAConebody(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 13 Gas-phase Photo-ionization and photo-dissociation.
    DO ic=nr%istart(13),nr%iend(13)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JAConebody(Neq,ic,JD,Y,PDJ2)
    END DO

    ! 14 Langmuir-Hinshelwood reactions for external surfaces
    DO ic=nr%istart(14),nr%iend(14)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JACsurf(Neq,ic,JD,Y,PDJ2)
    END DO

    ! 15 Langmuir-Hinshelwood reactions for internal surfaces
    DO ic=nr%istart(15),nr%iend(15)
    END DO
    
    ! 16 Langmuir-Hinshelwood reactions for external surfaces
    DO ic=nr%istart(16),nr%iend(16)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JACsurf(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 18 Acid-base reactions for external surfaces
    DO ic=nr%istart(18),nr%iend(18)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JACsurf(Neq,ic,JD,Y,PDJ2)
    END DO

    ! 19 Eley-Rideal reactions 
    ! JN + O -> NO
    DO ic=nr%istart(19),nr%iend(19)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JACtwobody1(Neq,ic,JD,Y,PDJ2)
    END DO

    ! 20 Gas-to-grain surface accretion
    DO ic=nr%istart(20),nr%iend(20)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JAConebody(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 21 Thermal evaporation. JX-->X
    DO ic=nr%istart(21),nr%iend(21)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JAConebody(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 22 Cosmic Ray induced general desorption following Hasegawa&Herbst 93
    DO ic=nr%istart(22),nr%iend(22)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JAConebody(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 23 Photo-desorption by induced-CR photons + background photons
    ! coming from wavelength-dependent exp studies (see Taquet et al. 2013)
    DO ic=nr%istart(23),nr%iend(23)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JAConebody(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 26 Photodissociation/photodesorption processes 
    ! coming from molecular dynamics simulations (see Taquet et al. 2013)
    DO ic=nr%istart(26),nr%iend(26)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JAConebody(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 28 Evaporation due to cloud-cloud collisions
    DO ic=nr%istart(28),nr%iend(28)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JAConebody(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 30 Non-porous surace->Pores rate exchange   JX -> QX
    DO ic=nr%istart(30),nr%iend(30)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JAConebody(Neq,ic,JD,Y,PDJ2)
    END DO
    
    ! 31 Pores->Non-porous rate    QX -> JX
    DO ic=nr%istart(31),nr%iend(31)
       IF (ph%Kreac(ic) .le. 1d-40) CYCLE
       CALL JAConebody(Neq,ic,JD,Y,PDJ2)
    END DO

    ! 32 Surface->Mantle exchange
!!$    IF (inp%chlayer .eq. 2 .and. nr%istart(32) .le. inp%Nreactions) THEN
!!$       CALL Freacs(Neq,Y,F)
!!$       CALL Fsurfmant(Neq,Y,F)
!!$       DO ic=nr%istart(32),nr%iend(32)
!!$          IF (abs(alphacc*F2tot/nstot) .le. 1d-40 .and. &
!!$               abs(alphades*F2tot/nmtot) .le. 1d-40) CYCLE
!!$          CALL JACsurfmant(Neq,ic,JD,Y,PDJ2)!,nstot,nmtot,alphacc,alphades)
!!$       END DO
!!$    END IF
    
    PDJ(1:Neq) = PDJ2(1:Neq)

  END SUBROUTINE JACode


  SUBROUTINE JAConebody(Neq,ic,JD,Y,PDJ2)

    USE SHARED_VARIABLES, only: ph, re, inp 
    INTEGER(KIND=LONG) :: i, j, j2, k, ir1, ir2, ir3, ip1, ip2, ip3, ip4
    !TYPE(physparam) :: ph
    !TYPE(reactions) :: re
    !TYPE(inputparam) :: inp
    INTEGER(KIND=LONG), INTENT(in) :: ic, JD, Neq
    REAL(KIND=DP), DIMENSION(Neq), INTENT(in) :: Y
    REAL(KIND=DP), DIMENSION(Neq+1), INTENT(inout) :: PDJ2
    REAL(KIND=DP) :: Rreactot

    Rreactot = 0.
    
    ! assign numbers to reactants and products
    ir1 = re%num_react1(ic)
    ir2 = re%num_react2(ic)
    ir3 = re%num_react3(ic)
    ip1 = re%num_prod1(ic)
    ip2 = re%num_prod2(ic)
    ip3 = re%num_prod3(ic)
    ip4 = re%num_prod4(ic)

    IF (ir1 .eq. JD) THEN
       Rreactot = ph%Kreac(ic)
       PDJ2(ir1) = PDJ2(ir1) - Rreactot
       PDJ2(ip1) = PDJ2(ip1) + Rreactot
       PDJ2(ip2) = PDJ2(ip2) + Rreactot
       PDJ2(ip3) = PDJ2(ip3) + Rreactot
       PDJ2(ip4) = PDJ2(ip4) + Rreactot
    END IF
    
  END SUBROUTINE JAConebody

  
  SUBROUTINE JACtwobody1(Neq,ic,JD,Y,PDJ2)

    USE SHARED_VARIABLES, only: ph, re, inp 
    INTEGER(KIND=LONG) :: i, j, j2, k, ir1, ir2, ir3, ip1, ip2, ip3, ip4
    INTEGER(KIND=LONG), INTENT(in) :: ic, JD, Neq
    REAL(KIND=DP), DIMENSION(Neq), INTENT(in) :: Y
    REAL(KIND=DP), DIMENSION(Neq+1), INTENT(inout) :: PDJ2
    REAL(KIND=DP) :: Rreactot, Ytotreact
    
    Rreactot = 0.
    
    ! assign numbers to reactants and products
    ir1 = re%num_react1(ic)
    ir2 = re%num_react2(ic)
    ir3 = re%num_react3(ic)
    ip1 = re%num_prod1(ic)
    ip2 = re%num_prod2(ic)
    ip3 = re%num_prod3(ic)
    ip4 = re%num_prod4(ic)
    
    IF (ir1 .eq. JD) THEN
       Ytotreact = Y(ir2)*ph%nH
       Rreactot = ph%Kreac(ic)*Ytotreact
       PDJ2(ir1) = PDJ2(ir1) - Rreactot
       PDJ2(ir2) = PDJ2(ir2) - Rreactot
       PDJ2(ip1) = PDJ2(ip1) + Rreactot
       PDJ2(ip2) = PDJ2(ip2) + Rreactot
       PDJ2(ip3) = PDJ2(ip3) + Rreactot
       PDJ2(ip4) = PDJ2(ip4) + Rreactot
    END IF
    IF (ir2 .eq. JD .and. ir1 .ne. JD) THEN
       Ytotreact = Y(ir1)*ph%nH
       Rreactot = ph%Kreac(ic)*Ytotreact
       PDJ2(ir1) = PDJ2(ir1) - Rreactot
       PDJ2(ir2) = PDJ2(ir2) - Rreactot
       PDJ2(ip1) = PDJ2(ip1) + Rreactot
       PDJ2(ip2) = PDJ2(ip2) + Rreactot
       PDJ2(ip3) = PDJ2(ip3) + Rreactot
       PDJ2(ip4) = PDJ2(ip4) + Rreactot
    END IF
    
  END SUBROUTINE JACtwobody1

  
  SUBROUTINE JACtwobody2(Neq,ic,JD,Y,PDJ2)

    USE SHARED_VARIABLES, only: ph, re, inp 
    INTEGER(KIND=LONG) :: i, j, j2, k, ir1, ir2, ir3, ip1, ip2, ip3, ip4
    INTEGER(KIND=LONG), INTENT(in) :: ic, JD, Neq
    REAL(KIND=DP), DIMENSION(Neq), INTENT(in) :: Y
    REAL(KIND=DP), DIMENSION(Neq+1), INTENT(inout) :: PDJ2
    REAL(KIND=DP) :: Rreactot, Ytotreact
    
    Rreactot = 0.
    
    ! assign numbers to reactants and products
    ir1 = re%num_react1(ic)
    ir2 = re%num_react2(ic)
    ir3 = re%num_react3(ic)
    ip1 = re%num_prod1(ic)
    ip2 = re%num_prod2(ic)
    ip3 = re%num_prod3(ic)
    ip4 = re%num_prod4(ic)

    IF (ir1 .eq. JD) THEN
       Ytotreact = Y(ir2)*ph%nH
       Rreactot = ph%Kreac(ic)*Ytotreact
       PDJ2(ir1) = PDJ2(ir1) - Rreactot
       PDJ2(ir2) = PDJ2(ir2) - Rreactot
       PDJ2(ip1) = PDJ2(ip1) + Rreactot
       PDJ2(ip2) = PDJ2(ip2) + Rreactot
       PDJ2(ip3) = PDJ2(ip3) + Rreactot
       PDJ2(ip4) = PDJ2(ip4) + Rreactot
    END IF
    IF (ir2 .eq. JD) THEN
       Ytotreact = Y(ir1)*ph%nH
       Rreactot = ph%Kreac(ic)*Ytotreact
       PDJ2(ir1) = PDJ2(ir1) - Rreactot
       PDJ2(ir2) = PDJ2(ir2) - Rreactot
       PDJ2(ip1) = PDJ2(ip1) + Rreactot
       PDJ2(ip2) = PDJ2(ip2) + Rreactot
       PDJ2(ip3) = PDJ2(ip3) + Rreactot
       PDJ2(ip4) = PDJ2(ip4) + Rreactot
    END IF
    
  END SUBROUTINE JACtwobody2

  
  SUBROUTINE JACsurf(Neq,ic,JD,Y,PDJ2)

    USE SHARED_VARIABLES, only: ph, re, inp, ab
    USE SHARED_VARIABLES, only: ns, nr
    INTEGER(KIND=LONG) :: i, j, j2, k, ir1, ir2, ir3, ip1, ip2, ip3, ip4, checkCO2reac
    INTEGER(KIND=LONG), INTENT(in) :: ic, JD, Neq
    REAL(KIND=DP), DIMENSION(Neq), INTENT(in) :: Y
    REAL(KIND=DP), DIMENSION(Neq+1), INTENT(inout) :: PDJ2
    REAL(KIND=DP) :: Rreactot, Ytotreact, COfact
    
    Rreactot = 0.; COfact = 1; checkCO2reac = 0
    
    ! assign numbers to reactants and products
    ir1 = re%num_react1(ic)
    ir2 = re%num_react2(ic)
    ir3 = re%num_react3(ic)
    ip1 = re%num_prod1(ic)
    ip2 = re%num_prod2(ic)
    ip3 = re%num_prod3(ic)
    ip4 = re%num_prod4(ic)

    ! compute the COfact factor to compute the CO2 formation from Garrod11
    IF ((re%num_react1(ic) .eq. ns%JO .and. &
         re%num_react2(ic) .eq. ns%JH) .or. &
         (re%num_react1(ic) .eq. ns%JO .and. &
         re%num_react2(ic) .eq. ns%JD) .and. &
         nr%reacCO2 .ne. 0) THEN
       COfact = 1-ab%abuCO
       ir3 = inp%Nspecies+1
    ELSE IF ((re%num_react1(ic) .eq. ns%JCO .and. &
         re%num_react2(ic) .eq. ns%JO .and. &
         re%num_react3(ic) .eq. ns%JH) .or. &
         (re%num_react1(ic) .eq. ns%JCO .and. &
         re%num_react2(ic) .eq. ns%JO .and. &
         re%num_react3(ic) .eq. ns%JD) .and. &
         nr%reacCO2 .ne. 0) THEN
       COfact = ab%abuCO
       checkCO2reac = 1
       !ir1 = re%num_react2(ic)
       !ir2 = re%num_react3(ic)
       !ir3 = inp%Nspecies+1
    END IF

    IF (checkCO2reac .eq. 0) THEN
       IF (ir1 .eq. JD) THEN
          Ytotreact = Y(ir2)*ph%nH
          Rreactot = COfact*ph%Kreac(ic)*Ytotreact
          PDJ2(ir1) = PDJ2(ir1) - Rreactot
          PDJ2(ir2) = PDJ2(ir2) - Rreactot
          PDJ2(ip1) = PDJ2(ip1) + Rreactot
          PDJ2(ip2) = PDJ2(ip2) + Rreactot
          PDJ2(ip3) = PDJ2(ip3) + Rreactot
          PDJ2(ip4) = PDJ2(ip4) + Rreactot
       END IF
       IF (ir2 .eq. JD) THEN
          Ytotreact = Y(ir1)*ph%nH
          Rreactot = COfact*ph%Kreac(ic)*Ytotreact
          PDJ2(ir1) = PDJ2(ir1) - Rreactot
          PDJ2(ir2) = PDJ2(ir2) - Rreactot
          PDJ2(ip1) = PDJ2(ip1) + Rreactot
          PDJ2(ip2) = PDJ2(ip2) + Rreactot
          PDJ2(ip3) = PDJ2(ip3) + Rreactot
          PDJ2(ip4) = PDJ2(ip4) + Rreactot
       END IF
    ELSE IF (checkCO2reac .eq. 1) THEN
       IF (ir2 .eq. JD) THEN
          Ytotreact = Y(ir3)*ph%nH
          Rreactot = COfact*ph%Kreac(ic)*Ytotreact
          PDJ2(ir1) = PDJ2(ir1) - Rreactot
          PDJ2(ir2) = PDJ2(ir2) - Rreactot
          PDJ2(ir3) = PDJ2(ir3) - Rreactot
          PDJ2(ip1) = PDJ2(ip1) + Rreactot
          PDJ2(ip2) = PDJ2(ip2) + Rreactot
          PDJ2(ip3) = PDJ2(ip3) + Rreactot
          PDJ2(ip4) = PDJ2(ip4) + Rreactot
       END IF
       IF (ir3 .eq. JD) THEN
          Ytotreact = Y(ir2)*ph%nH
          Rreactot = COfact*ph%Kreac(ic)*Ytotreact
          PDJ2(ir1) = PDJ2(ir1) - Rreactot
          PDJ2(ir2) = PDJ2(ir2) - Rreactot
          PDJ2(ir3) = PDJ2(ir3) - Rreactot
          PDJ2(ip1) = PDJ2(ip1) + Rreactot
          PDJ2(ip2) = PDJ2(ip2) + Rreactot
          PDJ2(ip3) = PDJ2(ip3) + Rreactot
          PDJ2(ip4) = PDJ2(ip4) + Rreactot
       END IF
    END IF
    
  END SUBROUTINE JACsurf

  
  SUBROUTINE JACsurfmant(Neq,ic,JD,Y,PDJ2)

    !---------------------------------------------------------------
    ! compute the ODEs
    !---------------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, re, nr, ab, gr, &
         F2tot, nstot, nmtot, alphacc, alphades
    USE SHARED_VARIABLES2

    INTEGER(KIND=LONG), INTENT(in) :: ic, JD, Neq
    REAL(KIND=DP), DIMENSION(Neq), INTENT(in) :: Y
    REAL(KIND=DP), DIMENSION(Neq+1), INTENT(inout) :: PDJ2
    REAL(KIND=DP), DIMENSION(inp%Nreactions) :: Rreactot
    REAL(KIND=DP), DIMENSION(Neq+1) :: Rform, Rdest, F2
    REAL(KIND=DP) :: Ytotreact, Rsurfmant
    INTEGER(KIND=LONG) :: j, ir1, ir2, ir3, ip1, ip2, ip3, ip4

    F2 = 0d0
    
    ! write(6,*) alphacc, alphades,nstot,nmtot,F2tot
    ! assign numbers to reactants and products
    ir1 = re%num_react1(ic)
    ip1 = re%num_prod1(ic)

    ! net accretion: surface->bulk exchange
    IF (ir1 .eq. JD .and. F2tot .ge. 0d0 .and. nstot .gt. 0d0) THEN
       !Rsurfmant = alphacc*F2tot*Y(ir1)/nstot
       PDJ2(ir1) = PDJ2(ir1) - alphacc*F2tot/nstot
       PDJ2(ip1) = PDJ2(ip1) + alphacc*F2tot/nstot
       ! net desorption: bulk->surface exchange
    ELSE IF (ip1 .eq. JD .and. F2tot .lt. 0d0 .and. nmtot .gt. 0d0) THEN
       !Rsurfmant = alphades*F2tot*Y(ip1)/nmtot
       PDJ2(ir1) = PDJ2(ir1) - alphades*F2tot/nmtot
       PDJ2(ip1) = PDJ2(ip1) + alphades*F2tot/nmtot
    END IF
       
    
  END SUBROUTINE JACsurfmant

  
  SUBROUTINE JACthreebody(Neq,ic,JD,Y,PDJ2)

    USE SHARED_VARIABLES, only: ph, re, inp
    
    INTEGER(KIND=LONG) :: i, j, j2, k, ir1, ir2, ir3, ip1, ip2, ip3, ip4
    INTEGER(KIND=LONG), INTENT(in) :: ic, JD, Neq
    REAL(KIND=DP), DIMENSION(Neq), INTENT(in) :: Y
    REAL(KIND=DP), DIMENSION(Neq+1), INTENT(inout) :: PDJ2
    REAL(KIND=DP) :: Rreactot, COfact, Ytotreact

    COfact = 1.
    Rreactot = 0.
        
    ! assign numbers to reactants and products
    ir1 = re%num_react1(ic)
    ir2 = re%num_react2(ic)
    ir3 = re%num_react3(ic)
    ip1 = re%num_prod1(ic)
    ip2 = re%num_prod2(ic)
    ip3 = re%num_prod3(ic)
    ip4 = re%num_prod4(ic)

    IF (ir1 .eq. JD) THEN
       Ytotreact = Y(ir2)*Y(ir3)*ph%nH*ph%nH
       Rreactot = COfact*ph%Kreac(ic)*Ytotreact
       PDJ2(ir1) = PDJ2(ir1) - Rreactot
       PDJ2(ir2) = PDJ2(ir2) - Rreactot
       PDJ2(ir3) = PDJ2(ir3) - Rreactot
       PDJ2(ip1) = PDJ2(ip1) + Rreactot
       PDJ2(ip2) = PDJ2(ip2) + Rreactot
       PDJ2(ip3) = PDJ2(ip3) + Rreactot
       PDJ2(ip4) = PDJ2(ip4) + Rreactot
    END IF
    IF (ir2 .eq. JD) THEN
       Ytotreact = Y(ir1)*Y(ir3)*ph%nH*ph%nH
       Rreactot = COfact*ph%Kreac(ic)*Ytotreact
       PDJ2(ir1) = PDJ2(ir1) - Rreactot
       PDJ2(ir2) = PDJ2(ir2) - Rreactot
       PDJ2(ir3) = PDJ2(ir3) - Rreactot
       PDJ2(ip1) = PDJ2(ip1) + Rreactot
       PDJ2(ip2) = PDJ2(ip2) + Rreactot
       PDJ2(ip3) = PDJ2(ip3) + Rreactot
       PDJ2(ip4) = PDJ2(ip4) + Rreactot
    END IF
    IF (ir3 .eq. JD) THEN 
       Ytotreact = Y(ir1)*Y(ir2)*ph%nH*ph%nH
       Rreactot = COfact*ph%Kreac(ic)*Ytotreact
       PDJ2(ir1) = PDJ2(ir1) - Rreactot
       PDJ2(ir2) = PDJ2(ir2) - Rreactot
       PDJ2(ir3) = PDJ2(ir3) - Rreactot
       PDJ2(ip1) = PDJ2(ip1) + Rreactot
       PDJ2(ip2) = PDJ2(ip2) + Rreactot
       PDJ2(ip3) = PDJ2(ip3) + Rreactot
       PDJ2(ip4) = PDJ2(ip4) + Rreactot
    END IF
    
  END SUBROUTINE JACthreebody

END MODULE ODES

