MODULE CHEMPROP

  !--------------------------------------------------------
  ! This module contains several subroutines that: 
  ! compute the rates of all chemical processes
  ! order the formation/destruction rates for each species
  !--------------------------------------------------------

  USE VARIABLES
!  USE INTLIB
  USE INPUT
  USE PHYSPROP

CONTAINS

  
!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------
!! REACTION RATES
!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------

  SUBROUTINE REAC_RATES(il,Y)

    !----------------------------------------------
    ! compute the rates of all chemical processes
    !----------------------------------------------

    USE SHARED_VARIABLES, only: inp, ph, sp, ns, re, nr, gr, ss, ab, ode

    INTEGER(KIND=LONG), INTENT(in) :: il
    REAL(KIND=DP), DIMENSION(inp%Neq), INTENT(in) :: Y
    REAL(KIND=DP) :: yieldUV, yieldCR, freac
    REAL(KIND=DP) :: abuH2, abuH, abuCO, abusp, abubare, COfact, abuH2O, test
    REAL(KIND=DP) :: nprot, p2g, Nlay, checknan, abuH2_prev, Pxy, Knondiff, maxcompet
    INTEGER(KIND=LONG) :: react0, in, ic, is

    ! compute the binding energy as function of the chemical abundance
    IF (inp%chevolen .eq. 1) THEN

       ! compute coverage of bare grains
       IF (il .eq. 1) THEN
          abubare = 1 - ab%Ptot/gr%Nstot(1)
          IF (abubare .lt. 0) abubare = 0
       ELSE IF (il .gt. 1) THEN
          abubare = 0
       END IF

       ! compute coverage of H2 on surface
       IF (ns%JH2 .ne. 0) THEN
          IF (ns%JoH2 .ne. 0) THEN 
             ab%abuH2 = (Y((ns%JH2)) + Y((ns%JoH2)))/ph%Nsndtot
          ELSE 
             ab%abuH2 = Y(ns%JH2)/ph%Nsndtot
          END IF
          abuH2_prev = ab%abuH2
          CALL SSabuH2(il,ab%abuH2,abubare,Y)
          IF (ab%abuH2 .lt. 0) ab%abuH2 = abuH2_prev
          IF (abubare .lt. 0) abubare = 0.
          IF (abubare .lt. 0) abubare = 0.
          !write(6,*) il, ab%abuH2, abubare, Y((ns%JH))
       END IF

       ! compute covereage of CO on surface
       ab%abuCO = Y(ns%JCO)/ph%Nsndtot

       ! derive Eb and Ed on surface
       DO is=inp%Nspgas+1,inp%Nspgas+inp%Nspgr
          abusp = Y(is)/ph%Nsndtot
          IF (abusp .ge. 1) abusp = 1.
          IF (ns%JH2 .ne. 0 .and. ab%abuH2 .ge. 0 .and. ab%abuH2 .lt. 0.5) THEN
            sp%Eb(is) = Eb_evol(il,ab%abuH2,ab%abuH,abusp,abubare,is,sp,ns)
          END IF
          IF (is .eq. ns%JH) THEN 
             sp%Ed(is) = inp%REdH*sp%Eb(is)
          ELSE
             sp%Ed(is) = inp%REdoth*sp%Eb(is)
          END IF
       END DO

!!$       IF (inp%chlayer .eq. 2) THEN
!!$
!!$          ! compute coverage of H2 in bulk
!!$          IF (ns%JH2 .ne. 0) THEN
!!$             IF (ns%JoH2 .ne. 0) THEN 
!!$                ab%abuH2bulk = (Y((ns%JH2+inp%Nspgr)) + Y((ns%JoH2+inp%Nspgr)))/ph%Nsndtot
!!$             ELSE 
!!$                ab%abuH2bulk = Y(ns%JH2+inp%Nspgr)/ph%Nsndtot
!!$             END IF
!!$             !CALL SSabuH2(il,ab%abuH2,abubare,Y)
!!$          END IF
!!$          
!!$          ! derive Eb and Ed in bulk
!!$          DO is=inp%Nspgas+inp%Nspgr+1,inp%Nspecies
!!$             abusp = Y(is)/ph%Nsndtot
!!$             IF (ns%JH2 .ne. 0 .and. ab%abuH2bulk .ge. 0) &
!!$                  sp%Eb(is) = Eb_evol(il,ab%abuH2bulk,ab%abuH,abusp,abubare,is,sp,ns)
!!$             !sp%Eb(is) = sp%Eb(is-inp%Nspgr)
!!$             sp%Ed(is) = inp%EdEbbulk*sp%Eb(is)
!!$          END DO
!!$       END IF
!!$       
    END IF

    ! check if binding energies are not null
    DO is=+inp%Nspgas+1,inp%Nspecies
    !write(6,*) sp%name(is), sp%Eb(is)
       IF (sp%Eb(is) .le. 0) THEN
          WRITE(6,*) is, sp%name(is), 'has binding energy = ', &
               sp%Eb(is), sp%Ed(is), ab%abuH2, abuH2_prev, ab%abuH2bulk,ab%abuH,abusp,abubare, Y((ns%JH2)) 
       END IF
    END DO
    !stop
    ! compute Nprotons/Ngrains
    DO is=1,inp%Nspgas
       IF (sp%charge(is) .gt. 0 .and. is .ne. ns%GP) nprot = nprot + Y(is)
    END DO
    !p2g = nprot/(Y(ns%GP)+Y(ns%GM)+Y(ns%G0))

    ph%Kreac = 0.

    ! 0 Gas-grain interaction, electron-grain recombination (Flower-PdF 03)
    DO ic=nr%istart(0),nr%iend(0)
       IF (trim(re%react1(ic)) .eq. 'G0' .or. &
            trim(re%react1(ic)) .eq. 'G+' .or. &
            trim(re%react1(ic)) .eq. 'G-') react0 = re%num_react2(ic)
       IF (trim(re%react1(ic)) .ne. 'G0' .and. &
            trim(re%react1(ic)) .ne. 'G+' .and. &
            trim(re%react1(ic)) .ne. 'G-') react0 = re%num_react1(ic)
       !write(6,*) ic, react0, sp%name(re%num_react1(ic)), re%type(ic)
       ph%Kreac(ic) = re%A(ic)*Pi*(gr%a(il)/2d0*1d-4)**2&
            *ph%Vth(react0)*0.5
       IF (sp%charge(re%num_react1(ic))*sp%charge(re%num_react2(ic)) .eq. -1d0) THEN
          ph%Kreac(ic) = ph%Kreac(ic) * (1+450d0/ph%Tg)
       END IF
       !write(6,*) ic, ph%Kreac(ic), re%A(ic),Pi,(gr%a(il)/2d0*1d-4)**2,ph%Vth(react0),ph%Tg
    END DO
        
    ! 10 electron detachment
    DO ic=nr%istart(10),nr%iend(10)
       ph%Kreac(ic) = Rcr(re%A(ic),ph%zeta)
    END DO
    
    ! 1 Direct cosmic-ray processes: Dissociation or ionization of species 
    ! due to direct collision with cosmic-ray particles.
    IF (nr%istart(1) .le. inp%Nreactions) THEN
       DO ic=nr%istart(1),nr%iend(1)
          ph%Kreac(ic) = Rcr(re%A(ic),ph%zeta)
       END DO
    END IF

    ! 2 Photo-processes induced by cosmic-rays (secondary photons): 
    ! Dissociation or ionization of species due to UV photons emitted 
    IF (nr%istart(2) .le. inp%Nreactions) THEN
       ! following H2 excitation.
       DO ic=nr%istart(2),nr%iend(2)
          ph%Kreac(ic) = Rcr(re%A(ic),ph%zeta)
       END DO
    END IF

    ! 3 Photo-processes: Dissociation or ionization of neutral species by UV 
    ! photons with a standard interstellar UV field.
    IF (nr%istart(3) .le. inp%Nreactions) THEN
       DO ic=nr%istart(3),nr%iend(3)
          ph%Kreac(ic) = Rphoto(re%A(ic),re%C(ic),ph%Av)
          ! compute the self shielding for H2, HD, D2, or CO
          CALL SHIELDING(ic)
       END DO
    END IF

    ! 4 Bimolecular reactions: Neutral-neutral (A + B -> C + D)
    !                          ion-neutral A+ + B -> C+ +D, A- + B -> C- + D
    !                          anion-cation (A+ + B- -> C + D) reactions
    !                          associative ionization (A + B -> AB+ + e-)
    IF (nr%istart(4) .le. inp%Nreactions) THEN
       DO ic=nr%istart(4),nr%iend(4)
          ! Standard Kooij
          IF (re%formula(ic) .eq. 3) THEN 
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN 
                ph%Kreac(ic) = Rkooij(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
             ! ionpol1
          ELSE IF (re%formula(ic) .eq. 4) THEN
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN 
                ph%Kreac(ic) = Rionpol1(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
             ! ionpol2
          ELSE IF (re%formula(ic) .eq. 5) THEN
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN 
                ph%Kreac(ic) = Rionpol2(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
          END IF
       END DO
    END IF

    ! 5 Charge exchange reactions: A+ + B -> A + B+ and A+ + B- -> A + B
    IF (nr%istart(5) .le. inp%Nreactions) THEN
       DO ic=nr%istart(5),nr%iend(5)
          ! Standard Kooij
          IF (re%formula(ic) .eq. 3) THEN 
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN 
                ph%Kreac(ic) = Rkooij(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
             ! ionpol1
          ELSE IF (re%formula(ic) .eq. 4) THEN
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN 
                ph%Kreac(ic) = Rionpol1(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
             ! ionpol2
          ELSE IF (re%formula(ic) .eq. 5) THEN
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN 
                ph%Kreac(ic) = Rionpol2(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
          END IF
       END DO
    END IF

    ! 6 Radiative associations: Association reactions between two species
    ! (neutral or ionized) stabilized by the emission of a photon 
    ! (A + B -> AB + photon or A+ + B -> AB+ + photon).
    IF (nr%istart(6) .le. inp%Nreactions) THEN
       DO ic=nr%istart(6),nr%iend(6)
          ! Standard Kooij
          IF (re%formula(ic) .eq. 3) THEN 
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN 
                ph%Kreac(ic) = Rkooij(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
             ! ionpol1
          ELSE IF (re%formula(ic) .eq. 4) THEN
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN 
                ph%Kreac(ic) = Rionpol1(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
             ! ionpol2
          ELSE IF (re%formula(ic) .eq. 5) THEN
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN 
                ph%Kreac(ic) = Rionpol2(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
          END IF
       END DO
    END IF

    ! 7 Associative detachment: Association of a neutral species and an anion
    ! resulting in the ejection of the extra electron (A- + B -> AB + e-).
    IF (nr%istart(7) .le. inp%Nreactions) THEN
       DO ic=nr%istart(7),nr%iend(7)
          ! Standard Kooij
          IF (re%formula(ic) .eq. 3) THEN 
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN 
                ph%Kreac(ic) = Rkooij(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
             ! ionpol1
          ELSE IF (re%formula(ic) .eq. 4) THEN
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN  
                ph%Kreac(ic) = Rionpol1(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
             ! ionpol2
          ELSE IF (re%formula(ic) .eq. 5) THEN
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN 
                ph%Kreac(ic) = Rionpol2(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
          END IF
       END DO
    END IF

    ! 8 Electronic recombination and attachment: Recombination of a positive
    ! ion with an electron resulting in the dissociation of the molecule 
    ! (AB+ + e- -> A + B) or the emission of a photon (AB+ + e- -> AB+photon)
    ! or the attachment of the electron (A + e- -> A- + photon)
    IF (nr%istart(8) .le. inp%Nreactions) THEN
       DO ic=nr%istart(8),nr%iend(8)
          ! Standard Kooij
          IF (re%formula(ic) .eq. 3) THEN 
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN  
                ph%Kreac(ic) = Rkooij(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
             ! ionpol1
          ELSE IF (re%formula(ic) .eq. 4) THEN
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN  
                ph%Kreac(ic) = Rionpol1(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
             ! ionpol2
          ELSE IF (re%formula(ic) .eq. 5) THEN
             IF (ph%Tg .le. re%Tmax(ic) .and. ph%Tg .gt. re%Tmin(ic)) THEN 
                ph%Kreac(ic) = Rionpol2(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
             END IF
          END IF
       END DO
    END IF

    ! 9 Third-body assisted association: Association reactions between two 
    ! species, stabilized by collision with a third body.
    !IF (nr%istart(9) .le. inp%Nreactions) THEN
       DO ic=nr%istart(9),nr%iend(9)
          ph%Kreac(ic) = Rkooij(re%A(ic),re%B(ic),re%C(ic),ph%Tg)
       END DO
    !END IF

    ! 13 Gas-phase Photo-ionization and photo-dissociation.
    IF (nr%istart(13) .le. inp%Nreactions) THEN
       DO ic=nr%istart(13),nr%iend(13)
          ph%Kreac(ic) = Rphoto(re%A(ic),re%C(ic),ph%Av)
       END DO
    END IF

    ! 14 Langmuir-Hinshelwood surface reactions for external surfaces
    IF (nr%istart(14) .le. inp%Nreactions) THEN
       DO ic=nr%istart(14),nr%iend(14)
          IF ((re%num_react1(ic) .le. inp%Nspgas+inp%Nspgr .and.  inp%chsurfchem .eq. 1) .or. &
               (re%num_react1(ic) .gt. inp%Nspgas+inp%Nspgr .and.  inp%chbulkchem .eq. 1)) THEN

             ! turn off chemical desorption
             !IF (re%num_prod1(ic) .le. inp%Nspgas) THEN
             !   re%B(ic) = 0.
             !   re%B(ic-1) = 1.
             !END IF
             
             ph%Kreac(ic) = Rsurf(ic,il)
             
             ! adjust branching ratio of reaction if product is a reactant of one non-diffusion reaction
             DO is=1,inp%Nspnd
                IF (re%num_prod1(ic).eq.re%num_ndreact1(is,1).or.re%num_prod2(ic).eq.re%num_ndreact1(is,1)) THEN
                   IF (re%num_ndreact2(is).le.inp%Nspgas+inp%Nspgr) &
                        ph%Kreac(ic) = ph%Kreac(ic)*ab%Xsurml(re%num_ndreact2(is)-inp%Nspgas)*inp%Fnd
                   IF (re%num_ndreact2(is).gt.inp%Nspgas+inp%Nspgr) &
                        ph%Kreac(ic) = ph%Kreac(ic)*ab%Xbulkml(re%num_ndreact2(is)-inp%Nspgas-inp%Nspgr)*inp%Fnd
                   !WRITE(6,*) ic, is, ph%Kreac(ic), sp%name(re%num_ndreact2(is))
                END IF
                IF (re%num_prod1(ic).eq.re%num_ndreact1(is,2).or.re%num_prod2(ic).eq.re%num_ndreact1(is,2)) THEN
                   IF (re%num_ndreact2(is).le.inp%Nspgas+inp%Nspgr) &
                        ph%Kreac(ic) = ph%Kreac(ic)*(1-ab%Xsurml(re%num_ndreact2(is)-inp%Nspgas)*inp%Fnd)
                   IF (re%num_ndreact2(is).gt.inp%Nspgas+inp%Nspgr) &
                        ph%Kreac(ic) = ph%Kreac(ic)*(1-ab%Xbulkml(re%num_ndreact2(is)-inp%Nspgas-inp%Nspgr)*inp%Fnd)
                   !WRITE(6,*) ic, is,ph%Kreac(ic), sp%name(re%num_ndreact2(is))
                END IF
             END DO
             
             !write(6,*) ic, ph%Kreac(ic)

          END IF
       END DO
    END IF

    ! 18 Acid-base-like reactions: reactions including one abundant molecule and
    ! not governed by diffusion but by their abundances
    maxcompet = 0.
    IF (nr%istart(18) .le. inp%Nreactions) THEN
       DO ic=nr%istart(18),nr%iend(18)
          ph%Kreac(ic) = re%A(ic)*re%B(ic)*max(nu0(re%num_react1(ic)),nu0(re%num_react2(ic)))*&
               max(re%proba(ic),exp(-re%C(ic)/ph%Td)) !/1d12 !/ph%Nsndtot
       END DO
    END IF
    
   
    ! 19 Eley-Rideal reactions 
    ! JN + O -> NO
    IF (nr%istart(19) .le. inp%Nreactions) THEN
       DO ic=nr%istart(19),nr%iend(19)
          ph%Kreac(ic) = re%A(ic)*Racc(re%num_react1(ic))/(ph%Nsndtot)
       END DO
    END IF

    ! 20 Gas-to-grain surface accretion
    IF (nr%istart(20) .le. inp%Nreactions) THEN
       DO ic=nr%istart(20),nr%iend(20)
          ph%Kreac(ic) = re%A(ic)*Racc(re%num_react1(ic))
          !IF (ph%time .ge. 1e6) ph%Kreac(ic) = re%A(ic)*Racc2(re%num_react1(ic))
          DO is=1,inp%Nspnd
             IF (re%num_prod1(ic).eq.re%num_ndreact1(is,1)) THEN
                IF (re%num_ndreact2(is).le.inp%Nspgas+inp%Nspgr) &
                     ph%Kreac(ic) = ph%Kreac(ic)*ab%Xsurml(re%num_ndreact2(is)-inp%Nspgas)*inp%Fnd
                IF (re%num_ndreact2(is).gt.inp%Nspgas+inp%Nspgr) &
                     ph%Kreac(ic) = ph%Kreac(ic)*ab%Xbulkml(re%num_ndreact2(is)-inp%Nspgas-inp%Nspgr)*inp%Fnd
                !WRITE(6,*) ic, is, ph%Kreac(ic), sp%name(re%num_ndreact2(is))
             END IF
             IF (re%num_prod1(ic).eq.re%num_ndreact1(is,2)) THEN
                IF (re%num_ndreact2(is).le.inp%Nspgas+inp%Nspgr) &
                     ph%Kreac(ic) = ph%Kreac(ic)*(1-ab%Xsurml(re%num_ndreact2(is)-inp%Nspgas)*inp%Fnd)
                IF (re%num_ndreact2(is).gt.inp%Nspgas+inp%Nspgr) &
                     ph%Kreac(ic) = ph%Kreac(ic)*(1-ab%Xbulkml(re%num_ndreact2(is)-inp%Nspgas-inp%Nspgr)*inp%Fnd)
                !WRITE(6,*) ic, is,ph%Kreac(ic), sp%name(re%num_ndreact2(is))
             END IF
          END DO
          !write(6,*) ic, sp%name(re%num_react1(ic)), ph%Kreac(ic)
       END DO
    END IF

    ! 21 Thermal evaporation. JX-->X
    IF (nr%istart(21) .le. inp%Nreactions) THEN
       DO ic=nr%istart(21),nr%iend(21)
          ph%Kreac(ic) = re%A(ic)*Revth(ph%Td,re%num_react1(ic))
       END DO
    END IF

    ! 22 Cosmic Ray induced general desorption following Hasegawa&Herbst 93
    IF (nr%istart(22) .le. inp%Nreactions) THEN
       DO ic=nr%istart(22),nr%iend(22)
          ph%Kreac(ic) = re%A(ic)*RevthCR(re%num_react1(ic),ph%zeta)
       END DO
    END IF

    
    ! 23 Photo-desorption by induced-CR photons + background photons
    ! coming from wavelength-dependent exp studies (see Taquet et al. 2013)
    IF (nr%istart(23) .le. inp%Nreactions) THEN
       DO ic=nr%istart(23),nr%iend(23)
          ph%Kreac(ic) = re%A(ic)*RevUV(il,re%num_react1(ic))
       END DO
    END IF

    ! 26 Photodissociation/photodesorption processes 
    ! coming from molecular dynamics simulations 
    IF (nr%istart(26) .le. inp%Nreactions) THEN
       DO ic=nr%istart(26),nr%iend(26)
          IF ((re%num_react1(ic) .le. inp%Nspgas+inp%Nspgr .and.  inp%chsurfchem .eq. 1) .or. &
               (re%num_react1(ic) .gt. inp%Nspgas+inp%Nspgr .and.  inp%chbulkchem .eq. 1)) THEN
             ph%Kreac(ic) = re%A(ic)*re%B(ic)*re%C(ic)*RabsUV(il,is)
          END IF
          !if (trim(re%react1(ic)) .eq. 'JH2O') write(6,*) ic, re%react1(ic), ph%Kreac(ic)
          !ph%Kreac(ic) = 0.
          DO is=1,inp%Nspnd
             IF (re%num_prod1(ic).eq.re%num_ndreact1(is,1)) THEN
                IF (re%num_ndreact2(is).le.inp%Nspgas+inp%Nspgr) &
                     ph%Kreac(ic) = ph%Kreac(ic)*ab%Xsurml(re%num_ndreact2(is)-inp%Nspgas)*inp%Fnd
                IF (re%num_ndreact2(is).gt.inp%Nspgas+inp%Nspgr) &
                     ph%Kreac(ic) = ph%Kreac(ic)*ab%Xbulkml(re%num_ndreact2(is)-inp%Nspgas-inp%Nspgr)*inp%Fnd
                !WRITE(6,*) ic, is, ph%Kreac(ic), sp%name(re%num_ndreact2(is))
             END IF
             IF (re%num_prod1(ic).eq.re%num_ndreact1(is,2)) THEN
                IF (re%num_ndreact2(is).le.inp%Nspgas+inp%Nspgr) &
                     ph%Kreac(ic) = ph%Kreac(ic)*(1-ab%Xsurml(re%num_ndreact2(is)-inp%Nspgas)*inp%Fnd)
                IF (re%num_ndreact2(is).gt.inp%Nspgas+inp%Nspgr) &
                     ph%Kreac(ic) = ph%Kreac(ic)*(1-ab%Xbulkml(re%num_ndreact2(is)-inp%Nspgas-inp%Nspgr)*inp%Fnd)
                !WRITE(6,*) ic, is,ph%Kreac(ic), sp%name(re%num_ndreact2(is))
             END IF
          END DO
          ph%Kreac(ic) = 0.
       END DO
    END IF
    
    ! 27 Ice sputtering scaled on H2O experiments by Dartois et al. (2015)
    IF (nr%istart(27) .le. inp%Nreactions) THEN
       DO ic=nr%istart(27),nr%iend(27)
          ph%Kreac(ic) = re%A(ic)*RsputCR(il,re%num_react1(ic))
          !if (trim(re%react1(ic)) .eq. 'JH2O') write(6,*) ic, re%react1(ic), ph%Kreac(ic)
          ph%Kreac(ic) = 0.
       END DO
    END IF
    
    ! 28 Evaporation induced by grain-grain collision in the case of cloud-cloud interaction
!!$    IF (nr%istart(28) .le. inp%Nreactions) THEN
!!$       DO ic=nr%istart(28),nr%iend(28)
!!$          ph%Kreac(ic) = re%A(ic)*Revcol(re%num_react1(ic),il)
!!$          IF (ph%time.le.1e6) ph%Kreac(ic) = 0.
!!$          IF (ph%time.gt.1e6) write(6,*)  ic, ph%Kreac(ic)
!!$       END DO
!!$    END IF

    
    ! 30 Non-porous surface->Pores rate exchange   JX -> QX
    IF (nr%istart(30) .le. inp%Nreactions) THEN
       DO ic=nr%istart(30),nr%iend(30)
          IF (inp%Fporsur .eq. 0) THEN
             ph%Kreac(ic) = 0d0
          ELSE IF (inp%Fporsur .gt. 0)  THEN
             ph%Kreac(ic) = (gr%Fedgetot(il)/gr%Fnptot(il))*&
                  Rdiff(ph%Td,re%num_react1(ic))
          END IF
       END DO
    END IF

    ! 31 Pores->Non-porous rate    QX -> JX
    IF (nr%istart(31) .le. inp%Nreactions) THEN
       DO ic=nr%istart(31),nr%iend(31)
          IF (inp%Fporsur .eq. 0) THEN
             ph%Kreac(ic) = 0d0
          ELSE IF (inp%Fporsur .gt. 0)  THEN
             ph%Kreac(ic) = (gr%Fedgetot(il)/gr%Fportot(il))*&
                  Rdiff(ph%Td,re%num_react1(ic))
          END IF
       END DO
    END IF
    

    ! 80 Isomerisation
    IF (nr%istart(80) .le. inp%Nreactions) THEN
       DO ic=nr%istart(80),nr%iend(80)
          ph%Kreac(ic) = re%A(ic)*nu0(re%num_react1(ic))*max(re%proba(ic),exp(-re%C(ic)/ph%Td))
          !write(6,*) ic,ph%Kreac(ic),re%C(ic)
       END DO
    END IF
    
    ph%Kreac(inp%Nreactions+1) = 0.

    checknan = 0.
    
    DO ic=1,inp%Nreactions
       ! No desorption for H and D when grains are bare
       !IF (inp%chlayer .eq. 2 .and. il .eq. 1 .and. &
       !     re%type(ic) .ge. 21 .and. re%type(ic) .le. 23) THEN
       !   IF (re%num_react1(ic) .eq. ns%JH .or. re%num_react1(ic) .eq. ns%JD) THEN
       !      ph%Kreac(ic) = 0d0
       !   END IF
       !END IF
       !IF ((sp%num(re%num_react1(ic)).eq.ns%JCO.and.sp%num(re%num_react2(ic)).eq.ns%JO).or.&
       !     (sp%num(re%num_react1(ic)).eq.ns%JO.and.sp%num(re%num_react2(ic)).eq.ns%JCO)) THEN
       !   ph%Kreac(ic) = 0.
       !END IF
       !IF ((ic) .eq. 3374) write(6,*) ph%Kreac(ic)
       ! check if no NaN in rates
       IF (ISNAN(ph%Kreac(ic))) THEN
          checknan = 1.
          WRITE(6,*) 'NaN rate found for reaction', ic, sp%name(re%num_react1(ic))
          stop
       END IF
       !ph%Kreac(ic) = 0.
       !write(6,*) ic, re%type(ic), ph%Kreac(ic)
    END DO
!stop
    
  END SUBROUTINE REAC_RATES

  
  SUBROUTINE SSabuH2(il,PH2,Pbare,Y)
    
    USE SHARED_VARIABLES2

    ! --------------------------------------------------------------------
    ! Compute the H2 abundance and binding energies in ices
    ! --------------------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, sp, ns, ph
    
    REAL(KIND=DP) :: nu, Kacc
    REAL(KIND=DP) :: Eb,EbH2O,EbH2
    INTEGER(KIND=LONG) :: k, Niter
    REAL(KIND=DP) :: PHatom, PH2_prev, PH2_prev2, nH2gas, KaccH2, KevH2, nH2ice
    REAL(KIND=DP), INTENT(inout) :: PH2
    INTEGER(KIND=LONG), INTENT(in) :: il
    REAL(KIND=DP), INTENT(inout) :: Pbare
    REAL(KIND=DP), DIMENSION(inp%Neq), INTENT(in) :: Y
    REAL(KIND=DP) :: KevH2a, KevH2b, KevH2c, EQa, EQb, EQc, PH2a, PH2b, PH2c
    REAL(KIND=DP) :: PH2a_prev, PH2b_prev, PH2c_prev, Pbarea, Pbareb, Pbarec
    
    k = 0
    nH2gas = Y(ns%H2)
    IF (ns%oH2 .ne. 0) nH2gas = Y((ns%H2)) + Y((ns%oH2))
    PH2_prev = PH2
    PH2 = 0.5
    PH2_prev2 = 0
    PHatom = Y((ns%JH))/ph%Nsndtot

    PH2a = 0
    PH2b = 1
    PH2c = 0.5
    PH2cprev = 0
    Pbarea = Pbare
    Pbareb = 0d0
    Pbarec = Pbare-PH2c
    IF (Pbarec .lt. 0d0) Pbarec = 0d0

    KaccH2 = Racc(ns%H2)

    sp%Eb(ns%JH2) = Eb_evol(il,PH2a,PHatom,PH2a,Pbarea,ns%JH2,sp,ns)
    KevH2a = Revth(ph%Td,ns%JH2) + RevthCR(ns%JH2,ph%zeta)
    sp%Eb(ns%JH2) = Eb_evol(il,PH2b,PHatom,PH2b,Pbareb,ns%JH2,sp,ns)
    KevH2b = Revth(ph%Td,ns%JH2) + RevthCR(ns%JH2,ph%zeta)
    sp%Eb(ns%JH2) = Eb_evol(il,PH2c,PHatom,PH2c,Pbarec,ns%JH2,sp,ns)
    KevH2c = Revth(ph%Td,ns%JH2) + RevthCR(ns%JH2,ph%zeta)

    !Ebb = Eb_evol(il,PH2b,PHatom,ns%JH2,sp)
    !Ebc = Eb_evol(il,PH2c,PHatom,ns%JH2,sp)

    EQa = KaccH2*nH2gas - KevH2a*PH2a/(ph%nH/ph%Nsndtot)
    EQb = KaccH2*nH2gas - KevH2b*PH2b/(ph%nH/ph%Nsndtot)
    EQc = KaccH2*nH2gas - KevH2c*PH2c/(ph%nH/ph%Nsndtot)

    IF (EQa*EQb .lt. 0d0) THEN

       DO WHILE (abs(EQc) .gt. 1d-15 .and. PH2c .ne. PH2cprev)
          
          PH2a_prev = PH2a
          PH2b_prev = PH2b
          PH2c_prev = PH2c

          IF (EQa*EQc .lt. 0) THEN
             PH2b = PH2c
          ELSE IF (EQa*EQc .gt. 0) THEN
             PH2a = PH2c
          END IF
          PH2c = (PH2a+PH2b)/2d0

          Pbarea = Pbarea - (PH2a-PH2a_prev)
          Pbareb = Pbareb - (PH2b-PH2b_prev)
          Pbarec = Pbarec - (PH2c-PH2c_prev)

          sp%Eb(ns%JH2) = Eb_evol(il,PH2a,PHatom,PH2a,Pbarea,ns%JH2,sp,ns)
          KevH2a = Revth(ph%Td,ns%JH2) + RevthCR(ns%JH2,ph%zeta)
          sp%Eb(ns%JH2) = Eb_evol(il,PH2b,PHatom,PH2b,Pbareb,ns%JH2,sp,ns)
          KevH2b = Revth(ph%Td,ns%JH2) + RevthCR(ns%JH2,ph%zeta)
          sp%Eb(ns%JH2) = Eb_evol(il,PH2c,PHatom,PH2c,Pbarec,ns%JH2,sp,ns)
          KevH2c = Revth(ph%Td,ns%JH2) + RevthCR(ns%JH2,ph%zeta)

          EQa = KaccH2*nH2gas - KevH2a*PH2a/(ph%nH/ph%Nsndtot)
          EQb = KaccH2*nH2gas - KevH2b*PH2b/(ph%nH/ph%Nsndtot)
          EQc = KaccH2*nH2gas - KevH2c*PH2c/(ph%nH/ph%Nsndtot)

          Niter = Niter + 1

       END DO

    END IF

    PH2 = PH2c 
    Pbare = Pbare - (PH2-PH2_prev)

  END SUBROUTINE SSabuH2

  
  ! self shielding
  SUBROUTINE SHIELDING(ic)

    USE SHARED_VARIABLES, only: inp, ph, sp, ns, re, nr, ss, ab
    
    INTEGER(KIND=LONG), INTENT(in) :: ic
    REAL(KIND=DP) :: XH2, XHD, XD2, XCO, NCOLLH2, NCOLLCO
    REAL(KIND=DP) :: teta, teta1, teta2, teta3
    
    ! get X(H2) for the self-shielding calculation
    IF (ns%oH2 .lt. inp%Nspecies+1 .and. ns%oH2 .gt. 0) THEN
       XH2 = ab%Xgas(ns%oH2) + ab%Xgas(ns%pH2)
    ELSE IF (ns%oH2 .eq. inp%Nspecies+1 .or. ns%oH2 .eq. 0) THEN
       XH2 = ab%Xgas(ns%H2)
    END IF
    
    ! computation of the H2 self-shielding by linear extrapolation
    ! from Lee et al. (1996) or from the Meudon PDR code
    IF (re%num_react1(ic) .eq. ns%oH2 .or. &
         re%num_react1(ic) .eq. ns%pH2 .or. &
         re%num_react1(ic) .eq. ns%H2) THEN   
       
       NCOLLH2=(ph%Av/5.34d-22*2)*XH2 
       teta=1
       DO in=1,202                                        
          IF ((ss%NH2_3(in).LE.NCOLLH2) .and. &
               (ss%NH2_3(in+1).GE.NCOLLH2)) THEN
             teta=ss%T_H2_3(in)+(NCOLLH2-ss%NH2_3(in))* &
                  (ss%T_H2_3(in+1)-ss%T_H2_3(in))/&
                  (ss%NH2_3(in+1)-ss%NH2_3(in))
          END IF
       END DO
       IF (NCOLLH2.GT.ss%NH2_3(202)) teta=ss%T_H2_2(105)
       
       ph%Kreac(ic) = re%A(ic)*teta
       
    ENDIF
    
    ! computation of the HD and D2 self-shielding by linear extrapolation
    ! from the Meudon PDR code
    IF(re%num_react1(ic) .eq. ns%HD .or. &
         re%num_react1(ic) .eq. ns%oD2 .or. &
         re%num_react1(ic) .eq. ns%pD2 .or. &
         re%num_react1(ic) .eq. ns%D2) THEN
       
       IF (re%num_react1(ic) .eq. ns%HD) THEN
          XHD = ab%Xgas(ns%HD)
          NCOLLH2=(ph%Av/5.34d-22*2)*XH2 
       ELSE IF (re%num_react1(ic) .eq. ns%oD2 .or. &
            re%num_react1(ic) .eq. ns%pD2 .or. &
            re%num_react1(ic) .eq. ns%D2) THEN
          XD2 = ab%Xgas(ns%D2)  + ab%Xgas(ns%oD2)
          NCOLLH2=(ph%Av/5.34d-22*2)*XH2 * XD2/ab%Xgas(ns%HD)
       END IF
       teta=1                                            
       DO in=1,202                                        
          IF ((ss%NH2_3(in).LE.NCOLLH2) .and. &
               (ss%NH2_3(in+1).GE.NCOLLH2)) THEN
             teta=ss%T_HD(in)+(NCOLLH2-ss%NH2_3(in))* &
                  (ss%T_HD(in+1)-ss%T_HD(in))/&
                  (ss%NH2_3(in+1)-ss%NH2_3(in))
          END IF
       END DO
       IF (NCOLLH2.GT.ss%NH2_3(202)) TETA=ss%T_H2_2(105)
       
       ph%Kreac(ic) = re%A(ic)*teta
       
    END IF
    
    ! computation of the CO self-shielding by linear extrapolation
    ! from Lee et al. (1996)
    IF (re%num_react1(ic) .eq. ns%CO) THEN   
       
       NCOLLH2=(ph%Av/5.34d-22)*XH2                    
       NCOLLCO=(ph%Av/5.34D-22)*ab%Xgas(ns%CO)
       
       teta1=1                                           
       teta2=1                                           
       teta3=1                                           
       DO in=1,51                                         
          IF ((ss%NCO(in).LE.NCOLLCO) .and. &
               (ss%NCO(in+1).GE.NCOLLH2)) THEN
             teta2=ss%T_CO(in)+(NCOLLCO-ss%NCO(in))*&
                  (ss%T_CO(in+1)-ss%T_CO(in))/(ss%NCO(in+1)-ss%NCO(in))    
          END IF
       END DO
       DO in=1,42                                         
          IF ((ss%NH2_1(in).LE.NCOLLH2) .and. &
               (ss%NH2_1(in+1).GE.NCOLLH2))  THEN
             teta1=ss%T_H2_1(in)+(NCOLLH2-ss%NH2_1(in))*&
                  (ss%T_H2_1(in+1)-ss%T_H2_1(in))&
                  /(ss%NH2_1(in+1)-ss%NH2_1(in)) 
          END IF
          IF ((ss%AV2(in).LE.ph%Av) .and. &
               (ss%AV2(in+1).GE.ph%Av)) THEN
             teta3 = ss%T_AV(in) + &
                  (ph%Av-ss%AV2(in))*(ss%T_AV(in+1)-ss%T_AV(in))&
                  /(ss%AV2(in+1)-ss%AV2(in))                     
          END IF
       END DO
       IF (NCOLLH2 .gt. ss%NH2_1(43)) teta1=ss%T_H2_1(43)
       IF (NCOLLCO .gt. ss%NCO(52)) teta2=ss%T_CO(52)
       IF (ph%Av .gt. ss%AV2(43)) TETA3=ss%T_AV(43)
       
       ph%Kreac(ic) = re%A(ic)*teta1*teta2*teta3    		      
       
    END IF
    
    
  END SUBROUTINE SHIELDING

  
  ! reaction probabilities
  SUBROUTINE REAC_PROBAS

    USE SHARED_VARIABLES, only: inp, re, sp
    
    REAL(KIND=DP) :: XH2, XHD, XD2, XCO, NCOLLH2, NCOLLCO
    REAL(KIND=DP) :: teta, teta1, teta2, teta3, red_mass
    REAL :: proba2

     re%proba = 0.
     DO ic=1,inp%Nreactions
        IF (re%type(ic) .eq. 14 .or. re%type(ic) .eq. 15 .or. &
             re%type(ic) .eq. 16 .or. re%type(ic) .eq. 18 .or. re%type(ic) .eq. 80) THEN
           ! 1- compute the reac proba using reacs.in and a square barrier
           re%proba(ic) = proba_square(ic)!/1d10
           ! write(6,*) ic,re%proba(ic)
           ! 2- compute the reac proba using eckart model and eckart_data.in
           IF (inp%checkart .eq. 2 .and. re%C(ic) .gt. 0d0) THEN
              re%proba(ic) = proba_eckart(ic)
              red_mass = sp%mass(re%num_react1(ic))*&
                   sp%mass(re%num_react2(ic))/&
                   (sp%mass(re%num_react1(ic)) + &
                   sp%mass(re%num_react2(ic)))*Mpr
              proba2 = re%proba(ic) !proba_eckart(ic)
              !write(6,*) ic, re%react1(ic), re%react2(ic), re%prod1(ic), re%prod2(ic), re%proba(ic)
              ! proba_square(ic),proba_eckart(ic),re%proba(ic),red_mass/Mpr
           END IF
           ! 3- compute relative rates and energies from Watanabe/Caselli
           IF (inp%chDnetwork .ne. 1) THEN
              re%proba(ic) = proba_COnetwork(ic)
           END IF
        END IF
     END DO

  END SUBROUTINE REAC_PROBAS

  ! chemical desorption
  SUBROUTINE CHEM_DESORPTION

    USE SHARED_VARIABLES, only: inp, re, sp, ent
    
    REAL(KIND=DP) :: XH2, XHD, XD2, XCO, NCOLLH2, NCOLLCO
    REAL(KIND=DP) :: mm, m, cd, exoreac, eps, Ndeg
    REAL :: proba2

    mm = 120.
    
    DO ic=1,inp%Nreactions

       cd = 0.
       
       IF (re%type(ic) .eq. 14 .or. re%type(ic) .eq. 15 .or. &
            re%type(ic) .eq. 16 .or. re%type(ic) .eq. 80) THEN
          
          IF (sp%dHf(re%num_react1(ic)).ge.-9999..and.sp%dHf(re%num_react2(ic)).ge.-9999..and.&
               sp%dHf(re%num_prod1(ic)).ge.-9999.) THEN
             
             exoreac = sp%dHf(re%num_react1(ic))+sp%dHf(re%num_react2(ic))-&
                  sp%dHf(re%num_prod1(ic))+sp%dHf(re%num_prod2(ic))
             exoreac = exoreac*120.
             
             eps = (mm-sp%mass(re%num_prod1(ic)))**2/(mm+sp%mass(re%num_prod1(ic)))**2
             !eps = 0.4
             
             Ndeg = 3*(sp%element(1,re%num_prod1(ic)) + &
                  sp%element(2,re%num_prod1(ic)) + &
                  sp%element(3,re%num_prod1(ic)) + &
                  sp%element(4,re%num_prod1(ic)) + &
                  sp%element(5,re%num_prod1(ic)) + &
                  sp%element(6,re%num_prod1(ic)) + &
                  sp%element(7,re%num_prod1(ic)) + &
                  sp%element(8,re%num_prod1(ic)) + &
                  sp%element(9,re%num_prod1(ic)) + &
                  sp%element(10,re%num_prod1(ic)) + &
                  sp%element(11,re%num_prod1(ic)) + &
                  sp%element(12,re%num_prod1(ic)) + &
                  sp%element(13,re%num_prod1(ic)) + &
                  sp%element(14,re%num_prod1(ic)))

             IF (exoreac .ge. 0) THEN
                cd = exp(-sp%Eb(re%num_prod1(ic))/(eps*exoreac/Ndeg))
                re%B(ic) = 1-cd
                re%B(ic+1) = cd
             END IF
             
             !write(6,*) ic, re%react1(ic),re%react2(ic),re%prod1(ic), cd, exoreac/12000., eps, Ndeg, sp%Eb(re%num_prod1(ic))
             
          END IF
       END IF
    END DO
    
   END SUBROUTINE CHEM_DESORPTION


!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------
!! FORMATION AND DESTRUCTION RATES
!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------


  SUBROUTINE COMPUTE_REACRATES(inp,ph,sp,ns,re,nr,Y)
    
    ! -----------------------------------------------------------------------
    ! order the importance of all reactions in species formation/destruction
    ! -----------------------------------------------------------------------

    USE SHARED_VARIABLES2

    IMPLICIT  NONE
    
    TYPE(inputparam), INTENT(in) :: inp
    TYPE(physparam), INTENT(inout) :: ph
    TYPE(species), INTENT(inout) :: sp
    TYPE(specnum), INTENT(in) :: ns
    TYPE(reactions), INTENT(in) :: re
    TYPE(reactionnumber), INTENT(in) :: nr
    REAL(KIND=DP), DIMENSION(inp%Neq), INTENT(in) :: Y

    REAL(KIND=DP) :: abuCO, COfact, Ytotreact
    INTEGER(KIND=LONG) :: ic, is

    !CALL REAC_RATES(chevolen,Tgas,Tdust,nH,Vth,stick,ph%Kreac)
    abuCO = (Y(ns%JCO))/ph%Nsndtot

    ! formation
    DO is=1,inp%Nsprate
       
       DO ic=1,sp%Nreacform(is)

          COfact = 1
          IF (re%num_react3(sp%numform(is,ic)) .ne. inp%Nspecies+1) THEN
             Ytotreact = Y(re%num_react1(sp%numform(is,ic)))*&
                  Y(re%num_react2(sp%numform(is,ic)))*&
                  Y(re%num_react3(sp%numform(is,ic)))*ph%nH*ph%nH
          ELSE IF (re%num_react3(sp%numform(is,ic)) .eq. inp%Nspecies+1) THEN
             IF (re%num_react2(sp%numform(is,ic)) .eq. inp%Nspecies+1) &
                  Ytotreact = Y(re%num_react1(sp%numform(is,ic)))
             IF (re%num_react2(sp%numform(is,ic)) .ne. inp%Nspecies+1) &
                  Ytotreact = Y(re%num_react1(sp%numform(is,ic)))*&
                  Y(re%num_react2(sp%numform(is,ic)))*ph%nH
          END IF
          sp%rateform(is,ic) = sp%degform(is,ic)*COfact*&
               ph%Kreac(sp%numform(is,ic))*Ytotreact


          IF ((re%num_react1(sp%numform(is,ic)) .eq. ns%JO .and. &
               re%num_react2(sp%numform(is,ic)) .eq. ns%JH) .or. &
               (re%num_react1(sp%numform(is,ic)) .eq. ns%JO .and. &
               re%num_react2(sp%numform(is,ic)) .eq. ns%JD) .or. &
               (re%num_react1(sp%numform(is,ic)) .eq. ns%JCO .and. &
               re%num_react2(sp%numform(is,ic)) .eq. ns%JO .and. &
               re%num_react3(sp%numform(is,ic)) .eq. ns%JH) .or. &
               (re%num_react1(sp%numform(is,ic)) .eq. ns%JCO .and. &
               re%num_react2(sp%numform(is,ic)) .eq. ns%JO .and. &
               re%num_react3(sp%numform(is,ic)) .eq. ns%JD) .and. &
               nr%reacCO2 .ne. 0) THEN
             
             IF (re%num_react1(sp%numform(is,ic)) .eq. ns%JO) THEN
                COfact = 1-abuCO
                sp%rateform(is,ic) = &
                     sp%degform(is,ic)*COfact*ph%Kreac(sp%numform(is,ic))*&
                     Y((re%num_react1(sp%numform(is,ic))))*&
                     Y((re%num_react2(sp%numform(is,ic))))*ph%nH
             ELSE IF (re%num_react1(sp%numform(is,ic)) .eq. ns%JCO) THEN
                COfact = abuCO
                sp%rateform(is,ic) = &
                     sp%degform(is,ic)*COfact* ph%Kreac(sp%numform(is,ic))*&
                     Y((re%num_react2(sp%numform(is,ic))))*&
                     Y((re%num_react3(sp%numform(is,ic))))*ph%nH
             END IF
          END IF

       END DO

       ! destruction
       DO ic=1,sp%Nreacdest(is)

          COfact = 1
          IF (re%num_react3(sp%numdest(is,ic)) .ne. inp%Nspecies+1) THEN
             Ytotreact = Y(re%num_react1(sp%numdest(is,ic)))*&
                  Y(re%num_react2(sp%numdest(is,ic)))*&
                  Y(re%num_react3(sp%numdest(is,ic)))*ph%nH*ph%nH
          ELSE IF (re%num_react3(sp%numdest(is,ic)) .eq. inp%Nspecies+1) THEN
             IF (re%num_react2(sp%numdest(is,ic)) .eq. inp%Nspecies+1) &
                  Ytotreact = Y(re%num_react1(sp%numdest(is,ic)))
             IF (re%num_react2(sp%numdest(is,ic)) .ne. inp%Nspecies+1) &
                  Ytotreact = Y(re%num_react1(sp%numdest(is,ic)))*&
                  Y(re%num_react2(sp%numdest(is,ic)))*ph%nH
          END IF
          sp%ratedest(is,ic) = sp%degdest(is,ic)*COfact*&
               ph%Kreac(sp%numdest(is,ic))*Ytotreact
               
          IF ((re%num_react1(sp%numdest(is,ic)) .eq. ns%JO .and. &
               re%num_react2(sp%numdest(is,ic)) .eq. ns%JH) .or. &
               (re%num_react1(sp%numdest(is,ic)) .eq. ns%JO .and. &
               re%num_react2(sp%numdest(is,ic)) .eq. ns%JD) .or. &
               (re%num_react1(sp%numdest(is,ic)) .eq. ns%JCO .and. &
               re%num_react2(sp%numdest(is,ic)) .eq. ns%JO .and. &
               re%num_react3(sp%numdest(is,ic)) .eq. ns%JH) .or. &
               (re%num_react1(sp%numdest(is,ic)) .eq. ns%JCO .and. &
               re%num_react2(sp%numdest(is,ic)) .eq. ns%JO .and. &
               re%num_react3(sp%numdest(is,ic)) .eq. ns%JD) .and. &
               nr%reacCO2 .ne. 0) THEN
             
             IF (re%num_react1(sp%numdest(is,ic)) .eq. &
                  ns%JO) THEN
                COfact = 1-abuCO
                sp%ratedest(is,ic) = &
                     sp%degdest(is,ic)*COfact*ph%Kreac(sp%numdest(is,ic))*&
                     Y((re%num_react1(sp%numdest(is,ic))))*&
                     Y((re%num_react2(sp%numdest(is,ic))))*ph%nH
             ELSE IF (re%num_react1(sp%numdest(is,ic)) .eq. &
                  ns%JCO) THEN
                COfact = abuCO
                sp%ratedest(is,ic) = &
                     sp%degdest(is,ic)*COfact*ph%Kreac(sp%numdest(is,ic))*&
                     Y((re%num_react2(sp%numdest(is,ic))))*&
                     Y((re%num_react3(sp%numdest(is,ic))))*ph%nH
             END IF
          END IF

       END DO
       
!!$       CALL SORT(sp%rateform(is,1:Nreacform(is)),&
!!$            sp%numform(is,1:Nreacform(is)),Nreacform(is))
!!$       CALL SORT(sp%ratedest(is,1:Nreacdest(is)),&
!!$            sp%numdest(is,1:Nreacdest(is)),Nreacdest(is))

    END DO

  END SUBROUTINE COMPUTE_REACRATES


  SUBROUTINE ORDER_RATES(Y)

    ! -----------------------------------------------------------------------
    ! compute the integrated contribution of every reaction for each species
    ! -----------------------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, ph, sp, ns, re, nr, out
    USE SHARED_VARIABLES2

    IMPLICIT NONE

    REAL(KIND=DP), DIMENSION(inp%Neq), INTENT(in) :: Y
    
    INTEGER(KIND=LONG), DIMENSION(out%Ninisteps) :: timecontrib
    REAL(KIND=DP) :: contribformtot, contribdesttot
    REAL(KIND=DP), DIMENSION(inp%Nsprate,inp%Nreactions) :: contribform0
    REAL(KIND=DP), DIMENSION(inp%Nsprate,inp%Nreactions) :: contribdest0
    REAL(KIND=DP), DIMENSION(inp%Nsprate,inp%Nreactions) :: contribform02
    REAL(KIND=DP), DIMENSION(inp%Nsprate,inp%Nreactions) :: contribdest02
    REAL(KIND=DP), DIMENSION(inp%Nsprate,inp%Nreactions) :: contribform2
    REAL(KIND=DP), DIMENSION(inp%Nsprate,inp%Nreactions) :: contribdest2
    INTEGER(KIND=LONG), DIMENSION(inp%Nsprate+1,inp%Nreactions) :: specnumform0
    INTEGER(KIND=LONG), DIMENSION(inp%Nsprate+1,inp%Nreactions) :: specnumdest0
    INTEGER(KIND=LONG) :: i2, is, ic, ic2

     i2 = 0

     specnumform0 = sp%numform
     specnumdest0 = sp%numdest
     
     DO is=1,inp%Nsprate

        ! compute the integrated contribution of each reaction
        contribformtot = 0d0
        contribdesttot = 0d0

        DO ic=1,sp%Nreacform(is)

           contribform0(is,ic) = (out%rateform(is,ic,1)+out%rateform(is,ic,2))/2d0*&
                (out%time(2)-out%time(1))*Rys
           DO i2=3,out%Nrealsteps
              !CALL avint(out%rateform(is,ic,1:i2),&
              !     out%time(1:i2)*Rys,i2,&
              !     out%time(1)*Rys,out%time(i2)*Rys,contribform0(is,ic))

              contribform0(is,ic) = contribform0(is,ic) + &
                   (out%rateform(is,ic,i2)+out%rateform(is,ic,i2-1))/2d0*&
                   (out%time(i2)-out%time(i2-1))*Rys
              !write(6,*) is,ic,i2,out%time(i2),out%rateform(is,ic,i2),&
              !contribform0(is,ic),contribform02(is,ic)
           END DO
           contribform2(is,ic) = contribform0(is,ic)
           contribformtot = contribformtot+contribform0(is,ic)

        END DO

        !contribform2(is,1:sp%Nreacform(is)) = contribform0(is,1:sp%Nreacform(is))
        contribform0(is,1:sp%Nreacform(is)) = &
             contribform0(is,1:sp%Nreacform(is))/contribformtot
        

        DO ic=1,sp%Nreacdest(is)

           contribdest0(is,ic) = (out%ratedest(is,ic,1)+out%ratedest(is,ic,2))/2d0*&
                (out%time(2)-out%time(1))*Rys
           DO i2=3,out%Nrealsteps
              !CALL avint(out%ratedest(is,ic,1:out%Nrealsteps),&
              !     out%time(1:out%Nrealsteps)*Rys,out%Nrealsteps,&
              !     out%time(1)*Rys,out%time(out%Nrealsteps)*Rys,contribdest0(is,ic))

              contribdest0(is,ic) = contribdest0(is,ic) + &
                   (out%ratedest(is,ic,i2)+out%ratedest(is,ic,i2-1))/2d0*&
                   (out%time(i2)-out%time(i2-1))*Rys
           END DO
           contribdest2(is,ic) = contribdest0(is,ic)
           contribdesttot = contribdesttot+contribdest0(is,ic)
        
          ! write(6,*) is,ic,contribdest0(is,ic),contribdesttot
        END DO


        !contribdest2(is,1:sp%Nreacdest(is)) = contribdest0(is,1:sp%Nreacdest(is))
        contribdest0(is,1:sp%Nreacdest(is)) = &
             contribdest0(is,1:sp%Nreacdest(is))/contribdesttot

        ! sort the reactions following their respective contributions
        CALL SORT(contribform0(is,1:sp%Nreacform(is)),&
             sp%numform(is,1:sp%Nreacform(is)),sp%Nreacform(is))
        CALL SORT(contribdest0(is,1:sp%Nreacdest(is)),&
             sp%numdest(is,1:sp%Nreacdest(is)),sp%Nreacdest(is))
     
        ! rearrange the reaction rates with the new order
        DO ic=1,sp%Nreacform(is) 
           DO ic2=1,sp%Nreacform(is)
              IF(sp%numform(is,ic) .eq. specnumform0(is,ic2)) THEN
                 out%rateform2(is,ic,1:out%Ninisteps) = &
                      out%rateform(is,ic2,1:out%Ninisteps)
                 out%numform2(is,ic,1:out%Ninisteps) = &
                      out%numform(is,ic2,1:out%Ninisteps)
                 out%contribform(is,ic) = contribform0(is,ic)
                 out%rateformtot(is,ic) = contribform2(is,ic2)
              END IF
           END DO
        END DO

        DO ic=1,sp%Nreacdest(is) 
           DO ic2=1,sp%Nreacdest(is)
              IF(sp%numdest(is,ic) .eq. specnumdest0(is,ic2)) THEN
                 out%ratedest2(is,ic,1:out%Ninisteps) = &
                      out%ratedest(is,ic2,1:out%Ninisteps)
                 out%numdest2(is,ic,1:out%Ninisteps) = &
                      out%numdest(is,ic2,1:out%Ninisteps)
                 out%contribdest(is,ic) = contribdest0(is,ic)
                 out%ratedesttot(is,ic) = contribdest2(is,ic2)
              END IF
           END DO
        END DO

        !DO ic=1,sp%Nreacform(is)
           !WRITE(6,*) is, ic, re%numform2(is,ic,1),re%contribform(is,ic)
!!$           specnumform0(is,ic), re%rateform(is,ic,out%Ninisteps-1),&
!!$                sp%numform(is,ic), re%rateform2(is,ic,out%Ninisteps-1)
        !END DO

     END DO

     
   END SUBROUTINE ORDER_RATES

   
  SUBROUTINE SORT(x, indent, Size)

    ! --------------------------------------------------------------------
    ! This subroutine receives an array x() and sorts it into ascending
    ! order.
    ! --------------------------------------------------------------------
    
    IMPLICIT  NONE
    
    INTEGER(KIND=LONG), INTENT(IN) :: Size
    REAL(KIND=DP), DIMENSION(Size), INTENT(INOUT) :: x
    INTEGER(KIND=LONG), DIMENSION(Size), INTENT(INOUT) :: indent
    INTEGER(KIND=LONG) :: i
    INTEGER(KIND=LONG) :: Location
    
    DO i = 1, Size-1			! except for the last
       Location = FINDMAXIMUM(x, Size, i, Size)	! find min from this to last
       CALL SWAP(x(i), x(Location), indent(i), indent(Location)) 
    END DO
    
  END SUBROUTINE SORT
  

  SUBROUTINE SWAP(a, b, c, d)

    ! --------------------------------------------------------------------
    ! SUBROUTINE  Swap():
    !    This subroutine swaps the values of its two formal arguments.
    ! --------------------------------------------------------------------
    
    IMPLICIT  NONE
    
    REAL(KIND=DP), INTENT(INOUT) :: a, b
    INTEGER(KIND=LONG), INTENT(INOUT) :: c, d
    REAL(KIND=DP) :: Temp
    INTEGER(KIND=LONG) :: Temp2
    
    Temp = a
    a    = b
    b    = Temp
    Temp2 = c
    c = d
    d = Temp2

  END SUBROUTINE SWAP

  
!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------
!! CHECKING ELEMENTS
!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------
   
   SUBROUTINE CHECK_ELEMS(it,checkcons)

     !---------------------------------------------------------------
     ! check if the system is conservative by looking 
     ! the elemental abundance and exit the layer loop if not
     !---------------------------------------------------------------

     USE SHARED_VARIABLES, only: ph, ab, sp 
     
     INTEGER(KIND=LONG), INTENT(in) :: it
     INTEGER(KIND=LONG), INTENT(out) :: checkcons
     REAL(KIND=DP), DIMENSION(Nelements) :: comp_abu
     INTEGER(KIND=LONG) :: ie

     checkcons = 1
     
     ! compare the elemental abundance at t and t=0
     DO ie=1,Nelements
        IF (ab%Xelem_ini(ie) .ne. 0d0) THEN
           IF (ie .ne. 1) THEN
              comp_abu(ie) = abs((ab%Xelem_ini(ie)-ab%Xelem(ie))/&
                   (ab%Xelem_ini(ie)+ab%Xelem(ie)))
              ! evolution of element abundance must remain lower than 0.1%
              IF (comp_abu(ie)*100 .gt. 0.1)  THEN 
                 checkcons = 2
              END IF
           END IF
        END IF
     END DO

     ! if not conservative: exit the loop
     IF (checkcons .eq. 2) THEN
        WRITE(6,*) "THE SYSTEM ISNT CONSERVATIVE ANYMORE, ",&
             "THE PRECISION HAS TO BE INCREASED"
        WRITE(6,*) "Step = ",it
        WRITE(6,*) "Time = ",ph%time," years"
        WRITE(6,*) "Elem            Xelem(ini)                Xelem(t)"
        DO ie=1,Nelements
           IF (ab%Xelem(ie) .ne. 0d0 .and. ie .ne. 1) THEN
              WRITE(6,*) sp%name_elem(ie), ab%Xelem_ini(ie), ab%Xelem(ie)
           END IF
        END DO
     END IF
     
   END SUBROUTINE CHECK_ELEMS




  !!-----------!!
  !!-----------!!
  !! FUNCTIONS !!
  !!-----------!!
  !!-----------!!

  FUNCTION FINDMAXIMUM(x, Size, Start, End)
    
    ! --------------------------------------------------------------------
    ! INTEGER FUNCTION  FindMinimum():
    !    This function returns the location of the minimum in the section
    ! between Start and End.
    ! --------------------------------------------------------------------
    
    IMPLICIT  NONE
    
    INTEGER(KIND=LONG) :: FindMaximum
    INTEGER(KIND=LONG), INTENT(IN) :: Size
    REAL(KIND=DP), DIMENSION(Size), INTENT(IN) :: x
    INTEGER(KIND=LONG), INTENT(IN) :: Start, End
    REAL(KIND=DP) :: Maximum
    INTEGER(KIND=LONG) :: Location
    INTEGER(KIND=LONG) :: i
    
    Maximum  = x(Start)		! assume the first is the min
    Location = Start		! record its position
    DO i = Start+1, End		! start with next elements
       IF (x(i) > Maximum) THEN	!   if x(i) less than the min?
          Maximum  = x(i)	!      Yes, a new minimum found
          Location = i          !      record its position
       END IF
    END DO
    FindMaximum = Location      ! return the position
    
  END FUNCTION FINDMAXIMUM
  
  
  ! Eb evolution with the coverage of H2 on the external layer
  REAL(KIND=DP) FUNCTION Eb_evol(il,PH2ini,PHini,Pspini,Pbareini,num_spec,sp,ns)
    
    REAL(KIND=DP), INTENT(in) :: PH2ini, PHini, Pspini, Pbareini
    TYPE(species), INTENT(in) :: sp
    TYPE(specnum), INTENT(in) :: ns
    INTEGER(KIND=LONG), INTENT(in) :: num_spec, il
    REAL(KIND=DP) :: EbH2O, EbH2, Ebpure, Ebcarb, PH2, PH, Psp, Pbare
    
    Ebcarb = sp%Eb_carb(num_spec)
    EbH2O = sp%Eb_wat(num_spec)
    EbH2 = sp%Eb_H2(num_spec)
    Ebpure = sp%Eb_pure(num_spec)

    PH2 = PH2ini
    PH = PHini
    Psp = Pspini
    Pbare = Pbareini

    IF ((PH2ini+PHini) .gt. 1) THEN 
       PH2 = 1d0
       PH = 0d0
    END IF
    IF (num_spec .eq. ns%JpH2 .or. num_spec .eq. ns%JoH2 .or. num_spec .eq. ns%JH) THEN
       IF (il .eq. 1) THEN
          Eb_evol = (1-PH2-PH)*Ebcarb + (PH2+PH)*EbH2
       ELSE IF (il .gt. 1) THEN
          Eb_evol = (1-PH2-PH)*EbH2O + (PH2+PH)*EbH2
       END IF
       !IF ((1-PH2-PH-Pbare) .lt. 0d0) Pbare = 1-PH2-PH
       !Eb_evol = (1-PH2-PH-Pbare)*EbH2O + (PH2+PH)*EbH2 + Pbare*Ebcarb
       !Eb_evol = 300.
       Eb_evol = (1-PH2-PH)*EbH2O + (PH2+PH)*EbH2
    ELSE
       Eb_evol = (1-PH2-PH-Pbare-Psp)*EbH2O + (PH2+PH)*EbH2 + Pbare*Ebcarb + Psp*Ebpure
       !if (trim(sp%name(num_spec)).eq.'JCO2') write(6,*) num_spec, Eb_evol,EbH2O, PH2, PH, Pbare, Psp
    END IF
    
    !Eb_evol = (1-PH2-PH-Pbare-Psp)*EbH2O + (PH2+PH)*EbH2 + Pbare*Ebcarb + Psp*Ebpure

    !write(6,*) num_spec, ns%JpH2, ns%JoH2, ns%JH


  END FUNCTION Eb_evol

  
  ! attempt rate nu0
  REAL(KIND=DP) FUNCTION nu0(num_spec)

    USE SHARED_VARIABLES, only: inp, sp
    
    INTEGER(KIND=LONG), INTENT(in) :: num_spec
    
    nu0 = sqrt(2*1/(inp%ds*1d-10)**2*sp%Eb(num_spec)*kb&
         /(Pi**2*sp%mass(num_spec)*Mpr))
    
  END FUNCTION nu0


  ! accretion rate: S*v(i)*sig*nd (s-1)
  REAL(KIND=DP) FUNCTION Racc(num_spec) 

    USE SHARED_VARIABLES, only: ph
    
    INTEGER(KIND=LONG), INTENT(in) :: num_spec

    Racc = ph%stick(num_spec)*ph%Vth(num_spec)*ph%signdtot

  END FUNCTION Racc

  ! accretion rate: S*v(i)*sig*nd (s-1)
  REAL(KIND=DP) FUNCTION Racc2(num_spec) 

    USE SHARED_VARIABLES, only: ph
    REAL(KIND=DP) :: Vcol, mgr, Tcrit, dEcrit, Rcol, Rev2
    
    INTEGER(KIND=LONG), INTENT(in) :: num_spec

    Vcol = 1. * 1e5 ! cm/s
    Racc2 = ph%stick(num_spec)*(ph%Vth(num_spec)+Vcol)*ph%signdtot

  END FUNCTION Racc2

  ! evaporation due to grain-grain collision
  REAL(KIND=DP) FUNCTION Revcol(num_spec,il) 

    USE SHARED_VARIABLES, only: ph, gr, inp, sp, ns
    
    INTEGER(KIND=LONG), INTENT(in) :: num_spec
    REAL(KIND=DP) :: Vcol, mgr, Tcrit, dEcrit, Rcol, Rev2

    Vcol = 1. * 1e3 ! collision velocity (m/s)
    mgr = inp%rhod*1e3*gr%Vd(1)*1e-18 ! kg m-3
    !mgr = mgr * 1e3 ! kg m-3

    dEcrit = 0.5 * 0.5 * mgr * Vcol**2 * 1e7 ! ergs
    !Tcrit = 30. * (dEcrit/5.4e-8)**(1/3.) * (gr%a(il)/2./0.1)
    Tcrit = 30. * (dEcrit/1.04e-7)**(1/2.3) * (gr%a(il)/2./0.1)**(-3/2.3)
    !IF (Tcrit.ge.100.) Tcrit = 100.
    !write(6,*) ' '
    !write(6,*) mgr, dEcrit, Vcol, Tcrit !, gr%a(il)
 
    Rcol = Vcol*1e2*ph%signdtot
    Rev2 = Revth(Tcrit,num_spec)
    IF (Rev2.ge.1d5) Rev2 = 1d5
    Revcol = Rcol !1d-5*Rcol*Rev2
    !IF (num_spec .eq. ns%JCH3OH) &
    !     write(6,*) Revcol, Rcol, ph%signdtot, Tcrit, Rev2, Revth(100d0,num_spec)
    
  END FUNCTION Revcol

  ! accretion rate: 1/4*Pi*S*v(i)*(dsite)^2 (cm3.s-1) to get in MLs/sec
  REAL(KIND=DP) FUNCTION RaccML(num_spec) 

    USE SHARED_VARIABLES, only: inp, ph
    
    INTEGER(KIND=LONG), INTENT(in) :: num_spec

    RaccML = 1/4d0*ph%stick(num_spec)*(inp%ds*1E-8)**2*ph%Vth(num_spec)

  END FUNCTION RaccML


  ! thermal evaporation rate 
  REAL(KIND=DP) FUNCTION Revth(temp,num_spec)

    USE SHARED_VARIABLES, only: inp, sp
    
    REAL(KIND=DP), INTENT(in) :: temp
    INTEGER(KIND=LONG), INTENT(in) :: num_spec
    REAL(KIND=DP) :: attempt

    Revth = nu0(num_spec)*exp(-sp%Eb(num_spec)/temp)

  END FUNCTION Revth

  
  ! evaporation rate induced by grain-heating from CRs from HH93
  REAL(KIND=DP) FUNCTION RevthCR(num_spec,zeta)

    USE SHARED_VARIABLES, only: inp, sp
    
    INTEGER(KIND=LONG), INTENT(in) :: num_spec
    REAL(KIND=DP), INTENT(in) :: zeta

    RevthCR = 3.6d-19*Revth(7d1,num_spec) * zeta/(3d-17)
    
  END FUNCTION RevthCR


  ! absorption rate of UV photons by water ice from MD simulations
  REAL(KIND=DP) FUNCTION RabsUV(il,num_spec)

    USE SHARED_VARIABLES, only: inp, ph, gr
    
    REAL(KIND=DP) :: UVflux
    INTEGER(KIND=LONG), INTENT(in) :: num_spec, il
    
    !yieldUV = (5.2d-3*ph%UVfluxCR+7.5d-3*FUV(inp%G0ext,ph%Av))
    UVflux = ph%UVfluxCR+FUV(inp%G0ext,ph%Av)
    RabsUV = UVflux*(inp%ds*1e-8)**2 !(4*Pi*gr%a(il)*1d-4/2d0)**2/gr%Nstot(il)
    !write(6,*) RabsUV, ph%UVfluxCR+FUV(inp%G0ext,ph%Av)
    !RabsUV = UVflux*4*Pi*(gr%a(il)*1d-4/2d0)**2/gr%Nstot(il)
    !write(6,*) 4*Pi*(gr%a(il)*1d-4/2d0)**2/gr%Nstot(il),(inp%ds*1e-8)**2
  END FUNCTION RabsUV

  ! sputtering by CRs scaled on H2O following Dartois et al. (2015) experiments
  REAL(KIND=DP) FUNCTION RsputCR(il,num_spec)

    USE SHARED_VARIABLES, only: inp, ph, gr, sp, ns
    
    REAL(KIND=DP) :: H2Osput
    INTEGER(KIND=LONG), INTENT(in) :: num_spec, il

    H2Osput = 8.3 + (4.1-8.3)*(ph%zeta-5.89e-17)/(2.12e-17-5.89e-17)
    
    RsputCR = H2Osput*(inp%ds*1e-8)**2 !*sp%Eb(ns%JH2O)/sp%Eb(num_spec) !(4*Pi*gr%a(il)*1d-4/2d0)**2/gr%Nstot(il)
    !write(6,*) sp%name(num_spec), RsputCR
    !write(6,*) RabsUV, ph%UVfluxCR+FUV(inp%G0ext,ph%Av)
    !RabsUV = UVflux*4*Pi*(gr%a(il)*1d-4/2d0)**2/gr%Nstot(il)
    !write(6,*) 4*Pi*(gr%a(il)*1d-4/2d0)**2/gr%Nstot(il),(inp%ds*1e-8)**2
  END FUNCTION RsputCR

  
  ! absorption rate of UV photons by water ice from MD simulations
  REAL(KIND=DP) FUNCTION RabsUVold(il,num_spec)
    
    USE SHARED_VARIABLES, only: inp, ph, gr
    
    REAL(KIND=DP) :: yieldUV
    INTEGER(KIND=LONG), INTENT(in) :: num_spec, il
    
    yieldUV = (5.2d-3*ph%UVfluxCR+7.5d-3*FUV(inp%G0ext,ph%Av))

    RabsUVold = 5*yieldUV*Pi*(gr%a(il)*1d-4/2d0)**2/gr%Nstot(il)
    
  END FUNCTION RabsUVold
  

    ! photodesorption rate from experiments or MD simulations
  REAL(KIND=DP) FUNCTION RevUV(il,num_spec)

    USE SHARED_VARIABLES, only: inp, ph, sp, ns, re, nr, gr
    
    REAL(KIND=DP) :: yieldUV, yieldUVH2O, Acore, Aedge
    INTEGER(KIND=LONG), INTENT(in) :: num_spec, il
    
    IF (trim(sp%name(num_spec)) .eq. 'JCO' .or. &
         trim(sp%name(num_spec)) .eq. 'JO2' .or. &
         trim(sp%name(num_spec)) .eq. 'JO3' .or. &
         trim(sp%name(num_spec)) .eq. 'JN2') THEN
        IF (trim(sp%name(num_spec)) .eq. 'JCO') THEN
           Acore = 1d-2
           Aedge = 1.3d-2
        ELSE IF (trim(sp%name(num_spec)) .eq. 'JO2' .or. &
             trim(sp%name(num_spec)) .eq. 'JO3') THEN
           Acore = 2.6d-3
           Aedge = 3.3d-3
        ELSE IF (trim(sp%name(num_spec)) .eq. 'JN2') THEN
           Acore = 2.2d-3
           Aedge = 2.6d-3
        END IF
       yieldUV = Acore*ph%UVfluxCR + Aedge*FUV(inp%G0ext,ph%Av)
    ELSE IF (trim(sp%name(num_spec)) .eq. 'JCO2' .or. &
         trim(sp%name(num_spec)) .eq. 'JtHOCO' .or. &
         trim(sp%name(num_spec)) .eq. 'JcHOCO' .or. &
         trim(sp%name(num_spec)) .eq. 'JtDOCO' .or. &
         trim(sp%name(num_spec)) .eq. 'JcDOCO' .or. &
         trim(sp%name(num_spec)) .eq. 'JHCOOH' .or. &
         trim(sp%name(num_spec)) .eq. 'JHCOOD' .or. &
         trim(sp%name(num_spec)) .eq. 'JDCOOH' .or. &
         trim(sp%name(num_spec)) .eq. 'JDCOOD') THEN
       yieldUV = 1.2*1d-3*(1-exp(-il/2.9d0))+1.1d-3*(1-exp(-il/4.6d0))*&
            (ph%UVfluxCR + FUV(inp%G0ext,ph%Av))
    ELSE 
       yieldUVH2O = re%C(nr%ref_photo)*&
            (5.2d-3*ph%UVfluxCR+7.5d-3*FUV(inp%G0ext,ph%Av))
       yieldUV = yieldUVH2O * sp%Eb(ns%JH2O)/sp%Eb(num_spec)
    END IF

    RevUV = yieldUV*Pi*(4*Pi*gr%a(il)*1d-4/2d0)**2/gr%Nstot(il)
    
  END FUNCTION RevUV

  
    ! photodesorption rate from experiments or MD simulations
  REAL(KIND=DP) FUNCTION Rsurfphoto(il,ic,num_spec)

    USE SHARED_VARIABLES, only: inp, ph, sp, ns, re, nr, gr
    
    REAL(KIND=DP) :: yieldUV, yieldUVH2O, Acore, Aedge, UVflux
    INTEGER(KIND=LONG), INTENT(in) :: num_spec, il, ic
    
!!$    IF (trim(sp%name(num_spec)) .eq. 'JCO' .or. &
!!$         trim(sp%name(num_spec)) .eq. 'JO2' .or. &
!!$         trim(sp%name(num_spec)) .eq. 'JO3' .or. &
!!$         trim(sp%name(num_spec)) .eq. 'JN2') THEN
!!$        IF (trim(sp%name(num_spec)) .eq. 'JCO') THEN
!!$           Acore = 1d-2
!!$           Aedge = 1.3d-2
!!$        ELSE IF (trim(sp%name(num_spec)) .eq. 'JO2' .or. &
!!$             trim(sp%name(num_spec)) .eq. 'JO3') THEN
!!$           Acore = 2.6d-3
!!$           Aedge = 3.3d-3
!!$        ELSE IF (trim(sp%name(num_spec)) .eq. 'JN2') THEN
!!$           Acore = 2.2d-3
!!$           Aedge = 2.6d-3
!!$        END IF
!!$       yieldUV = Acore*ph%UVfluxCR + Aedge*FUV(inp%G0ext,ph%Av)
!!$    ELSE IF (trim(sp%name(num_spec)) .eq. 'JCO2' .or. &
!!$         trim(sp%name(num_spec)) .eq. 'JtHOCO' .or. &
!!$         trim(sp%name(num_spec)) .eq. 'JcHOCO' .or. &
!!$         trim(sp%name(num_spec)) .eq. 'JtDOCO' .or. &
!!$         trim(sp%name(num_spec)) .eq. 'JcDOCO' .or. &
!!$         trim(sp%name(num_spec)) .eq. 'JHCOOH' .or. &
!!$         trim(sp%name(num_spec)) .eq. 'JHCOOD' .or. &
!!$         trim(sp%name(num_spec)) .eq. 'JDCOOH' .or. &
!!$         trim(sp%name(num_spec)) .eq. 'JDCOOD') THEN
!!$       yieldUV = 1.2*1d-3*(1-exp(-il/2.9d0))+1.1d-3*(1-exp(-il/4.6d0))*&
!!$            (ph%UVfluxCR + FUV(inp%G0ext,ph%Av))
!!$    ELSE 
!!$       yieldUVH2O = re%C(nr%ref_photo)*&
!!$            (5.2d-3*ph%UVfluxCR+7.5d-3*FUV(inp%G0ext,ph%Av))
!!$       yieldUV = yieldUVH2O * sp%Eb(ns%JH2O)/sp%Eb(num_spec)
!!$    END IF


    IF (re%A(ic) .ne. 1 .or. re%B(ic) .ne. 1 .or. re%C(ic) .ne. 1) THEN
       UVflux = ph%UVfluxCR+FUV(inp%G0ext,ph%Av)
       Rsurfphoto = UVflux*Pi*(4*Pi*gr%a(il)*1d-4/2d0)**2/gr%Nstot(il)
    ELSE IF (re%A(ic) .eq. 1 .and. re%B(ic) .eq. 1 .and. re%C(ic) .eq. 1) THEN
       yieldUVH2O = re%C(nr%ref_photo)*&
            (5.2d-3*ph%UVfluxCR+7.5d-3*FUV(inp%G0ext,ph%Av))
       yieldUV = yieldUVH2O * sp%Eb(ns%JH2O)/sp%Eb(num_spec)
       Rsurfphoto = yieldUV*Pi*(4*Pi*gr%a(il)*1d-4/2d0)**2/gr%Nstot(il)
    END IF
    
  END FUNCTION Rsurfphoto

  

  ! cosmic-ray induced rate
  REAL(KIND=DP) FUNCTION Rcr(A,zeta)
    
    REAL(KIND=DP), INTENT(in) :: A, zeta
    
    Rcr = A*zeta
    
  END FUNCTION Rcr

  ! photo-process rate
  REAL(KIND=DP) FUNCTION Rphoto(A,C,Av)
    
    REAL(KIND=DP), INTENT(in) :: A, C, Av
    
    Rphoto = A*exp(-C*Av)
    
  END FUNCTION Rphoto

  ! Kooij rate
  REAL(KIND=DP) FUNCTION Rkooij(A,B,C,Tg)
    
    REAL(KIND=DP), INTENT(in) :: A, B, C, Tg
    
    Rkooij = A*(Tg/3d2)**B*exp(-C/Tg)
    
  END FUNCTION Rkooij

  ! ionpol1 rate
  REAL(KIND=DP) FUNCTION Rionpol1(A,B,C,Tg)
    
    REAL(KIND=DP), INTENT(in) :: A, B, C, Tg
    
    Rionpol1 = A*B*(0.62+0.4767*C*sqrt(3d2/Tg))
    
  END FUNCTION Rionpol1

  ! ionpol2 rate
  REAL(KIND=DP) FUNCTION Rionpol2(A,B,C,Tg)
    
    REAL(KIND=DP), INTENT(in) :: A, B, C, Tg
    
    Rionpol2 = A*B*(1+0.0967*C*sqrt(3d2/Tg)+C**2/10.526*3d2/Tg)
    
  END FUNCTION Rionpol2

  
  ! LH reaction rate
  REAL(KIND=DP) FUNCTION Rsurf(ic,il)

    USE SHARED_VARIABLES, only: inp, ph, sp, ns, re, nr, gr, ss, ab
    
   INTEGER(KIND=LONG), INTENT(in) :: ic, il
    REAL(KIND=DP) :: freac
    
    freac=0d0
    
    ! 3-body reaction for the formation of CO2 (see Taquet et al. 13)
    IF ((re%num_react1(ic) .eq. ns%JCO .and. &
         re%num_react2(ic) .eq. ns%JO .and. &
         re%num_react3(ic) .eq. ns%JH) .or. &
         (re%num_react1(ic) .eq. ns%JCO .and. &
         re%num_react2(ic) .eq. ns%JO .and. &
         re%num_react3(ic) .eq. ns%JD) .and. &
         nr%reacCO2 .ne. 0) THEN
       ! reaction probability
       IF (inp%chreacproba .eq. 1) THEN
          ! reaction proba directly deduced from transmission proba
          freac = max(re%proba(ic),exp(-re%C(ic)/ph%Td))
       ELSE IF (inp%chreacproba .eq. 2) THEN
          ! reaction proba deduced from competition process
          freac = max(nu0(re%num_react2(ic)),&
               nu0(re%num_react3(ic)))*re%proba(ic)/&
               (max(nu0(re%num_react2(ic)),&
               nu0(re%num_react3(ic)))*re%proba(ic)+&
               Rdiff(ph%Td,re%num_react2(ic))+&
               Rdiff(ph%Td,re%num_react3(ic)))
       END IF
       Rsurf = re%A(ic)*re%B(ic)*freac*&
            (Rdiff(ph%Td,re%num_react2(ic))+&
            Rdiff(ph%Td,re%num_react3(ic)))/&
            (gr%Fnptot(il)*ph%Nsndtot)
       ! Standard 2-body Langmuir-Hinshelwood surface reaction
    ELSE
       ! reaction probability
       IF (inp%chreacproba .eq. 1) THEN
          ! reaction proba directly deduced from transmission proba
          freac = max(re%proba(ic),exp(-re%C(ic)/ph%Td))
       ELSE IF (inp%chreacproba .eq. 2) THEN
          ! reaction proba deduced from competition process
          freac = max(nu0(re%num_react1(ic)),&
               nu0(re%num_react2(ic)))*re%proba(ic)/&
               (max(nu0(re%num_react1(ic)),&
               nu0(re%num_react2(ic)))*re%proba(ic)+&
               Rdiff(ph%Td,re%num_react1(ic))+&
               Rdiff(ph%Td,re%num_react2(ic)))
       END IF
       Rsurf = re%A(ic)*re%B(ic)*freac*&
            (Rdiff(ph%Td,re%num_react1(ic))+&
            Rdiff(ph%Td,re%num_react2(ic)))/&
            (gr%Fnptot(il)*ph%Nsndtot)
!!$       if (ic .eq. 8170) then &
!!$            re%A(ic)*re%B(ic)*freac*&
!!$            (Rdiff(ph%Td,re%num_react1(ic))+&
!!$            Rdiff(ph%Td,re%num_react2(ic)))/&
!!$            (gr%Fnptot(il)*ph%Nsndtot)
    END IF

    !IF (ISNAN(Rsurf)) write(6,*) Rsurf, freac,Rdiff(ph%Td,re%num_react1(ic)),&
    !     Rdiff(ph%Td,re%num_react2(ic)),gr%Fnptot(il),ph%Nsndtot
    
  END FUNCTION Rsurf

  
  ! probability of chemical desorption upon formation
  REAL(KIND=DP) FUNCTION Pevreac(num_reac)
    
    REAL(KIND=DP) :: P, a
    INTEGER(KIND=LONG), INTENT(in) :: num_reac

    a = 0.012
    
    P = 1
    Pevreac = a*P/(1+a*P)
    
  END FUNCTION Pevreac



  ! diffusion rate by thermal hopping (+ tunnelling for H if needed)
  REAL(KIND=DP) FUNCTION Rdiff(temp,num_spec)

    USE SHARED_VARIABLES, only: inp, sp, ns
    
    REAL(KIND=DP), INTENT(in) :: temp
    INTEGER(KIND=LONG), INTENT(in) :: num_spec
    REAL(KIND=DP) :: attempt
    REAL(KIND=DP) :: diffth_rate
    REAL(KIND=DP) :: difftun_rate


    IF ((num_spec == ns%JH .or. num_spec .eq. ns%JD .or. &
         num_spec == ns%JHD .or. num_spec .eq. ns%JH2 .or. & 
         num_spec == ns%JD2) .and. inp%chHtunnel == 1) THEN
       diffth_rate =  nu0(num_spec)*exp(-sp%Ed(num_spec)/temp)
       Rdiff = diffth_rate
    ELSE IF ((num_spec == ns%JH .or. num_spec .eq. ns%JD .or. &
         num_spec == ns%JHD .or. num_spec .eq. ns%JH2 .or. & 
         num_spec == ns%JD2) .and. inp%chHtunnel == 2) THEN 
       difftun_rate =  nu0(num_spec)*exp(-(2*inp%diffwid*1D-10/hb)*&
            sqrt(2*sp%mass(num_spec)*Mpr*sp%Ed(num_spec)*kb))
       Rdiff = difftun_rate
    ELSE IF ((num_spec == ns%JH .or. num_spec .eq. ns%JD .or. &
         num_spec == ns%JHD .or. num_spec .eq. ns%JH2 .or. & 
         num_spec == ns%JD2) .and. inp%chHtunnel == 3) THEN 
       diffth_rate =  nu0(num_spec)*exp(-sp%Ed(num_spec)/temp)
       difftun_rate =  nu0(num_spec)*exp(-(2*inp%diffwid*1D-10/hb)*&
            sqrt(2*sp%mass(num_spec)*Mpr*sp%Ed(num_spec)*kb))
       Rdiff = diffth_rate + difftun_rate
    ELSE IF ((num_spec == ns%JC .or. num_spec .eq. ns%JN .or. &
         num_spec == ns%JO) .and. inp%chCtunnel == 2) THEN 
       difftun_rate =  nu0(num_spec)*exp(-(2*inp%diffwid*1D-10/hb)*&
            sqrt(2*sp%mass(num_spec)*Mpr*sp%Ed(num_spec)*kb))
       Rdiff = difftun_rate
    ELSE IF ((num_spec == ns%JC .or. num_spec .eq. ns%JN .or. &
         num_spec == ns%JO) .and. inp%chCtunnel == 3) THEN 
       diffth_rate =  nu0(num_spec)*exp(-sp%Ed(num_spec)/temp)
       difftun_rate =  nu0(num_spec)*exp(-(2*inp%diffwid*1D-10/hb)*&
            sqrt(2*sp%mass(num_spec)*Mpr*sp%Ed(num_spec)*kb))
       Rdiff = diffth_rate + difftun_rate
    ELSE 
       diffth_rate =  nu0(num_spec)*exp(-sp%Ed(num_spec)/temp)
       Rdiff = diffth_rate
    END IF
    
  END FUNCTION Rdiff

  
  ! transmission probability through rectangular barrier
  REAL(KIND=DP) FUNCTION proba_square(num_reac)

    USE SHARED_VARIABLES, only: inp, sp, re
    
    INTEGER(KIND=LONG), INTENT(in) :: num_reac
    REAL(KIND=DP) :: red_mass, proba_ref,red_mass_ref, Ea_ref

    red_mass = sp%mass(re%num_react1(num_reac))*&
         sp%mass(re%num_react2(num_reac))/&
         (sp%mass(re%num_react1(num_reac)) + &
         sp%mass(re%num_react2(num_reac)))*Mpr

    proba_square = exp(-(2*inp%rectwid*1E-10/hb)*&
         sqrt(2*red_mass*re%C(num_reac)*kb))
    
  END FUNCTION proba_square
       

  ! reaction rates of the methanol network from Caselli+02 or Watanabe+08
  REAL(KIND=DP) FUNCTION proba_COnetwork(num_reac)

    USE SHARED_VARIABLES, only: inp, sp, ns, re, nr, con
    
    INTEGER(KIND=LONG), INTENT(in) :: num_reac
    REAL(KIND=DP) :: proba_ref, Ea_ref, Ea, red_mass
    INTEGER(KIND=LONG) :: i

    proba_COnetwork = re%proba(num_reac)

    ! probabilities given by the Watanabe's relative rates
    IF (inp%chDnetwork .eq. 2) THEN 

       proba_ref = re%proba(nr%COH)

       DO i=1,NCOnetreacs

          IF (trim(re%react1(num_reac)) .eq. trim(con%react1(i)) .and. &
            trim(re%react2(num_reac)) .eq. trim(con%react2(i)) .and. &
            trim(re%prod1(num_reac)) .eq. trim(con%prod1(i))) THEN

             proba_COnetwork = proba_ref*con%A(i)

          END IF
       END DO


    ! probabilities given by the Caselli's energies
    ELSE IF (inp%chDnetwork .eq. 3) THEN

       Ea_ref = re%C(nr%COH)

       DO i=1,NCOnetreacs
          IF (trim(re%react1(num_reac)) .eq. trim(con%react1(i)) .and. &
            trim(re%react2(num_reac)) .eq. trim(con%react2(i)) .and. &
            trim(re%prod1(num_reac)) .eq. trim(con%prod1(i))) THEN

             red_mass = sp%mass(re%num_react1(num_reac))*&
                  sp%mass(re%num_react2(num_reac))/&
                  (sp%mass(re%num_react1(num_reac)) + &
                  sp%mass(re%num_react2(num_reac)))*Mpr

             Ea = Ea_ref*con%B(i)

             proba_COnetwork = exp(-(2*inp%rectwid*1E-10/hb)*&
                  sqrt(2*red_mass*Ea*kb))

          END IF
       END DO

    END IF


  END FUNCTION proba_COnetwork


  ! transmission probability from the Eckart model
  REAL(KIND=DP) FUNCTION proba_eckart(num_reac)

    USE SHARED_VARIABLES, only: inp, sp, ns, re, ec
    
    INTEGER(KIND=LONG), INTENT(in) :: num_reac
    REAL(KIND=DP) :: red_mass, dv1, dv2, nus, mu, fs, A, B, L, x, y, vx
    REAL(KIND=DP) :: En, alpha1, alpha2, ksi, frac, aa, bb, dd, nu00
    REAL(KIND=DP) :: Ktoj
    REAL(KIND=DP), DIMENSION(51) :: kapa
    INTEGER(KIND=LONG) :: ic, i, num_react1, num_react2
    CHARACTER(len=15) :: react1, react2, prod1, prod2

    proba_eckart = re%proba(num_reac)

    react1 = re%react1(num_reac)
    react2 = re%react2(num_reac)
    prod1 = re%prod1(num_reac)
    prod2 = re%prod2(num_reac)


!!$    IF (re%num_prod1(num_reac) .gt. inp%Nspgas .and. re%num_prod1(num_reac) .le. inp%Nspgas+inp%Nspgr) THEN 
!!$       prod1 = re%prod1(num_reac)
!!$       prod2 = re%prod2(num_reac)
!!$       ! if product is gaseous (via chemical desorption), re-adjust to match
!!$    ELSE IF (re%num_prod1(num_reac) .le. inp%Nspgas) THEN 
!!$       prod1 = 'J'//re%prod1(num_reac)
!!$       IF (re%prod2(num_reac) .ne. '          ') prod2 = 'J'//re%prod2(num_reac)
!!$    ELSE IF (re%num_prod1(num_reac) .ge. inp%Nspgas+inp%Nspgr) THEN 
!!$       prod1 = re%prod1(num_reac-inp%Nspgr) !'J'//re%prod1(num_reac)
!!$       IF (re%prod2(num_reac) .ne. '          ') &
!!$            prod2 = re%prod2(num_reac-inp%Nspgr) !'J'//re%prod2(num_reac)
!!$    END IF

    !write(6,*) num_reac, (react1),(react2),(prod1),(prod2)
    
    IF (react1(1:1) .eq. 'J' .and. prod1(1:1) .eq. 'J') THEN 
       prod1 = re%prod1(num_reac)
       prod2 = re%prod2(num_reac)
       ! if product is gaseous (via chemical desorption), re-adjust to match
    ELSE IF (react1(1:1) .eq. 'J' .and. prod1(1:1) .ne. 'J') THEN 
       prod1 = re%prod1(num_reac-1)
       prod2 = re%prod2(num_reac-1)
    ELSE IF (react1(1:1) .eq. 'Q') THEN 
       react1 = re%react1(num_reac-2)
       react2 = re%react2(num_reac-2)
       prod1 = re%prod1(num_reac-2)
       prod2 = re%prod2(num_reac-2)
    END IF

    DO is=1,inp%Nspnd
       IF (trim(react1).eq.trim(re%ndreact1(is,1))) THEN
          !write(6,*) trim(react1),' ',trim(re%ndreact1(is,1)),' ',trim(re%ndreact1(is,2))
          react1 = re%ndreact1(is,2)
       END IF
    END DO

    !write(6,*) num_reac, (react1),(react2),(prod1),(prod2)
    
    DO ic=1,Neckartreacs

       !write(6,*) ic, ec%react1(ic), ec%react2(ic), ec%prod1(ic), ec%prod2(ic)
       
       ! find the good reaction from data file
       IF (trim(react1) .eq. trim(ec%react1(ic)) .and. &
            trim(react2) .eq. trim(ec%react2(ic)) .and. &
            trim(prod1) .eq. trim(ec%prod1(ic)) .and. &
            trim(prod2) .eq. trim(ec%prod2(ic)) .or. &
            trim(react1) .eq. trim(ec%react2(ic)) .and. &
            trim(react2) .eq. trim(ec%react1(ic)) .and. &
            trim(prod1) .eq. trim(ec%prod1(ic)) .and. &
            trim(prod2) .eq. trim(ec%prod2(ic)) .or. &
            trim(react1) .eq. trim(ec%react1(ic)) .and. &
            trim(react2) .eq. trim(ec%react2(ic)) .and. &
            trim(prod1) .eq. trim(ec%prod2(ic)) .and. &
            trim(prod2) .eq. trim(ec%prod1(ic)) .or. &
            trim(react1) .eq. trim(ec%react1(ic)) .and. &
            trim(react2) .eq. trim(ec%react2(ic)) .and. &
            trim(prod1) .eq. trim(ec%prod2(ic)) .and. &
            trim(prod2) .eq. trim(ec%prod1(ic)) .and. &
            ec%Vf(ic) .gt. 0) THEN

          num_react1 = re%num_react1(num_reac)
          num_react2 = re%num_react2(num_reac)

          ! specific case: O..CO + H -> CO2 + H 
          IF (trim(react1) .eq. 'JOCO') THEN
             num_react1 = ns%JCO
             IF (trim(react2) .eq. 'JH') THEN
                num_react2 = ns%JOH
             ELSE IF (trim(react2) .eq. 'JD') THEN
                num_react2 = ns%JOD                
             END IF
          END IF

          ! compute the reduced mass
          IF (num_react2 .ne. inp%Nspecies+1) THEN
             red_mass = sp%mass(num_react1)*sp%mass(num_react2)/&
                  (sp%mass(num_react1) + sp%mass(num_react2))
          ELSE IF (num_react2 .eq. inp%Nspecies+1 .or. num_react2 .eq. 0) THEN
             red_mass = 1
          END IF
          
          ! convert the input parameters in the right unity
          dv1 = ec%Vf(ic)*kb
          dv2 = ec%Vr(ic)*kb
          mu = red_mass

          IF (re%type(num_reac) .eq. 14 .or. re%type(num_reac) .eq. 15) THEN
             nu00 = max(nu0(re%num_react1(num_reac)),&
                  nu0(re%num_react2(num_reac)))
          ELSE IF (re%type(num_reac) .eq. 80) THEN
             nu00 = nu0(re%num_react1(num_reac))
             IF (trim(re%react1(num_reac)) .eq. 'JHCOH') nu00=3.75340e+13
             IF (trim(re%react1(num_reac)) .eq. 'JCH3O') nu00=1.49896e+13
          END IF
          nus = ec%wf(ic)*cmtos
          fs = mu*((nus*nus)/conv)
          mu = mu*Mpr
          fs = fs*mdyneatonm

          ! compute constants
          A=dv1-dv2
          B=(sqrt(dv2)+sqrt(dv1))**2
          L=(((pi*pi)*((A**2-B**2)**2))/(2*fs*B**3))
       
          i=0
          ! vibrational energy within harmonic approximation
          En =hb*(2*Pi)*nu00*(i+0.5d0) 
          
          alpha1 = 2.d0*pi*dv1/(hb*(2*Pi)*nus)
          alpha2 = 2.d0*pi*dv2/(hb*(2*Pi)*nus)
          ksi = En/dv1
          
          frac = 1.d0/(alpha1**(-0.5)+alpha2**(-0.5))
          
          aa = (2.d0*sqrt(alpha1*ksi)*frac)/(2.d0*pi)
          bb = (2.d0*sqrt((1.d0-ksi)*(-alpha1)+alpha2)*frac)/(2.d0*pi)
          dd = 2.d0*sqrt(alpha1*alpha2-4*pi**2/16.d0)/(2.d0*pi)
          
          ! compute the transmission proba for fundamental energy
          kapa(1)=2.d0*dsinh(2.d0*pi*aa)*dsinh(2.d0*pi*bb)/&
               (dcosh(2.d0*pi*(aa+bb))+dcosh(2.d0*pi*dd))

          proba_eckart = kapa(1)

          IF (dv1 .gt. dv2) proba_eckart = 0d0


          !write(6,*) ic,react1, react2, prod1, prod2, red_mass, proba_eckart!, dv1,dv2,nus,mu,nu00,En,ksi
          
       ! probability = 1 when Ea = 0
       ELSE IF (trim(react1) .eq. trim(ec%react1(ic)) .and. &
            trim(react2) .eq. trim(ec%react2(ic)) .and. &
            trim(prod1) .eq. trim(ec%prod1(ic)) .and. &
            trim(prod2) .eq. trim(ec%prod2(ic)) .and. &
            ec%Vf(ic) .eq. 0) THEN
          proba_eckart = 1
       END IF

    END DO

    !write(6,*) 'end', num_reac

  END FUNCTION proba_eckart



END MODULE CHEMPROP
