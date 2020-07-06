MODULE PHYSPROP

  !----------------------------------------------------------------------------
  ! This module contains several subroutines that: 
  ! compute the grain properties: size, abundance, number of sites, volume, etc
  ! compute the evolution of grain properties with the mantle growth
  ! compute the evolution of physical properties (density, temperature, etc) 
  !----------------------------------------------------------------------------

  USE VARIABLES
  USE INPUT

CONTAINS

  SUBROUTINE GRAIN_PROP(grinput,grab)

    !---------------------------------------------------------------
    ! compute the abundance and other properties of each grain size 
    !---------------------------------------------------------------
  
    USE SHARED_VARIABLES, only: inp, gr, ph, out

    INTEGER(KIND=LONG), INTENT(in) :: grinput
    REAL(KIND=DP), INTENT(out) :: grab
    REAL(KIND=DP) :: Npore2, diam, abund
    REAL(KIND=DP) :: distrib_graintot
    REAL(KIND=DP) :: a, ai, da
    INTEGER(KIND=LONG) :: i


    ph%Nsndtot = 0d0
    ph%signdtot = 0d0
       
    !mass_ratio = 0
    diam = inp%acst
    abund = 6*inp%Rdg*Mpr*1d3/(Pi*(diam*1E-4)**3*inp%rhod)
    gr%a(1) = diam
    gr%Nssur(1) = 4*Pi*((diam)/2d0)**2/(inp%ds*1E-4)**2
    IF (out%Ninilayer .gt. 0) THEN
       DO i=1,out%Ninilayer
          IF (inp%chgrainevol .eq. 1) THEN
             gr%a(i) = gr%a(1) + 2*inp%laywid*1d-4*i
          ELSE IF (inp%chgrainevol .eq. 2) THEN
             gr%a(i) = gr%a(1)
          END IF
          gr%Nssur(i) = 4*Pi*(gr%a(i)/2d0)**2/(inp%ds*1E-4)**2
       END DO
    END IF
    gr%Vd(1) = 4*Pi/3d0*(diam/2d0)**3
    ! compute the porosity parameters and their influence on the grain surface
    IF (inp%Fporsur .gt. 0) THEN
       Npore2 = inp%Fporsur*gr%Nssur(1)/(inp%Slat**2)
       gr%Sdep(1) = gr%Vd(1)*inp%Fvacuum/(Npore2*inp%Slat**2*(inp%ds*1d-4)**3)
       gr%Vvacuum(1) = Npore2 * inp%Slat**2*gr%Sdep(1)*(inp%ds*1d-4)**3
       gr%Spore(1) = inp%Slat**2 + 4*inp%Slat*gr%Sdep(1)
       gr%Nstot(1) = gr%Nssur(1)*(inp%Fnpsur + &
            gr%Spore(1)*inp%Fporsur/inp%Slat**2)
       gr%Fnptot(1) = inp%Fnpsur * gr%Nssur(1)/gr%Nstot(1)
       gr%Fportot(1) = 1d0 - gr%Fnptot(1)
       gr%Fedgetot(1) = inp%Fedgesur * gr%Nssur(1)/gr%Nstot(1)
    ! no porosity -> no modification
    ELSE IF (inp%Fporsur .eq. 0) THEN
       Npore2 = 0d0
       gr%Sdep(1) = 0d0
       gr%Vvacuum(1) = 0d0
       gr%Spore(1) = 0d0
       gr%Nstot(1) = gr%Nssur(1)
       gr%Fnptot(1) = 1d0
       gr%Fportot(1) = 0d0
       gr%Fedgetot(1) = 0d0
    END IF
    gr%X = abund
    !ph%DTOGN = abund
    !ph%Nsndtot = gr%Nstot(1)*gr%X*nH
    !ph%signdtot = Pi*(gr%a(1)*1d-4/2d0)**2*gr%X*nH
    
    IF (grinput .eq. 1) THEN
       grab = abund
    END IF

  END SUBROUTINE GRAIN_PROP


  SUBROUTINE GRAIN_EVOLUTION(il)

    !----------------------------------------------------------------------
    ! compute the evolution of the grain properties with the mantle growth
    !----------------------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, gr, ph,ab

    INTEGER(KIND=LONG) :: k
    REAL(KIND=DP) :: Npore2, aprev
    INTEGER(KIND=LONG), INTENT(in) :: il

    IF (inp%chgrainevol .eq. 1) THEN
       IF (il .eq. 1) aprev = inp%acst
       IF (il .gt. 1) aprev = gr%a(il-1)
       gr%ad = gr%a(1) + 2*inp%laywid*1d-4*ab%Nlay
    ELSE IF (inp%chgrainevol .eq. 2) THEN
       gr%ad = gr%a(1)
    END IF
    !gr%Nssur(il) = 4*Pi*(gr%a(il)/2d0)**2/(inp%ds*1E-4)**2
    gr%Vd(il) = 4*Pi/3d0*(gr%a(il)/2d0)**3
    ! compute the porosity parameters and their influence on the surface
    IF (inp%Fporsur .gt. 0) THEN
       Npore2 = inp%Fporsur*gr%Nssur(il-1)/(inp%Slat**2)
       gr%Sdep(il) = gr%Vd(il)*inp%Fvacuum/&
            (Npore2*inp%Slat**2*(inp%ds*1d-4)**3)
       gr%Spore(il) = inp%Slat**2 + 4*inp%Slat*gr%Sdep(il)
       gr%Nstot(il) = gr%Nssur(il)*&
            (inp%Fnpsur + gr%Spore(il)*inp%Fporsur/inp%Slat**2)
       gr%Fnptot(il) = inp%Fnpsur * gr%Nssur(il)/gr%Nstot(il)
       gr%Fportot(il) = 1d0 - gr%Fnptot(il)
       gr%Fedgetot(il) = inp%Fedgesur*gr%Nssur(il)/gr%Nstot(il)
       ! no porosity -> no modification
    ELSE IF (inp%Fporsur .eq. 0) THEN
       Npore2 = 0d0
       gr%Sdep(il) = 0d0
       gr%Spore(il) = 0d0
       gr%Nstot(il) = gr%Nssur(il)
       gr%Fnptot(il) = 1d0
       gr%Fportot(il) = 0d0
       gr%Fedgetot(il) = 0d0
    END IF

  END SUBROUTINE GRAIN_EVOLUTION


  SUBROUTINE PHYS_PROP(il,t)

    !----------------------------------------------------------------
    ! compute the physical properties (nH, Tgas, Tdust, velocity, etc)
    !----------------------------------------------------------------
  
    USE SHARED_VARIABLES, only: inp, iph, sp, ns, gr, ph
    
    INTEGER(KIND=LONG), INTENT(in) :: il
    REAL(KIND=DP), INTENT(in) :: t
    REAL(KIND=DP) :: dens0
    INTEGER(KIND=LONG) :: k, is

    ! constant physical conditions
    IF (inp%chphys .eq. 'no') THEN

       ph%nH = inp%nHini
       ph%Tg = inp%Tgini
       ph%Td = inp%Tdini
       ph%dnH = 0d0
       ph%Av = inp%Avini
       ph%zeta = inp%zetaini
       ph%UVfluxCR = inp%UVfluxCRini

    ! evolving physical conditions
    ELSE IF (inp%chphys .eq. 'yes') THEN
       
       IF (inp%chgrid .eq. 3 .or. (inp%chphys .eq. 'yes' .and. inp%chinphys .eq. 'yes')) ph%rad = radius_time(t*Rys)
       IF (inp%chgrid .eq. 3 .or. (inp%chphys .eq. 'yes' .and. inp%chinphys .eq. 'yes')) ph%v = veloc_time(t*Rys)
       dens0 = inp%nHini
       ph%nH = nH_time(t*Rys)
       ph%dnH = d_density(t*Rys,ph%nH)
       ph%Av = visext(t*Rys,ph%nH)
       ph%Tg = Tg_time(t*Rys,ph%Av)
       ph%Td = ph%Tg !Td_time(t*Rys,ph%Av)
       !ph%Tg = ph%Td
       ph%zeta = inp%zetaini
       ph%UVfluxCR = inp%UVfluxCRini

    END IF
    
    ph%DTOGN = gr%X
    ph%Nsndtot = gr%Nstot(il)*gr%X*ph%nH
    ph%signdtot = Pi*(gr%ad*1d-4/2d0)**2*gr%X*ph%nH
    !write(6,*) ph%Nsndtot,gr%Nstot(il), il, gr%Nssur(il)
    DO is=1,inp%Nspgas
       ph%Vth(is) = velocity(ph%Tg,is)
       ph%stick(is) = stickcoeff(ph%Tg,ph%Td,is)
    END DO
    
  END SUBROUTINE PHYS_PROP


  
  !------------------------!
  !------------------------!
  ! ENVIRONMENT PARAMETERS !
  !------------------------!
  !------------------------!

  ! Radius evolution if spatial grid
  FUNCTION radius_time(t)

    USE SHARED_VARIABLES, only: inp, ph, iph
    
    REAL(KIND=DP) :: radius_time
    REAL(KIND=DP), INTENT(in) :: t
    INTEGER(KIND=LONG) :: itp

    loop : DO itp=1,iph%Nphysteps-1
       IF (t/Rys .ge. iph%time(itp) .and. t/Rys .le. iph%time(itp+1)) THEN 
          radius_time = iph%rad(itp) + &
               (iph%rad(itp+1)-iph%rad(itp))*&
               (t/Rys-iph%time(itp))/&
               (iph%time(itp+1)-iph%time(itp))
          EXIT loop
       END IF
    END DO loop

  END FUNCTION radius_time

  ! Velocity evolution if spatial grid
  FUNCTION veloc_time(t)

    USE SHARED_VARIABLES, only: inp, ph, iph
    
    REAL(KIND=DP) :: veloc_time
    REAL(KIND=DP), INTENT(in) :: t
    INTEGER(KIND=LONG) :: itp

    loop : DO itp=1,iph%Nphysteps-1
       IF (t/Rys .ge. iph%time(itp) .and. t/Rys .le. iph%time(itp+1)) THEN 
          veloc_time = iph%veloc(itp) + &
               (iph%veloc(itp+1)-iph%veloc(itp))*&
               (t/Rys-iph%time(itp))/&
               (iph%time(itp+1)-iph%time(itp))
          EXIT loop
       END IF
    END DO loop

  END FUNCTION veloc_time

  ! Tg evolution from power-law function or input file
  FUNCTION Tg_time(t,Av2) ! in Kelvins

    USE SHARED_VARIABLES, only: inp, ph, iph
    
    REAL(KIND=DP), INTENT(in) :: t, Av2
    REAL:: t_max
    REAL(KIND=DP) :: Tg_time
    INTEGER(KIND=LONG) :: itp

    t_max = inp%tmax*Rys
       
    IF (inp%chinphys .ne. 'yes') THEN

       IF (inp%temp_evol .gt. 0) THEN 
          Tg_time = (inp%Tfin-inp%Tini)*(t/t_max)**inp%temp_evol + inp%Tini
       ELSE IF (inp%temp_evol .eq. 0) THEN
          IF (Av2 .le. 15) THEN
             Tg_time = inp%Tini + &
                  (inp%Tfin-inp%Tini)*(Av2-inp%Avini)/(15-inp%Avini)
          ELSE IF (Av2 .gt. 15) THEN
             Tg_time = inp%Tfin
          END IF
       END IF

    ELSE IF (inp%chinphys .eq. 'yes') THEN       

       loop : DO itp=1,iph%Nphysteps-1
          IF (t/Rys .ge. iph%time(itp) .and. t/Rys .le. iph%time(itp+1)) THEN 
             Tg_time = iph%Tg(itp) + &
                  (iph%Tg(itp+1)-iph%Tg(itp))*&
                  (t/Rys-iph%time(itp))/&
                  (iph%time(itp+1)-iph%time(itp))
             EXIT loop
          END IF
             
       END DO loop

    END IF

  END FUNCTION Tg_time


  ! Td evolution from power-law function or input file
  FUNCTION Td_time(t,Av2) ! in Kelvins

    USE SHARED_VARIABLES, only: inp, ph, iph
    
    REAL(KIND=DP), INTENT(in) :: t, Av2
    REAL:: t_max
    REAL(KIND=DP) :: Td_time
    INTEGER(KIND=LONG) :: itp

    t_max = inp%tmax*Rys
       
    IF (inp%chinphys .ne. 'yes') THEN

       IF (inp%temp_evol .gt. 0) THEN 
          Td_time = (inp%Tfin-inp%Tini)*(t/t_max)**inp%temp_evol + inp%Tini
       ELSE IF (inp%temp_evol .eq. 0) THEN
          IF (Av2 .le. 15) THEN
             Td_time = inp%Tini + &
                  (inp%Tfin-inp%Tini)*(Av2-inp%Avini)/(15-inp%Avini)
          ELSE IF (Av2 .gt. 15) THEN
             Td_time = inp%Tfin
          END IF
       END IF

    ELSE IF (inp%chinphys .eq. 'yes') THEN       

       loop : DO itp=1,iph%Nphysteps-1
          IF (t/Rys .ge. iph%time(itp) .and. t/Rys .le. iph%time(itp+1)) THEN
             Td_time = iph%Td(itp) + &
                  (iph%Td(itp+1)-iph%Td(itp))*&
                  (t/Rys-iph%time(itp))/&
                  (iph%time(itp+1)-iph%time(itp))
             EXIT loop
          END IF
       END DO loop

    END IF

  END FUNCTION Td_time

  ! density evolution from power-law function or input file
  FUNCTION nH_time(t) ! in cm-3

    USE SHARED_VARIABLES, only: inp, ph, iph
    
    REAL(KIND=DP), INTENT(in) :: t
    REAL :: t_max
    REAL(KIND=DP) :: nH_time
    INTEGER(KIND=LONG) :: itp

    t_max = inp%tmax*Rys
    
    IF (inp%chinphys .ne. 'yes') THEN
       nH_time = (inp%nHfin-inp%nHini)*(t/t_max)**inp%dens_evol + &
            inp%nHini

    ELSE IF (inp%chinphys .eq. 'yes') THEN       

       loop : DO itp=1,iph%Nphysteps-1
          IF (t/Rys .ge. iph%time(itp) .and. t/Rys .le. iph%time(itp+1)) THEN 
             nH_time = iph%nH(itp) + &
                  (iph%nH(itp+1)-iph%nH(itp))*&
                  (t/Rys-iph%time(itp))/&
                  (iph%time(itp+1)-iph%time(itp))
             EXIT loop
          END IF
       END DO loop
       
    END IF

  END FUNCTION nH_time


  ! density derivative from power-law function or free-fall collapse
  FUNCTION d_density(t,dens)

    USE SHARED_VARIABLES, only: inp, ph, iph
    
    REAL(KIND=DP), INTENT(in) :: t, dens
    REAL :: t_max
    REAL(KIND=DP) :: rad_ini2
    REAL(KIND=DP) :: d_density
    INTEGER(KIND=LONG) :: itp
    !REAL(KIND=DP) :: dens_evol

    t_max = inp%tmax*Rys
       
    IF (inp%chinphys .ne. 'yes') THEN

       IF (inp%dens_evol .gt. 0) THEN 
          d_density = inp%dens_evol*(inp%nHfin-inp%nHini)*&
               t**(inp%dens_evol-1)/t_max**inp%dens_evol
       ELSE IF (inp%dens_evol .eq. 0) THEN
          d_density = (dens**4/inp%nHini*1d18)**(1/3d0)*&
               (24*Pi*G*Mpr*inp%nHini*1d6*((dens/inp%nHini)**(1/3d0)-1))**0.5&
               * 1d-6
       END IF

    ELSE IF (inp%chinphys .eq. 'yes') THEN     
  
       loop : DO itp=1,iph%Nphysteps-1
          IF (t/Rys .ge. iph%time(iph%Nphysteps-itp)) THEN 
             d_density = (iph%nH(iph%Nphysteps-itp+1)-iph%nH(iph%Nphysteps-itp))/&
                  ((iph%time(iph%Nphysteps-itp+1)-iph%time(iph%Nphysteps-itp))*Rys)
             EXIT loop
          END IF
       END DO loop

    END IF

  END FUNCTION d_density

  ! visual extinction evolution with density
  REAL(KIND=DP) FUNCTION visext(t,nH)

    USE SHARED_VARIABLES, only: inp, ph, iph
    
    REAL(KIND=DP), INTENT(in) :: t
    REAL(KIND=DP), INTENT(in) :: nH
    INTEGER(KIND=LONG) :: itp

    IF (inp%chinphys .ne. 'yes') THEN

       visext = inp%Avini*(nH/inp%nHini)**(2/3d0)

    ELSE IF (inp%chinphys .eq. 'yes') THEN     
  
       loop : DO itp=1,iph%Nphysteps-1
          IF (t/Rys .ge. iph%time(iph%Nphysteps-itp)) THEN 
             visext = iph%Av(iph%Nphysteps-itp) + &
                  (iph%Av(iph%Nphysteps-itp+1)-iph%Av(iph%Nphysteps-itp))*&
                  (t/Rys-iph%time(iph%Nphysteps-itp))/&
                  (iph%time(iph%Nphysteps-itp+1)-iph%time(iph%Nphysteps-itp))
             EXIT loop
          END IF
       END DO loop

    END IF

  END FUNCTION visext


  ! thermal velocity
  FUNCTION velocity(temp,num_spec) ! in cm.s-1

    USE SHARED_VARIABLES, only: inp, ph, iph, sp
    
    REAL(KIND=DP), INTENT(in) :: temp
    INTEGER(KIND=LONG), INTENT(in) :: num_spec
    REAL(KIND=DP) :: velocity

    IF (trim(sp%name(num_spec)) .ne. 'E') THEN 
       velocity = sqrt(8*kb*temp/(Pi*sp%mass(num_spec)*Mpr))*1d2
    ELSE IF (trim(sp%name(num_spec)) .eq. 'E') THEN
       velocity = sqrt(8*kb*temp/(Pi*Me))*1d2
    END IF

  END FUNCTION velocity

  
  ! ISRF as function of the visual extinction
  REAL(KIND=DP) FUNCTION FUV(G0,Av2)
        
    REAL(KIND=DP), INTENT(in) :: G0, Av2

    FUV = G0*2d8*exp(-1.8*Av2)

  END FUNCTION FUV

  ! sticking coefficient from Tielens (2005)
  REAL(KIND=DP) FUNCTION stickcoeff(Tg,Td,num_spec)
    
    USE SHARED_VARIABLES, only: ns, inp
    
    REAL(KIND=DP), INTENT(in) :: Tg, Td 
    INTEGER(KIND=LONG), INTENT(in) :: num_spec

    IF (num_spec .eq. ns%H .or. num_spec .eq. ns%D) THEN
       stickcoeff = inp%stick_coef
       stickcoeff = 1/(1+4d-2*(Tg+Td)**(5d-1) &
            + 2d-3*Tg + 8d-6*Tg**2)
    ELSE
       stickcoeff = inp%stick_coef
    END IF

  END FUNCTION stickcoeff


END MODULE PHYSPROP
