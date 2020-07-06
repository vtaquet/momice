MODULE NUMPROP

  USE VARIABLES
  USE ALLOCATION
  USE CHEMPROP
  USE PHYSPROP
  USE ODES

  IMPLICIT NONE

  !----------------------------------------------------------------------------
  ! This module contains several subroutines that: 
  ! read the various input files
  !----------------------------------------------------------------------------

CONTAINS 

  
!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------
!! TIME PARAMETERS
!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------
  
  ! compute the initial number of timesteps
  SUBROUTINE COMPUTE_NSTEPS

    USE SHARED_VARIABLES, only: inp, ph, sp, ns, iph, out

    INTEGER(KIND=LONG) :: is
    REAL :: Nlay_max, tmin, tmax
    
    out%Ninisteps = 0d0; out%Nplotsteps = 0d0
    out%Ninilayer = 0d0; Nlay_max = 0d0

    tmin = inp%tmax/1d4
    tmax = inp%tmax
    IF (trim(inp%chphys) .ne. 'yes' .or. trim(inp%chinphys) .ne. 'yes') THEN
       out%Nplotsteps = (alog10(tmax)-alog10(tmin))*30+1
    ELSE IF (trim(inp%chphys) .eq. 'yes' .and. trim(inp%chinphys) .eq. 'yes') THEN
       out%Nplotsteps = iph%Nphysteps
       !write(6,*) out%Nplotsteps, iph%Nphysteps
    END IF

    Nlay_max = 0d0
    DO is=1,inp%Nspgas
       IF (sp%num(is) .ne. ns%H2 .and. sp%num(is) .ne. ns%oH2 .and. &
            sp%num(is) .ne. ns%He) THEN
          Nlay_max = Nlay_max + sp%Xini(is)*ph%nH/(ph%Nsndtot*inp%Fnpsur)
       END IF
    END DO
    out%Ninilayer = IFIX(Nlay_max+1+1000)
    IF (out%Ninilayer == 0) out%Ninilayer = 1
    
    out%Ninisteps = out%Nplotsteps + 2*out%Ninilayer
    out%Ninisteps = out%Ninisteps * 2.
    
  END SUBROUTINE COMPUTE_NSTEPS

  ! compute the output timesteps
  SUBROUTINE COMPUTE_TIME

    USE SHARED_VARIABLES, only: inp, ph, sp, ns, iph, out

    INTEGER(KIND=LONG) :: it
    REAL(KIND=DP) :: time_min
    
    ! compute the time at which the abundances are plotted
    IF (trim(inp%chphys) .ne. 'yes' .or. trim(inp%chinphys) .ne. 'yes') THEN
       time_min = inp%tmax/1d7
       DO it=1,out%Nplotsteps
          out%timefix(it) = 10**(log10(time_min) + &
               (it-1)*(log10(inp%tmax)-log10(time_min))/(out%Nplotsteps-1))
       END DO
       
    ELSE IF (trim(inp%chphys) .eq. 'yes' .and. trim(inp%chinphys) .eq. 'yes') THEN
       out%Nplotsteps = iph%Nphysteps
       DO it=1,out%Nplotsteps
          out%timefix(it) = iph%time(it)
       END DO
    END IF


  END SUBROUTINE COMPUTE_TIME


  SUBROUTINE COMPUTE_DTIME(it,iplot,il,timeout,timein)

    USE SHARED_VARIABLES, only: inp, out, ph, ode, ns, sp, F2tot
    USE SHARED_VARIABLES2

    IMPLICIT NONE

    REAL(KIND=DP), INTENT(out) :: timeout
    REAL(KIND=DP), INTENT(in) :: timein
    INTEGER(KIND=LONG), INTENT(in) :: iplot, il, it
    REAL(KIND=DP) :: lay_full, divdNgr, t_min, popsurf
    REAL(KIND=DP), DIMENSION(inp%Neq) :: F
    INTEGER(KIND=LONG) :: ic, is
    REAL(KIND=DP) :: dNgr, dNacc, dNev, dNgas, Ftot
    REAL(KIND=DP) :: pop_int_prev, dtime1, dtime2, dNgr2, Ptot

    ! accretion rate for grains -> needed for timestep
    dNgr = 0d0; dNacc = 0d0; dNev = 0d0; dNgr2 = 0d0; dNgas = 0d0; Ptot = 0d0
    dtime1 = 0d0; dtime2 = 0d0

    ! 3-phase: timestep depends on population evolution
    IF (inp%chlayer .eq. 2) THEN
       CALL Fode0(inp%Neq,timein,ode%Y,F)
       DO is=inp%Nspgas+1,inp%Nspgas+inp%Nspgr
          IF (is .ne. ns%JH2 .and. is .ne. ns%JoH2 .and. &
               is .ne. ns%He) Ftot = Ftot + F(is)
          !WRITE(6,*) sp%name(is), Ftot,F(is), ODE%Y(IS)
       END DO
       dNgr = Ftot*ph%nH/ph%Nsndtot
       dtime1 = abs(1/dNgr)
       
       ! net accretion
       IF (Ftot.gt.0) THEN

          ! low timestep to follow gradual formation of ices
          ph%dtime = dtime1 !max(dtime1,dtime2)

          IF (ph%dtime/Rys .lt. 1d-2) THEN
             ph%dtime = 1d-2*Rys
          END IF
       
          ! modify the timestep to get the timestep for output vectors
          IF (ph%dtime .gt. out%timefix(iplot)*Rys-timein) THEN
             ph%dtime = out%timefix(iplot)*Rys-timein
          END IF

       ! net evaporation
       ELSE IF (Ftot.le.0) THEN
          
          ph%dtime = out%timefix(iplot)*Rys-timein
          
       END IF
       
       timeout = timein + ph%dtime
       
    END IF

    ! 2-phase: timestep depends on initially fixed timestep
    IF (inp%chlayer .eq. 1) THEN
       timeout = out%timefix(iplot)*Rys
       ph%dtime = timeout-timein
    END IF
    timeout = out%timefix(iplot)*Rys
    ph%dtime = timeout-timein

    !timeout = out%timefix(iplot)*Rys
    !ph%dtime = timeout-timein

    ! evolving physical conditions: limit the timestep from nH evolution
    !IF (inp%chphys .eq. 'yes') THEN
    !   IF (ph%dtime*ph%dnH .gt. ph%nH*1d-2) THEN
    !   write(6,*) 'ooh'
    !      dtime2 = ph%nH*1d-2/ph%dnH
    !      ph%dtime = max(ph%dtime,dtime2)
    !      timeout = timein + ph%dtime
    !   END IF
    !END IF


  END SUBROUTINE COMPUTE_DTIME


!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------
!! PHYSICAL/CHEMICAL OUTPUT PARAMETERS
!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------

  SUBROUTINE INI_PHYSCHEMPARAMS(iplot,ilout)

    !----------------------------------------------------------------------
    ! initialize the vectors for each model
    !----------------------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, out, ph, sp, ns, ab, ode, gr 
    USE SHARED_VARIABLES2

    IMPLICIT NONE
    
    INTEGER(KIND=LONG), INTENT(inout) :: iplot, ilout
    REAL(KIND=DP), DIMENSION(inp%Nspgr) :: Pbulk
    INTEGER(KIND=LONG) :: is, ie, il3, ilint, i
    REAL(KIND=DP) :: il2
    REAL(KIND=DP) :: Npart, Npart2, Nsurf

   
    ! timestep
    ph%dtime = 0d0; ph%time = 0d0
    
    ! initialize output vectors
    out%time = 0d0; ph%time_prev = 0d0; ph%rad = 0d0; ph%rad_prev = 0d0
    out%nH = 0d0; out%Tg = 0d0; out%Td = 0d0; out%Av = 0d0; ab%Nlay = 0d0
    out%Xgas = 0d0; out%Xice = 0d0; out%Xsur = 0d0  
    ph%nH_prev = 0d0; ph%Av_prev = 0d0; ph%Tg_prev = 0d0; ph%Td_prev = 0d0
    ab%Xgas = 0d0; ab%Xgas_prev = 0d0; ab%Xice = 0d0; ab%Xice_prev = 0d0
    ab%Pice = 0d0; ab%Pice_prev = 0d0; ab%Xsurml = 0d0;  ab%Xsurml_prev = 0d0
    Pbulk = 0d0; ab%Psur = 0d0;  ab%Psur_prev = 0d0
    ode%Y = 0d0; Npart = 0d0; Npart2 = 0d0; Nsurf = 0d0; ilout = 1; ilint = 1

    ! if Xini(H) ne 0 then Xini(H) = 1/nH
    IF (sp%Xini(ns%H) .ne. 0) sp%Xini(ns%H) = (2.3*ph%zeta)/(2*ph%Vth(ns%H)*Pi*(inp%acst*1d-4/2d0)**2*gr%X*ph%nH)
    !write(6,*) sp%Xini(ns%H),sp%Xini(ns%H)/sp%Xini(ns%O)!,(2.3*ph%zeta),ph%Vth(ns%H),Pi*(inp%acst*1d-4/2d0)**2,gr%X,ph%nH
        
    ! 2-phase method
    IF (inp%chlayer .eq. 1) THEN

       DO is=1,inp%Nspecies
          ode%Y(is) = sp%Xini(is)
          ode%Yini(is) = ode%Y(is) 
       END DO

       ab%Xgas(1:inp%Nspgas) = sp%Xini(1:inp%Nspgas)
       ab%Xice(1:inp%Nspgr) = sp%Xini(inp%Nspgas+1:inp%Nspgas+inp%Nspgr)

       ! 3-phase method
    ELSE IF (inp%chlayer .eq. 2 .and. inp%Nspgr .gt. 0) THEN

       ! gas phase species
       DO is=1,inp%Nspgas
          ode%Y(is) = sp%Xini(is)
          ode%Yini(is) = ode%Y(is) 
       END DO
       ab%Xgas(1:inp%Nspgas) = sp%Xini(1:inp%Nspgas)

       ! count the initial number of solid particules
       DO is=inp%Nspgas+1,inp%Nspecies
          Npart = Npart + sp%Xini(is)/gr%X
       END DO

       ! if ice is initially present
       IF (Npart .gt. 0d0) THEN
          
          ! count the exact number of corresponding layers (due to grain growth)
          il2 = 1; Npart2 = 0d0
          DO WHILE (Npart2 .lt. Npart) 
             ilint = il2
             Npart2 = Npart2 + gr%Nstot(ilint)
             IF (Npart2 .lt. Npart) THEN
                il2 = il2+1
             ELSE IF (Npart2 .ge. Npart) THEN
                Npart2 = Npart2 - gr%Nstot(ilint)
                il2 = il2 + (Npart-Npart2)/gr%Nstot(ilint)-1
                Npart2 = Npart
             END IF
          END DO

          ! attribute the composition to each layer in the surface and bulk
          ilint = floor(il2)+1
          ilout = floor(il2)
          IF (ilint .gt. 1) THEN
             DO is=inp%Nspgas+1,inp%Nspecies
                ode%Y(is) = sp%Xini(is)*gr%Nstot(ilint)/Npart
                ode%Yini(is) = ode%Y(is) 
                Nsurf = Nsurf + ode%Y(is)/gr%X
                ode%Y(inp%Nspgr+is) = sp%Xini(is)-ode%Y(is)
                ab%Xsur(is) = ode%Y(is)
                ab%Xbulk(is) = ode%Y(is+inp%Nspgr)
             END DO
          ELSE IF (ilint .eq. 1) THEN
             DO is=inp%Nspgas+1,inp%Nspecies
                ode%Y(is) = sp%Xini(is)*gr%Nstot(ilint)/Npart
                ode%Yini(is) = ode%Y(is) 
                ode%Y(inp%Nspgr+is) = 0.
                ab%Xsur(is) = ode%Y(is)
                ab%Xbulk(is) = ode%Y(is+inp%Nspgr)
             END DO
          END IF

       END IF

    END IF
    
    ! initial abundances
    ab%Xtot = 0d0; ab%Xelem_ini = 0d0; ab%Xelem = 0d0
    ab%Xcharge = 0d0; ab%Xelem_lay = 0d0
     
    ! total abundance
    !gaseous species > icy species and no porosity
    IF (inp%Nspgas .gt. inp%Nspgr .and. inp%Nspgr .ne. 2*inp%Nspgr2) THEN
       ! species in the gas and in ices (basically, neutrals)
       DO is=1,inp%Nspgr
          ab%Xtot(is) = ab%Xgas(is) + ab%Xice(is)
       END DO
       DO is=inp%Nspgr+1,inp%Nspgas
          ab%Xtot(is) = ab%Xgas(is)
       END DO
       !gaseous species > icy species + porosity
    ELSE IF (inp%Nspgas .lt. inp%Nspgr .and. inp%Nspgr .eq. 2*inp%Nspgr2) THEN
        ! species in the gas and in ices (basically, neutrals)
        DO is=1,inp%Nspgr2
          ab%Xtot(is) = ab%Xgas(is) + ab%Xice(is) + ab%Xice(is+inp%Nspgr2)
        END DO
        ! species only in the gas (basically ions)
        DO is=inp%Nspgr2+1,inp%Nspgas
           ab%Xtot(is) = ab%Xgas(is)
        END DO

     ! gaseous species = icy species (-> no porosity)
     ELSE IF (inp%Nspgas .eq. inp%Nspgr) THEN
        DO is=1,inp%Nspgr
           ab%Xtot(is) = ab%Xgas(is) + ab%Xice(is)
        END DO
     END IF

     ! elemental abundances and charges
     DO is=1,inp%Nspgas
        ab%Xcharge(1) = ab%Xcharge(1) + ab%Xtot(is)*sp%charge(is)
        DO ie=1,Nelements
           ab%Xelem(ie) = ab%Xelem(ie) + ab%Xtot(is)*sp%element(ie,is)
           ab%Xelem_lay(ie) = ab%Xelem(ie)
           ab%Xelem_ini(ie) = ab%Xelem(ie)
        END DO
     END DO

     IF (out%timefix(1) .eq. 0d0) THEN 
        out%nH(1) = ph%nH; out%Tg(1) = ph%Tg; out%Td(1) = ph%Td; out%Av(1) = ph%Av
        out%Xgas(1:inp%Nspgas,1) = ab%Xgas(1:inp%Nspgas)
        out%Xice(1:inp%Nspgr,1) = ab%Xice(1:inp%Nspgr)
        iplot = 2
     END IF

   END SUBROUTINE INI_PHYSCHEMPARAMS


  SUBROUTINE SAVE_PREVPARAMS(iprev,il)

    !----------------------------------------------------
    ! save the parameters at the previous time
    !----------------------------------------------------

    USE SHARED_VARIABLES, only: inp, ph, ab, ode, out
    
    INTEGER(KIND=LONG), INTENT(in) :: iprev
    INTEGER(KIND=LONG), INTENT(inout) :: il
    
    IF (iprev .eq. 1) THEN

       ph%rad_prev = ph%rad
       ph%v_prev = ph%v
       ph%time_prev = ph%time
       ph%Tg_prev = ph%Tg
       ph%Td_prev = ph%Td
       ph%nH_prev = ph%nH
       ph%dnH_prev = ph%dnH
       ph%Av_prev = ph%Av
       ph%zeta_prev = ph%zeta
       ph%UVfluxCR_prev = ph%UVfluxCR

       ab%Xgas_prev(1:inp%Nspgas) = ab%Xgas(1:inp%Nspgas)
       ab%Xice_prev(1:inp%Nspgr) = ab%Xice(1:inp%Nspgr)
       ab%Xsur_prev(1:inp%Nspgr) = ab%Xsur(1:inp%Nspgr)
       ab%Xbulk_prev(1:inp%Nspgr) = ab%Xbulk(1:inp%Nspgr)
       ab%Xtot_prev(1:inp%Nspgas) = ab%Xtot(1:inp%Nspgas)
       ab%Xelem_prev(1:Nelements) = ab%Xelem(1:Nelements)
       ab%Xcharge_prev = ab%Xcharge
       ab%Psur_prev(1:inp%Nspgr) = ab%Psur(1:inp%Nspgr)
       ab%Xsurml_prev(1:inp%Nspgr) = ab%Xsurml(1:inp%Nspgr)
       ab%Pice_prev(1:inp%Nspgr) = ab%Pice(1:inp%Nspgr)
       ab%Nlay_prev = ab%Nlay

       il = int(ab%Nlay)+1
       
    ELSE IF (iprev .eq. 2) THEN

       ph%rad = ph%rad_prev
       ph%v = ph%v_prev
       ph%time = ph%time_prev
       ph%Tg = ph%Tg_prev
       ph%Td = ph%Td_prev
       ph%nH = ph%nH_prev
       ph%dnH = ph%dnH_prev
       ph%Av = ph%Av_prev
       ph%zeta = ph%zeta_prev
       ph%UVfluxCR = ph%UVfluxCR_prev

       ab%Xgas(1:inp%Nspgas) = ab%Xgas_prev(1:inp%Nspgas)
       ab%Xice(1:inp%Nspgr) = ab%Xice_prev(1:inp%Nspgr)
       ab%Xsur(1:inp%Nspgr) = ab%Xsur_prev(1:inp%Nspgr)
       ab%Xbulk(1:inp%Nspgr) = ab%Xbulk_prev(1:inp%Nspgr)
       ab%Xtot(1:inp%Nspgr) = ab%Xtot_prev(1:inp%Nspgr)
       ab%Xelem(1:Nelements) = ab%Xelem_prev(1:Nelements)
       ab%Xcharge = ab%Xcharge_prev
       ab%Psur(1:inp%Nspgr) = ab%Psur_prev(1:inp%Nspgr) 
       ab%Xsurml(1:inp%Nspgr) = ab%Xsurml_prev(1:inp%Nspgr)
       ab%Pice(1:inp%Nspgr) = ab%Pice_prev(1:inp%Nspgr)
       ab%Nlay = ab%Nlay_prev
    END IF

  END SUBROUTINE SAVE_PREVPARAMS


  SUBROUTINE CONVERT_CHEMPARAMS(il,layloop)

    !---------------------------------------------
    ! convert and save chemical parameters
    !---------------------------------------------

    USE SHARED_VARIABLES, only: inp, ph, sp, ab, ode, gr, ns
    USE SHARED_VARIABLES2
    
    INTEGER(KIND=LONG), INTENT(in) :: il
    REAL(KIND=DP), INTENT(out) :: layloop
    INTEGER(KIND=LONG) :: is, ie, i, il2
    REAL(KIND=DP), DIMENSION(inp%Nspgr) :: dPsp, Pbulk
    REAL(KIND=DP), DIMENSION(il) :: MLfrac
    REAL(KIND=DP) :: Prem, MLact, MLact2, Ptotfracml, Pbulktot, Pbulkacttot, Pbulkacttot2, ML


    ! initialize vectors
    ab%Xgas = 0d0; ab%Xtot = 0d0; ab%Xelem = 0d0; ab%Xcharge = 0d0
    ab%Ptotml(il) = 0d0; ab%pop_int = 0d0; ab%Xice = 0d0;  ab%Ptot = 0d0; ab%Psurtot = 0d0
    ab%Pice = 0d0; ab%Xelem = 0d0; ab%Psur = 0d0; Pbulktot = 0d0
    dPsp = 0d0; Pbulk = 0d0; ab%Xbulk = 0d0; ab%Nlay = 0d0
    
    ! limit low values to error tolerances
    DO is=1,inp%Nspecies
       IF (ode%Y(is) .le. 0d0) ode%Y(is) = 0d0
    END DO
    
    ! convert gas phase abundances
    DO is=1,inp%Nspgas
       ab%Xgas(is) = ode%Y(is)
    END DO

    DO is=1,inp%Nspgr

       ! ice abundances
       ab%Xice(is) = ode%Y(inp%Nspgas+is)
       ab%Xsur(is) = ode%Y(inp%Nspgas+is)
       IF (inp%chmulti .eq. 2 .and. inp%chlayer .eq. 2) THEN
          ab%Xbulk(is) = ode%Y(inp%Nspgas+inp%Nspgr+is)
          ab%Xice(is) = ab%Xice(is) + ode%Y(inp%Nspgas+inp%Nspgr+is)
       END IF

       ! number of ice particles
       ab%Psur(is) = ab%Xsur(is)/gr%X
       ab%Pbulk(is) = ab%Xbulk(is)/gr%X
       ab%Xsurml(is) = ab%Xsur(is)/gr%X
       ab%Xbulkml(is) = ab%Xbulk(is)/gr%X
       ab%Pice(is) = ab%Xice(is)/gr%X
       ab%Ptot = ab%Ptot + ab%Pice(is)
       ab%Psurtot = ab%Psurtot + ab%Psur(is)
       ab%Pbulktot = ab%Pbulktot + ab%Pbulk(is)

       !write(6,*) sp%name(is), ab%Xsurml(is), ab%Xbulkml(is)
       
    END DO

    IF (ab%Psurtot.ne.0) THEN
       ab%Xsurml(1:inp%Nspgr) = ab%Xsurml(1:inp%Nspgr)/ab%Psurtot
    ELSE
       ab%Xsurml(1:inp%Nspgr) = 0.
    END IF
    IF (ab%Pbulktot.ne.0) THEN
       ab%Xbulkml(1:inp%Nspgr) = ab%Xbulkml(1:inp%Nspgr)/ab%Pbulktot
    ELSE
       ab%Xbulkml(1:inp%Nspgr) = 0.
    END IF

    !write(6,*) ab%Xsurml(ns%H2O), ab%Xbulkml(ns%H2O)
    
    ! compute the number of layers
    Prem = ab%Ptot; i = 1; ML = 0
    DO WHILE (Prem .gt. 0)
       Prem = Prem - gr%Nssur(i)
       IF (Prem .gt. 0) THEN
          ML = ML + 1
       ELSE IF (Prem .le. 0) THEN
          ML = ML + 1+Prem/gr%Nssur(i)
       END IF
       i=i+1
    END DO

    !write(6,*) gr%Nssur(i), gr%X, ML, ab%Ptot
    
    ab%Nlay = ML

    ! compute the total abundance (grain+gas) of each species
    IF (inp%Nspgas .gt. inp%Nspgr) THEN
       DO is=1,inp%Nspgr
          ab%Xtot(is) = ab%Xgas(is) + ab%Xice(is)
       END DO
       DO is=inp%Nspgr+1,inp%Nspgas
          ab%Xtot(is) = ab%Xgas(is)
       END DO
    ELSE IF (inp%Nspgas .lt. inp%Nspgr .and. &
         inp%Nspgr .eq. 2*inp%Nspgr2) THEN
       DO is=1,inp%Nspgr2
          ab%Xtot(is) = ab%Xgas(is) + ab%Xice(is) + &
               ab%Xice(is+inp%Nspgr2)
       END DO
       DO is=inp%Nspgr2+1,inp%Nspgas
          ab%Xtot(is) = ab%Xgas(is)
       END DO
    ELSE IF (inp%Nspgas .eq. inp%Nspgr) THEN
       DO is=1,inp%Nspgr
          ab%Xtot(is) = ab%Xgas(is) + ab%Xice(is)
       END DO
    END IF
    
    ! compute the abundance of elements/charge to check the conservation
    DO is=1,inp%Nspgas
       ab%Xcharge(1) = ab%Xcharge(1) + ab%Xtot(is)*sp%charge(is)
       DO ie=1,Nelements
          ab%Xelem(ie) = ab%Xelem(ie) + ab%Xtot(is)*sp%element(ie,is)
       END DO
    END DO
    
    layloop = ab%Ptotml(il)


  END SUBROUTINE CONVERT_CHEMPARAMS


  SUBROUTINE SAVE_OUTPARAMS(isave,il,iplot,iplot2)

    !-----------------------------------------------------------------
    ! save output parameters at specific timesteps
    !-----------------------------------------------------------------

    USE SHARED_VARIABLES, only: inp, ph, ab, out, sp, ns, re, nr, ode, gr 
    
    INTEGER(KIND=LONG), INTENT(in) :: isave, il
    INTEGER(KIND=LONG), INTENT(inout) :: iplot, iplot2
    INTEGER(KIND=LONG) :: ic, is

    IF (isave .eq. 1) THEN 
       !IF (inp%chgrid .eq. 3 .or. (inp%chphys .eq. 'yes' .and. inp%chinphys .eq. 'yes')) 
       out%rad(iplot2) = ph%rad
       !IF (inp%chgrid .eq. 3 .or. (inp%chphys .eq. 'yes' .and. inp%chinphys .eq. 'yes')) 
       out%v(iplot2) = ph%v
       out%time(iplot2) = ph%time
       out%nH(iplot2) = ph%nH
       out%Av(iplot2) = ph%Av
       out%Tg(iplot2) = ph%Tg
       out%Td(iplot2) = ph%Td
       out%Nlay(iplot2) = ab%Nlay
       out%Xgas(1:inp%Nspgas,iplot2) = ab%Xgas(1:inp%Nspgas)
       out%Xice(1:inp%Nspgr,iplot2) = ab%Xice(1:inp%Nspgr)
       out%Xsur(1:inp%Nspgr,iplot2) = ab%Xsur(1:inp%Nspgr)
       out%Xsurml(1:inp%Nspgr,iplot2) = ab%Xsurml(1:inp%Nspgr)
    ELSE IF (isave .eq. 2) THEN
       !IF (inp%chgrid .eq. 3 .or. (inp%chphys .eq. 'yes' .and. inp%chinphys .eq. 'yes')) 
       out%rad(iplot2) = ph%rad_prev
       !IF (inp%chgrid .eq. 3 .or. (inp%chphys .eq. 'yes' .and. inp%chinphys .eq. 'yes')) 
       out%v(iplot2) = ph%v_prev
       out%time(iplot2) = ph%time_prev
       out%nH(iplot2) = ph%nH_prev
       out%Av(iplot2) = ph%Av_prev
       out%Tg(iplot2) = ph%Tg_prev
       out%Td(iplot2) = ph%Td_prev
       out%Nlay(iplot2) = ab%Nlay_prev
       out%Xgas(1:inp%Nspgas,iplot2) = ab%Xgas_prev(1:inp%Nspgas)
       out%Xice(1:inp%Nspgr,iplot2) = ab%Xice_prev(1:inp%Nspgr)
       out%Xsur(1:inp%Nspgr,iplot2) = ab%Xsur_prev(1:inp%Nspgr)
       out%Xsurml(1:inp%Nspgr,iplot2) = ab%Xsurml_prev(1:inp%Nspgr)
    END IF

    ! order the importance of reactions for each species
    IF (inp%chreacrates .eq. 1) THEN
       
       CALL COMPUTE_REACRATES(inp,ph,sp,ns,re,nr,ode%Y)
       
       DO is=1,inp%Nsprate
          DO ic=1,sp%maxNreac
             out%numform(is,ic,iplot2) = sp%numform(is,ic)
             out%numdest(is,ic,iplot2) = sp%numdest(is,ic)
             out%rateform(is,ic,iplot2) = sp%rateform(is,ic)
             out%ratedest(is,ic,iplot2) = sp%ratedest(is,ic)
             !write(6,*) is, 'Form:',sp%numform(is,ic), sp%rateform(is,ic)
             !write(6,*) is, 'Dest:',sp%numdest(is,ic), sp%ratedest(is,ic)
          END DO
       END DO
    END IF

    IF (isave .eq. 1 .and. abs(ph%time-out%timefix(iplot)) .lt. 1/Rys) &
         iplot = iplot + 1
    iplot2 = iplot2 + 1
    
  END SUBROUTINE SAVE_OUTPARAMS


!!$  SUBROUTINE SAVE_POPLAY(ifin,checklay,il,it)
!!$
!!$    !--------------------------------------------------------------------------
!!$    ! save fractional composition of the layer and adjust the ODE variables
!!$    !--------------------------------------------------------------------------
!!$
!!$    USE SHARED_VARIABLES, only: inp, ph, sp, ab, gr, ode 
!!$    USE SHARED_VARIABLES2
!!$
!!$    INTEGER(KIND=LONG), INTENT(in) :: il, it, ifin, checklay
!!$    INTEGER(KIND=LONG) :: is, i
!!$    REAL(KIND=DP) :: poplay_last
!!$
!!$    ! save the number of species in the full ML 
!!$    !DO is=1,inp%Nspecies
!!$    !   sp%Eb_lay(is,il) = sp%Eb(is)
!!$    !END DO
!!$
!!$    ! Multilayer method from Taquet et al. (2012): 
!!$    ! re-initialize Y and save the surface compo from surface population
!!$    IF (inp%chmulti .eq. 1) THEN 
!!$
!!$       IF (ifin .eq. 1) THEN
!!$          DO is=1,inp%Nspgas
!!$             ode%Y(is) = ab%Xgas(is)
!!$          END DO
!!$          DO is=inp%Nspgas+1,inp%Nspecies
!!$             ode%Y(is) = 0d0
!!$          END DO
!!$          ode%Ysave(1:inp%Neq) = ode%Y(1:inp%Neq) 
!!$          ab%pop_int = 0d0
!!$       END IF
!!$
!!$    ! Multilayer method from Hasegawa & Herbst (1993):
!!$    ! save the surface compo from external bulk population
!!$    ELSE IF (inp%chlayer .eq. 2 .and. inp%chmulti .eq. 2 .and. ifin .eq. 1) THEN
!!$
!!$       DO is=1,inp%Nspgr
!!$          poplay_last = poplay_last + ab%Pfracml(il,is)
!!$       END DO
!!$
!!$       ! A REVOIR !
!!$    ELSE IF (inp%chlayer .eq. 2 .and. inp%chmulti .eq. 2 .and. ifin .eq. 2) THEN
!!$
!!$       CALL GRAIN_EVOLUTION(il+1,1d0)
!!$
!!$       ! at the end of the calculation, compute the compo of the last two MLs
!!$       DO is=1,inp%Nspgr
!!$          poplay_last = poplay_last + ab%Pfracml(il,is)
!!$       END DO
!!$
!!$       DO is=1,inp%Nspgr
!!$          ab%Pfracml(il,is) = ab%Pfracml(il,is) + &
!!$               (1-poplay_last)*ode%Y(inp%Nspgas+is)/(gr%X)
!!$          ab%Pfracml(il+1,is) = poplay_last*ode%Y(inp%Nspgas+is)/(gr%X)
!!$       END DO
!!$
!!$    END IF
!!$
!!$    ! save the time and the timestep 
!!$    ab%Pfracml(il,inp%Nspgr+1) = ph%time_prev
!!$    ab%Pfracml(il,inp%Nspgr+2) = it-1
!!$    ab%Xelem_lay(1:Nelements) = ab%Xelem_prev(1:Nelements)
!!$
!!$
!!$  END SUBROUTINE SAVE_POPLAY


!!$  SUBROUTINE RESTART_ML(il,it)
!!$
!!$    !-----------------------------------
!!$    ! return to the beginning of the ML
!!$    !-----------------------------------
!!$
!!$    USE SHARED_VARIABLES, only: inp, ph, ode, ab
!!$    
!!$    INTEGER(KIND=LONG), INTENT(inout) :: il,it
!!$    
!!$    IF (il == 1) THEN
!!$       it = 1
!!$       ph%time = 0d0
!!$    ELSE
!!$       it = int(ab%Pfracml(il-1,inp%Nspgr+2))+1
!!$       IF (il .eq. 1) THEN
!!$          ph%time = 0d0
!!$       ELSE IF (il .gt. 1) THEN
!!$          ph%time = ab%Pfracml(il-1,inp%Nspgr+1)
!!$       END IF
!!$    END IF
!!$
!!$    ode%Y(1:inp%Neq) = ode%Ysave(1:inp%Neq)
!!$    ab%Xelem(1:Nelements) = ab%Xelem_lay(1:Nelements)
!!$    ab%Ptotml(il) = 0d0
!!$
!!$
!!$  END SUBROUTINE RESTART_ML


END MODULE NUMPROP
    
