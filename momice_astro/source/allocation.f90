MODULE ALLOCATION

  USE VARIABLES

  IMPLICIT NONE


  !----------------------------------------------------------------------------
  ! This module contains several subroutines which
  ! allocate/deallocate memory for tables 
  !----------------------------------------------------------------------------


CONTAINS 

  SUBROUTINE ALLO_SPECIES

    USE SHARED_VARIABLES, only: inp, sp
    
    ALLOCATE(sp%num(inp%Neq+1), sp%name(inp%Neq+1), stat=ok)
    ALLOCATE(sp%charge(inp%Neq+1), sp%mass(inp%Neq+1), stat=ok)
    ALLOCATE(sp%element(Nelements,inp%Neq+1))
    ALLOCATE(sp%name_elem(Nelements), stat=ok)
    ALLOCATE(sp%Eb_carb(inp%Neq+1), sp%Eb_wat(inp%Neq+1), stat=ok)
    ALLOCATE(sp%Eb_pure(inp%Neq+1), sp%Ed_pure(inp%Neq+1), stat=ok)
    ALLOCATE(sp%Eb_HH93(inp%Neq+1), sp%Eb_H2(inp%Neq+1), stat=ok)
    ALLOCATE(sp%Eb_ini(inp%Neq+1), sp%Eb(inp%Neq+1), stat=ok)
    ALLOCATE(sp%Ed_carb(inp%Neq+1), sp%Ed_wat(inp%Neq+1), stat=ok)
    ALLOCATE(sp%Ed(inp%Neq+1), sp%Ed_H2(inp%Neq+1), sp%dHf(inp%Neq+1), stat=ok)
    ALLOCATE(sp%grproc(inp%Neq,5), stat=ok) 
    ALLOCATE(sp%formreac(inp%Neq,inp%Nreactions), stat=ok)
    ALLOCATE(sp%destreac(inp%Neq,inp%Nreactions), stat=ok)
    ALLOCATE(sp%maxreac(2,inp%Neq))
    ALLOCATE(sp%cr_coef(inp%Neq+1), sp%Xini(inp%Neq+1), stat=ok)

    IF (inp%chreacrates .eq. 1) THEN
       ALLOCATE(sp%rate_name(inp%Nsprate), stat=ok)
       ALLOCATE(sp%rate_num(inp%Nsprate), stat=ok)
       ALLOCATE(sp%rateform(inp%Nsprate+1,inp%Nreactions), stat=ok)
       ALLOCATE(sp%ratedest(inp%Nsprate+1,inp%Nreactions), stat=ok)
       ALLOCATE(sp%numform(inp%Nsprate+1,inp%Nreactions), stat=ok)
       ALLOCATE(sp%numdest(inp%Nsprate+1,inp%Nreactions), stat=ok)
       ALLOCATE(sp%degform(inp%Nsprate+1,inp%Nreactions), stat=ok)
       ALLOCATE(sp%degdest(inp%Nsprate+1,inp%Nreactions), stat=ok)
       ALLOCATE(sp%Nreacform(inp%Nspecies+1),sp%Nreacdest(inp%Nspecies+1), stat=ok)
    END IF

  END SUBROUTINE ALLO_SPECIES


  SUBROUTINE DEALLO_SPECIES
    
    USE SHARED_VARIABLES, only: inp, sp

    DEALLOCATE(sp%num, sp%name)
    DEALLOCATE(sp%charge, sp%mass)
    DEALLOCATE(sp%element, sp%name_elem)
    DEALLOCATE(sp%Eb_carb, sp%Eb_wat, sp%Eb_HH93, sp%dHf)
    DEALLOCATE(sp%Eb_H2, sp%Eb_pure, sp%Eb_ini, sp%Eb, sp%grproc)
    DEALLOCATE(sp%Ed_carb, sp%Ed_wat, sp%Ed_H2, sp%Ed_pure, sp%Ed)
    DEALLOCATE(sp%cr_coef, sp%Xini, sp%formreac, sp%maxreac) 

    IF (inp%chreacrates .eq. 1) THEN
       DEALLOCATE(sp%rate_name, sp%rate_num)
       DEALLOCATE(sp%rateform, sp%ratedest)
       DEALLOCATE(sp%numform, sp%numdest)
       DEALLOCATE(sp%degform, sp%degdest)
       DEALLOCATE(sp%Nreacform, sp%Nreacdest)
    END IF

  END SUBROUTINE DEALLO_SPECIES



  SUBROUTINE ALLO_REACTIONS

    USE SHARED_VARIABLES, only: inp, re

    ALLOCATE(re%num_react1(inp%Nreactions+1),  stat=ok)
    ALLOCATE(re%react1(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%num_react2(inp%Nreactions+1),  stat=ok)
    ALLOCATE(re%react2(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%num_react3(inp%Nreactions+1),  stat=ok)
    ALLOCATE(re%react3(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%num_prod1(inp%Nreactions+1),  stat=ok)
    ALLOCATE(re%prod1(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%num_prod2(inp%Nreactions+1),  stat=ok)
    ALLOCATE(re%prod2(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%num_prod3(inp%Nreactions+1),  stat=ok)
    ALLOCATE(re%prod3(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%num_prod4(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%prod4(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%A(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%B(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%C(inp%Nreactions+1),stat=ok)
    ALLOCATE(re%type(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%formula(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%Nrate(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%sim(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%Tmin(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%Tmax(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%num(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%act_en(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%proba(inp%Nreactions+1), stat=ok)
    ALLOCATE(re%ndreact1(inp%Nspgr,2), stat=ok)
    ALLOCATE(re%num_ndreact1(inp%Nspgr,2), stat=ok)
    ALLOCATE(re%ndreact2(inp%Nspgr), stat=ok)
    ALLOCATE(re%num_ndreact2(inp%Nspgr), stat=ok)

  END SUBROUTINE ALLO_REACTIONS


  SUBROUTINE DEALLO_REACTIONS

    USE SHARED_VARIABLES, only: re

    DEALLOCATE(re%num_react1, re%react1)
    DEALLOCATE(re%num_react2, re%react2)
    DEALLOCATE(re%num_react3, re%react3)
    DEALLOCATE(re%num_prod1, re%prod1)
    DEALLOCATE(re%num_prod2, re%prod2)
    DEALLOCATE(re%num_prod3, re%prod3)
    DEALLOCATE(re%num_prod4, re%prod4)
    DEALLOCATE(re%A,re%B,re%C, re%formula, re%sim)
    DEALLOCATE(re%type, re%num, re%Tmin, re%Tmax, re%Nrate)
    DEALLOCATE(re%act_en, re%proba)
    DEALLOCATE(re%ndreact1, re%num_ndreact1)
    DEALLOCATE(re%ndreact2, re%num_ndreact2)

  END SUBROUTINE DEALLO_REACTIONS


  SUBROUTINE ALLO_GRID

    USE SHARED_VARIABLES, only: gi, gp
    
    ALLOCATE(gp%nH(gi%Nmodels), gp%T(gi%Nmodels), stat=ok)
    ALLOCATE(gp%XO(gi%Nmodels), stat=ok)
    ALLOCATE(gp%acst(gi%Nmodels), gp%direct(gi%Nmodels), stat=ok)
    ALLOCATE(gp%treatment(gi%Nmodels), gp%enratio(gi%Nmodels), stat=ok)
    ALLOCATE(gp%Fpor(gi%Nmodels), gp%ds(gi%Nmodels), stat=ok)
    ALLOCATE(gp%Ea(gi%Nmodels),gp%treatment2(gi%Nmodels), stat=ok)
    ALLOCATE(gp%EbH(gi%Nmodels), gp%UV(gi%Nmodels), stat=ok)
    ALLOCATE(gp%op(gi%Nmodels),gp%zeta(gi%Nmodels),gp%Av(gi%Nmodels), stat=ok)
    ALLOCATE(gp%Xini_file(gi%Nmodels), stat=ok)

  END SUBROUTINE ALLO_GRID


  SUBROUTINE DEALLO_GRID

    USE SHARED_VARIABLES, only: gp

    DEALLOCATE (gp%nH, gp%T)
    DEALLOCATE (gp%acst, gp%direct)
    DEALLOCATE (gp%treatment, gp%enratio, gp%XO)
    DEALLOCATE (gp%Fpor, gp%ds, gp%Ea)
    DEALLOCATE (gp%EbH, gp%UV, gp%op, gp%zeta, gp%Av)
    DEALLOCATE (gp%Xini_file)

  END SUBROUTINE DEALLO_GRID


  SUBROUTINE ALLO_GRID_INI

    USE SHARED_VARIABLES, only: gi, gp, inp
    
    ALLOCATE(gp%Xgasini(gi%Nmodels,inp%Nspgas), stat=ok)
    ALLOCATE(gp%Xiceini(gi%Nmodels,inp%Nspgr), stat=ok)

   END SUBROUTINE ALLO_GRID_INI


  SUBROUTINE DEALLO_GRID_INI

    USE SHARED_VARIABLES, only: gp

    DEALLOCATE(gp%Xgasini, gp%Xiceini)

  END SUBROUTINE DEALLO_GRID_INI


  SUBROUTINE ALLO_PHYS
 
    USE SHARED_VARIABLES, only: inp, ph

    ALLOCATE(ph%Vth(inp%Nspecies), ph%stick(inp%Nspecies), stat=ok)
    ALLOCATE(ph%Kreac(inp%Nreactions+1), stat=ok)

  END SUBROUTINE ALLO_PHYS


  SUBROUTINE DEALLO_PHYS

    USE SHARED_VARIABLES, only: ph
    
    DEALLOCATE(pH%stick, ph%Vth, ph%Kreac)

  END SUBROUTINE DEALLO_PHYS


  SUBROUTINE ALLO_INPHYS

    USE SHARED_VARIABLES, only: iph

    ALLOCATE(iph%time(iph%Nphysteps), iph%rad(iph%Nphysteps), stat=ok)
    ALLOCATE(iph%veloc(iph%Nphysteps), iph%nH(iph%Nphysteps), iph%Tg(iph%Nphysteps), stat=ok)
    ALLOCATE(iph%Td(iph%Nphysteps), iph%Av(iph%Nphysteps), stat=ok)
    
  END SUBROUTINE ALLO_INPHYS


  SUBROUTINE DEALLO_INPHYS

    USE SHARED_VARIABLES, only: iph

    DEALLOCATE(iph%time, iph%rad, iph%veloc)
    DEALLOCATE(iph%nH, iph%Tg)
    DEALLOCATE(iph%Td, iph%Av)
    
  END SUBROUTINE DEALLO_INPHYS


  SUBROUTINE ALLO_GRAIN(Nlayer) 

    USE SHARED_VARIABLES, only: inp, gr
    
    INTEGER(KIND=LONG), INTENT(in) :: Nlayer

    ALLOCATE(gr%a(Nlayer), stat=ok)
    ALLOCATE(gr%Nssur(Nlayer), stat=ok)
    ALLOCATE(gr%Nstot(Nlayer), stat=ok)
    ALLOCATE(gr%Vd(Nlayer), gr%Vvacuum(Nlayer), stat=ok)
    ALLOCATE(gr%Sdep(Nlayer), gr%Spore(Nlayer), stat=ok)
    ALLOCATE(gr%Fnptot(Nlayer), stat=ok)
    ALLOCATE(gr%Fportot(Nlayer), stat=ok)
    ALLOCATE(gr%Fedgetot(Nlayer), stat=ok)
    ALLOCATE(gr%Npore(Nlayer), stat=ok)

   END SUBROUTINE ALLO_GRAIN


   SUBROUTINE DEALLO_GRAIN

     USE SHARED_VARIABLES, only: gr
     
     DEALLOCATE (gr%a)
     DEALLOCATE (gr%Nssur, gr%Nstot)
     DEALLOCATE (gr%Vd, gr%Vvacuum)
     DEALLOCATE (gr%Sdep, gr%Spore)
     DEALLOCATE (gr%Fnptot, gr%Fportot)
     DEALLOCATE (gr%Fedgetot, gr%Npore)
     
   END SUBROUTINE DEALLO_GRAIN


   SUBROUTINE ALLO_PLOT

     USE SHARED_VARIABLES, only: inp, sp, out

     ALLOCATE(out%rad(out%Ninisteps+1), stat=ok)
     ALLOCATE(out%v(out%Ninisteps+1), stat=ok)
     ALLOCATE(out%time(out%Ninisteps+1), stat=ok)
     ALLOCATE(out%timefix(out%Nplotsteps+1), stat=ok)
     ALLOCATE(out%nH(out%Ninisteps+1), stat=ok)
     ALLOCATE(out%Tg(out%Ninisteps+1), stat=ok)
     ALLOCATE(out%Td(out%Ninisteps+1), stat=ok)
     ALLOCATE(out%Av(out%Ninisteps+1), stat=ok)
     ALLOCATE(out%Nlay(out%Ninisteps+1), stat=ok)
     ALLOCATE(out%Xgas(inp%Nspgas,out%Ninisteps+1), stat=ok)
     ALLOCATE(out%Xice(inp%Nspgr,out%Ninisteps+1), stat=ok)
     ALLOCATE(out%Xsur(inp%Nspgr,out%Ninisteps+1), stat=ok)
     ALLOCATE(out%Xsurml(inp%Nspgr,out%Ninisteps+1), stat=ok)
     IF (inp%chreacrates .eq. 1) THEN
        ALLOCATE (out%rateform(inp%Nsprate,sp%maxNreac,out%Ninisteps),stat=ok)
        ALLOCATE (out%ratedest(inp%Nsprate,sp%maxNreac,out%Ninisteps),stat=ok)
        ALLOCATE (out%rateformtot(inp%Nsprate,sp%maxNreac),stat=ok)
        ALLOCATE (out%ratedesttot(inp%Nsprate,sp%maxNreac),stat=ok)
        ALLOCATE (out%rateform2(inp%Nsprate,sp%maxNreac,out%Ninisteps),stat=ok)
        ALLOCATE (out%ratedest2(inp%Nsprate,sp%maxNreac,out%Ninisteps),stat=ok)
        ALLOCATE (out%numform(inp%Nsprate,sp%maxNreac,out%Ninisteps),stat=ok)
        ALLOCATE (out%numdest(inp%Nsprate,sp%maxNreac,out%Ninisteps),stat=ok)
        ALLOCATE (out%numform2(inp%Nsprate,sp%maxNreac,out%Ninisteps),stat=ok)
        ALLOCATE (out%numdest2(inp%Nsprate,sp%maxNreac,out%Ninisteps),stat=ok)
        ALLOCATE (out%contribform(inp%Nsprate,sp%maxNreac), stat=ok)
        ALLOCATE (out%contribdest(inp%Nsprate,sp%maxNreac), stat=ok)
     END IF

   END SUBROUTINE ALLO_PLOT


   SUBROUTINE DEALLO_PLOT

     USE SHARED_VARIABLES, only: inp, out

     DEALLOCATE (out%nH)
     DEALLOCATE (out%Tg)
     DEALLOCATE (out%Td)
     DEALLOCATE (out%Av)
     DEALLOCATE (out%Xgas)
     DEALLOCATE (out%Xice)
     DEALLOCATE (out%Xsur, out%Xsurml)
     DEALLOCATE (out%rad,out%time,out%timefix)
     IF (inp%chreacrates .eq. 1) THEN
        DEALLOCATE(out%rateform)
        DEALLOCATE(out%ratedest)
        DEALLOCATE(out%rateformtot)
        DEALLOCATE(out%ratedesttot)
        DEALLOCATE(out%rateform2)
        DEALLOCATE(out%ratedest2)
        DEALLOCATE(out%numform)
        DEALLOCATE(out%numdest)
        DEALLOCATE(out%numform2)
        DEALLOCATE(out%numdest2)
        DEALLOCATE(out%contribdest)
        DEALLOCATE(out%contribform)
     END IF

   END SUBROUTINE DEALLO_PLOT


   SUBROUTINE ALLO_ABU

     USE SHARED_VARIABLES, only: inp, out, ab

     ALLOCATE(ab%Xgas(inp%Nspgas), ab%Xgas_prev(inp%Nspgas), stat=ok)
     ALLOCATE(ab%Xtot(inp%Nspecies), ab%Xtot_prev(inp%Nspecies), stat=ok)
     ALLOCATE(ab%Xelem(Nelements), ab%Xelem_prev(Nelements), stat=ok)
     ALLOCATE(ab%Xcharge(2), ab%Xcharge_prev(2), stat=ok)
     ALLOCATE(ab%Xelem_ini(Nelements), ab%Xelem_lay(Nelements), stat=ok)
     ALLOCATE(ab%Xice(inp%Nspgr), ab%Xice_prev(inp%Nspgr), stat=ok)
     ALLOCATE(ab%Xsur(inp%Nspgr), ab%Xsur_prev(inp%Nspgr), stat=ok)
     ALLOCATE(ab%Xbulk(inp%Nspgr), ab%Xbulk_prev(inp%Nspgr), stat=ok)
     ALLOCATE(ab%Psur(inp%Nspgr), ab%Psur_prev(inp%Nspgr), stat=ok)
     ALLOCATE(ab%Pbulk(inp%Nspgr), ab%Pbulk_prev(inp%Nspgr), stat=ok)
     ALLOCATE(ab%Xsurml(inp%Nspgr), ab%Xsurml_prev(inp%Nspgr), stat=ok)
     ALLOCATE(ab%Xbulkml(inp%Nspgr), ab%Xbulkml_prev(inp%Nspgr), stat=ok)
     ALLOCATE(ab%Ptotml(out%Ninilayer), ab%Ptotml_prev(out%Ninilayer), stat=ok)
     ALLOCATE(ab%Pice(inp%Nspgr), ab%Pice_prev(inp%Nspgr), stat=ok)

   END SUBROUTINE ALLO_ABU


   SUBROUTINE DEALLO_ABU

     USE SHARED_VARIABLES, only: ab

     DEALLOCATE(ab%Xgas, ab%Xtot,ab%Xelem,ab%Xcharge,ab%Xice,ab%Xsur,ab%Xbulk)
     DEALLOCATE(ab%Xgas_prev,ab%Xtot_prev,ab%Xelem_prev,ab%Xcharge_prev)
     DEALLOCATE(ab%Xice_prev,ab%Xsur_prev,ab%Xbulk_prev)
     DEALLOCATE(ab%Xelem_lay,ab%Xelem_ini) 
     DEALLOCATE(ab%Psur,ab%Psur_prev,ab%Pice,ab%Pice_prev)
     DEALLOCATE(ab%Pbulk, ab%Pbulk_prev, stat=ok)
     DEALLOCATE(ab%Xbulkml, ab%Xbulkml_prev, stat=ok)

   END SUBROUTINE DEALLO_ABU


   SUBROUTINE ALLO_VARED(Neq)

     USE SHARED_VARIABLES, only: ode
     USE SHARED_VARIABLES2

     INTEGER(KIND=LONG), INTENT(in) :: Neq
     
     ALLOCATE(ode%Y(Neq), ode%Ysave(Neq), ode%Yini(Neq))
     ALLOCATE(ode%ATOL(Neq),ode%ATOLini(Neq),ode%RTOL(Neq), stat=ok)
     ALLOCATE(ode%RPAR(Neq),ode%IPAR(Neq), stat=ok)
     ALLOCATE(ode%RSTATS(22), ode%ISTATS(31), stat=ok)
     ALLOCATE(ode%RWORK(ode%LRW), ode%IWK(ode%LRW), stat=ok)
     ALLOCATE(ode%IWORK(ode%LIW), stat=ok)
     ALLOCATE(ode%IOUT(21), ode%ROUT(6), stat=ok)
     ALLOCATE(ode%nonzero(Neq+1,Neq+1),ode%IA(Neq+1),ode%JA(Neq*Neq), stat=ok)
     ALLOCATE(inzIord(Neq), inzIord2(Neq), inzNZord(Neq), stat=ok)
     ALLOCATE(inzIordini(Neq), inzIord2ini(Neq), Pbulkact(Neq), Pbulkact0(Neq), stat=ok)

   END SUBROUTINE ALLO_VARED


   SUBROUTINE DEALLO_VARED

     USE SHARED_VARIABLES, only: ode
     USE SHARED_VARIABLES2

     DEALLOCATE (ode%Y,ode%Ysave, ode%Yini)
     DEALLOCATE (ode%RWORK,ode%IWORK,ode%nonzero,ode%IA,ode%JA)
     DEALLOCATE (ode%ATOL,ode%RTOL,ode%ATOLini,ode%RPAR,ode%IPAR)
     DEALLOCATE (ode%RSTATS,ode%ISTATS,ode%IOUT,ode%ROUT)
     DEALLOCATE (inzIord,inzIord2,inzIordini,inzIord2ini,inzNZord, Pbulkact,Pbulkact0)

   END SUBROUTINE DEALLO_VARED



END MODULE
