MODULE FUNCS

  USE VARIABLES

  IMPLICIT NONE

  !----------------------------------------------------------------------------
  ! This module contains various functions that compute: 
  ! the physical conditions, energies, rates, transmission probabilities
  !----------------------------------------------------------------------------


CONTAINS 





  !----------!
  !----------!
  ! ENERGIES !
  !----------!
  !----------!


  ! derivative of the Eb evolution with the coverage of H2 on the external layer
!!$  REAL(KIND=DP) FUNCTION Eb_evol_prim(PH2,num_spec)
!!$    
!!$    REAL(KIND=DP), INTENT(in) :: PH2
!!$    INTEGER(KIND=LONG), INTENT(in) :: num_spec
!!$    REAL(KIND=DP) :: EbH2O, EbH2
!!$    
!!$    EbH2O = spec_Eb_ini(num_spec)
!!$    EbH2 = spec_Eb_H2(num_spec)
!!$
!!$    Eb_evol_prim = -EbH2O + EbH2
!!$
!!$  END FUNCTION Eb_evol_prim



  !-------!
  !-------!
  ! RATES !
  !-------!
  !-------!





  !----------------------------------------------------!
  !----------------------------------------------------!
  ! ACTIVATION ENERGIES AND TRANSMISSION PROBABILITIES !
  !----------------------------------------------------!
  !----------------------------------------------------!



END MODULE FUNCS

