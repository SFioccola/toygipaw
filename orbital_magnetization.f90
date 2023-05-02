!-----------------------------------------------------------------------
MODULE orbital_magnetization
  !-----------------------------------------------------------------------
  !
  ! ... the title is self explaining
  !
  USE kinds, ONLY : DP
  USE gipaw_module, ONLY: dudk_method
  IMPLICIT NONE
  SAVE
  REAL(DP) :: orb_magn_LC(3)
  REAL(DP) :: orb_magn_IC(3)
  REAL(DP) :: berry_curvature(3)
  REAL(DP) :: orb_magn_tot(3)
  REAL(DP) :: delta_M_bare(3)
  REAL(DP) :: delta_M_para(3)
  REAL(DP) :: delta_M_dia(3)
  COMPLEX(dp), ALLOCATABLE :: dbecp(:,:,:)
  COMPLEX(dp), ALLOCATABLE :: paw_dbecp(:,:,:)
  INTEGER, ALLOCATABLE :: igk(:,:)
  real(dp) :: delta_k


 !-----------------------------------------------------------------------
END MODULE orbital_magnetization
!-----------------------------------------------------------------------


