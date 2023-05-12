!-----------------------------------------------------------------------
SUBROUTINE newscf
!-----------------------------------------------------------------------
!! Recalculate self-consistent
  USE constants, ONLY : bohr_radius_angs, fpi
  USE io_global, ONLY : stdout
  USE cell_base, ONLY : omega
  USE ions_base, ONLY : nat, tau
  USE fft_base,  ONLY : dfftp
  USE scf,       ONLY : rho, rho_core
  USE xc_lib,    ONLY : xclib_set_threshold, dmxc
  USE basis, ONLY: starting_wfc, starting_pot, natomwfc
  USE cellmd,ONLY: lmovecell
  USE gvecs, ONLY: doublegrid
  USE gvect, ONLY: gstart
  USE wvfct, ONLY: btype
  USE klist, ONLY: nkstot
  USE wvfct, ONLY: nbnd, nbndx
  USE noncollin_module, ONLY: report
  USE check_stop,    ONLY : check_stop_init
  USE fft_base,      ONLY : dfftp
  USE symm_base,     ONLY : nsym
  USE io_files,      ONLY : iunwfc, prefix, tmp_dir, postfix
  USE ldaU,          ONLY : lda_plus_u
  USE control_flags, ONLY : restart, io_level, lscf, iprint, &
                            david, max_cg_iter, nexxiter, &
                            isolve, tr2, ethr, mixing_beta, nmix, niter
  USE extrapolation, ONLY : extrapolate_charge
  USE kinds,         ONLY : dp
  USE input_parameters, ONLY : startingpot, startingwfc, restart_mode
  USE dfunct,        ONLY : newd

  !
  IMPLICIT NONE
  !
  INTEGER :: i, j
  CHARACTER(LEN=256) :: dirname
  INTEGER :: iter
  LOGICAL               :: exst
  REAL(DP), DIMENSION(3,3) :: chi(3,3)
  REAL(DP), DIMENSION(dfftp%nnr,1) ::  rhotot, sign_r
  REAL(DP) :: exxen

  !  dft='Same as Before'
!  restart  =.false.
  io_level = 0
  lscf=.true.
  lda_plus_u=.false.
  doublegrid=.false.
  lmovecell=.false.
  iprint=10000
!  starting_wfc='file'
  starting_wfc='atomic'
  starting_pot='file'
  report=0
  CALL check_stop_init()
  CALL setup_para ( dfftp%nr3, 1, nbnd )
  CALL export_gstart_2_solvers(gstart)
  if ( .not. allocated (btype) ) then
     allocate( btype( nbnd, nkstot ) )
     btype(:,:) = 1
  end if
  !
  !  since we use only Gamma we don't need symmetries
  !
  nsym=1
  !
  ! these must be tuned for fast convergence
  !
  david = 4
  nbndx = max (nbndx, david*nbnd)
  max_cg_iter=20
  isolve=0
  tr2 =1.d-8
  ethr=1.d-8
  mixing_beta=0.7d0
  nmix=4
  niter=50
  nexxiter=100
  !
  iunwfc  = 9 ! change the number of unit iunwfc
  call openfil
  call hinit0 ( )
  call potinit ( )
  call newd ( )
  call wfcinit ( )
!  call openfil 

!  dirname = TRIM(tmp_dir) //TRIM(prefix) // postfix
!  CALL extrapolate_charge( dirname, 1 )
!  CALL hinit1
  CALL electrons_gipaw ( )
!  exxen = 0.d0
!  call electrons_scf (2, exxen )
  !
  CLOSE(unit=iunwfc, status='keep')
  !
  !
  RETURN
END SUBROUTINE newscf




