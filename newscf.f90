!-----------------------------------------------------------------------
SUBROUTINE newscf
!-----------------------------------------------------------------------
!! Recalculate self-consistent
  USE constants, ONLY : bohr_radius_angs, fpi
  USE io_global, ONLY : stdout
  USE cell_base, ONLY : omega
  USE ions_base, ONLY : nat, tau, ityp, atm
  USE fft_base,  ONLY : dfftp, dffts
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
  USE fft_base,      ONLY : dfftp, dffts
  USE symm_base,     ONLY : nsym
  USE io_files,      ONLY : iunwfc, prefix, tmp_dir, postfix
  USE ldaU,          ONLY : lda_plus_u
  USE control_flags, ONLY : restart, io_level, lscf, iprint, &
                            david, max_cg_iter, nexxiter, &
                            isolve, tr2, ethr, mixing_beta, nmix, niter, &
                            iverbosity
  USE extrapolation, ONLY : extrapolate_charge
  USE kinds,         ONLY : dp
  USE input_parameters, ONLY : startingpot, startingwfc, restart_mode
  USE dfunct,        ONLY : newd
  USE scf,        ONLY: rho, rho_core, v, vltot, vrs, kedtau
  USE lsda_mod,   ONLY: nspin, current_spin
  USE orbital_magnetization, ONLY : dvrs
  USE wvfct,    ONLY : g2kin, nbndx, nbnd
  USE gipaw_module, ONLY : conv_threshold
  USE mp_pools,        ONLY : intra_pool_comm, inter_pool_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: i, j, nrxxs
  CHARACTER(LEN=256) :: dirname
  INTEGER :: iter
  LOGICAL               :: exst
  REAL(DP), DIMENSION(3,3) :: chi(3,3)
  REAL(DP), DIMENSION(dfftp%nnr,1) ::  rhotot, sign_r
  REAL(DP) :: exxen
  
  
  !  dft='starting from scratch'
  restart  =.false.
  io_level = 0
  lscf=.true.
  lda_plus_u=.false.
  doublegrid=.false.
  lmovecell=.false.
  iprint=10000
!  starting_wfc='file'  ! read wfc from preav. scf
!  starting_pot='file'  ! read potential from preav scf
  starting_wfc='atomic' ! Initialization of wfc
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
  isolve=0
  ethr = conv_threshold
  !
  nmix=8
  niter=100
  nexxiter=100
  !
  iunwfc  = 9 ! change the number of unit iunwfc
  call openfil
  call hinit0 ( )
  call potinit ( )
  call newd ( )
  call wfcinit_gipaw ( )
  CALL electrons_gipaw ( )

!!! Alternative routine to perform SCF:
!  call openfil 

!  dirname = TRIM(tmp_dir) //TRIM(prefix) // postfix
!  CALL extrapolate_charge( dirname, 1 )
!  CALL hinit1
!  CALL electrons
!  exxen = 0.d0
  !
  CLOSE(unit=iunwfc, status='keep')
  !
  !
  RETURN
END SUBROUTINE newscf




