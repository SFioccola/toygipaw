!-----------------------------------------------------------------------
SUBROUTINE newscf
!-----------------------------------------------------------------------
!! Recalculate self-consistent
  USE constants, ONLY : bohr_radius_angs, fpi, rytoev
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
  USE klist, ONLY: nkstot, lgauss, ltetra
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
                            iverbosity, do_makov_payne, ldftd3, noinv
  USE extrapolation, ONLY : extrapolate_charge
  USE kinds,         ONLY : dp
  USE input_parameters, ONLY : startingpot, startingwfc, restart_mode
  USE dfunct,        ONLY : newd
  USE scf,        ONLY: rho, rho_core, v, vltot, vrs, kedtau
  USE lsda_mod,   ONLY: nspin, current_spin
  USE orbital_magnetization, ONLY : dvrs
  USE wvfct,    ONLY : g2kin, nbndx, nbnd
  USE gipaw_module, ONLY : conv_threshold, assume_isolated
  USE mp_pools,        ONLY : intra_pool_comm, inter_pool_comm
  USE ener,                 ONLY : ef, ef_up, ef_dw, ef_cond
  USE martyna_tuckerman, ONLY : do_comp_mt
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
  REAL(dp) :: ehomo, elumo ! highest occupied and lowest unoccupied levels
  allocate (dvrs (dfftp%nnr, nspin, 3))

  
  
  !  dft='starting from scratch'
  restart  =.false.
  io_level = 0
  lscf=.true.
  lda_plus_u=.false.
  doublegrid=.false.
  lmovecell=.false.
  iprint=10000
  starting_wfc='file'  ! read wfc from preav. scf
  starting_pot='file'  ! read potential from preav scf
!  starting_wfc='atomic' ! Initialization of wfc
  report=0
  CALL check_stop_init()
  CALL setup_para ( dfftp%nr3, 1, nbnd )
  CALL export_gstart_2_solvers(gstart)
  if ( .not. allocated (btype) ) then
     allocate( btype( nbnd, nkstot ) )
     btype(:,:) = 1
  end if
  !
  !
  SELECT CASE( trim( assume_isolated ) )
      !
    CASE( 'makov-payne', 'm-p', 'mp' )
      !
      do_makov_payne = .true.
      !
    CASE( 'martyna-tuckerman', 'm-t', 'mt' )
      !
      do_comp_mt     = .true.
      !
      CONTINUE
      !
    CASE DEFAULT
      !
      do_makov_payne = .false.
      do_comp_mt     = .false.
      !
  END SELECT
  !  since we use only Gamma we don't need symmetries
  !
  nsym=1
  noinv=.true.
  !
  ! these must be tuned for fast convergence
  !
  david = 6
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
  call summary ( )
  call hinit0 ( )
  call potinit ( )
  
  CALL set_dvrs( dvrs, vrs, dfftp%nnr, nspin )

  call newd ( )
  call wfcinit_gipaw ( )
  CALL electrons_gipaw ( )
  
  ! compute the Fermi Energy for not a metal 
  IF ( .NOT. lgauss .OR. ltetra ) THEN
          ! ... presumably not a metal
          CALL get_homo_lumo (ehomo, elumo)
          ef = (ehomo + elumo)/2.d0
  ENDIF

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




