! toygipaw
PROGRAM toygipaw

  USE constants, ONLY : eps12
  USE control_flags, ONLY : gamma_only
  USE environment, ONLY : environment_start, environment_end
  USE io_files, ONLY : prefix, tmp_dir
  USE io_global, ONLY : ionode, ionode_id
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin, starting_magnetization !FZ: added starting_magnetization for spinors
  USE mp, ONLY : mp_bcast
  USE mp_world, ONLY : world_comm
  USE mp_global, ONLY : mp_startup
  USE paw_variables, ONLY : okpaw
  USE scf, ONLY : rho_core, rhog_core
  USE uspp, ONLY : okvan
  !USE funct, ONLY : dft_is_hybrid, dft_is_meta, dft_is_gradient !FZ: for hybrid functional, metaGGA, add dft_is_gradient for spinors
  USE xc_lib, ONLY : xclib_dft_is  !FZ: for qe6.7, for hybrid functional, metaGGA, add dft_is_gradient for spinors
  USE ldaU, ONLY : lda_plus_u
  USE scf,   ONLY : rho, rho_core, rhog_core, v, vrs!FZ: test spinors added vrs 
  USE fft_rho,  ONLY : rho_g2r  
  USE fft_base,  ONLY : dfftp    
  USE noncollin_module,   ONLY : noncolin, npol, m_loc, ux, angle1, angle2  !FZ: test for spinors
  USE ions_base,          ONLY : nat, tau, ntyp => nsp, ityp, zv  !FZ: test for spinors

  IMPLICIT NONE

  character(len=6) :: codename = 'TOYGIPAW'
  character(len=80) :: diagonalization, dudk_method
  real(dp) :: q_gipaw, diago_thr_init, conv_thr
  character ( len = 256 ) :: outdir
  character (len=256), external :: trimcheck
  integer :: ios
  

  NAMELIST / input_toygipaw / prefix, outdir, &
                        diagonalization, q_gipaw, dudk_method, &
                        diago_thr_init, conv_thr

#if defined(__MPI)
  CALL mp_startup ( )
#endif

  CALL environment_start ( codename )

! read input file
  prefix = 'prefix'
  CALL get_environment_variable ( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM ( outdir ) == ' ' ) outdir = './'
  diagonalization = 'david'
  q_gipaw = 0.01d0
  dudk_method = 'covariant'
  diago_thr_init = 1d-4
  conv_thr = 1e-8

  IF ( ionode ) THEN
    CALL input_from_file ( )
    READ ( 5, input_toygipaw, iostat = ios )
    IF ( ios /= 0 ) CALL errore ( codename, 'toygipaw', abs ( ios ) )
  ENDIF

  tmp_dir = trimcheck ( outdir )
  CALL mp_bcast ( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast ( prefix, ionode_id, world_comm )
  CALL mp_bcast ( diagonalization, ionode_id, world_comm )
  CALL mp_bcast ( q_gipaw, ionode_id, world_comm )
  CALL mp_bcast ( dudk_method, ionode_id, world_comm )
  CALL mp_bcast ( diago_thr_init, ionode_id, world_comm )
  CALL mp_bcast ( conv_thr, ionode_id, world_comm )

!read ground state wavefunctions  
  CALL read_file ( )

!-------------------------------------------------------------------------------

END PROGRAM toygipaw
