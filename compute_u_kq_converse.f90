
! Copyright (C) 2001-2007 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE compute_u_kq_converse( ik, q, iter, ethr_q )
  !----------------------------------------------------------------------------
  !
  ! ... diagonalize at k+q
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : nwordatwfc, iunsat, iunwfc, &
                                   nwordwfc
  USE klist,                ONLY : nkstot, nks, xk, ngk, igk_k
  USE uspp,                 ONLY : vkb, nkb
  USE gvect,                ONLY : g, ngm, ngl
  USE gvecw,                ONLY : gcutw
  USE wvfct,                ONLY : et, nbnd, npwx, npw, current_k, &
                                   nbndx, g2kin
  USE control_flags,        ONLY : ethr, isolve, io_level, lscf, istep
!  USE ldaU,                 ONLY : lda_plus_u, swfcatom
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions,        ONLY : evc
  USE bp,                   ONLY : lelfield
  USE random_numbers,       ONLY : randy
  USE constants,            ONLY : tpi, RytoeV
  USE cell_base,            ONLY : tpiba2
  USE buffers,              ONLY : get_buffer
  USE input_parameters,     ONLY : etot_conv_thr
  USE uspp_init,            ONLY : init_us_2
  USE gipaw_module,         ONLY : q_gipaw, evq, conv_threshold, iverbosity
  IMPLICIT NONE
  INTEGER :: ik, iter       ! k-point, current iterations
  REAL(DP) :: q(3)          ! q-vector
  REAL(DP) :: avg_iter, ethr_q
  INTEGER :: i, ig
  REAL(DP) :: xkold(3)
  REAL(DP), allocatable :: et_old(:,:)
  REAL(DP) :: rr, arg

  CALL start_clock( 'c_bands' )
  print*, 'entra in u k+q coverse'

  iter = 2
  istep = 0
  lscf = .false.
  ethr = conv_threshold

  npw = ngk(ik)

  ! save eigenvalues
  allocate( et_old(nbnd,nkstot) )
  et_old = et

  !! debug
    write(stdout, '(5X,"compute_u_kq: q = (",F10.4,",",F10.4,",",F10.4,")")') q

  IF (iverbosity > 1) THEN
    IF ( isolve == 0 ) THEN
       WRITE( stdout, '(5X,"Davidson diagonalization with overlap")' )
    ELSE IF ( isolve == 1 ) THEN
       WRITE( stdout, '(5X,"CG style diagonalization")')
    ELSE
       CALL errore ( 'c_bands', 'invalid type of diagonalization', isolve)
       !!! WRITE( stdout, '(5X,"DIIS style diagonalization")')
    END IF
  ENDIF

  avg_iter = 0.D0

  current_k = ik
  IF ( lsda ) current_spin = isk(ik)
  CALL gk_sort( xk(1,ik), ngm, g, gcutw, npw, igk_k(1,ik), g2kin )

  ! set the k-point
  xkold(:) = xk(:,ik)
  xk(:,ik) = xk(:,ik) + q(:)

  ! various initializations
  IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )

  ! kinetic energy
  call g2_kin(ik)

  ! read in wavefunctions from the previous iteration
  CALL get_buffer(evc, nwordwfc, iunwfc, ik)

  ! Needed for LDA+U
!  IF ( lda_plus_u ) CALL davcio( swfcatom, nwordatwfc, iunsat, ik, -1 )

  ! diagonalization of bands for k-point ik
  call diag_bands_gipaw ( iter, ik, avg_iter )

!  call poolreduce( 1, avg_iter )
  avg_iter = avg_iter / nkstot
  if (iverbosity > 0) &
    WRITE( stdout, &
         '( 5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1 )' ) &
         ethr, avg_iter

  ! check if diagonalization was ok
  write(stdout,'(8F9.4)') et_old(1:nbnd,ik)*RytoeV
  write(stdout,'(8F9.4)') et(1:nbnd,ik)*RytoeV
  do i = 1, nbnd
    if (abs(et(i,ik) - et_old(i,ik))*RytoeV > 0.2d0) then
      write(stdout,'(5X,''ATTENTION: ik='',I4,''  ibnd='',I3,$)') ik, i
      write(stdout,'(2X,''eigenvalues differ too much!'')')
      write(stdout,'(5X,2(F10.4,2X))') et_old(i,ik)*RytoeV, et(i,ik)*RytoeV
    endif
  enddo


  ! restore the k-point and eigenvalues
  xk(:,ik) = xkold(:)
  et = et_old
  deallocate(et_old)

  CALL stop_clock( 'c_bands' )  
  RETURN

END SUBROUTINE compute_u_kq_converse
