!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!-----------------------------------------------------------------------
! Compute the covariant derivative
!-----------------------------------------------------------------------
SUBROUTINE compute_dudk_converse(dudk_method)
  USE kinds,     ONLY : dp
  USE io_files,  ONLY : nwordwfc, iunwfc, diropn  
  USE wvfct,     ONLY : npw, g2kin, nbnd, npwx
  USE klist,     ONLY : xk, nks, lgauss, igk_k
  USE io_global, ONLY : stdout
  USE mp_global, ONLY : my_pool_id, me_pool, root_pool
  USE mp,        ONLY : mp_barrier, mp_sum
  USE mp_pools,  ONLY : inter_pool_comm
  USE gipaw_module,  ONLY : q_gipaw, alpha_pv
  implicit none
  integer, parameter :: iundudk1 = 75, iundudk2 = 76, iundudk3 = 77
  character(80), intent(in) :: dudk_method
  complex(dp), allocatable :: dudk(:,:,:)
  real(dp):: emin, emax
  logical :: exst
  integer :: ik, ipol, occ
  
  print*, 'entra in dudk converse'
  ! open files
#ifdef __MPI
  call mp_barrier (inter_pool_comm)
  call diropn(iundudk1, 'dudk1_', 2*nwordwfc, exst)
  call diropn(iundudk2, 'dudk2_', 2*nwordwfc, exst)
  call diropn(iundudk3, 'dudk3_', 2*nwordwfc, exst)
#else
  call diropn(iundudk1, 'dudk1', 2*nwordwfc, exst)
  call diropn(iundudk2, 'dudk2', 2*nwordwfc, exst)
  call diropn(iundudk3, 'dudk3', 2*nwordwfc, exst)
#endif

  ! allocate covariant derivatives
  allocate (dudk(npwx,nbnd,3) )

  write(stdout,*)
  write(stdout,'(5X,''Computing du/dk '',$)')

  if (trim(dudk_method) == 'read') then
    write(stdout,'(''(read from disk)'')')

  elseif (trim(dudk_method) == 'covariant') then
    write(stdout,'(''(covariant derivative)'')') 
    do ik = 1, nks
      call find_nbnd_occ_converse(ik, occ, emin, emax)
      write(stdout,'(5X,''k-point:'',I4,4X,''occ='',I3)') ik, occ
      call dudk_covariant_converse(ik, occ, dudk)
      do ipol = 1, 3
        call davcio( dudk(1,1,ipol), 2*nwordwfc, iundudk1 + ipol - 1, ik, +1 )
      enddo
    enddo

!  elseif (trim(dudk_method) == 'singlepoint') then
!    write(stdout,'(''(single point)'')') 
!    do ik = 1, nks
!      call find_nbnd_occ(ik, occ, emin, emax)
!      write(stdout,'(5X,''k-point:'',I4,4X,''occ='',I3)') ik, occ
!      call dudk_covariant_single_point(ik, occ, dudk)
!      do ipol = 1, 3
!        call davcio( dudk(1,1,ipol), 2*nwordwfc, iundudk1 + ipol - 1, ik, +1 )
!      enddo
!    enddo
    
!  elseif (trim(method) == 'kdotp') then
!    write(stdout,'(''(k \dot p perturbation)'')') 
!    do ik = 1, nks
!      call find_nbnd_occ(ik, occ, emin, emax)
!      ! determine alpha_pv
!      alpha_pv = emax - emin
!      if (.not. lgauss) alpha_pv = 2.d0 * alpha_pv
!      alpha_pv = max(alpha_pv, 1.d-2)
!      if (me_pool == root_pool) &
!        write(*,'(5X,''k-point:'',I4,4X,''pool:'',I4,4X,''occ='',I3,4X,''alpha_pv='',F10.4)') ik, my_pool_id+1, occ, alpha_pv
!      call dudk_kdotp(ik, occ, alpha_pv, dudk(1,1,1))
!      do ipol = 1, 3
!        call davcio( dudk(1,1,ipol), 2*nwordwfc, iundudk1 + ipol - 1, ik, +1 )
!      enddo
!    enddo

  elseif (trim(dudk_method) == 'null') then
    write(stdout,'(''Skipping du/dk calculation!!!!'')')

  else
    write(stdout,*)
    call errore('compute_dudk', 'unknown du/dk method: '//trim(dudk_method), 1)
  endif

!  write(stdout,'(5X,''done with du/dk'')')
  write(stdout,'(5X,''done with du/dk'',3X,''q_gipaw = '',F8.6)') q_gipaw
  deallocate(dudk)

END SUBROUTINE compute_dudk_converse



!-----------------------------------------------------------------------
! covariant derivative
!-----------------------------------------------------------------------
  SUBROUTINE dudk_covariant_converse(ik, occ, dudk)
  USE kinds,                ONLY : dp  
  USE cell_base,            ONLY : bg, tpiba2, tpiba
  USE wvfct,                ONLY : npw, nbnd, npwx
  USE wavefunctions,        ONLY : evc
  USE klist,                ONLY : xk
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE io_global,            ONLY : stdout
  USE buffers,              ONLY : get_buffer, save_buffer
  USE gipaw_module,         ONLY : q_gipaw, evq, conv_threshold
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp,                   ONLY : mp_barrier, mp_sum                         
  IMPLICIT NONE
  complex(dp), external :: zdotc
  complex(dp), allocatable :: overlap(:,:), evc0(:,:)
  complex(dp) :: dudk(npwx,nbnd,3)
  real(dp) :: q_gipaw3(3), ethr_q, delta_k
  integer :: ik, i, occ, sig, iter
  integer :: ibnd, jbnd, ipol

  if (occ == 0) return

  ! allocate overlap matrix and evc0
  allocate(overlap(occ,occ), evc0(npwx,nbnd))

  ! read the wavefunction
  call get_buffer(evc, nwordwfc, iunwfc, ik)
  q_gipaw3(:) = 0.d0
  ethr_q = conv_threshold

  call compute_u_kq_converse(ik, q_gipaw3, iter, ethr_q)
  call save_buffer(evc, nwordwfc, iunwfc, ik)
  evc0(1:npwx,1:nbnd) = evc(1:npwx,1:nbnd)
!  delta_k = q_gipaw/2.d0
  delta_k = q_gipaw/tpiba

  ! loop over crystal directions
  do ipol = 1, 3
    dudk(1:npwx,1:occ,ipol) = (0.d0,0.d0)
    
    ! loop over +/-1
    do sig = -1, 1, 2

      ! set the k-point
      q_gipaw3(:) = 0.d0
      q_gipaw3(ipol) = delta_k * sig
      !!write(*,'(5X,''overlap: ik='',I4,'' ipol='',I1,'' sig='',I2)') &
      !!      ik, ipol, sig
      evc(1:npwx,1:nbnd) = evc0(1:npwx,1:nbnd)
      ethr_q = conv_threshold
      call compute_u_kq_converse(ik, q_gipaw3, iter, ethr_q)

      ! compute overlaps
      do ibnd = 1, occ
        do jbnd = 1, occ
          overlap(ibnd,jbnd) = zdotc(npw, evc0(1,ibnd), 1, evc(1,jbnd), 1)
       enddo
      enddo
#if defined(__MPI)
  CALL mp_sum( overlap, intra_bgrp_comm )
#endif


      !!write(*,'(5X,''inverse: ik='',I4,'' ipol='',I1,'' sig='',I2)') &
      !!      ik, ipol, sig
      call invert_matrix_converse(occ, overlap)

      ! compute the covariant derivative
      !!write(*,'(5X,''dual   : ik='',I4,'' ipol='',I1,'' sig='',I2)') &
      !!      ik, ipol, sig
      do ibnd = 1, occ
        do jbnd = 1, occ
          dudk(1:npw,ibnd,ipol) = dudk(1:npw,ibnd,ipol) + &
                        sig * 0.5d0/(delta_k*tpiba) * &
                        overlap(jbnd,ibnd) * evc(1:npw,jbnd)
        enddo
      enddo

    enddo ! sig
  enddo ! ipol

  deallocate(overlap, evc0)
END SUBROUTINE dudk_covariant_converse





!-----------------------------------------------------------------------
! (k dot p) perturbation
!-----------------------------------------------------------------------
!SUBROUTINE dudk_kdotp(ik, occ, alpha_pv, dudk)
!  USE kinds,                ONLY : dp  
!  USE cell_base,            ONLY : bg, tpiba2, tpiba
!  USE wvfct,                ONLY : npw, nbnd, npwx, current_k, igk, g2kin
!  USE lsda_mod,             ONLY : lsda, current_spin, isk
!  USE wavefunctions_module, ONLY : evc
!  USE klist,                ONLY : xk
!  USE io_files,             ONLY : nwordwfc, iunwfc, nwordatwfc, iunsat
!  USE io_global,            ONLY : stdout
!  USE buffers,              ONLY : get_buffer, save_buffer
!  USE ldaU,                 ONLY : swfcatom, lda_plus_u
 ! USE paw,                  ONLY : paw_vkb, paw_becp, paw_nkb
!  USE uspp,                 ONLY : nkb, vkb
!  USE gvect,                ONLY : ngm, gcutm, ecutwfc, g
!  USE becmod,               ONLY : becp
!  USE g_tensor_module,      ONLY : q_conv_thr, q_gipaw
!  IMPLICIT NONE
!  real(dp) :: alpha_pv
!  complex(dp), external :: zdotc
!  complex(dp), allocatable :: overlap(:,:), evc0(:,:), v_evc(:,:)
!  complex(dp) :: dudk(npwx,nbnd,3)
!  real(dp) :: ethr_q
!  integer :: ik, i, occ, sig, iter
!  integer :: ibnd, jbnd, ipol

!  allocate(v_evc(npwx,nbnd))

  !!if (occ == 0) return
!  dudk(:,:,:) = (0.d0,0.d0)

  ! read the wavefunction
!  call get_buffer(evc, nwordwfc, iunwfc, ik)

!  ethr_q = q_conv_thr
!  current_k = ik
!  if (lsda) current_spin = isk(ik)
!  if (lda_plus_u) call davcio(swfcatom, nwordatwfc, iunsat, ik, -1)

  ! loop over crystal directions
!  do ipol = 1, 3
!    dudk(1:npwx,1:occ,ipol) = (0.d0,0.d0)

    ! apply the velocity operator
!    call apply_vel(evc, v_evc, ik, ipol)
    !!v_evc(:,:) = (0.d0,1.d0)*v_evc(:,:)

    ! solve (k dot p) system
!    call greenfunction(ik, v_evc, occ, alpha_pv, dudk(1,1,ipol), ethr_q)

!  enddo

!END SUBROUTINE dudk_kdotp





!-----------------------------------------------------------------------
SUBROUTINE invert_matrix_converse(n, a)
!-----------------------------------------------------------------------
  USE kinds, ONLY : dp
  implicit none
  integer :: n, ipiv(n), info, lwork
  complex(dp) :: a(n,n)  
  complex(dp), allocatable, dimension(:) :: work

  if (n == 0) return
  if (n == 1) then
    a(1,1) = 1.d0/a(1,1)
    return
  endif

  call ZGETRF(n, n, a, n, ipiv, info)
  if (info /= 0) call errore('invert_matrix', 'ZGETRF info=', info)

  allocate(work(1))
  lwork = -1
  call ZGETRI(n, a, n, ipiv, work, lwork, info)
  if (info /= 0) call errore('invert_matrix', 'ZGETRI(1) info=', info)

  lwork = nint(real(work(1)))
  deallocate(work)
  allocate(work(lwork))
  call ZGETRI(n, a, n, ipiv, work, lwork, info)
  if (info /= 0) call errore('invert_matrix', 'ZGETRI(2) info=', info)
END SUBROUTINE invert_matrix_converse





!-----------------------------------------------------------------------
! apply e^(i G r) to a wavefunction in real space
!-----------------------------------------------------------------------
!SUBROUTINE apply_eigr(sig, gv, psi, res)
!  USE kinds,     ONLY : dp
!  USE constants, ONLY : tpi
!  USE gvect,     ONLY : nr1, nr2, nr3
!  USE gsmooth,   ONLY : nr1s, nr2s, nr3s, nrxxs, nrx1s, nrx2s, nrx3s
!  USE mp_global, ONLY : me_pool
!  USE pffts,     ONLY : npps
!  implicit none
!  real(dp) :: gv(3), phase
!  complex(dp) :: psi(nrxxs), res(nrxxs)
!  integer :: sig, i1, i2, i3, itmp, izub, izlb, ir

!#ifdef __PARA
!     izub=0
!     do itmp=1,me_pool + 1
!        izlb=izub+1
!        izub=izub+npps(itmp)
!     enddo
!#else
!     izlb=1
!     izub=nr3s
!#endif

!  do i1 = 1, nr1s
!    do i2 = 1, nr2s
!      do i3 = izlb, izub
!        phase = sig* tpi * ( gv(1) * dble(i1)/dble(nr1s) + &
!                             gv(2) * dble(i2)/dble(nr2s) + &
!                             gv(3) * dble(i3)/dble(nr3s) )
!        itmp = i3 - izlb + 1
!        ir = i1 + (i2-1)*nrx1s + (itmp-1)*nrx1s*nrx2s
!        res(ir) = cmplx(cos(phase),sin(phase)) * psi(ir)
!      enddo
!    enddo
!  enddo
!END SUBROUTINE apply_eigr




!-----------------------------------------------------------------------
SUBROUTINE find_nbnd_occ_converse(ik, nbnd_occ, emin, emax)
!-----------------------------------------------------------------------
  USE kinds,     ONLY : dp
  USE wvfct,     ONLY : wg, nbnd
  USE klist,     ONLY : wk, lgauss, degauss, ngauss
  USE constants, ONLY : pi
  USE pwcom
  implicit none
  integer, intent(in) :: ik
  integer, intent(out):: nbnd_occ
  real(dp), intent(out) :: emin, emax
  real(dp) :: small, xmax, fac, e_target
  integer :: ibnd

  if (lgauss) then ! metallic case
    small = 6.9626525973374d-5
    xmax = sqrt(-log(sqrt(pi)*small))
    if (ngauss == -99) then
      fac = 1.d0 / sqrt(small)
      xmax = 2.d0 * log(0.5d0*(fac+sqrt(fac*fac-4.d0)))
    endif
    e_target = ef + xmax * degauss
  endif

  nbnd_occ = 0
  emin = 1d6
  emax = -1d6
  do ibnd = 1, nbnd
    if (lgauss) then ! metallic case
      emin = min(emin, et(ibnd,ik))
      emax = e_target
      if (et(ibnd,ik) < e_target) nbnd_occ = ibnd
    elseif (wg(ibnd,ik) > 1d-4*wk(ik)) then ! insulator
      emin = min(emin, et(ibnd,ik))
      emax = max(emax, et(ibnd,ik))
      nbnd_occ = nbnd_occ + 1
    endif
  enddo

END SUBROUTINE find_nbnd_occ_converse


