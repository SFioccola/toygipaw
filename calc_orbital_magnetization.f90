!-----------------------------------------------------------------------
  SUBROUTINE calc_orbital_magnetization
  !-----------------------------------------------------------------------
  USE kinds,                 ONLY : dp
  USE io_files,              ONLY : nwordwfc, iunwfc
  USE wvfct,                 ONLY : npw, g2kin, nbnd, npwx, wg, et, &
                                   current_k
  USE klist,                 ONLY : xk, nks, lgauss, igk_k, ngk
  USE fft_base,              ONLY : dffts, dfftp
  USE cell_base,             ONLY : tpiba2, tpiba
  USE gvecw,                 ONLY : ecutwfc
  USE gvect,                 ONLY : ngm, g
  USE io_global,             ONLY : stdout
  USE io_files,              ONLY : iunsat, nwordatwfc
  USE uspp,                  ONLY : nkb, vkb
  USE becmod,                ONLY : becp
  USE wavefunctions,         ONLY : evc
  USE lsda_mod,              ONLY : current_spin, lsda, isk, nspin
  USE gipaw_module,          ONLY : lambda_so, dudk_method
  USE nmr_mod,               ONLY : m_0_atom, m_0, calc_chemical_shift
  USE paw_gipaw,             ONLY : paw_vkb, paw_nkb, paw_becp
  USE ener,                  ONLY : ef
!  USE ldaU,                 ONLY : swfcatom, lda_plus_u
  USE control_flags,         ONLY : iverbosity
  USE buffers,               ONLY : get_buffer
  USE scf,                   ONLY : vrs
  USE mp_global,             ONLY : my_pool_id, me_pool, root_pool
  USE mp,                    ONLY : mp_barrier
  USE gipaw_module,          ONLY : q_gipaw, nmr_shift_core
  USE orbital_magnetization, ONLY : delta_k
  implicit none
  complex(dp), external :: ZDOTC
  integer, parameter :: iundudk1 = 75, iundudk2 = 76, iundudk3 = 77
  real(dp), parameter :: rydtohar = 0.5d0
  complex(dp), allocatable :: dudk_bra(:,:), dudk_ket(:,:), hpsi(:)
  complex(dp), allocatable :: vkb_save(:,:), aux(:,:)
  complex(dp) :: braket
  real(dp) :: kp_berry(3), kp_M_LC(3), kp_M_IC(3), tmp1(3), tmp2(3)
  integer :: ik, ibnd, jbnd, kk, ii, jj, occ, nrxxs, nr1, nr2, nr3
  real(dp) :: tmp(3), emin, emax
  ! index for the cross product
  integer :: ind(2,3)
  ind(:,1) = (/ 2, 3 /)
  ind(:,2) = (/ 3, 1 /)
  ind(:,3) = (/ 1, 2 /)

  ! various initializations:
  nrxxs = dffts%nnr
  nr1=dfftp%nr1
  nr2=dfftp%nr2
  nr3=dfftp%nr3

  do ik = 1, nks
    npw = ngk (ik) !number of plane waves for each k point
  enddo
  delta_k = q_gipaw/tpiba
  
  ! compute the covariant derivatives
!  call compute_dudk(dudk_method)
  ! allocate memory
!  allocate (dudk_bra(npwx,nbnd), dudk_ket(npwx,nbnd), hpsi(npwx))


  END SUBROUTINE calc_orbital_magnetization



