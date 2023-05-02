!-----------------------------------------------------------------------
  SUBROUTINE calc_orbital_magnetization
  !-----------------------------------------------------------------------
  USE kinds,                ONLY : dp
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE wvfct,                ONLY : npw, g2kin, nbnd, npwx, wg, et, &
                                   current_k
  USE orbital_magnetization, ONLY : igk
  USE klist,                ONLY : xk, nks, lgauss
!  USE gsmooth,              ONLY : nrxxs
  USE cell_base,            ONLY : tpiba2, tpiba
!  USE gvect,                ONLY : nr1, nr2, nr3, ngm, ecutwfc, g
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunsat, nwordatwfc
  USE uspp,                 ONLY : nkb, vkb
  USE becmod,               ONLY : becp
!  USE wavefunctions_module, ONLY : evc
  USE lsda_mod,             ONLY : current_spin, lsda, isk, nspin
  USE gipaw_module,         ONLY : lambda_so, dudk_method
!  USE nmr_mod,              ONLY : m_0_atom, m_0, calc_chemical_shift
!  USE paw,                  ONLY : paw_vkb
  USE ener,                 ONLY : ef
!  USE ldaU,                 ONLY : swfcatom, lda_plus_u
  USE control_flags,        ONLY : iverbosity
!  USE paw,                  ONLY : paw_nkb, paw_becp
  USE buffers,              ONLY : get_buffer
!  USE scf,                  ONLY : dvrs, vrs
!  USE mp_global,            ONLY : my_pool_id, me_pool, root_pool, mpime
  USE mp,                   ONLY : mp_barrier
!  USE g_tensor_module,      ONLY : init_paw_2_no_phase, init_us_2_no_phase
!  USE g_tensor_module,      ONLY : q_gipaw, nmr_shift_core


 END SUBROUTINE calc_orbital_magnetization



