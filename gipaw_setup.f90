!-----------------------------------------------------------------------
SUBROUTINE gipaw_setup
  !-----------------------------------------------------------------------
  !
  ! ... GIPAW setup
  !
  USE kinds,         ONLY : dp
  USE io_global,     ONLY : stdout
  USE wvfct,         ONLY : nbnd, et, wg
  USE lsda_mod,      ONLY : nspin
  USE scf,           ONLY : v, vrs, vltot, kedtau, rho
  USE fft_base,      ONLY : dfftp
  USE gvecs,         ONLY : doublegrid
  USE gvect,         ONLY : ecutrho, ngm, g, gg, eigts1, eigts2, eigts3
  USE klist,         ONLY : degauss, ngauss, nks, lgauss, wk, two_fermi_energies, ltetra
  USE ions_base,     ONLY : nat, nsp, ityp, tau
  USE noncollin_module,  ONLY : noncolin
  USE constants,     ONLY : degspin, pi
  USE mp_pools,      ONLY : inter_pool_comm
  USE mp,            ONLY : mp_max, mp_min
  USE dfunct,        ONLY : newd
  USE pwcom,         ONLY : ef
  USE constants,     ONLY : rytoev
  USE cell_base,     ONLY : alat, at, bg, omega
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE gipaw_module
  USE ions_base, only: tau, ityp

  implicit none
  integer :: ik, ibnd
  real(dp) :: emin, emax, xmax, small, fac, target
  
  write(*,*) 'VEB', iverbosity
  
  ! initialize pseudopotentials
  call init_us_1(nat, ityp, omega, ngm, g, gg, intra_bgrp_comm)

! setup GIPAW operators
  call gipaw_setup_integrals


END SUBROUTINE gipaw_setup

!-----------------------------------------------------------------------
SUBROUTINE gipaw_setup_integrals
  !-----------------------------------------------------------------------
  !
  ! ... Setup the GIPAW integrals: NMR core contribution, diamagnetic 'E'
  ! ... and paramagnetic 'F' terms, relativistic mass corrections
  !
  !
  USE gipaw_module
  USE kinds,         ONLY : dp
  USE ions_base,     ONLY : ntyp => nsp, atm
  USE atom,          ONLY : rgrid
  USE paw_gipaw,     ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp, set_paw_upf
  USE uspp_param,    ONLY : upf
  USE io_global,     ONLY : stdout
  USE wvfct,         ONLY : nbnd, npwx

  implicit none

  real(dp), allocatable :: work(:), kinetic_aephi(:), kinetic_psphi(:)
  real(dp), allocatable :: aephi_dvloc_dr(:), psphi_dvloc_dr(:)
  integer :: nt, il, il1, il2, l1, l2, j, kkpsi, nrc
  integer :: core_orb
  real(dp) :: integral, occupation
  

  ! initialize data, also for the case that no GIPAW is present
  if ( .not. allocated(paw_recon) ) allocate(paw_recon(ntyp))
  paw_recon(:)%gipaw_data_in_upf_file = .false.
  paw_recon(:)%paw_nbeta = 0
  paw_recon(:)%paw_nh = 0
  paw_recon(:)%paw_kkbeta = 0
  paw_recon(:)%gipaw_ncore_orbital = 0
  paw_recon(:)%vloc_present = .false.
  do nt = 1, ntyp
    paw_recon(nt)%paw_lll(:) = 0
  end do
  
  ! setup GIPAW projectors
  do nt = 1, ntyp
      call set_paw_upf(nt, upf(nt))
!     call read_recon(file_reconstruction(nt), nt, paw_recon(nt))
  enddo

   do nt = 1, ntyp
     do il = 1, paw_recon(nt)%paw_nbeta
       if ( paw_recon(nt)%psphi(il)%label%rc < -0.99d0 ) then
        rc(nt,paw_recon(nt)%psphi(il)%label%l) = 1.6d0
        rc(nt,paw_recon(nt)%aephi(il)%label%l) = 1.6d0
        paw_recon(nt)%psphi(il)%label%rc = rc(nt,paw_recon(nt)%psphi(il)%label%l)
        paw_recon(nt)%aephi(il)%label%rc = rc(nt,paw_recon(nt)%aephi(il)%label%l)
       else
        rc(nt,paw_recon(nt)%psphi(il)%label%l) = paw_recon(nt)%psphi(il)%label%rc
        rc(nt,paw_recon(nt)%aephi(il)%label%l) = paw_recon(nt)%aephi(il)%label%rc
       endif
     enddo
    enddo

    call init_gipaw_1()

    ! allocate GIPAW projectors
  allocate ( paw_vkb(npwx,paw_nkb) )
  allocate ( paw_becp(paw_nkb,nbnd) )
  allocate ( paw_becp2(paw_nkb,nbnd) )
  allocate ( paw_becp3(paw_nkb,nbnd) )

  ! allocate GIPAW integrals
  allocate ( radial_integral_diamagnetic(nbrx,nbrx,ntypx) )
  allocate ( radial_integral_paramagnetic(nbrx,nbrx,ntypx) )
  allocate ( radial_integral_diamagnetic_so(nbrx,nbrx,ntypx) )
  allocate ( radial_integral_paramagnetic_so(nbrx,nbrx,ntypx) )
  allocate ( radial_integral_rmc(nbrx,nbrx,ntypx) )
  radial_integral_diamagnetic = 0.d0
  radial_integral_paramagnetic = 0.d0
  radial_integral_diamagnetic_so = 0.d0
  radial_integral_paramagnetic_so = 0.d0
  radial_integral_rmc = 0.d0


     ! calculate GIPAW integrals
  do nt = 1, ntyp

    do il1 = 1, paw_recon(nt)%paw_nbeta
      l1 = paw_recon(nt)%psphi(il1)%label%l

      kkpsi = paw_recon(nt)%aephi(il1)%kkpsi
      nrc = paw_recon(nt)%psphi(il1)%label%nrc

      allocate ( work(kkpsi) )

      do il2 = 1, paw_recon(nt)%paw_nbeta
        l2 = paw_recon(nt)%psphi(il2)%label%l

        if ( l1 /= l2 ) cycle

        ! NMR diamagnetic: (1/r)
        do j = 1, nrc
           work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * paw_recon(nt)%aephi(il2)%psi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     / rgrid(nt)%r(j)
        enddo
        call simpson( nrc, work, rgrid(nt)%rab(:nrc), radial_integral_diamagnetic(il1,il2,nt) )

        
        ! NMR paramagnetic: (1/r^3)
        do j = 1, nrc
           work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * paw_recon(nt)%aephi(il2)%psi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     / rgrid(nt)%r(j)
        enddo
        call simpson( nrc, work, rgrid(nt)%rab(:nrc), radial_integral_diamagnetic(il1,il2,nt) )

        
        ! NMR paramagnetic: (1/r^3)
        do j = 1, nrc
           work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * paw_recon(nt)%aephi(il2)%psi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     / rgrid(nt)%r(j) ** 3
        enddo
        call simpson( nrc, work, rgrid(nt)%rab(:nrc), radial_integral_paramagnetic(il1,il2,nt) )

        ! calculate the radial integral only if the radial potential is present
        if ( .not. paw_recon(nt)%vloc_present ) cycle

        ! g-tensor relativistic mass correction: (-nabla^2)
        allocate ( kinetic_aephi ( kkpsi ), kinetic_psphi ( kkpsi ) )
        call radial_kinetic_energy (nrc, l2, rgrid(nt)%r(:nrc), paw_recon(nt)%aephi(il2)%psi(:nrc), &
                                    kinetic_aephi(:nrc))
        call radial_kinetic_energy (nrc, l2, rgrid(nt)%r(:nrc), paw_recon(nt)%psphi(il2)%psi(:nrc), &
                                    kinetic_psphi(:nrc))
        do j = 1, nrc
           work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * kinetic_aephi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * kinetic_psphi(j) )
        enddo
        deallocate ( kinetic_aephi, kinetic_psphi )
        call simpson ( nrc, work, rgrid(nt)%rab(:nrc), radial_integral_rmc(il1,il2,nt) )

        ! calculate dV/dr
        allocate ( aephi_dvloc_dr(nrc), psphi_dvloc_dr(nrc) )
        call radial_derivative (nrc, rgrid(nt)%r(:nrc), paw_recon(nt)%gipaw_ae_vloc(:nrc), aephi_dvloc_dr(:nrc))
        call radial_derivative (nrc, rgrid(nt)%r(:nrc), paw_recon(nt)%gipaw_ps_vloc(:nrc), psphi_dvloc_dr(:nrc))
         
                ! g tensor diamagnetic: (r*dV/dr)
        do j = 1, nrc
           work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * aephi_dvloc_dr(j) * paw_recon(nt)%aephi(il2)%psi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * psphi_dvloc_dr(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     * rgrid(nt)%r(j)
        enddo
        call simpson( nrc, work, rgrid(nt)%rab(:nrc), radial_integral_diamagnetic_so(il1,il2,nt) )

        ! g tensor paramagnetic: (1/r*dV/dr)
        do j = 1, nrc
          work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * aephi_dvloc_dr(j) * paw_recon(nt)%aephi(il2)%psi(j) &
                    - paw_recon(nt)%psphi(il1)%psi(j) * psphi_dvloc_dr(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                    / rgrid(nt)%r(j)
        enddo
        call simpson( nrc,work,rgrid(nt)%rab(:nrc), radial_integral_paramagnetic_so(il1,il2,nt) )


        deallocate ( aephi_dvloc_dr, psphi_dvloc_dr )


        enddo  ! l2

      deallocate ( work )

     enddo  ! l1
  enddo  ! nt

   ! Compute the shift due to core orbitals (purely diamagnetic)
  do nt = 1, ntyp
    if ( paw_recon(nt)%gipaw_ncore_orbital == 0 ) cycle

    allocate ( work(rgrid(nt)%mesh) )
    nmr_shift_core(nt) = 0.0

    do core_orb = 1, paw_recon(nt)%gipaw_ncore_orbital
      do j = 1, size(work)
        work(j) = paw_recon(nt)%gipaw_core_orbital(j,core_orb) ** 2 / rgrid(nt)%r(j)
      end do
      call simpson( size(work), work, rgrid(nt)%rab(:), integral )
      occupation = 2 * ( 2 * paw_recon(nt)%gipaw_core_orbital_l(core_orb) + 1 )
      nmr_shift_core(nt) = nmr_shift_core(nt) + occupation * integral
    enddo
    deallocate ( work )

    nmr_shift_core(nt) = nmr_shift_core(nt) * 17.75045395 * 1e-6
  enddo

 
        
  ! print integrals
!  if (iverbosity > 10) then
    write(stdout,'(5X,''GIPAW integrals: -------------------------------------------'')')
    write(stdout,'(5X,''Atom  i/j   nmr_para    nmr_dia     epr_rmc    epr_para     epr_dia'')')
    do nt = 1, ntyp
      do il1 = 1, paw_recon(nt)%paw_nbeta
        l1 = paw_recon(nt)%psphi(il1)%label%l
        do il2 = 1, paw_recon(nt)%paw_nbeta
          l2 = paw_recon(nt)%psphi(il2)%label%l

          if (l1 /= l2) cycle
          if (il1 < il2) cycle

          write(stdout,1000) atm(nt), il1, il2, &
             radial_integral_paramagnetic(il1,il2,nt), &
             radial_integral_diamagnetic(il1,il2,nt), &
             radial_integral_rmc(il1,il2,nt), &
             radial_integral_paramagnetic_so(il1,il2,nt), &
             radial_integral_diamagnetic_so(il1,il2,nt)

        enddo
      enddo
    enddo
    write(stdout,'(5X,''------------------------------------------------------------'')')
    write(stdout,*)
!  endif
1000 format(5X,A5,1X,I1,1X,I1,2X,5(E10.4,2x))

  
  END SUBROUTINE gipaw_setup_integrals 

