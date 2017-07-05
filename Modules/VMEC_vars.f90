!------------------------------------------------------------------------------!
!   Variables that concern the output of VMEC                                  !
!------------------------------------------------------------------------------!
module VMEC_vars
#include <PB3D_macros.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: &
        &dp, max_str_ln, pi
    use read_wout_mod, only: &                                                  ! from LIBSTELL
        &is_asym_V => lasym, VMEC_version => version_, is_freeb_V => lfreeb, &  ! stellerator symmetry, version number, free boundary or not
        &mpol_V => mpol, ntor_V => ntor, &                                      ! mpol, ntor
        &mnmax_V => mnmax, nfp_V => nfp, &                                      ! mnmax, nfp
        &aspr_V => aspect, &                                                    ! aspect ratio
        &beta_V => betaxis, &                                                   ! beta on axis
        &gam_V => gamma, &                                                      ! gamma in adiabatic law (not important here, incompressibility)
        &rmax_surf, rmin_surf, zmax_surf                                        ! max and min values of R, Z
    
    implicit none
    private
    public dealloc_VMEC, &
        &R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, jac_V_c, jac_V_s, mnmax_V, &
        &mpol_V, ntor_V, mn_V, pres_V, rot_t_V, q_saf_V, is_asym_V, flux_t_V, &
        &flux_p_V, VMEC_version, gam_V, is_freeb_V, nfp_V, B_0_V, rmin_surf, &
        &rmax_surf, aspr_V, beta_V
#if ldebug
    public B_V_sub_s, B_V_sub_c, B_V_c, B_V_s
#endif
    
    ! local variables
    integer, allocatable :: mn_V(:,:)                                           ! m and n of modes
    real(dp) :: B_0_V                                                           ! the magnitude of B at the magnetic axis, theta = zeta = 0
    real(dp), allocatable :: flux_t_V(:,:)                                      ! toroidal flux
    real(dp), allocatable :: flux_p_V(:,:)                                      ! poloidal flux
    real(dp), allocatable :: pres_V(:,:)                                        ! pressure
    real(dp), allocatable :: rot_t_V(:,:)                                       ! rotational transform
    real(dp), allocatable :: q_saf_V(:,:)                                       ! safety factor
    real(dp), allocatable :: R_V_c(:,:,:), R_V_s(:,:,:)                         ! Coeff. of R in (co)sine series (FM) and norm. deriv.
    real(dp), allocatable :: Z_V_c(:,:,:), Z_V_s(:,:,:)                         ! Coeff. of Z in (co)sine series (FM) and norm. deriv.
    real(dp), allocatable :: L_V_c(:,:,:), L_V_s(:,:,:)                         ! Coeff. of lambda in (co)sine series (HM) and norm. deriv.
    real(dp), allocatable :: jac_V_c(:,:,:), jac_V_s(:,:,:)                     ! Jacobian in VMEC coordinates (HM and FM) and norm. deriv.
#if ldebug
    real(dp), allocatable :: B_V_sub_c(:,:,:), B_V_sub_s(:,:,:)                 ! Coeff. of B_i in (co)sine series (r,theta,phi) (FM)
    real(dp), allocatable :: B_V_c(:,:), B_V_s(:,:)                             ! Coeff. of magnitude of B (HM and FM)
#endif
contains
    ! deallocates VMEC quantities that are not used anymore
    subroutine dealloc_VMEC
#if ldebug
        use num_vars, only: rank, print_mem_usage
        
        ! local variables
        integer :: mem_diff                                                     ! difference in memory
        
        ! memory usage before deallocation
        if (print_mem_usage) mem_diff = get_mem_usage()
#endif
        
        ! deallocate
        deallocate(rot_t_V)
        deallocate(q_saf_V)
        deallocate(pres_V)
        deallocate(flux_t_V)
        deallocate(flux_p_V)
        deallocate(mn_V)
        deallocate(R_V_c)
        deallocate(R_V_s)
        deallocate(Z_V_c)
        deallocate(Z_V_s)
        deallocate(L_V_c)
        deallocate(L_V_s)
        deallocate(jac_V_c)
        deallocate(jac_V_s)
#if ldebug
        deallocate(B_V_sub_c)
        deallocate(B_V_sub_s)
        deallocate(B_V_c)
        deallocate(B_V_s)
        
        ! memory usage difference after deallocation
        if (print_mem_usage) then
            mem_diff = mem_diff - get_mem_usage()
            call writo('Rank '//trim(i2str(rank))//' liberated '//&
                &trim(i2str(mem_diff))//'kB deallocating VMEC')
        end if
#endif
    end subroutine dealloc_VMEC
end module VMEC_vars
