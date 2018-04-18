!------------------------------------------------------------------------------!
!> Variables that have to do with HELENA quantities.
!------------------------------------------------------------------------------!
module HELENA_vars
#include <PB3D_macros.h>
    use num_vars, only: pi
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, max_str_ln
    use grid_vars, only: grid_type
    use X_vars, only: X_1_type, X_2_type
    
    implicit none
    private
    public dealloc_HEL, &
        &pres_H, q_saf_H, rot_t_H, flux_p_H, flux_t_H, nchi, chi_H, ias, &
        &RBphi_H, R_H, Z_H, h_H_11, h_H_12, h_H_33, RMtoG_H, BMtoG_H
    
    ! global variables
    real(dp) :: RMtoG_H                                                         !< R_geo/R_mag
    real(dp) :: BMtoG_H                                                         !< B_geo/B_mag
    real(dp), allocatable :: chi_H(:)                                           !< poloidal angle
    real(dp), allocatable :: flux_p_H(:,:)                                      !< poloidal flux
    real(dp), allocatable :: flux_t_H(:,:)                                      !< toroidal flux
    real(dp), allocatable :: pres_H(:,:)                                        !< pressure profile
    real(dp), allocatable :: q_saf_H(:,:)                                       !< safety factor
    real(dp), allocatable :: rot_t_H(:,:)                                       !< rotational transform
    real(dp), allocatable :: RBphi_H(:,:)                                       !< \f$R B_\phi (= F) \f$
    real(dp), allocatable :: h_H_11(:,:)                                        !< upper metric factor \f$h_{11}\f$ (\c gem11)
    real(dp), allocatable :: h_H_12(:,:)                                        !< upper metric factor \f$h_{12}\f$ (\c gem12)
    real(dp), allocatable :: h_H_33(:,:)                                        !< upper metric factor \f$h_{32}\f$ (1 / \c gem12)
    real(dp), allocatable :: R_H(:,:)                                           !< major radius \f$R\f$ (xout)
    real(dp), allocatable :: Z_H(:,:)                                           !< height \f$Z\f$ (yout)
    integer :: nchi                                                             !< nr. of poloidal points
    integer :: ias                                                              !< 0 if top-bottom symmetric, 1 if not

contains
    !> Deallocates HELENA quantities that are not used any more.
    subroutine dealloc_HEL
#if ldebug
        use num_vars, only: rank, print_mem_usage
        
        ! local variables
        integer :: mem_diff                                                     ! difference in memory
        
        ! memory usage before deallocation
        if (print_mem_usage) mem_diff = get_mem_usage()
#endif
        
        ! deallocate
        deallocate(chi_H)
        deallocate(flux_p_H)
        deallocate(flux_t_H)
        deallocate(pres_H)
        deallocate(q_saf_H)
        deallocate(rot_t_H)
        deallocate(RBphi_H)
        deallocate(R_H)
        deallocate(Z_H)
        deallocate(h_H_11)
        deallocate(h_H_12)
        deallocate(h_H_33)
        
#if ldebug
        ! memory usage difference after deallocation
        if (print_mem_usage) then
            mem_diff = mem_diff - get_mem_usage()
            call writo('Rank '//trim(i2str(rank))//' liberated '//&
                &trim(i2str(mem_diff))//'kB deallocating HELENA')
        end if
#endif
    end subroutine dealloc_HEL
end module HELENA_vars
