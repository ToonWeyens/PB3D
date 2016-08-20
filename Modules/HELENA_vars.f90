!------------------------------------------------------------------------------!
!   Variables that have to do with HELENA quantities                           !
!------------------------------------------------------------------------------!
module HELENA_vars
#include <PB3D_macros.h>
    use num_vars, only: pi
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, max_str_ln
    use grid_vars, only: grid_type, disc_type
    use X_vars, only: X_1_type, X_2_type
    
    implicit none
    private
    public dealloc_HEL, &
        &pres_H, qs_H, flux_p_H, flux_t_H, Dflux_p_H, Dflux_t_H, nchi, chi_H, &
        &ias, RBphi_H, R_H, Z_H, h_H_11, h_H_12, h_H_33
    
    ! global variables
    real(dp), allocatable :: chi_H(:)                                           ! poloidal angle
    real(dp), allocatable :: flux_p_H(:)                                        ! poloidal flux
    real(dp), allocatable :: flux_t_H(:)                                        ! toroidal flux
    real(dp), allocatable :: Dflux_p_H(:)                                       ! normal derivative of poloidal flux
    real(dp), allocatable :: Dflux_t_H(:)                                       ! normal derivative of toroidal flux
    real(dp), allocatable :: pres_H(:)                                          ! pressure profile
    real(dp), allocatable :: qs_H(:)                                            ! safety factor
    real(dp), allocatable :: RBphi_H(:)                                         ! R B_phi (= F)
    real(dp), allocatable :: h_H_11(:,:)                                        ! adapted upper metric factor 11 (gem11)
    real(dp), allocatable :: h_H_12(:,:)                                        ! adapted upper metric factor 12 (gem12)
    real(dp), allocatable :: h_H_33(:,:)                                        ! adapted upper metric factor 33 (1/gem33)
    real(dp), allocatable :: R_H(:,:)                                           ! major radius R (xout)
    real(dp), allocatable :: Z_H(:,:)                                           ! height Z (yout)
    integer :: nchi                                                             ! nr. of poloidal points (nchi)
    integer :: ias                                                              ! 0 if top-bottom symmetric, 1 if not

contains
    ! deallocates HELENA quantities that are not used any more
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
        deallocate(Dflux_p_H)
        deallocate(Dflux_t_H)
        deallocate(pres_H)
        deallocate(qs_H)
        deallocate(RBphi_H)
        deallocate(R_H)
        deallocate(Z_H)
        
#if ldebug
        deallocate(h_H_11)
        deallocate(h_H_12)
        deallocate(h_H_33)
        
        ! memory usage difference after deallocation
        if (print_mem_usage) then
            mem_diff = mem_diff - get_mem_usage()
            call writo('Rank '//trim(i2str(rank))//' liberated '//&
                &trim(i2str(mem_diff))//'kB deallocating HELENA')
        end if
#endif
    end subroutine dealloc_HEL
end module HELENA_vars
