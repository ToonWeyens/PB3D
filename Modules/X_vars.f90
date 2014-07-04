!------------------------------------------------------------------------------!
!   Holds variables pertaining to the perturbation quantities
!------------------------------------------------------------------------------!
module X_vars
    use num_vars, only: dp
    use output_ops, only: lvl_ud, writo

    implicit none
    private
    public calc_rho, &
        &rho, n_x, m_X, n_r_X

    ! global variables
    real(dp), allocatable :: rho(:,:)                                           ! density
    integer :: n_X                                                              ! toroidal mode number
    integer, allocatable :: m_X(:)                                              ! vector of poloidal mode numbers
    integer :: n_r_X                                                            ! number of normal points for perturbations
    
contains
    ! calculate rho from user input
    ! TO BE IMPLEMENTED. TEMPORARILY SET TO 1 EVERYWHERE !!!
    subroutine calc_rho(n_par,n_r)
        ! input variables
        integer, intent(in) :: n_par, n_r                                       ! dimensions of the size to be allocated for the rho matrix
        
        call writo('calc_rho NOT YET IMPLEMENTED!!!')
        
        ! allocate rho
        if (allocated(rho)) deallocate(rho)
        allocate(rho(n_par,n_r))
        
        ! TEMPORAL REPLACEMENT !!!
        rho = 1
    end subroutine
end module
