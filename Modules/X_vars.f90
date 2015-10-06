!------------------------------------------------------------------------------!
!   Variables pertaining to the perturbation quantities                        !
!------------------------------------------------------------------------------!
module X_vars
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, max_str_ln, iu
    use grid_vars, only: grid_type
    
    implicit none
    
    private
    public dealloc_X, create_X, &
        &X_type, &
        &min_n_r_X, min_n_X, max_n_X, min_m_X, max_m_X, min_r_X, &
        &max_r_X
    
    ! global variables
    ! (These  variables  should  be   used   only  until  the  grids  have  been
    ! established. They are put here for lack of a better place.)
    real(dp) :: min_r_X, max_r_X                                                ! min. and max. normal range for pert. (either pol. or tor., depending on use_pol_flux_F)
    integer :: min_n_r_X                                                        ! min. of n_r_X (e.g. first value in Richardson loop)
    integer :: min_n_X                                                          ! lowest poloidal mode number m_X
    integer :: max_n_X                                                          ! highest poloidal mode number m_X
    integer :: min_m_X                                                          ! lowest poloidal mode number m_X
    integer :: max_m_X                                                          ! highest poloidal mode number m_X
    
    ! perturbation type
    ! The arrays here are of the form:
    !   - (angle_1,angle_2,r,n_mod)         for U_X_i, DU_X_i
    !   - (angle_1,angle_2,r,nn_mod)        for PVi, KVi, J_exp_ang_par_F
    !   - (nn_mod,angle_2,r,3)              for KV_int, PV_int
    ! where it is refered to the discussion  of the grid type for an explanation
    ! of the angles angle_1 and angle_2.
    ! nn_mod refers to the number of nonzero entries, as discussed in create_X.
    type :: X_type
        integer :: n_mod                                                        ! size of n and m (nr. of modes)
        integer :: min_n                                                        ! lowest poloidal mode number m
        integer :: max_n                                                        ! highest poloidal mode number m
        integer :: min_m                                                        ! lowest poloidal mode number m
        integer :: max_m                                                        ! highest poloidal mode number m
        integer, allocatable :: n(:)                                            ! vector of poloidal mode numbers
        integer, allocatable :: m(:)                                            ! vector of poloidal mode numbers
        complex(dp), allocatable :: vec(:,:,:)                                  ! Eigenvector solution
        complex(dp), allocatable :: val(:)                                      ! Eigenvalue solution
        complex(dp), allocatable :: U_0(:,:,:,:)                                ! U_m(X_m) = [ U_m^0 + U_m^1 i/n d/dx] (X_m)
        complex(dp), allocatable :: U_1(:,:,:,:)                                ! U_m(X_m) = [ U_m^0 + U_m^1 i/n d/dx] (X_m)
        complex(dp), allocatable :: DU_0(:,:,:,:)                               ! d(U_m(X_m))/dtheta = [ DU_m^0 + DU_m^1 i/n d/dx] (X_m)
        complex(dp), allocatable :: DU_1(:,:,:,:)                               ! d(U_m(X_m))/dtheta = [ DU_m^0 + DU_m^1 i/n d/dx] (X_m)
        complex(dp), allocatable :: PV_0(:,:,:,:)                               ! ~PV^0 coefficient
        complex(dp), allocatable :: PV_1(:,:,:,:)                               ! ~PV^1 coefficient
        complex(dp), allocatable :: PV_2(:,:,:,:)                               ! ~PV^2 coefficient
        complex(dp), allocatable :: KV_0(:,:,:,:)                               ! ~KV^0 coefficient
        complex(dp), allocatable :: KV_1(:,:,:,:)                               ! ~KV^1 coefficient
        complex(dp), allocatable :: KV_2(:,:,:,:)                               ! ~KV^2 coefficient
        complex(dp), allocatable :: J_exp_ang_par_F(:,:,:,:)                    ! J_F exp(i (k-m) ang_par_F)
        complex(dp), allocatable :: PV_int_0(:,:,:)                             ! <~PV^0 e^i(k-m)ang_par_F> coefficient
        complex(dp), allocatable :: PV_int_1(:,:,:)                             ! <~PV^1 e^i(k-m)ang_par_F> coefficient
        complex(dp), allocatable :: PV_int_2(:,:,:)                             ! <~PV^2 e^i(k-m)ang_par_F> coefficient
        complex(dp), allocatable :: KV_int_0(:,:,:)                             ! <~KV^0 e^i(k-m)ang_par_F> coefficient
        complex(dp), allocatable :: KV_int_1(:,:,:)                             ! <~KV^1 e^i(k-m)ang_par_F> coefficient
        complex(dp), allocatable :: KV_int_2(:,:,:)                             ! <~KV^2 e^i(k-m)ang_par_F> coefficient
        complex(dp), allocatable :: vac_res(:,:)                                ! vacuum response
    end type
    
contains
    ! initialize the variable m and check and/or plot it
    ! Optionally, the limits  for n_X and/or m_X can be  specified. If not, they
    ! are taken from the global variables in this module.
    ! Note: intent(out) automatically deallocates the variable
    subroutine create_X(grid,X,lim_n_X,lim_m_X)
        use num_vars, only: use_pol_flux_F
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid
        type(X_type), intent(out) :: X                                          ! perturbation variables
        integer, intent(in), optional :: lim_n_X(2)                             ! min. and max. of n_X
        integer, intent(in), optional :: lim_m_X(2)                             ! min. and max. of m_X
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: loc_n_r                                                      ! local nr. of normal points
        integer :: n_par, n_geo                                                 ! tot. nr. of angular points in parallel and geodesic direction
        integer :: nn_mod_1, nn_mod_2                                           ! number of indices for a quantity that is symmetric or not
        integer :: min_n_X_loc, max_n_X_loc, min_m_X_loc, max_m_X_loc           ! local versions of min_n_X, max_n_X, min_m_X, max_m_X
        
        ! set local variables
        loc_n_r = grid%loc_n_r
        n_par = grid%n(1)
        n_geo = grid%n(2)
        min_n_X_loc = min_n_X
        if (present(lim_n_X)) min_n_X_loc = lim_n_X(1)
        max_n_X_loc = max_n_X
        if (present(lim_n_X)) max_n_X_loc = lim_n_X(2)
        min_m_X_loc = min_m_X
        if (present(lim_m_X)) min_m_X_loc = lim_m_X(1)
        max_m_X_loc = max_m_X
        if (present(lim_m_X)) max_m_X_loc = lim_m_X(2)
        
        ! set n_mod
        X%n_Mod = (max_m_X_loc-min_m_X_loc+1)*(max_n_X_loc-min_n_X_loc+1)
        
        ! set nn_mod_1 and nn_mod_2
        nn_mod_1 = X%n_mod**2
        nn_mod_2 = X%n_mod*(X%n_mod+1)/2
        
        ! n and m
        allocate(X%n(X%n_mod),X%m(X%n_mod))
        if (use_pol_flux_F) then
            X%n = min_n_X_loc
            X%m = [(id, id = min_m_X_loc, max_m_X_loc)]
        else
            X%n = [(id, id = min_n_X_loc, max_n_X_loc)]
            X%m = min_m_X_loc
        end if
        
        ! U_0
        allocate(X%U_0(n_par,n_geo,loc_n_r,X%n_mod))
        
        ! U_1
        allocate(X%U_1(n_par,n_geo,loc_n_r,X%n_mod))
        
        ! DU_0
        allocate(X%DU_0(n_par,n_geo,loc_n_r,X%n_mod))
        
        ! DU_1
        allocate(X%DU_1(n_par,n_geo,loc_n_r,X%n_mod))
        
        ! PV_i
        allocate(X%PV_0(n_par,n_geo,loc_n_r,nn_mod_2))                          ! symmetric
        allocate(X%PV_1(n_par,n_geo,loc_n_r,nn_mod_1))                          ! not symmetric
        allocate(X%PV_2(n_par,n_geo,loc_n_r,nn_mod_2))                          ! symmetric
        
        ! KV_i
        allocate(X%KV_0(n_par,n_geo,loc_n_r,nn_mod_2))                          ! symmetric
        allocate(X%KV_1(n_par,n_geo,loc_n_r,nn_mod_1))                          ! not symmetric
        allocate(X%KV_2(n_par,n_geo,loc_n_r,nn_mod_2))                          ! symmetric
        
        ! J_exp_ang_par_F
        allocate(X%J_exp_ang_par_F(n_par,n_geo,loc_n_r,nn_mod_2))               ! symmetric
        
        ! PV_int_i
        allocate(X%PV_int_0(nn_mod_2,n_geo,loc_n_r))                            ! symmetric
        allocate(X%PV_int_1(nn_mod_1,n_geo,loc_n_r))                            ! not symmetric
        allocate(X%PV_int_2(nn_mod_2,n_geo,loc_n_r))                            ! symmetric
        
        ! KV_int_i
        allocate(X%KV_int_0(nn_mod_2,n_geo,loc_n_r))                            ! symmetric
        allocate(X%KV_int_1(nn_mod_1,n_geo,loc_n_r))                            ! not symmetric
        allocate(X%KV_int_2(nn_mod_2,n_geo,loc_n_r))                            ! symmetric
        
        ! vacuum response
        allocate(X%vac_res(X%n_mod,X%n_mod))
    end subroutine create_X
    
    ! deallocates perturbation variables
    subroutine dealloc_X(X)
        ! input / output
        type(X_type), intent(inout) :: X                                        ! perturbation variables to be deallocated
        
        ! deallocate allocatable variables
        call dealloc_X_final(X)
    contains
        ! Note: intent(out) automatically deallocates the variable
        subroutine dealloc_X_final(X)
            ! input / output
            type(X_type), intent(out) :: X                                      ! perturbation variables to be deallocated
        end subroutine dealloc_X_final
    end subroutine dealloc_X
end module
