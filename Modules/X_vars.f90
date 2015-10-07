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
        &X_1_type, X_2_type, &
        &min_n_r_X, min_n_X, max_n_X, min_m_X, max_m_X, min_r_X, max_r_X
    
    ! global variables
    ! (These  variables  should  be   used   only  until  the  grids  have  been
    ! established. They are put here for lack of a better place.)
    real(dp) :: min_r_X, max_r_X                                                ! min. and max. normal range for pert. (either pol. or tor., depending on use_pol_flux_F)
    integer :: min_n_r_X                                                        ! min. of n_r_X (e.g. first value in Richardson loop)
    integer :: min_n_X                                                          ! lowest poloidal mode number m_X
    integer :: max_n_X                                                          ! highest poloidal mode number m_X
    integer :: min_m_X                                                          ! lowest poloidal mode number m_X
    integer :: max_m_X                                                          ! highest poloidal mode number m_X
    
    ! vectorial perturbation type with arrays of the form:
    !   - (angle_1,angle_2,r,n_mod)         for U_X_i, DU_X_i
    ! where it is refered to the discussion  of the grid type for an explanation
    ! of the angles angle_1 and angle_2.
    type :: X_1_type
        integer :: n_mod                                                        ! size of n and m (nr. of modes)
        integer, allocatable :: n(:)                                            ! vector of poloidal mode numbers
        integer, allocatable :: m(:)                                            ! vector of poloidal mode numbers
        complex(dp), allocatable :: U_0(:,:,:,:)                                ! U_m(X_m) = [ U_m^0 + U_m^1 i/n d/dx] (X_m)
        complex(dp), allocatable :: U_1(:,:,:,:)                                ! U_m(X_m) = [ U_m^0 + U_m^1 i/n d/dx] (X_m)
        complex(dp), allocatable :: DU_0(:,:,:,:)                               ! d(U_m(X_m))/dtheta = [ DU_m^0 + DU_m^1 i/n d/dx] (X_m)
        complex(dp), allocatable :: DU_1(:,:,:,:)                               ! d(U_m(X_m))/dtheta = [ DU_m^0 + DU_m^1 i/n d/dx] (X_m)
    end type
    
    ! tensorial perturbation type with arrays of the form:
    !   - (angle_1,angle_2,r,n_mod^2)       for PVi, KVi, J_exp_ang_par_F
    !   - (n_mod^2,angle_2,r,3)             for KV_int, PV_int
    ! where it is refered to the discussion  of the grid type for an explanation
    ! of the angles angle_1 and angle_2.
    type :: X_2_type
        integer :: n_mod                                                        ! size of n and m (nr. of modes)
        integer, allocatable :: n(:,:)                                          ! vector of poloidal mode numbers
        integer, allocatable :: m(:,:)                                          ! vector of poloidal mode numbers
        complex(dp), allocatable :: PV_0(:,:,:,:)                               ! ~PV^0 coefficient
        complex(dp), allocatable :: PV_1(:,:,:,:)                               ! ~PV^1 coefficient
        complex(dp), allocatable :: PV_2(:,:,:,:)                               ! ~PV^2 coefficient
        complex(dp), allocatable :: KV_0(:,:,:,:)                               ! ~KV^0 coefficient
        complex(dp), allocatable :: KV_1(:,:,:,:)                               ! ~KV^1 coefficient
        complex(dp), allocatable :: KV_2(:,:,:,:)                               ! ~KV^2 coefficient
        complex(dp), allocatable :: PV_int_0(:,:,:)                             ! <~PV^0 e^i(k-m)ang_par_F> coefficient
        complex(dp), allocatable :: PV_int_1(:,:,:)                             ! <~PV^1 e^i(k-m)ang_par_F> coefficient
        complex(dp), allocatable :: PV_int_2(:,:,:)                             ! <~PV^2 e^i(k-m)ang_par_F> coefficient
        complex(dp), allocatable :: KV_int_0(:,:,:)                             ! <~KV^0 e^i(k-m)ang_par_F> coefficient
        complex(dp), allocatable :: KV_int_1(:,:,:)                             ! <~KV^1 e^i(k-m)ang_par_F> coefficient
        complex(dp), allocatable :: KV_int_2(:,:,:)                             ! <~KV^2 e^i(k-m)ang_par_F> coefficient
        complex(dp), allocatable :: vac_res(:,:)                                ! vacuum response
    end type
    
    ! interfaces
    interface create_X
        module procedure create_X_1, create_X_2
    end interface
    interface dealloc_X
        module procedure dealloc_X_1, dealloc_X_2
    end interface
    
contains
    ! Create  a  vectorial  or  tensorial perturbation  type  and  allocate  the
    ! variables, the number of modes, as well as n and m
    ! Optionally, the  secondary mode  numbers can be  specified (m  if poloidal
    ! flux is used and n if toroidal  flux). By default, they are taken from the
    ! global X_vars variables.
    integer function create_X_1(grid_eq,X,lim_sec_X) result(ierr)               ! vectorial version
        use num_vars, only: use_pol_flux_F
        
        character(*), parameter :: rout_name = 'create_X_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(X_1_type), intent(inout) :: X                                      ! vectorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2)                           ! limits of m_X (pol. flux) or n_X (tor. flux)
        
        ! local variables
        integer :: id                                                           ! counters
        integer :: loc_n_r                                                      ! local nr. of normal points
        integer :: n_par, n_geo                                                 ! tot. nr. of angular points in parallel and geodesic direction
        integer :: lim_n_X(2)                                                   ! min. and max. of n_X
        integer :: lim_m_X(2)                                                   ! min. and max. of m_X
        
        ! initialize ierr
        ierr = 0
        
        ! set lim_n_X and lim_m_X
        if (use_pol_flux_F) then
            lim_n_X = [min_n_X,max_n_X]
            if (present(lim_sec_X)) then
                lim_m_X = lim_sec_X
            else
                lim_m_X = [min_m_X,max_m_X]
            end if
        else
            if (present(lim_sec_X)) then
                lim_n_X = lim_sec_X
            else
                lim_n_X = [min_n_X,max_n_X]
            end if
            lim_m_X = [min_m_X,max_m_X]
        end if
        
        ! set local variables
        loc_n_r = grid_eq%loc_n_r
        n_par = grid_eq%n(1)
        n_geo = grid_eq%n(2)
        
        ! set n_mod
        X%n_mod = (lim_m_X(2)-lim_m_X(1)+1)*(lim_n_X(2)-lim_n_X(1)+1)
        
        ! set n and m
        allocate(X%n(X%n_mod),X%m(X%n_mod))
        if (use_pol_flux_F) then
            X%n = lim_n_X(1)
            X%m = [(id, id = lim_m_X(1), lim_m_X(2))]
        else
            X%n = [(id, id = lim_n_X(1), lim_n_X(2))]
            X%m = lim_m_X(1)
        end if
        
        ! allocate U_0
        allocate(X%U_0(n_par,n_geo,loc_n_r,X%n_mod))
        
        ! allocate U_1
        allocate(X%U_1(n_par,n_geo,loc_n_r,X%n_mod))
        
        ! allocate DU_0
        allocate(X%DU_0(n_par,n_geo,loc_n_r,X%n_mod))
        
        ! allocate DU_1
        allocate(X%DU_1(n_par,n_geo,loc_n_r,X%n_mod))
    end function create_X_1
    integer function create_X_2(grid_eq,X,lim_sec_X) result(ierr)               ! tensorial version
        use num_vars, only: use_pol_flux_F
        
        character(*), parameter :: rout_name = 'create_X_2'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(X_2_type), intent(inout) :: X                                      ! tensorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol. flux) or n_X (tor. flux) for both dimensions
        
        ! local variables
        integer :: id, jd                                                       ! counters
        integer :: loc_n_r                                                      ! local nr. of normal points
        integer :: n_par, n_geo                                                 ! tot. nr. of angular points in parallel and geodesic direction
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: lim_n_X(2,2)                                                 ! min. and max. of n_X for both dimensions
        integer :: lim_m_X(2,2)                                                 ! min. and max. of m_X for both dimensions
        integer :: n_mod_loc(2)                                                 ! local n_mod for both dimensions
        integer :: nn_mod                                                       ! either n_mod^2, n_mod*(n_mod+1)/2 or 0
        
        ! initialize ierr
        ierr = 0
        
        ! set local variables
        loc_n_r = grid_eq%loc_n_r
        n_par = grid_eq%n(1)
        n_geo = grid_eq%n(2)
        
        ! set lim_n_X and lim_m_X
        if (use_pol_flux_F) then
            do jd = 1,2
                lim_n_X(:,jd) = [min_n_X,max_n_X]
            end do
            if (present(lim_sec_X)) then
                lim_m_X = lim_sec_X
            else
                do jd = 1,2
                    lim_m_X(:,jd) = [min_m_X,max_m_X]
                end do
            end if
        else
            if (present(lim_sec_X)) then
                lim_n_X = lim_sec_X
            else
                do jd = 1,2
                    lim_n_X(:,jd) = [min_n_X,max_n_X]
                end do
            end if
            do jd = 1,2
                lim_m_X(:,jd) = [min_m_X,max_m_X]
            end do
        end if
        
        ! set local variables
        loc_n_r = grid_eq%loc_n_r
        n_par = grid_eq%n(1)
        n_geo = grid_eq%n(2)
        
        ! set n_mod
        n_mod_loc = (lim_m_X(2,:)-lim_m_X(1,:)+1)*(lim_n_X(2,:)-lim_n_X(1,:)+1)
        if (n_mod_loc(1).eq.n_mod_loc(2)) then
            X%n_mod = n_mod_loc(1)
        else
            ierr = 1
            err_msg = 'Nr. of modes has to be equal for both dimensions'
            CHCKERR(err_msg)
        end if
        
        ! set nnmod for symmetric quantities
        if (use_pol_flux_F) then
            if (lim_m_X(1,1).gt.lim_m_X(1,2)) then                              ! below diagonal: dense
                nn_mod = X%n_mod**2
            else if (lim_m_X(1,1).eq.lim_m_X(1,2)) then                         ! on diagonal: partly filled
                nn_mod = X%n_mod*(X%n_mod+1)/2
            else                                                                ! above diagonal: not filled
                nn_mod = 0
            end if
        else
            if (lim_n_X(1,1).gt.lim_n_X(1,2)) then                              ! below diagonal: dense
                nn_mod = X%n_mod**2
            else if (lim_n_X(1,1).eq.lim_n_X(1,2)) then                         ! on diagonal: party filled
                nn_mod = X%n_mod*(X%n_mod+1)/2
            else                                                                ! above diagonal :: not filled
                nn_mod = 0
            end if
        end if
        
        ! set n and m
        allocate(X%n(X%n_mod,2),X%m(X%n_mod,2))
        if (use_pol_flux_F) then
            X%n = lim_n_X(1,1)
            do jd = 1,2
                X%m(:,jd) = [(id, id = lim_m_X(1,jd), lim_m_X(2,jd))]
            end do
        else
            do jd = 1,2
                X%n(:,jd) = [(id, id = lim_n_X(1,jd), lim_n_X(2,jd))]
            end do
            X%m = lim_m_X(1,1)
        end if
        
        ! allocate PV_i
        allocate(X%PV_0(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        allocate(X%PV_1(n_par,n_geo,loc_n_r,X%n_mod**2))                        ! not symmetric
        allocate(X%PV_2(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        
        ! allocate KV_i
        allocate(X%KV_0(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        allocate(X%KV_1(n_par,n_geo,loc_n_r,X%n_mod**2))                        ! not symmetric
        allocate(X%KV_2(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        
        ! allocate PV_int_i
        allocate(X%PV_int_0(nn_mod,n_geo,loc_n_r))                              ! symmetric
        allocate(X%PV_int_1(X%n_mod**2,n_geo,loc_n_r))                          ! not symmetric
        allocate(X%PV_int_2(nn_mod,n_geo,loc_n_r))                              ! symmetric
        
        ! allocate KV_int_i
        allocate(X%KV_int_0(nn_mod,n_geo,loc_n_r))                              ! symmetric
        allocate(X%KV_int_1(X%n_mod**2,n_geo,loc_n_r))                          ! not symmetric
        allocate(X%KV_int_2(nn_mod,n_geo,loc_n_r))                              ! symmetric
        
        ! allocate vacuum response
        allocate(X%vac_res(X%n_mod,X%n_mod))
    end function create_X_2
    
    ! deallocates perturbation variables
    subroutine dealloc_X_1(X)                                                   ! vectorial version
        ! input / output
        type(X_1_type), intent(out) :: X                                        ! perturbation variables to be deallocated
    end subroutine dealloc_X_1
    subroutine dealloc_X_2(X)                                                   ! tensorial version
        ! input / output
        type(X_2_type), intent(out) :: X                                        ! perturbation variables to be deallocated
    end subroutine dealloc_X_2
end module X_vars
