!------------------------------------------------------------------------------!
!   Variables pertaining to the perturbation quantities                        !
!------------------------------------------------------------------------------!
module X_vars
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, max_name_ln, iu
    use grid_vars, only: grid_type
    
    implicit none
    
    private
    public init_X_vars, create_X, dealloc_X, set_nm_X, set_nn_mod, &
        &X_1_type, X_2_type, &
        &n_mod_X, prim_X, min_sec_X, max_sec_X, min_nm_X, min_n_X, max_n_X, &
        &min_m_X, max_m_X, min_r_sol, max_r_sol, X_1_var_names, X_2_var_names, &
        &n_X, m_X, sec_X_ind
    
    ! global variables
    integer :: prim_X                                                           ! n_X (pol. flux) or m_X (tor. flux)
    integer :: min_sec_X, max_sec_X                                             ! m_X (pol. flux) or n_X (tor. flux) (only for X style 1)
    integer :: n_mod_X                                                          ! size of m_X (pol. flux) or n_X (tor. flux)
    integer :: min_nm_X = 5                                                     ! minimum for the high-n theory (arbitrary and probably too low)
    integer, allocatable :: min_n_X(:)                                          ! lowest poloidal mode number m_X, in total eq grid
    integer, allocatable :: max_n_X(:)                                          ! highest poloidal mode number m_X, in total eq grid
    integer, allocatable :: min_m_X(:)                                          ! lowest poloidal mode number m_X, in total eq grid
    integer, allocatable :: max_m_X(:)                                          ! highest poloidal mode number m_X, in total eq grid
    integer, allocatable :: n_X(:,:)                                            ! n for all modes, in total X grid
    integer, allocatable :: m_X(:,:)                                            ! m for all modes, in total X grid
    integer, allocatable :: sec_X_ind(:,:)                                      ! index of m or n for all possible modes, in total X grid
    real(dp) :: min_r_sol, max_r_sol                                            ! min. and max. normal range for pert.
    character(len=max_name_ln), allocatable :: X_1_var_names(:)                 ! internal vectorial perturbation variables names
    character(len=max_name_ln), allocatable :: X_2_var_names(:)                 ! internal tensorial perturbation variables names
    
    ! vectorial perturbation type with arrays of the form:
    !   - (angle_1,angle_2,r,n_mod)         for U_X_i, DU_X_i
    ! where it is refered to the discussion  of the grid type for an explanation
    ! of the angles angle_1 and angle_2.
    type :: X_1_type
        integer :: n_mod                                                        ! size of n and m (nr. of modes)
        integer, allocatable :: n(:,:)                                          ! vector of poloidal mode numbers
        integer, allocatable :: m(:,:)                                          ! vector of poloidal mode numbers
        complex(dp), allocatable :: U_0(:,:,:,:)                                ! U_m(X_m) = [ U_m^0 + U_m^1 i/n d/dx] (X_m)
        complex(dp), allocatable :: U_1(:,:,:,:)                                ! U_m(X_m) = [ U_m^0 + U_m^1 i/n d/dx] (X_m)
        complex(dp), allocatable :: DU_0(:,:,:,:)                               ! d(U_m(X_m))/dtheta = [ DU_m^0 + DU_m^1 i/n d/dx] (X_m)
        complex(dp), allocatable :: DU_1(:,:,:,:)                               ! d(U_m(X_m))/dtheta = [ DU_m^0 + DU_m^1 i/n d/dx] (X_m)
    end type
    
    ! tensorial perturbation type with arrays of the form:
    !   - (angle_1,angle_2,r,n_mod^2)       for PVi, KVi
    !   - (n_mod^2,angle_2,r,3)             for KV_int, PV_int
    ! where it is refered to the discussion  of the grid type for an explanation
    ! of the angles angle_1 and angle_2.
    type :: X_2_type
        integer :: n_mod(2)                                                     ! size of n and m (nr. of modes)
        integer, allocatable :: n_1(:,:)                                        ! vector of toroidal mode numbers of dimension 1
        integer, allocatable :: n_2(:,:)                                        ! vector of toroidal mode numbers of dimension 2
        integer, allocatable :: m_1(:,:)                                        ! vector of poloidal mode numbers of dimension 1
        integer, allocatable :: m_2(:,:)                                        ! vector of poloidal mode numbers of dimension 2
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
    interface set_nm_X
        module procedure set_nm_X_1, set_nm_X_2
    end interface
    interface dealloc_X
        module procedure dealloc_X_1, dealloc_X_2
    end interface
    
contains
    ! Initalizes some of the common variables
    subroutine init_X_vars()
        allocate(X_1_var_names(8))
        X_1_var_names(1) = 'RE_U_0'
        X_1_var_names(2) = 'IM_U_0'
        X_1_var_names(3) = 'RE_U_1'
        X_1_var_names(4) = 'IM_U_1'
        X_1_var_names(5) = 'RE_DU_0'
        X_1_var_names(6) = 'IM_DU_0'
        X_1_var_names(7) = 'RE_DU_1'
        X_1_var_names(8) = 'IM_DU_1'
        
        allocate(X_2_var_names(12))
        X_2_var_names(1) = 'RE_PV_int_0'
        X_2_var_names(2) = 'IM_PV_int_0'
        X_2_var_names(3) = 'RE_PV_int_1'
        X_2_var_names(4) = 'IM_PV_int_1'
        X_2_var_names(5) = 'RE_PV_int_2'
        X_2_var_names(6) = 'IM_PV_int_2'
        X_2_var_names(7) = 'RE_KV_int_0'
        X_2_var_names(8) = 'IM_KV_int_0'
        X_2_var_names(9) = 'RE_KV_int_1'
        X_2_var_names(10) = 'IM_KV_int_1'
        X_2_var_names(11) = 'RE_KV_int_2'
        X_2_var_names(12) = 'IM_KV_int_2'
    end subroutine init_X_vars
    
    ! Create  a  vectorial  or  tensorial perturbation  type  and  allocate  the
    ! variables, the number of modes, as well as n and m
    ! Optionally, the  secondary mode  numbers can be  specified (m  if poloidal
    ! flux is used and n if toroidal  flux). By default, they are taken from the
    ! global X_vars variables.
    ! Note: The  lowest limits of the grid  need to be 1; e.g.  grid_X%i_min = 1
    ! for first process.
    subroutine create_X_1(grid_X,X,lim_sec_X)                                   ! vectorial version
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_1_type), intent(inout) :: X                                      ! vectorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2)                           ! limits of m_X (pol. flux) or n_X (tor. flux)
        
        ! local variables
        integer :: loc_n_r                                                      ! local nr. of normal points
        integer :: n_par, n_geo                                                 ! tot. nr. of angular points in parallel and geodesic direction
        
        ! set local variables
        loc_n_r = grid_X%loc_n_r
        n_par = grid_X%n(1)
        n_geo = grid_X%n(2)
        
        ! set mode numbers
        call set_nm_X(grid_X,X%n,X%m,lim_sec_X)
        
        ! set n_mod
        X%n_mod = size(X%n,2)
        
        ! allocate U_0
        allocate(X%U_0(n_par,n_geo,loc_n_r,X%n_mod))
        
        ! allocate U_1
        allocate(X%U_1(n_par,n_geo,loc_n_r,X%n_mod))
        
        ! allocate DU_0
        allocate(X%DU_0(n_par,n_geo,loc_n_r,X%n_mod))
        
        ! allocate DU_1
        allocate(X%DU_1(n_par,n_geo,loc_n_r,X%n_mod))
    end subroutine create_X_1
    subroutine create_X_2(grid_X,X,lim_sec_X)                                   ! tensorial version
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_2_type), intent(inout) :: X                                      ! tensorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol. flux) or n_X (tor. flux) for both dimensions
        
        ! local variables
        integer :: loc_n_r                                                      ! local nr. of normal points
        integer :: n_par, n_geo                                                 ! tot. nr. of angular points in parallel and geodesic direction
        integer :: nn_mod                                                       ! either n_mod^2, n_mod*(n_mod+1)/2 or 0
        
        ! set local variables
        loc_n_r = grid_X%loc_n_r
        n_par = grid_X%n(1)
        n_geo = grid_X%n(2)
        
        ! set mode numbers
        call set_nm_X(grid_X,X%n_1,X%m_1,X%n_2,X%m_2,lim_sec_X)
        
        ! set n_mod
        X%n_mod(1) = size(X%n_1,2)
        X%n_mod(2) = size(X%n_2,2)
        
        ! set nnmod for symmetric quantities
        nn_mod = set_nn_mod(lim_sec_X)
        
        ! allocate PV_i
        allocate(X%PV_0(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        allocate(X%PV_1(n_par,n_geo,loc_n_r,product(X%n_mod)))                  ! not symmetric
        allocate(X%PV_2(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        
        ! allocate KV_i
        allocate(X%KV_0(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        allocate(X%KV_1(n_par,n_geo,loc_n_r,product(X%n_mod)))                  ! not symmetric
        allocate(X%KV_2(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        
        ! allocate PV_int_i
        allocate(X%PV_int_0(n_geo,loc_n_r,nn_mod))                              ! symmetric
        allocate(X%PV_int_1(n_geo,loc_n_r,product(X%n_mod)))                    ! not symmetric
        allocate(X%PV_int_2(n_geo,loc_n_r,nn_mod))                              ! symmetric
        
        ! allocate KV_int_i
        allocate(X%KV_int_0(n_geo,loc_n_r,nn_mod))                              ! symmetric
        allocate(X%KV_int_1(n_geo,loc_n_r,product(X%n_mod)))                    ! not symmetric
        allocate(X%KV_int_2(n_geo,loc_n_r,nn_mod))                              ! symmetric
        
        ! allocate vacuum response
        allocate(X%vac_res(X%n_mod(1),X%n_mod(2)))
    end subroutine create_X_2
    
    ! Sets number of entries for symmetric tensorial perturbation variables.
    integer function set_nn_mod(lim_sec_X) result(nn_mod)
        use num_vars, only: use_pol_flux_F
        
        ! input / output
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        integer :: id, jd                                                       ! counters
        integer :: lim_sec_X_loc(2,2)                                           ! local version of lim_sec_X
        
        ! set local lim_sec_X
        lim_sec_X_loc(:,1) = [1,n_mod_X]
        lim_sec_X_loc(:,2) = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! set nnmod for symmetric quantities: discard indices above diagonal
        nn_mod = 0
        do jd = lim_sec_X_loc(1,2),lim_sec_X_loc(2,2)
            do id = lim_sec_X_loc(1,1),lim_sec_X_loc(2,1)
                if (use_pol_flux_F) then
                    if (id.ge.jd) nn_mod = nn_mod + 1
                else
                    if (id.ge.jd) nn_mod = nn_mod + 1
                end if
            end do
        end do
    end function set_nn_mod
    
    ! Sets  n_X  and  m_X  using  by default  global  variables  but  optionally
    ! different limits for the secundary mode  numbers (m_X for poloidal flux or
    ! n_X for toroidal flux).
    ! Note: The  lowest limits of the grid  need to be 1; e.g.  grid_X%i_min = 1
    ! for first process.
    subroutine set_nm_X_1(grid_X,n_X_loc,m_X_loc,lim_sec_X)                     ! vectorial version
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        integer, intent(inout), allocatable :: n_X_loc(:,:), m_X_loc(:,:)       ! toroidal and poloidal mode numbers
        integer, intent(in), optional :: lim_sec_X(2)
        
        ! local variables
        integer :: lim_sec_X_loc(2)                                             ! local version of lim_sec_X
        
        ! set local lim_sec_X
        lim_sec_X_loc = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! set n and m
        allocate(n_X_loc(grid_X%loc_n_r,lim_sec_X_loc(2)-lim_sec_X_loc(1)+1))
        allocate(m_X_loc(grid_X%loc_n_r,lim_sec_X_loc(2)-lim_sec_X_loc(1)+1))
        n_X_loc = n_X(grid_X%i_min:grid_X%i_max,&
            &lim_sec_X_loc(1):lim_sec_X_loc(2))
        m_X_loc = m_X(grid_X%i_min:grid_X%i_max,&
            &lim_sec_X_loc(1):lim_sec_X_loc(2))
    end subroutine set_nm_X_1
    subroutine set_nm_X_2(grid_X,n_X_1,m_X_1,n_X_2,m_X_2,lim_sec_X)             ! tensorial version
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        integer, intent(inout), allocatable :: n_X_1(:,:), m_X_1(:,:)           ! toroidal and poloidal mode numbers for dimension 1
        integer, intent(inout), allocatable :: n_X_2(:,:), m_X_2(:,:)           ! toroidal and poloidal mode numbers for dimension 2
        integer, intent(in), optional :: lim_sec_X(2,2)
        
        ! call vectorial version
        if (present(lim_sec_X)) then
            call set_nm_X(grid_X,n_X_1,m_X_1,lim_sec_X(:,1))
            call set_nm_X(grid_X,n_X_2,m_X_2,lim_sec_X(:,2))
        else
            call set_nm_X(grid_X,n_X_1,m_X_1)
            call set_nm_X(grid_X,n_X_2,m_X_2)
        end if
    end subroutine set_nm_X_2
    
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
