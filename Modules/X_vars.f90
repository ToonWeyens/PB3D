!------------------------------------------------------------------------------!
!   Variables pertaining to the perturbation quantities                        !
!------------------------------------------------------------------------------!
module X_vars
#include <PB3D_macros.h>
    use str_utilities
    use messages
    use num_vars, only: dp, max_name_ln, iu, weight_dp
    use grid_vars, only: grid_type
    
    implicit none
    
    private
    public set_nm_X, set_nn_mod, &
        &X_1_type, X_2_type, &
        &n_mod_X, prim_X, min_sec_X, max_sec_X, min_nm_X, min_n_X, max_n_X, &
        &min_m_X, max_m_X, min_r_sol, max_r_sol, n_X, m_X, sec_X_ind
#if ldebug
    public n_alloc_X_1s, n_alloc_X_2s
#endif
    
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
#if ldebug
    integer :: n_alloc_X_1s                                                     ! nr. of allocated X_1's
    integer :: n_alloc_X_2s                                                     ! nr. of allocated X_2's
#endif
    
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
#if ldebug
        real(dp) :: estim_mem_usage                                             ! estimated memory usage
#endif
    contains
        procedure :: init => init_X_1
        procedure :: dealloc => dealloc_X_1
    end type
    
    ! tensorial perturbation type with arrays of the form:
    !   - (angle_1,angle_2,r,n_mod^2)       for PVi, KVi
    !   - (n_mod,n_mod)                     for vac_res
    ! where it is refered to the discussion  of the grid type for an explanation
    ! of the angles angle_1 and angle_2.
    ! Note that this type is also used for field-averaged tensorial perturbation
    ! variables, with angle_1 = 1.
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
        complex(dp), allocatable :: vac_res(:,:)                                ! vacuum response
#if ldebug
        real(dp) :: estim_mem_usage                                             ! estimated memory usage
#endif
    contains
        procedure :: init => init_X_2
        procedure :: dealloc => dealloc_X_2
    end type
    
    ! interfaces
    interface set_nm_X
        module procedure set_nm_X_1, set_nm_X_2
    end interface
    
contains
    ! Initializes a  vectorial or tensorial  perturbation type and  allocate the
    ! variables, the number of modes, as well as n and m
    ! Optionally, the  secondary mode  numbers can be  specified (m  if poloidal
    ! flux is used and n if toroidal  flux). By default, they are taken from the
    ! global X_vars variables.
    ! Furthermore, the tensorial  perturbation type can also be  used for field-
    ! aligned  variables, in  which  case the  first index  is  assumed to  have
    ! dimension 1 only. This can be triggered using "is_field_averaged".
    ! Note: The  lowest limits of the grid  need to be 1; e.g.  grid_X%i_min = 1
    ! for first process.
    subroutine init_X_1(X,grid_X,lim_sec_X)                                     ! vectorial version
#if ldebug
        use num_vars, only: print_mem_usage, rank
#endif
        
        ! input / output
        class(X_1_type), intent(inout) :: X                                     ! vectorial perturbation variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        integer, intent(in), optional :: lim_sec_X(2)                           ! limits of m_X (pol. flux) or n_X (tor. flux)
        
        ! local variables
        integer :: loc_n_r                                                      ! local nr. of normal points
        integer :: n_par, n_geo                                                 ! tot. nr. of angular points in parallel and geodesic direction
        
        ! set local variables
        loc_n_r = grid_X%loc_n_r
        n_par = grid_X%n(1)
        n_geo = grid_X%n(2)
        
#if ldebug
        ! initialize memory usage
        if (print_mem_usage) X%estim_mem_usage = 0._dp
#endif
        
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
        
#if ldebug
        ! set estimated memory usage
        if (print_mem_usage) X%estim_mem_usage = X%estim_mem_usage + &
            &(n_par*n_geo*loc_n_r)*(X%n_mod*4)
        
        ! increment n_alloc_X_1s
        n_alloc_X_1s = n_alloc_X_1s + 1
        
        ! print memory usage
        if (print_mem_usage) call writo('[rank '//trim(i2str(rank))//&
            &' - Expected memory usage of X_1: '//&
            &trim(r2strt(X%estim_mem_usage*weight_dp*2))//' kB]',alert=.true.)
#endif
    end subroutine init_X_1
    subroutine init_X_2(X,grid_X,lim_sec_X,is_field_averaged)                   ! tensorial version
#if ldebug
        use num_vars, only: print_mem_usage, rank
#endif
        
        ! input / output
        class(X_2_type), intent(inout) :: X                                     ! tensorial perturbation variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol. flux) or n_X (tor. flux) for both dimensions
        logical, intent(in), optional :: is_field_averaged                      ! if field-aligned, only one dimension for first index
        
        ! local variables
        integer :: loc_n_r                                                      ! local nr. of normal points
        integer :: n_par, n_geo                                                 ! tot. nr. of angular points in parallel and geodesic direction
        integer :: nn_mod                                                       ! either n_mod^2, n_mod*(n_mod+1)/2 or 0
        
        ! set local variables
        n_par = grid_X%n(1)
        if (present(is_field_averaged)) then
            if (is_field_averaged) n_par = 1
        end if
        n_geo = grid_X%n(2)
        loc_n_r = grid_X%loc_n_r
        
#if ldebug
        ! initialize memory usage
        if (print_mem_usage) X%estim_mem_usage = 0._dp
#endif
        
        ! set mode numbers
        call set_nm_X(grid_X,X%n_1,X%m_1,X%n_2,X%m_2,lim_sec_X)
        
        ! set n_mod
        X%n_mod(1) = size(X%n_1,2)
        X%n_mod(2) = size(X%n_2,2)
        
        ! set nnmod for symmetric quantities
        nn_mod = set_nn_mod(.true.,lim_sec_X)
        
        ! allocate PV_i
        allocate(X%PV_0(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        allocate(X%PV_1(n_par,n_geo,loc_n_r,product(X%n_mod)))                  ! not symmetric
        allocate(X%PV_2(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        
        ! allocate KV_i
        allocate(X%KV_0(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        allocate(X%KV_1(n_par,n_geo,loc_n_r,product(X%n_mod)))                  ! not symmetric
        allocate(X%KV_2(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        
        ! allocate vacuum response
        allocate(X%vac_res(X%n_mod(1),X%n_mod(2)))
        
#if ldebug
        ! set estimated memory usage
        if (print_mem_usage) X%estim_mem_usage = &
            &X%estim_mem_usage + (n_par*n_geo*loc_n_r)*&
            &(nn_mod*4+product(X%n_mod)*2) + product(X%n_mod)
        
        ! increment n_alloc_X_2s
        n_alloc_X_2s = n_alloc_X_2s + 1
        
        ! print memory usage
        if (print_mem_usage) call writo('[rank '//trim(i2str(rank))//&
            &' - Expected memory usage of X_2: '//&
            &trim(r2strt(X%estim_mem_usage*weight_dp*2))//' kB]',alert=.true.)
#endif
    end subroutine init_X_2
    
    ! Sets number of entries for tensorial perturbation variables.
    integer function set_nn_mod(sym,lim_sec_X) result(nn_mod)
        ! input / output
        logical, intent(in) :: sym                                              ! whether the variable is symmetric
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        integer :: id, jd                                                       ! counters
        integer :: lim_sec_X_loc(2,2)                                           ! local version of lim_sec_X
        
        ! set local lim_sec_X
        lim_sec_X_loc(:,1) = [1,n_mod_X]
        lim_sec_X_loc(:,2) = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        if (sym) then
            ! set nnmod for symmetric quantities: discard indices above diagonal
            nn_mod = 0
            do jd = lim_sec_X_loc(1,2),lim_sec_X_loc(2,2)
                do id = lim_sec_X_loc(1,1),lim_sec_X_loc(2,1)
                    if (id.ge.jd) nn_mod = nn_mod + 1
                end do
            end do
        else
            ! set nn_mod for asymmetric quantities: don't discard anything
            nn_mod = product(lim_sec_X_loc(2,:)-lim_sec_X_loc(1,:)+1)
        end if
    end function set_nn_mod
    
    ! Sets  n_X  and  m_X  using  by default  global  variables  but  optionally
    ! different limits for the secundary mode  numbers (m_X for poloidal flux or
    ! n_X for toroidal flux).
    ! Note: n_X  and m_X need to  have been set up  with the same limits  as the
    ! grid used here.
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
#if ldebug
        use num_vars, only: rank, print_mem_usage
#endif
        
        ! input / output
        class(X_1_type), intent(inout) :: X                                     ! perturbation variables to be deallocated
        
#if ldebug
        ! local variables
        integer :: mem_diff                                                     ! difference in memory
        real(dp) :: estim_mem_usage                                             ! estimated memory usage
        
        ! memory usage before deallocation
        if (print_mem_usage) then
            mem_diff = get_mem_usage()
            estim_mem_usage = X%estim_mem_usage
        end if
#endif
        
        ! deallocate allocatable variables
        call dealloc_X_1_final(X)
        
#if ldebug
        ! decrement n_alloc_X_1s
        n_alloc_X_1s = n_alloc_X_1s - 1
        
        ! memory usage difference after deallocation
        if (print_mem_usage) then
            mem_diff = mem_diff - get_mem_usage()
            call writo('[Rank '//trim(i2str(rank))//' - liberated '//&
                &trim(i2str(mem_diff))//'kB deallocating X_1 ('//&
                &trim(i2str(nint(100*mem_diff/&
                &(estim_mem_usage*weight_dp*2))))//&
                &'% of estimated)]',alert=.true.)
        end if
#endif
    contains
        ! Note: intent(out) automatically deallocates the variable
        subroutine dealloc_X_1_final(X)
            ! input / output
            type(X_1_type), intent(out) :: X                                    ! equilibrium to be deallocated
        end subroutine dealloc_X_1_final
    end subroutine dealloc_X_1
    subroutine dealloc_X_2(X)                                                   ! tensorial version
#if ldebug
        use num_vars, only: rank, print_mem_usage
#endif
        
        ! input / output
        class(X_2_type), intent(inout) :: X                                     ! perturbation variables to be deallocated
        
#if ldebug
        ! local variables
        integer :: mem_diff                                                     ! difference in memory
        real(dp) :: estim_mem_usage                                             ! estimated memory usage
        
        ! memory usage before deallocation
        if (print_mem_usage) then
            mem_diff = get_mem_usage()
            estim_mem_usage = X%estim_mem_usage
        end if
#endif
        
        ! deallocate allocatable variables
        call dealloc_X_2_final(X)
        
#if ldebug
        ! decrement n_alloc_X_2s
        n_alloc_X_2s = n_alloc_X_2s - 1
        
        ! memory usage difference after deallocation
        if (print_mem_usage) then
            mem_diff = mem_diff - get_mem_usage()
            call writo('[Rank '//trim(i2str(rank))//' - liberated '//&
                &trim(i2str(mem_diff))//'kB deallocating X_2 ('//&
                &trim(i2str(nint(100*mem_diff/&
                &(estim_mem_usage*weight_dp*2))))//&
                &'% of estimated)]',alert=.true.)
        end if
#endif
    contains
        ! Note: intent(out) automatically deallocates the variable
        subroutine dealloc_X_2_final(X)
            ! input / output
            type(X_2_type), intent(out) :: X                                    ! equilibrium to be deallocated
        end subroutine dealloc_X_2_final
    end subroutine dealloc_X_2
end module X_vars
