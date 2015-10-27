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
    public init_X_vars, dealloc_X, create_X, get_suffix, set_nm_X, &
        &is_necessary_X, set_nn_mod, &
        &X_1_type, X_2_type, &
        &min_n_X, max_n_X, min_m_X, max_m_X, min_n_r_sol, min_r_sol, &
        &max_r_sol, X_1_var_names, X_2_var_names
    
    ! global variables
    real(dp) :: min_r_sol, max_r_sol                                            ! min. and max. normal range for pert. (either pol. or tor., depending on use_pol_flux_F)
    integer :: min_n_r_sol                                                      ! min. of n_r_sol (e.g. first value in Richardson loop)
    integer :: min_n_X                                                          ! lowest poloidal mode number m_X
    integer :: max_n_X                                                          ! highest poloidal mode number m_X
    integer :: min_m_X                                                          ! lowest poloidal mode number m_X
    integer :: max_m_X                                                          ! highest poloidal mode number m_X
    character(len=max_str_ln), allocatable :: X_1_var_names(:)                  ! internal vectorial perturbation variables names
    character(len=max_str_ln), allocatable :: X_2_var_names(:)                  ! internal tensorial perturbation variables names
    
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
    !   - (angle_1,angle_2,r,n_mod^2)       for PVi, KVi
    !   - (n_mod^2,angle_2,r,3)             for KV_int, PV_int
    ! where it is refered to the discussion  of the grid type for an explanation
    ! of the angles angle_1 and angle_2.
    type :: X_2_type
        integer :: n_mod(2)                                                     ! size of n and m (nr. of modes)
        integer, allocatable :: n_1(:)                                          ! vector of toroidal mode numbers of dimension 1
        integer, allocatable :: n_2(:)                                          ! vector of toroidal mode numbers of dimension 2
        integer, allocatable :: m_1(:)                                          ! vector of poloidal mode numbers of dimension 1
        integer, allocatable :: m_2(:)                                          ! vector of poloidal mode numbers of dimension 2
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
    interface set_nm_X
        module procedure set_nm_X_1, set_nm_X_2
    end interface
    interface get_suffix
        module procedure get_suffix_1, get_suffix_2, get_suffix_1_lim, &
            &get_suffix_2_lim
    end interface
    interface set_nn_mod
        module procedure set_nn_mod_X, set_nn_mod_lim
    end interface
    interface is_necessary_X
        module procedure is_necessary_X_X, is_necessary_X_lim
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
    subroutine create_X_1(grid_eq,X,lim_sec_X)                                  ! vectorial version
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(X_1_type), intent(inout) :: X                                      ! vectorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2)                           ! limits of m_X (pol. flux) or n_X (tor. flux)
        
        ! local variables
        integer :: loc_n_r                                                      ! local nr. of normal points
        integer :: n_par, n_geo                                                 ! tot. nr. of angular points in parallel and geodesic direction
        
        ! set local variables
        loc_n_r = grid_eq%loc_n_r
        n_par = grid_eq%n(1)
        n_geo = grid_eq%n(2)
        
        ! set mode numbers
        call set_nm_X(X%n,X%m,lim_sec_X)
        
        ! set n_mod
        X%n_mod = size(X%n)
        
        ! allocate U_0
        allocate(X%U_0(n_par,n_geo,loc_n_r,X%n_mod))
        
        ! allocate U_1
        allocate(X%U_1(n_par,n_geo,loc_n_r,X%n_mod))
        
        ! allocate DU_0
        allocate(X%DU_0(n_par,n_geo,loc_n_r,X%n_mod))
        
        ! allocate DU_1
        allocate(X%DU_1(n_par,n_geo,loc_n_r,X%n_mod))
    end subroutine create_X_1
    subroutine create_X_2(grid_eq,X,lim_sec_X)                                  ! tensorial version
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(X_2_type), intent(inout) :: X                                      ! tensorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol. flux) or n_X (tor. flux) for both dimensions
        
        ! local variables
        integer :: loc_n_r                                                      ! local nr. of normal points
        integer :: n_par, n_geo                                                 ! tot. nr. of angular points in parallel and geodesic direction
        integer :: nn_mod                                                       ! either n_mod^2, n_mod*(n_mod+1)/2 or 0
        
        ! set local variables
        loc_n_r = grid_eq%loc_n_r
        n_par = grid_eq%n(1)
        n_geo = grid_eq%n(2)
        
        ! set mode numbers
        call set_nm_X(X%n_1,X%m_1,X%n_2,X%m_2,lim_sec_X)
        
        ! set n_mod
        X%n_mod(1) = size(X%n_1)
        X%n_mod(2) = size(X%n_2)
        
        ! set nnmod for symmetric quantities
        nn_mod = set_nn_mod(X)
        
        ! allocate PV_i
        allocate(X%PV_0(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        allocate(X%PV_1(n_par,n_geo,loc_n_r,product(X%n_mod)))                  ! not symmetric
        allocate(X%PV_2(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        
        ! allocate KV_i
        allocate(X%KV_0(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        allocate(X%KV_1(n_par,n_geo,loc_n_r,product(X%n_mod)))                  ! not symmetric
        allocate(X%KV_2(n_par,n_geo,loc_n_r,nn_mod))                            ! symmetric
        
        ! allocate PV_int_i
        allocate(X%PV_int_0(nn_mod,n_geo,loc_n_r))                              ! symmetric
        allocate(X%PV_int_1(product(X%n_mod),n_geo,loc_n_r))                    ! not symmetric
        allocate(X%PV_int_2(nn_mod,n_geo,loc_n_r))                              ! symmetric
        
        ! allocate KV_int_i
        allocate(X%KV_int_0(nn_mod,n_geo,loc_n_r))                              ! symmetric
        allocate(X%KV_int_1(product(X%n_mod),n_geo,loc_n_r))                    ! not symmetric
        allocate(X%KV_int_2(nn_mod,n_geo,loc_n_r))                              ! symmetric
        
        ! allocate vacuum response
        allocate(X%vac_res(X%n_mod(1),X%n_mod(2)))
    end subroutine create_X_2
    
    ! Sets number of entries for symmetric tensorial perturbation variables.
    integer function set_nn_mod_X(X) result(nn_mod)                             ! version using X
        use num_vars, only: use_pol_flux_F
        
        ! input / output
        type(X_2_type), intent(inout) :: X                                      ! tensorial perturbation variables
        
        ! local variables
        integer :: id, jd                                                       ! counters
        
        ! set nnmod for symmetric quantities: discard indices above diagonal
        nn_mod = 0
        do jd = 1,X%n_mod(2)
            do id = 1,X%n_mod(1)
                if (use_pol_flux_F) then
                    if (X%m_1(id).ge.X%m_2(jd)) nn_mod = nn_mod + 1
                else
                    if (X%n_1(id).ge.X%n_2(jd)) nn_mod = nn_mod + 1
                end if
            end do
        end do
    end function set_nn_mod_X
    integer function set_nn_mod_lim(lim_sec_X) result(nn_mod)                   ! version using limits
        use num_vars, only: use_pol_flux_F
        
        ! input / output
        integer, intent(in) :: lim_sec_X(2,2)                                   ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        integer :: id, jd                                                       ! counters
        
        ! set nnmod for symmetric quantities: discard indices above diagonal
        nn_mod = 0
        do jd = lim_sec_X(1,2),lim_sec_X(2,2)
            do id = lim_sec_X(1,1),lim_sec_X(2,1)
                if (use_pol_flux_F) then
                    if (id.ge.jd) nn_mod = nn_mod + 1
                else
                    if (id.ge.jd) nn_mod = nn_mod + 1
                end if
            end do
        end do
    end function set_nn_mod_lim
    
    ! Sets  n_X  and  m_X  using  by default  global  variables  but  optionally
    ! different limits for the secundary mode  numbers (m_X for poloidal flux or
    ! n_X for toroidal flux)
    subroutine set_nm_X_1(n_X,m_X,lim_sec_X)                                    ! vectorial version
        use num_vars, only: use_pol_flux_F
        
        ! input / output
        integer, intent(inout), allocatable :: n_X(:), m_X(:)                   ! toroidal and poloidal mode numbers
        integer, intent(in), optional :: lim_sec_X(2)
        
        ! local variables
        integer :: n_mod                                                        ! nr. of modes
        integer :: id                                                           ! counter
        
        ! set nr. of modes
        if (present(lim_sec_X)) then
            n_mod = lim_sec_X(2)-lim_sec_X(1)+1
        else
            n_mod = (max_n_X-min_n_X+1)*(max_m_X-min_m_X+1)
        end if
        
        ! set n and m
        allocate(n_X(n_mod),m_X(n_mod))
        if (use_pol_flux_F) then
            n_X = min_n_X
            if (present(lim_sec_X)) then
                m_X = [(id, id = lim_sec_X(1),lim_sec_X(2))]
            else
                m_X = [(id, id = min_m_X,max_m_X)]
            end if
        else
            if (present(lim_sec_X)) then
                n_X = [(id, id = lim_sec_X(1),lim_sec_X(2))]
            else
                n_X = [(id, id = min_n_X,max_n_X)]
            end if
            m_X = min_m_X
        end if
    end subroutine set_nm_X_1
    subroutine set_nm_X_2(n_X_1,m_X_1,n_X_2,m_X_2,lim_sec_X)                    ! tensorial version
        ! input / output
        integer, intent(inout), allocatable :: n_X_1(:), m_X_1(:)               ! toroidal and poloidal mode numbers for dimension 1
        integer, intent(inout), allocatable :: n_X_2(:), m_X_2(:)               ! toroidal and poloidal mode numbers for dimension 2
        integer, intent(in), optional :: lim_sec_X(2,2)
        
        ! call vectorial version
        if (present(lim_sec_X)) then
            call set_nm_X(n_X_1,m_X_1,lim_sec_X(:,1))
            call set_nm_X(n_X_2,m_X_2,lim_sec_X(:,2))
        else
            call set_nm_X(n_X_1,m_X_1)
            call set_nm_X(n_X_2,m_X_2)
        end if
    end subroutine set_nm_X_2
    
    ! Sets the suffix used to refer to a perturbation quantity.
    ! Note  that for  the limit  variants  the global  variables 'min_n_X',  and
    ! 'min_m_X' have to be correct.
    character(len=max_str_ln) function get_suffix_1(X,id) result(res)           ! vectorial version
        ! input / output
        type(X_1_type), intent(in) :: X                                         ! vectorial perturbation variables
        integer, intent(in) :: id                                               ! mode index
        
        ! set suffix
        res = trim(i2str(X%n(id)))//'_'//trim(i2str(X%m(id)))
    end function get_suffix_1
    character(len=max_str_ln) function get_suffix_2(X,id,jd) result(res)        ! tensorial version
        ! input / output
        type(X_2_type), intent(in) :: X                                         ! tensorial perturbation variables
        integer, intent(in) :: id, jd                                           ! mode indices
        
        ! set suffix
        res = trim(i2str(X%n_1(id)))//'_'//trim(i2str(X%m_1(id)))//&
            &'_'//trim(i2str(X%n_2(jd)))//'_'//trim(i2str(X%m_2(jd)))
    end function get_suffix_2
    character(len=max_str_ln) function get_suffix_1_lim(lim_sec_X,id) &
        &result(res)                                                            ! vectorial version, using limits
        use num_vars, only: use_pol_flux_F
        
        ! input / output
        integer, intent(in) :: lim_sec_X(2)                                     ! limits of m_X (pol. flux) or n_X (tor. flux)
        integer, intent(in) :: id                                               ! mode index
        
        ! set suffix
        if (use_pol_flux_F) then
            res = trim(i2str(min_n_X))//'_'//trim(i2str(lim_sec_X(1)-1+id))
        else
            res = trim(i2str(lim_sec_X(1)-1+id))//'_'//trim(i2str(min_m_X))
        end if
    end function get_suffix_1_lim
    character(len=max_str_ln) function get_suffix_2_lim(lim_sec_X,id,jd) &
        &result(res)                                                            ! tensorial version, using limits
        use num_vars, only: use_pol_flux_F
        
        ! input / output
        integer, intent(in) :: lim_sec_X(2,2)                                   ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        integer, intent(in) :: id, jd                                           ! mode indices
        
        ! set suffix
        if (use_pol_flux_F) then
            res = trim(i2str(min_n_X))//'_'//&
                &trim(i2str(lim_sec_X(1,1)-1+id))//'_'//&
                &trim(i2str(min_n_X))//'_'//&
                &trim(i2str(lim_sec_X(1,2)-1+jd))
        else
            res = trim(i2str(lim_sec_X(1,1)-1+id))//'_'//&
                &trim(i2str(min_m_X))//'_'//&
                &trim(i2str(lim_sec_X(1,2)-1+jd))//'_'//&
                &trim(i2str(min_m_X))
        end if
    end function get_suffix_2_lim
    
    ! deallocates perturbation variables
    subroutine dealloc_X_1(X)                                                   ! vectorial version
        ! input / output
        type(X_1_type), intent(out) :: X                                        ! perturbation variables to be deallocated
    end subroutine dealloc_X_1
    subroutine dealloc_X_2(X)                                                   ! tensorial version
        ! input / output
        type(X_2_type), intent(out) :: X                                        ! perturbation variables to be deallocated
    end subroutine dealloc_X_2
    
    ! Determines whether a variable needs to be  considered: Only if it is on or
    ! below the diagonal for symmetric quantities.
    logical function is_necessary_X_X(X,sym,sec_X_id) result(res)               ! version using X
        use num_vars, only: use_pol_flux_F
        
        ! input / output
        type(X_2_type), intent(in) :: X                                         ! tensorial perturbation variables
        logical, intent(in) :: sym                                              ! whether the variable is symmetric
        integer, intent(in) :: sec_X_id(2)                                      ! mode indices
        
        ! initialize res
        res = .true.
        
        ! modify res depending on symmetry
        if (sym) then
            if (use_pol_flux_F) then
                if (X%m_1(sec_X_id(1)).lt.X%m_2(sec_X_id(2))) res = .false.
            else
                if (X%n_1(sec_X_id(1)).lt.X%n_2(sec_X_id(2))) res = .false.
            end if
        end if
    end function is_necessary_X_X
    logical function is_necessary_X_lim(lim_sec_X,sym,sec_X_id) result(res)     ! version using limits
        ! input / output
        integer, intent(in) :: lim_sec_X(2,2)                                   ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        logical, intent(in) :: sym                                              ! whether the variable is symmetric
        integer, intent(in) :: sec_X_id(2)                                      ! mode indices
        
        ! initialize res
        res = .true.
        
        ! modify res depending on symmetry
        if (sym) then
            if (lim_sec_X(1,1)+sec_X_id(1).lt.lim_sec_X(1,2)+sec_X_id(2)) &
                &res = .false.
        end if
    end function is_necessary_X_lim
end module X_vars
