!------------------------------------------------------------------------------!
!   Numerical utilities related to perturbation operations                     !
!------------------------------------------------------------------------------!
module X_utilities
#include <PB3D_macros.h>
    use str_utilities
    use messages
    use num_vars, only: dp, max_name_ln, iu, max_str_ln
    use grid_vars, only: grid_type
    use X_vars, only: n_mod_X
    
    implicit none
    
    private
    public sec_ind_loc2tot, is_necessary_X, get_sec_X_range, calc_memory_X
    
    ! interfaces
    interface sec_ind_loc2tot
        module procedure sec_ind_loc2tot_1, sec_ind_loc2tot_2
    end interface
    
contains
    ! Sets the sec_ind_tot used to refer to a perturbation quantity.
    function sec_ind_loc2tot_1(id,lim_sec_X) result(res)                        ! vectorial version
        ! input / output
        integer, intent(in) :: id                                               ! mode index
        integer, intent(in), optional :: lim_sec_X(2)                           ! limits of m_X (pol. flux) or n_X (tor. flux)
        integer :: res                                                          ! output
        
        ! local variables
        integer :: lim_sec_X_loc(2)                                             ! local lim_sec_X
        
        ! set local lim_sec_X
        lim_sec_X_loc = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! set sec_ind_tot
        res = lim_sec_X_loc(1)-1+id
    end function sec_ind_loc2tot_1
    function sec_ind_loc2tot_2(id,jd,lim_sec_X) result(res)                     ! tensorial version
        ! input / output
        integer, intent(in) :: id, jd                                           ! mode indices
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        integer :: res(2)                                                       ! output
        
        ! local variables
        integer :: lim_sec_X_loc(2,2)                                           ! local lim_sec_X
        
        ! set local lim_sec_X
        lim_sec_X_loc(:,1) = [1,n_mod_X]
        lim_sec_X_loc(:,2) = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! set sec_ind_tot
        res = [lim_sec_X_loc(1,1)-1+id,lim_sec_X_loc(1,2)-1+jd]
    end function sec_ind_loc2tot_2
    
    ! Gets  one of  the the  local ranges  of contiguous  tensorial perturbation
    ! variables to be printed or read  during one call of the corresponding HDF5
    ! variables. 
    ! More specifically,  a range of indices  k in the first  dimension is given
    ! for every value  of the indices m  in the second dimension.  An example is
    ! now given for the subrange [2:3,2:5] of  a the total range [1:5, 1:5]. For
    ! asymmetric variables the situation is simple: The k range is [2:3] for all
    ! 5 values of m. However, for symmetric variables, the upper diagonal values
    ! are not stored, which  gives k ranges [2:3], [3:3] and no range  for m = 4
    ! and 5.
    ! This routine then translates these  ranges to the corresponding 1-D ranges
    ! that  are used  in  the actual  variables. For  above  example, the  total
    ! indices are
    !   [  1  6 11 16 21 ]
    !   [  2  7 12 17 22 ]
    !   [  3  8 13 18 23 ] -> [7:8], [12:13], [17:18] and [22:23],
    !   [  4  9 14 19 24 ]
    !   [  5 10 15 20 25 ]
    ! for asymmetric variables and
    !   [  1  2  3  4  5 ]
    !   [  2  6  7  8  9 ]
    !   [  3  7 10 11 12 ] -> [6:7], [10:10], [:] and [:],
    !   [  4  8 11 13 14 ]
    !   [  5  9 12 14 15 ]
    ! for symmetric variables.
    ! These can  then related  to the  local indices for  the variables  in this
    ! perturbation job. For above example, the results are:
    !   [1:2], [3:4], [5:6] and [7:8],
    ! for asymmetric variables and
    !   [1:2], [3:3], [:] and [:], 
    ! for symmetric variables.
    ! As can be seen, the local ranges of the variables in the submatrix of this
    ! perturbation job are (designed to be)  contiguous, but the total ranges of
    ! the variables in the submatrix are clearly not in general.
    ! The procedure outputs both the local and total ranges.
    subroutine get_sec_X_range(sec_X_range_loc,sec_X_range_tot,m,sym,lim_sec_X)
        use X_vars, only: n_mod_X
        use num_utilities, only: c
        
        ! input / output
        integer, intent(inout) :: sec_X_range_loc(2)                            ! start and end of local range in dimension 1 (vertical)
        integer, intent(inout) :: sec_X_range_tot(2)                            ! start and end of total range in dimension 1 (vertical)
        integer, intent(in) :: m                                                ! dimension 2 (horizontal)
        logical, intent(in) :: sym                                              ! whether the variable is symmetric
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol. flux) or n_X (tor. flux)
        
        ! local variables
        integer :: k_range_loc(2)                                               ! local range of k
        integer :: k                                                            ! counter
        integer :: n_mod                                                        ! local n_mod
        
        ! set number of modes of dimension 1
        if (present(lim_sec_X)) then
            n_mod = lim_sec_X(2,1)-lim_sec_X(1,1)+1
        else
            n_mod = n_mod_X
        endif
        
        ! initialize secondary X range
        k_range_loc = [n_mod+1,0]                                               ! initialize to inverted values, out of bounds
        
        ! find start and end of k range
        find_start: do k = 1,n_mod
            if (is_necessary_X(sym,[k,m],lim_sec_X)) then                       ! first necessary pair [k,m]
                k_range_loc(1) = k
                exit find_start
            end if
        end do find_start
        find_stop: do k = n_mod,k_range_loc(1),-1
            if (is_necessary_X(sym,[k,m],lim_sec_X)) then                       ! last necessary pair [k,m]
                k_range_loc(2) = k
                exit find_stop
            end if
        end do find_stop
        
        ! translate k range and m value to 1-D index
        do k = 1,2
            sec_X_range_loc(k) = c([k_range_loc(k),m],sym,n_mod_X,lim_sec_X)
            sec_X_range_tot(k) = c(sec_ind_loc2tot(k_range_loc(k),m,lim_sec_X),&
                &sym,n_mod_X)
        end do
    end subroutine get_sec_X_range
    
    ! Determines whether a variable needs to be  considered: Only if it is on or
    ! below the diagonal for symmetric quantities.
    logical function is_necessary_X(sym,sec_X_id,lim_sec_X) result(res)
        ! input / output
        logical, intent(in) :: sym                                              ! whether the variable is symmetric
        integer, intent(in) :: sec_X_id(2)                                      ! mode indices
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        integer :: lim_sec_X_loc(2,2)                                           ! local version of lim_sec_X
        
        ! initialize res
        res = .true.
        
        ! modify res depending on symmetry
        if (sym) then
            ! set local lim_sec_X
            lim_sec_X_loc(:,1) = [1,n_mod_X]
            lim_sec_X_loc(:,2) = [1,n_mod_X]
            if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
            
            ! set res
            if (lim_sec_X_loc(1,1)+sec_X_id(1).lt.&
                &lim_sec_X_loc(1,2)+sec_X_id(2)) res = .false.
        end if
    end function is_necessary_X
    
    ! Calculate memory in MB necessary for X variables of a certain order
    !   - order 1: 4x n_par_X x n_geo x loc_n_r x n_mod
    !   - order 2: 2x n_par_X x n_geo x loc_n_r x n_mod^2
    !              4x n_par_X x n_geo x loc_n_r x n_mod(n_mod+1)/2
    !   - higher order: not used
    ! where n_par_X x  n_geo x loc_n_r should be passed  as 'arr_size' and n_mod
    ! as well
    integer function calc_memory_X(ord,arr_size,n_mod,mem_size) result(ierr)
        use ISO_C_BINDING
        use num_vars, only: mem_scale_fac
        
        character(*), parameter :: rout_name = 'calc_memory_X'
        
        ! input / output
        integer, intent(in) :: ord                                              ! order of data
        integer, intent(in) :: arr_size                                         ! size of part of X array
        integer, intent(in) :: n_mod                                            ! number of modes
        real(dp), intent(inout) :: mem_size                                     ! total size
        
        ! local variables
        integer(C_SIZE_T) :: dp_size                                            ! size of dp
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! get size of complex variable
        dp_size = 2*sizeof(1._dp)                                               ! complex variable
        
        ! calculate memory depending on order
        select case(ord)
            case (1)                                                            ! vectorial data: U, DU
                ! set memory size
                mem_size = 4*arr_size
                mem_size = mem_size*n_mod**ord
                mem_size = mem_size*dp_size
            case (2)                                                            ! tensorial data: PV, KV
                ! set memory size
                mem_size = arr_size
                mem_size = mem_size*(2*n_mod**ord+4*n_mod*(n_mod+1)/2)
                mem_size = mem_size*dp_size
            case default
                ierr = 1
                err_msg = 'Orders > 2 are not implemented'
                CHCKERR(err_msg)
        end select
        
        ! convert B to MB
        mem_size = mem_size*1.E-6_dp
        
        ! scale memory to account for rough estimation
        mem_size = mem_size*mem_scale_fac
        
        ! test overflow
        if (mem_size.lt.0) then
            ierr = 1
            CHCKERR('Overflow occured')
        end if
        
        call lvl_ud(-1)
    end function calc_memory_X
end module X_utilities
