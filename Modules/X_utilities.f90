!------------------------------------------------------------------------------!
!> Numerical utilities related to perturbation operations.
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
    public sec_ind_loc2tot, is_necessary_X, get_sec_X_range, calc_memory_X, &
        &do_X, trim_modes
    
    ! interfaces
    
    !> \public  Returns the  \c  sec_ind_tot  used to  refer  to a  perturbation
    !! quantity.
    interface sec_ind_loc2tot
        !> \public
        module procedure sec_ind_loc2tot_1
        !> \public
        module procedure sec_ind_loc2tot_2
    end interface
    
contains
    !> \private vectorial version
    function sec_ind_loc2tot_1(id,lim_sec_X) result(res)
        ! input / output
        integer, intent(in) :: id                                               !< mode index
        integer, intent(in), optional :: lim_sec_X(2)                           !< limits of \c m_X (pol. flux) or \c n_X (tor. flux)
        integer :: res                                                          !< output
        
        ! local variables
        integer :: lim_sec_X_loc(2)                                             ! local lim_sec_X
        
        ! set local lim_sec_X
        lim_sec_X_loc = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! set sec_ind_tot
        res = lim_sec_X_loc(1)-1+id
    end function sec_ind_loc2tot_1
    !> \private tensorial version
    function sec_ind_loc2tot_2(id,jd,lim_sec_X) result(res)
        ! input / output
        integer, intent(in) :: id                                               !< mode index for dimension 1
        integer, intent(in) :: jd                                               !< mode index for dimension 1
        integer, intent(in), optional :: lim_sec_X(2,2)                         !< limits of \c m_X (pol flux) or \c n_X (tor flux) for both dimensions
        integer :: res(2)                                                       !< output
        
        ! local variables
        integer :: lim_sec_X_loc(2,2)                                           ! local lim_sec_X
        
        ! set local lim_sec_X
        lim_sec_X_loc(:,1) = [1,n_mod_X]
        lim_sec_X_loc(:,2) = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! set sec_ind_tot
        res = [lim_sec_X_loc(1,1)-1+id,lim_sec_X_loc(1,2)-1+jd]
    end function sec_ind_loc2tot_2
    
    !> Gets one  of the  the local ranges  of contiguous  tensorial perturbation
    !! variables to be printed or read during one call of the corresponding HDF5
    !! variables.
    !!
    !! More specifically,  a range  of indices  \c k in  the first  dimension is
    !! given for every value of the indices \c m in the second dimension.
    !!
    !! An example  is now  given for  the subrange  <tt>[2:3,2:5]</tt> of  a the
    !! total range  <tt>[1:5, 1:5]</tt>. For asymmetric  variables the situation
    !! is simple:  The \c k range  is <tt>[2:3]</tt> for  all 5 values of  \c m.
    !! However,  for symmetric  variables,  the upper  diagonal  values are  not
    !! stored, which  gives \c  k ranges  <tt>[2:3]</tt>, <tt>[3:3]</tt>  and no
    !! range for \c m = 4 and 5.
    !!
    !! This routine then translates these  ranges to the corresponding 1-D ranges
    !! that  are used  in  the actual  variables. For  above  example, the  total
    !! indices are
    !!  \f[
    !!      \left(\begin{array}{ccccc}
    !!          1 &  6 & 11 & 16 & 21 \\
    !!          2 &  7 & 12 & 17 & 22 \\
    !!          3 &  8 & 13 & 18 & 23 \\
    !!          4 &  9 & 14 & 19 & 24 \\
    !!          5 & 10 & 15 & 20 & 25
    !!      \end{array}\right)
    !!      \rightarrow \left[7:8\right], \left[12:13\right],
    !!      \left[17:18\right] \ \text{and} \ \left[22:23\right],
    !!  \f]
    !! for asymmetric variables and
    !!  \f[
    !!      \left(\begin{array}{ccccc}
    !!          1 & 2 &  3 &  4 &  5 \\
    !!          2 & 6 &  7 &  8 &  9 \\
    !!          3 & 7 & 10 & 11 & 12 \\
    !!          4 & 8 & 11 & 13 & 14 \\
    !!          5 & 9 & 12 & 14 & 15
    !!      \end{array}\right)
    !!      \rightarrow \left[6:7\right], \left[10:10\right],
    !!      \left[:\right] \ \text{and} \ \left[:\right],
    !!  \f]
    !! for symmetric variables.
    !!
    !! These can  then related  to the  local indices for  the variables  in this
    !! perturbation job.
    !!
    !! For above example, the results are:
    !!
    !!   <tt>[1:2], [3:4], [5:6] and [7:8]</tt>,
    !!
    !! for asymmetric variables and
    !!
    !!   <tt>[1:2], [3:3], [:] and [:]</tt>, 
    !!
    !! for symmetric variables.
    !!
    !! As can  be seen, the  local ranges of the  variables in the  submatrix of
    !! this  perturbation job  are (designed  to be)  contiguous, but  the total
    !! ranges of the variables in the submatrix are clearly not in general.
    !!
    !! The procedure outputs both the local and total ranges.
    subroutine get_sec_X_range(sec_X_range_loc,sec_X_range_tot,m,sym,lim_sec_X)
        use X_vars, only: n_mod_X
        use num_utilities, only: c
        
        ! input / output
        integer, intent(inout) :: sec_X_range_loc(2)                            !< start and end of local range in dimension 1 (vertical)
        integer, intent(inout) :: sec_X_range_tot(2)                            !< start and end of total range in dimension 1 (vertical)
        integer, intent(in) :: m                                                !< dimension 2 (horizontal)
        logical, intent(in) :: sym                                              !< whether the variable is symmetric
        integer, intent(in), optional :: lim_sec_X(2,2)                         !< limits of \c m_X (pol. flux) or \c n_X (tor. flux)
        
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
    
    !> Tests whether this perturbatino job should be done.
    !!
    !! Also increments \c X_job_nr
    logical function do_X()
        use num_vars, only: X_jobs_lims, X_job_nr
        
        ! increment perturbation job nr.
        X_job_nr = X_job_nr + 1
        
        if (X_job_nr.le.size(X_jobs_lims,2)) then
            do_X = .true.
        else
            do_X = .false.
        end if
    end function do_X
    
    !> Determines whether a variable needs to be  considered:
    !! 
    !! This depends on whether the quantity is symmetric or not:
    !!  - Only if it is on or below the diagonal for symmetric quantities.
    !!  - Always for asymmetric quantities
    logical function is_necessary_X(sym,sec_X_id,lim_sec_X) result(res)
        ! input / output
        logical, intent(in) :: sym                                              !< whether the variable is symmetric
        integer, intent(in) :: sec_X_id(2)                                      !< mode indices
        integer, intent(in), optional :: lim_sec_X(2,2)                         !< limits of \c m_X (pol flux) or \c n_X (tor flux) for both dimensions
        
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
    
    !> Calculate memory in MB necessary for X variables.
    !!
    !! This depends on the order:
    !!  - order 1:
    !!      - <tt>4x n_par_X x n_geo x loc_n_r x n_mod</tt>
    !!  - order 2:
    !!      - <tt>2x n_par_X x n_geo x loc_n_r x n_mod^2</tt>
    !!      - <tt>4x n_par_X x n_geo x loc_n_r x n_mod(n_mod+1)/2</tt>
    !!  - higher order: not used
    !!
    !! where <tt>n_par_X x n_geo x loc_n_r</tt>  should be passed as \c arr_size
    !! and \c n_mod as well.
    !!
    !! \return ierr
    integer function calc_memory_X(ord,arr_size,n_mod,mem_size) result(ierr)
        use ISO_C_BINDING
        
        character(*), parameter :: rout_name = 'calc_memory_X'
        
        ! input / output
        integer, intent(in) :: ord                                              !< order of data
        integer, intent(in) :: arr_size                                         !< size of part of X array
        integer, intent(in) :: n_mod                                            !< number of modes
        real(dp), intent(inout) :: mem_size                                     !< total size
        
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
        
        ! apply 50% safety factor (empirical)
        mem_size = mem_size*1.5_dp
        
        ! test overflow
        if (mem_size.lt.0) then
            ierr = 1
            CHCKERR('Overflow occured')
        end if
        
        call lvl_ud(-1)
    end function calc_memory_X

    !> Limit input mode range to output mode range.
    integer function trim_modes(mds_i,mds_o,id_lim_i,id_lim_o) result(ierr)
        use X_vars, only: modes_type
        
        character(*), parameter :: rout_name = 'trim_modes'
        
        ! input / output
        type(modes_type), intent(in) :: mds_i                                   !< general modes variables for input
        type(modes_type), intent(in) :: mds_o                                   !< general modes variables for output
        integer, intent(inout) :: id_lim_i(2)                                   !< limits on input modes
        integer, intent(inout) :: id_lim_o(2)                                   !< limits on output modes
        
        ! local variables
        integer :: m                                                            ! counter
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! find limits point where input modes and output modes coincide
        ! (input grid should comprise output grid)
        id_lim_o = [1,size(mds_o%sec,1)]
        id_lim_i = [-1,-1]
        do m = 1,size(mds_i%sec,1)
            if (mds_i%sec(m,1).eq.mds_o%sec(id_lim_o(1),1)) then
                id_lim_i(1) = m
                exit
            end if
        end do
        do m = size(mds_i%sec,1),1,-1
            if (mds_i%sec(m,1).eq.mds_o%sec(id_lim_o(2),1)) then
                id_lim_i(2) = m
                exit
            end if
        end do
        
        ! test
        if (any(id_lim_o.le.0)) then
            ierr = 1
            err_msg = 'Cannot find limits of input modes'
            CHCKERR(err_msg)
        end if
    end function trim_modes
end module X_utilities
