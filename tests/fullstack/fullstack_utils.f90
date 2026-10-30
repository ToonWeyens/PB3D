!------------------------------------------------------------------------------!
!> Rank-agnostic helpers for the full-stack tests.
!!
!! The vacuum matrices G and H live in the 2-D block-cyclic ScaLAPACK
!! distribution described in vac_vars.init_vac and the response vac%res only
!! on the last process, so the tests never index them as if they were
!! global: they assert on quantities evaluated identically on every process,
!! through vac_utilities (n_col1_loc, vec_glob2dis, dis_matvec) and
!! gather_res below. Because every process asserts on the same numbers,
!! test-drive's control flow cannot diverge between ranks, which would
!! otherwise deadlock collective calls in subsequent tests.
!------------------------------------------------------------------------------!
module fullstack_utils
    use num_vars, only: dp
    use vac_vars, only: vac_type
    implicit none
    private
    public gather_res

contains
    !> Broadcast the vacuum response from the last process (the only one
    !! where calc_vac_res leaves it) to all.
    integer function gather_res(vac,res) result(ierr)
        use MPI
        use num_vars, only: rank, n_procs

        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        complex(dp), intent(inout) :: res(:,:)                                  !< global response on all processes

        ierr = 0
        if (rank.eq.n_procs-1) res = vac%res
        call MPI_Bcast(res,size(res),MPI_DOUBLE_COMPLEX,n_procs-1,&
            &MPI_COMM_WORLD,ierr)
    end function gather_res
end module fullstack_utils
