!------------------------------------------------------------------------------!
!> Rank-agnostic helpers for the full-stack tests.
!!
!! The vacuum matrices G and H live in the 2-D block-cyclic ScaLAPACK
!! distribution described in vac_vars.init_vac (whole-matrix block for a
!! single process, blocksize-16 blocks otherwise), and the response vac%res
!! only on the last process. The tests must therefore never index them as if
!! they were global. Instead, these helpers evaluate everything the tests
!! assert on as *global* quantities, identically on every process:
!!
!!  - dis_matvec: y = A x for a distributed square matrix A and globally
!!    replicated x, through pdgemv on a properly distributed column vector,
!!    gathering y on all processes. This exercises the same distribution
!!    machinery (descriptors, lims_r/lims_c index bookkeeping, dgsum2d
!!    reductions) that the production solves use.
!!  - vec_glob2dis: the scatter half of the above, also used to build
!!    distributed right-hand sides for solve_Phi_BEM.
!!  - gather_res: broadcast vac%res from the last process to all.
!!
!! Because every process ends up asserting on the same numbers, test-drive's
!! per-test control flow cannot diverge between ranks, which would otherwise
!! deadlock collective calls in subsequent tests.
!------------------------------------------------------------------------------!
module fullstack_utils
    use num_vars, only: dp
    use vac_vars, only: vac_type, in_context, BLACSCTXTSIZE

    implicit none
    private
    public n_col1_loc, vec_glob2dis, dis_matvec, gather_res

contains
    !> Number of local columns this process owns of a global 1-column matrix:
    !! the single column sits in the first block, owned by process column 0.
    integer function n_col1_loc(vac) result(res)
        type(vac_type), intent(in) :: vac                                       !< vacuum variables

        if (in_context(vac%ctxt_HG) .and. vac%ind_p(2).eq.0) then
            res = 1
        else
            res = 0
        end if
    end function n_col1_loc

    !> Scatter a global column vector into its local block-cyclic rows.
    !!
    !! \c x_dis is the flattened local array of size n_loc(1) * n_col1_loc();
    !! processes that do not own the column receive nothing.
    subroutine vec_glob2dis(vac,x_glob,x_dis)
        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        real(dp), intent(in) :: x_glob(:)                                       !< global vector
        real(dp), intent(out) :: x_dis(:)                                       !< local part (flattened)

        integer :: i_rd, rd, rdl                                                ! subrow index, global and local row

        x_dis = 0._dp
        if (n_col1_loc(vac).eq.0) return
        do i_rd = 1,size(vac%lims_r,2)
            do rd = vac%lims_r(1,i_rd),vac%lims_r(2,i_rd)
                rdl = sum(vac%lims_r(2,1:i_rd-1)-vac%lims_r(1,1:i_rd-1)+1) + &
                    &rd-vac%lims_r(1,i_rd)+1
                x_dis(rdl) = x_glob(rd)
            end do
        end do
    end subroutine vec_glob2dis

    !> Distributed matrix-vector product y = A x.
    !!
    !! \c A is a local part of a distributed n_bnd x n_bnd matrix with
    !! descriptor \c desc_A (e.g. vac%H with vac%desc_H), \c x_glob is
    !! globally replicated, and \c y_glob is returned complete on all
    !! processes in the context (zero outside it).
    integer function dis_matvec(vac,A,desc_A,x_glob,y_glob) result(ierr)
        use vac_utilities, only: vec_dis2loc

        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        real(dp), intent(in) :: A(:,:)                                          !< local part of distributed matrix
        integer, intent(in) :: desc_A(BLACSCTXTSIZE)                            !< descriptor of A
        real(dp), intent(in) :: x_glob(:)                                       !< global input vector
        real(dp), intent(out) :: y_glob(:)                                      !< global result vector

        integer :: n_col1                                                       ! local number of columns of the vector
        integer :: n_flat                                                       ! true local size of the vector
        integer :: desc_x(BLACSCTXTSIZE)                                        ! descriptor for x and y
        real(dp), allocatable :: x_dis(:), y_dis(:)                             ! local vector parts

        ierr = 0
        y_glob = 0._dp
        if (.not.in_context(vac%ctxt_HG)) return

        n_col1 = n_col1_loc(vac)
        n_flat = vac%n_loc(1)*n_col1
        allocate(x_dis(max(1,n_flat)),y_dis(max(1,n_flat)))                     ! at least 1 element for ScaLAPACK dummies
        call vec_glob2dis(vac,x_glob,x_dis)
        y_dis = 0._dp

        call descinit(desc_x,vac%n_bnd,1,vac%bs,vac%bs,0,0,vac%ctxt_HG,&
            &max(1,vac%n_loc(1)),ierr)
        if (ierr.ne.0) return

        call pdgemv('N',vac%n_bnd,vac%n_bnd,1._dp,A,1,1,desc_A,&
            &x_dis,1,1,desc_x,1,0._dp,y_dis,1,1,desc_x,1)

        ierr = vec_dis2loc(vac%ctxt_HG,y_dis(1:n_flat),vac%lims_r,y_glob)       ! all processes receive
    end function dis_matvec

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
