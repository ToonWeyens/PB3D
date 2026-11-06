! gfortran14_bounds_repro.f90 -- standalone reproducer for a gfortran 14 bug
! that breaks PB3D's Debug build (issue #4, found on Snellius, 30 October 2026).
!
! With -fcheck=bounds (part of -fcheck=all), gfortran 14.2.0 and 14.3.0
! miscompile the bounds check of an array section whose subscript calls a
! function with an array-constructor argument, when that section is an actual
! argument:  eq_2%T_EF(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)).
! The check evaluates the subscript to garbage ("Index '2770510' of dimension
! 4 ... outside of expected range") or segfaults inside c, where the argument
! arrives at a wild address. PB3D's eq_ops.f90 has this pattern at every
! metric-tensor element; the first call, in calc_T_HF_ind, kills every Debug
! regression deck on Snellius (foss/2025a, GCC 14.2.0).
!
!   gfortran -O0 -g -fcheck=bounds gfortran14_bounds_repro.f90 && ./a.out
!
! 14.2.0, 14.3.0: runtime error or SIGSEGV.  13.3.0, 16.1.0: prints
!   ierr = 0  c([1,1]) = 1  sum = -126716.0
! Hoisting the call (id = c([1,1],.false.); ... T_EF(:,:,:,id,...)) avoids it;
! the generic interface, the negated section and the optional arguments of c
! are not needed to trigger it. CMakeLists.txt drops -fcheck=bounds from the
! Debug flags for gfortran 14.x because of this.
module m_idx
    implicit none
    private
    public c
contains
    integer function c(ij,sym,n,lim_n)
        integer, intent(in) :: ij(2)
        logical, intent(in) :: sym
        integer, intent(in), optional :: n
        integer, intent(in), optional :: lim_n(2,2)
        integer :: n_loc, js, ij_loc(2), kd
        n_loc = 3
        if (present(n)) n_loc = n
        if (sym) then
            ij_loc = ij
            if (present(lim_n)) ij_loc = ij + lim_n(1,:) - 1
            if (ij_loc(2).gt.ij_loc(1)) then
                kd = ij_loc(2); ij_loc(2) = ij_loc(1); ij_loc(1) = kd
            end if
            c = nint((ij_loc(2)-1)*n_loc + ij_loc(1) - (ij_loc(2)-1)*ij_loc(2)*0.5)
            if (present(lim_n)) then
                js = max(0,min(ij_loc(2)-lim_n(1,2)+1,lim_n(1,1)-lim_n(1,2)+1))
                c = c - nint((lim_n(1,2)-1)*(n_loc+1-lim_n(1,2)*0.5) &
                    &+js*(lim_n(1,1)-lim_n(1,2)+0.5-js*0.5) &
                    &+(n_loc-lim_n(2,1))*(ij_loc(2)-lim_n(1,2)))
            end if
        else
            if (present(lim_n)) then
                c = (ij(2)-1)*(lim_n(2,1)-lim_n(1,1)+1) + ij(1)
            else
                c = (ij(2)-1)*n_loc + ij(1)
            end if
        end if
    end function c
end module m_idx

module m_ops
    implicit none
    private
    public add_arr_mult
    interface add_arr_mult
        module procedure add_arr_mult_3_3
        module procedure add_arr_mult_3_1
    end interface
contains
    integer function add_arr_mult_3_3(arr_1,arr_2,arr_3,deriv) result(ierr)
        real(8), intent(in) :: arr_1(1:,1:,1:,0:,0:,0:)
        real(8), intent(in) :: arr_2(1:,1:,1:,0:,0:,0:)
        real(8), intent(out) :: arr_3(1:,1:,1:)
        integer, intent(in) :: deriv(3)
        ierr = 0
        arr_3 = arr_1(:,:,:,0,0,0)*arr_2(:,:,:,0,0,0) + deriv(1)
    end function add_arr_mult_3_3
    integer function add_arr_mult_3_1(arr_1,arr_2,arr_3,deriv) result(ierr)
        real(8), intent(in) :: arr_1(1:,1:,1:,0:,0:,0:)
        real(8), intent(in) :: arr_2(1:,0:)
        real(8), intent(out) :: arr_3(1:,1:,1:)
        integer, intent(in) :: deriv(3)
        integer :: kd
        ierr = 0
        do kd = 1,size(arr_3,3)
            arr_3(:,:,kd) = arr_1(:,:,kd,0,0,0)*arr_2(kd,0) + deriv(1)
        end do
    end function add_arr_mult_3_1
end module m_ops

module m_types
    implicit none
    type :: eq_1_type
        real(8), allocatable :: q_saf_E(:,:)
    end type
    type :: eq_2_type
        real(8), allocatable :: T_EF(:,:,:,:,:,:,:)
    end type
end module m_types

program repro
    use m_idx, only: c
    use m_ops, only: add_arr_mult
    use m_types
    implicit none
    type(eq_1_type) :: eq_1
    type(eq_2_type) :: eq_2
    real(8), allocatable :: theta_s(:,:,:,:,:,:)
    integer :: ierr, deriv(3), dims(3)
    deriv = 0
    dims = [401,1,158]
    allocate(eq_1%q_saf_E(dims(3),0:2)); eq_1%q_saf_E = 2.0d0
    allocate(eq_2%T_EF(dims(1),dims(2),dims(3),9,0:1,0:1,0:1)); eq_2%T_EF = 0.0d0
    allocate(theta_s(dims(1),dims(2),dims(3),0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1))
    theta_s = 1.0d0
    ierr = add_arr_mult(theta_s,-eq_1%q_saf_E(:,1:),&
        &eq_2%T_EF(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)),&
        &deriv)
    print *, 'ierr =', ierr, ' c([1,1]) =', c([1,1],.false.), ' sum =', sum(eq_2%T_EF)
end program repro
