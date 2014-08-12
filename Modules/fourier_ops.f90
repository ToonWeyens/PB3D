!------------------------------------------------------------------------------!
!   Variables, subroutines and  functions that  have  to do  with the  fourier !
!   representation used and its relations with the real space                  !
!------------------------------------------------------------------------------!
module fourier_ops
#include <PB3D_macros.h>
    use output_ops, only: writo, print_ar_1, print_ar_2
    use num_vars, only: dp, max_str_ln
    use str_ops, only: i2str

    implicit none
    private
    public repack, calc_mesh_cs, f2r

contains
    ! Inverse Fourier transformation, VMEC style Also calculates the poloidal or
    ! toroidal derivatives, as indicated by the variable deriv(2)
    ! (Normal derivative is done discretely, outside of this function)
    function f2r(fun_cos,fun_sin,ang_factor,mpol,ntor,nfp,deriv,ierr)
        character(*), parameter :: rout_name = 'f2r'
        
        ! input / output
        integer, intent(in) :: mpol, ntor
        integer, intent(in), optional :: deriv(2)
        real(dp) :: f2r
        real(dp), intent(in) :: fun_cos(0:mpol-1,-ntor:ntor)                    ! cos part of Fourier variables (coeff. of the sum)
        real(dp), intent(in) :: fun_sin(0:mpol-1,-ntor:ntor)                    ! sin part of Fourier variables (coeff. of the sum)
        real(dp), intent(in) :: ang_factor(0:mpol-1,-ntor:ntor,2)               ! (co)sine factors on mesh
        integer, intent(in) :: nfp                                              ! common denominator nfp in toroidal mode numbers
        integer, intent(inout) :: ierr                                          ! error
        
        ! local variables
        integer :: m,n                                                          ! counters for mode numbers
        integer :: id                                                           ! counters
        real(dp) :: fac_cos, fac_sin                                            ! factor in front of cos and sin, after taking derivatives
        real(dp) :: fac_cos_temp, fac_sin_temp                                  ! when calculating factors for angular derivatives
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! some tests
        if (mpol.lt.1 .and. ntor.lt. 1) then 
            err_msg = 'The number of modes has to be at least 1'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! initiate
        f2r = 0.0_dp
        
        ! sum over modes
        do n = -ntor,ntor
            do m = 0,mpol-1
                ! initialize factors in front of cos and sin
                fac_cos = fun_cos(m,n)
                fac_sin = fun_sin(m,n)
                
                ! angular derivatives
                if (present(deriv)) then 
                    ! apply possible poloidal derivatives
                    do id = 1,deriv(1)
                        fac_cos_temp = m * fac_sin
                        fac_sin_temp = - m * fac_cos
                        fac_cos = fac_cos_temp
                        fac_sin = fac_sin_temp
                    end do
                    ! apply possible toroidal derivatives
                    do id = 1,deriv(2)
                        fac_cos_temp = - n * nfp * fac_sin
                        fac_sin_temp = n * nfp * fac_cos
                        fac_cos = fac_cos_temp
                        fac_sin = fac_sin_temp
                    end do
                end if
                
                ! orthonormalization has been taken  care of by considering only
                ! half the summation in the poloidal variable
                f2r = f2r + fac_cos*ang_factor(m,n,1) &
                    &+ fac_sin*ang_factor(m,n,2)
            end do
        end do
    end function f2r
 
    ! Calculate the cosine and sine factors  on a mesh (0:mpol-1, -ntor:ntor) at
    ! a given poloidal and toroidal position (theta,zeta)
    ! The first index contains the cosine  factors and the second one the sines.
    ! CHANGE THIS USING THE GONIOMETRIC IDENTITIES TO SAVE COMPUTING TIME
    integer function calc_mesh_cs(mesh_cs,mpol,ntor,nfp,theta,zeta) result(ierr)
        character(*), parameter :: rout_name = 'mesh_cs'
        
        ! input / output
        integer, intent(in) :: mpol, ntor
        real(dp), intent(inout) :: mesh_cs(0:mpol-1,-ntor:ntor,2)
        real(dp), intent(in) :: theta, zeta
        integer, intent(in) :: nfp
        
        ! local variables
        integer :: m, n
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! test the given inputs
        if (mpol.lt.1) then
            err_msg = 'mpol has to be at least 1'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        mesh_cs = 0.0_dp
        
        do n = -ntor,ntor
            do m = 0,mpol-1
                ! cos factor
                mesh_cs(m,n,1) = cos(m*theta - n*nfp*zeta)
                ! sin factor
                mesh_cs(m,n,2) = sin(m*theta - n*nfp*zeta)
            end do
        end do
    end function calc_mesh_cs

    ! Repack  variables  representing the  Fourier  composition  such as  R,  Z,
    ! lambda, ...  In VMEC these  are stored as  (1:mnmax, 1:ns) with  mnmax the
    ! total  number of  all modes.  Here  they are  to be  stored as  (0:mpol-1,
    ! -ntor:ntor, 1:ns), which  is valid due to the symmetry  of the modes (only
    ! one of either theta  or zeta has to be able to change  sign because of the
    ! (anti)-symmetry of the (co)sine.
    ! It is  possible that the  input variable is  not allocated. In  this case,
    ! output zeros
    ! [MPI] only global master
    !       (this is a precaution: only the global master should use it)
    function repack(var_VMEC,mnmax,n_r,mpol,ntor,xm,xn)
        use num_vars, only: glob_rank
        
        integer, intent(in) :: mnmax, n_r, mpol, ntor
        real(dp), intent(in) :: xm(mnmax), xn(mnmax)
        real(dp), allocatable :: var_VMEC(:,:)
            
        integer :: mode, m, n
            
        real(dp) :: repack(0:mpol-1,-ntor:ntor,1:n_r)
            
        if (allocated(var_VMEC) .and. glob_rank.eq.0) then                      ! only global rank
            ! check if the  values in xm and xn don't  exceed the maximum number
            ! of poloidal and toroidal modes (xm  and xn are of length mnmax and
            ! contain the pol/tor mode number)
            if (maxval(xm).gt.mpol .or. maxval(abs(xn)).gt.ntor) then
                call writo('WARNING: In repack, less modes are used than in the&
                    & VMEC format')
            end if
                
            repack = 0.0_dp
            ! copy the VMEC modes using the PB3D format
            do mode = 1,mnmax
                m = nint(xm(mode))
                n = nint(xn(mode))
                ! if the modes don't fit, cycle
                if (m.gt.mpol .or. abs(n).gt.ntor) then
                    call writo('WARNING: In f2r, m > mpol or n > ntor!')
                    cycle
                end if
                
                repack(m,n,:) = var_VMEC(mode,:)
            end do
        else
            repack = 0.0_dp
        end if
    end function repack
end module fourier_ops
