!-------------------------------------------------------
!   Variables,  subroutines and  functions that  have  to do  with the  fourier
!   representation used and its relations with the real space
!-------------------------------------------------------
module fourier_ops
    use output_ops, only: writo, print_ar_1, print_ar_2
    use num_vars, only: dp

    implicit none
    private
    public repack, mesh_cs, f2r

contains
    ! Inverse Fourier transformation, VMEC style
    ! Also calculates  the poloidal and toroidal  derivatives. Normal derivative
    ! is done discretely, outside of this function
    function f2r(fun_cos,fun_sin,ang_factor,mpol,ntor)
        integer, intent(in) :: mpol, ntor
        real(dp), allocatable :: f2r(:)
        real(dp), intent(in) :: fun_cos(0:mpol-1,-ntor:ntor)                    ! cos part of Fourier variables (coeff. of the sum)
        real(dp), intent(in) :: fun_sin(0:mpol-1,-ntor:ntor)                    ! sin part of Fourier variables (coeff. of the sum)
        real(dp), intent(in) :: ang_factor(0:mpol-1,-ntor:ntor,2)               ! (co)sine factors on mesh
        
        integer :: m,n
        
        ! some tests
        if (mpol.lt.1 .and. ntor.lt. 1) then 
            call writo('ERROR: number of modes has to be at least 1')
            stop
        end if
        
        ! allocate and initiate
        allocate(f2r(4))
        f2r = 0.0_dp
        
        ! sum over all poloidal and toroidal modes
        do n = -ntor,ntor
            do m = 0,mpol-1
                f2r(1) = f2r(1) + fun_cos(m,n)*ang_factor(m,n,1) &              ! the variable itself
                    &+ fun_sin(m,n)*ang_factor(m,n,2)
                f2r(3) = f2r(3) + m * (-fun_cos(m,n)*ang_factor(m,n,2) &        ! theta derivative: m(-f_c*s+f_s*c)
                    &+ fun_sin(m,n)*ang_factor(m,n,1))                          
                f2r(4) = f2r(4) + n * (-fun_cos(m,n)*ang_factor(m,n,2) &        ! zeta derivative: n(-f_c*s+f_s*c)
                    &+ fun_sin(m,n)*ang_factor(m,n,1))
            end do
        end do
    end function f2r
 
    ! Calculate the cosine and sine factors  on a mesh (0:mpol-1, -ntor:ntor) at
    ! a given poloidal and toroidal position (theta,zeta)
    ! The first index contains the cosine factors and the second one the sines.
    function mesh_cs(mpol,ntor,theta,zeta)
        ! input / output
        integer, intent(in) :: mpol, ntor
        real(dp), allocatable :: mesh_cs(:,:,:)
        real(dp) :: theta, zeta
        
        ! local variables
        integer :: m, n
        
        ! test the given inputs
        if (mpol.lt.1) then
            call writo('ERROR: mpol has to be at least 1')
            stop
        end if
        
        allocate(mesh_cs(0:mpol-1,-ntor:ntor,2))
        mesh_cs = 0.0_dp
        
        do n = -ntor,ntor
            do m = 0,mpol-1
                ! cos factor
                mesh_cs(m,n,1) = cos(m*theta + n*zeta)
                ! sin factor
                mesh_cs(m,n,2) = sin(m*theta + n*zeta)
            end do
        end do
    end function mesh_cs

    ! Repack variables representing the Fourier  composition of R, Z and lambda.
    ! In VMEC these are stored as (1:mnmax, 1:ns) with mnmax the total number of
    ! all modes.  Here they are  to be  stored as (0:mpol-1,  -ntor:ntor, 1:ns),
    ! which is valid due to the symmetry  of the modes (only one of either theta
    ! or zeta has  to be able to  change sign because of  the (anti)-symmetry of
    ! the (co)sine.
    ! It is  possible that the  input variable is  not allocated. In  this case,
    ! output zero's
    function repack(var_VMEC,mnmax,ns,mpol,ntor,xm,xn)
        integer, intent(in) :: mnmax, ns, mpol, ntor
        real(dp), intent(in) :: xm(mnmax), xn(mnmax)
        real(dp), allocatable :: var_VMEC(:,:)
            
        integer :: mode, m, n
            
        real(dp) :: repack(0:mpol-1,-ntor:ntor,1:ns)
            
        if (allocated(var_VMEC)) then
            ! check if the  values in xm and  xn don't exceed the  maximum number of
            ! poloidal and toroidal modes (xm and xn are of length mnmax and contain
            ! the pol/tor mode number)
            if (maxval(xm).gt.mpol .or. maxval(abs(xn)).gt.ntor) then
                call writo('WARNING: In repack, less modes are used than in the &
                    &VMEC format')
            end if
                
            repack = 0.0_dp
            ! copy the VMEC modes using the PB3D format
            do mode = 1,mnmax
                m = nint(xm(mode))
                n = nint(xn(mode))
                ! if the modes don't fit, cycle
                if (m.gt.mpol .or. abs(n).gt.ntor) cycle
                
                repack(m,n,:) = var_VMEC(mode,:)
            end do
        else
            repack = 0.0_dp
        end if
    end function repack
end module fourier_ops
