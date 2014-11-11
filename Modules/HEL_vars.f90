!------------------------------------------------------------------------------!
!   functions and routines that concern variables given by HELENA              !
!------------------------------------------------------------------------------!
module HEL_vars
#include <PB3D_macros.h>
    use message_ops, only: writo, lvl_ud
    
    implicit none
    private
    public read_HEL, dealloc_HEL

contains
    ! Reads the HELENA equilibrium data
    ! (from HELENA routine IODSK)
    ! [MPI] only global master
    integer function read_HEL() result(ierr)
        use num_vars, only: glb_rank, eq_name
        
        character(*), parameter :: rout_name = 'read_HEL'
        
        ! local variables
        integer :: id                                                           ! counter
        real :: scaleq, rbphi02, dj0, dje, dp0, dpe, drbphi0, drbphie, wurzel, &
            &dqec, raxis
        real :: dummy(3)
        integer  js0
        
        ! initialize ierr
        ierr = 0
        
        if (glb_rank.eq.0) then                                                 ! only global master
            call writo('Reading data from HELENA output "' &
                &// trim(eq_name) // '"')
            call lvl_ud(1)
            
            !*CALL COMMAX
            !*CALL COMPIO
            !*CALL COMGEM
            !*CALL COMEQUI
            !*CALL COMIOD
            !*CALL COMSPL
            !*CALL COMBND
            !*CALL COMVAC
            
            ! rewind input file
            !rewind(eq_i)
            
            ! read mapped equilibrium from disk
            !read(eq_i,*) js0, (cs(js),js=1,js0+1), &
                !&(qs(js),js=1,js0+1), dqs(1), dqec, (dqs(js),js=2,js0+1), &
                !&(curj(js),js=1,js0+1), dj0, dje, &
                !&nchi, (chi(jc),jc=1,nchi), &
                !&(gem11(j), j=nchi+1, (js0+1)*nchi), &
                !&(gem12(j), j=nchi+1, (js0+1)*nchi), &
                !&cpsurf, radius
            
            !nvchi = nchi
            !npsi = js0 + 1
            !ng = npsi*nchi
            
            !do id = 1,ng
                !gem22(id) = 1.0
            !end do
            
            !! for toroidal geometry additional input gem22(j) = r**2
            !if (nltore) read(eq_i,*) (gem22(j),j=nchi+1,ng), raxis
            
            !read(eq_i,*) (p0(js),js=1,npsi), dp0, dpe, &
                !&(rbphi(js),js=1,npsi), drbphi0, drbphie
            
            !! read additional boundary data
            !if (rwall.gt.1.) then
                !read(eq_i,*) (vx(js),js=1,nchi)
                !read(eq_i,*) (vy(js),js=1,nchi)
                !read(eq_i,*) aspi
            !end if
            
            !do jc = 1,nchi
                !gem11(jc) = 0.0
            !end do
            
            !if (nltore) then
                !do jc = 1,nchi
                    !gem22(jc) = raxis**2
                !end do
            !end if
            
            !do jc = 1,nchi
                !call dcopy(npsi-1,gem12(nchi+jc),nchi,c1,1)
                !call spline(npsi-1,cs(2),c1,0.0,0.0,3,q1,q2,q3,q4)
                !gem12(jc) = spwert(npsi-1,0.0,q1,q2,q3,q4,cs(2),dummy)
            !end do
            
            !! scale to q on axis
            !if (q0zyl.gt.0.0) then
                !scaleq = qs(1)/q0zyl
                
                !cpsurf = cpsurf*scaleq
                !call dscal(npsi,scaleq,curj,1)
                !call dscal(npsi,scaleq**2,p0,1)
                !call dscal(npsi*nchi,-scaleq,gem12,1)
                !call dscal(npsi*nchi,scaleq**2,gem11,1)
                
                !rbphi02 = rbphi(1)**2
                
                !dj0 = dj0*scaleq
                !dje = dje*scaleq
                !dp0 = dp0*scaleq**2
                !dpe = dpe*scaleq**2
                !drbphi0 = drbphi0*scaleq**2
                !drbphie = drbphie*scaleq/sqrt(1.+rbphi02/rbphi(npsi)**2 *&
                    !&(1./scaleq**2-1.))
                
                !do id = 1,npsi
                    !wurzel   = sqrt(1.+rbphi02/rbphi(id)**2*(1./scaleq**2-1.))
                    !qs(id)    = wurzel * qs(id)
                    !rbphi(id) = wurzel * scaleq * rbphi(id)
                !end do
            !else
                !call dscal(npsi*nchi,-1.,gem12,1)
            !endif
            
            !51 FORMAT(//' AFTER Q ON AXIS SCALING : SCALE FACTOR =',1P,E12.4,0P)
            !52 FORMAT(/' CPSURF = ',1P,E12.4,0P)
            !53 FORMAT(/' QS'/(1X,1P,6E12.4,0P))
            !54 FORMAT(/' P0'/(1X,1P,6E12.4,0P))
            !55 FORMAT(/' RBPHI'/(1X,1P,6E12.4,0P))
            !56 FORMAT(/' DJ0,DJE',1P,2E12.4,0P/' CURJ'/(1X,1P,5E12.4,0P))
            
            call lvl_ud(-1)
            call writo('Grid parameters successfully read')
        end if
    end function read_HEL
    
    ! deallocates HELENA quantities that are not used anymore
    subroutine dealloc_HEL
    end subroutine dealloc_HEL
end module HEL_vars
