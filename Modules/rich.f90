!------------------------------------------------------------------------------!
!   Operations and variables concerning Richardson extrapolation               !
!------------------------------------------------------------------------------!
module rich
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln, max_it_rich, tol_rich

    implicit none
    private
    public init_rich, term_rich, start_rich_lvl, stop_rich_lvl, do_rich, &
        &rich_info, rich_info_short, calc_rich_ex, &
        &rich_lvl, no_guess, use_guess, n_r_sol, min_n_r_sol
    
    ! global variables
    integer :: rich_lvl                                                         ! current level of Richardson extrapolation
    integer :: n_r_sol                                                          ! nr. of normal points in sol grid for Richardson lvl
    integer :: min_n_r_sol                                                      ! min. of n_r_sol (e.g. first value in Richardson loop)
    logical :: use_guess                                                        ! whether a guess is formed from previous level of Richardson
    logical :: no_guess                                                         ! disable guessing Eigenfunction from previous level of Richardson
    logical :: rich_conv                                                        ! if Richarson extrapolation has converged
    real(dp), allocatable :: x_axis(:,:)                                        ! x axis for plot of Eigenvalues with Richardson level
    complex(dp), allocatable :: X_val_rich(:,:,:)                               ! Richardson array of eigenvalues
    real(dp), allocatable :: max_rel_error(:)                                   ! maximum relative error for all Richardson levels
    integer, allocatable :: loc_max_rel_err(:,:)                                ! location of maximum of relative error
    
contains
    subroutine init_rich
        use num_vars, only: n_sol_requested
        
        ! set variables
        rich_lvl = 1
        rich_conv = .false.
        no_guess = .true.
        
        if (max_it_rich.gt.1) then                                              ! only when more than one level
            ! allocate variables
            allocate(X_val_rich(1:max_it_rich,1:max_it_rich,1:n_sol_requested))
            X_val_rich = 0.0_dp
            allocate(x_axis(1:max_it_rich,1:n_sol_requested))
            allocate(max_rel_error(1:max_it_rich-1))
            allocate(loc_max_rel_err(1:max_it_rich-1,1))
            
            ! user output
            call writo('Starting Richardson extrapolation loop')
            call lvl_ud(1)
            
            call writo('Maximum number of iterations: '//&
                &trim(i2str(max_it_rich)))
            call writo('Tolerance requested: '//trim(r2str(tol_rich)))
            
            call writo('')
            call lvl_ud(-1)
        end if
    end subroutine init_rich
    
    subroutine term_rich
        use num_vars, only: rank,norm_disc_prec_sol
        
        ! local variables
        integer :: id                                                           ! counter
        
        if (max_it_rich.gt.1) then                                              ! only when more than one level
            ! user output
            call writo('Finishing Richardson extrapolation loop')
            call lvl_ud(1)
            
            ! user output
            if (rich_conv) then
                call writo('Convergence reached in '//trim(i2str(rich_lvl-1))//&
                    &' steps:')
            else
                call writo('After '//trim(i2str(rich_lvl-1))//&
                    &' steps, no convergence reached:')
            end if
            call lvl_ud(1)
            
            call writo('Maximum relative error '//&
                &trim(r2str(max_rel_error(rich_lvl-2)))//&
                &' for Eigenvalue '//trim(i2str(loc_max_rel_err(rich_lvl-2,1))))
            call writo('Tolerance requested: '//trim(r2str(tol_rich)))
            call writo('Resulting best guesses for Eigenvalues:')
            call lvl_ud(1)
            do id = 1,size(X_val_rich,3)
                call writo('For Eigenvalue '//trim(i2str(id))//': '//&
                    &trim(c2str(X_val_rich(rich_lvl-1,rich_lvl-1,id)))//',')
            end do
            call lvl_ud(-1)
            call writo('with, theoretically, an error of the order O(Δ^'//&
                &trim(i2str(2*norm_disc_prec_sol*rich_lvl))//'),')
            call writo('compared to O(Δ^'//trim(i2str(2*norm_disc_prec_sol))//&
                &') without Richardson extrapolation')
            call writo('Notice that this is ONLY valid if the problem in the &
                &different Richardson stages is similar:')
            call lvl_ud(1)
            call writo('- It has to be part of the point spectrum')
            call writo('- The Eigenvectors should look similar')
            call lvl_ud(-1)
            
            call lvl_ud(-1)
            
            ! draw X values
            if (rank.eq.0) call draw_X_val_rich()
            
            call writo('')
            call lvl_ud(-1)
        end if
    contains
        ! Draws  the   Eigenvalues  for  the  different   levels  of  Richardson
        ! extrapolation as a function of the number of normal points.
        subroutine draw_X_val_rich()
            use num_vars, only: n_sol_requested
            
            ! local variables
            character(len=max_str_ln) :: plot_title                             ! title for plots
            character(len=max_str_ln) :: plot_name                              ! name of plot
            
            ! user output
            call writo('Plotting Eigenvalues as function of nr. of normal &
                &points')
            call lvl_ud(1)
            
            ! print, using rich_lvl-1 as it has been aumented
            plot_title = 'Eigenvalues as function of nr. of normal points'
            plot_name = 'Eigenvalues_richardson'
            call print_GP_2D(plot_title,plot_name,&
                &realpart(X_val_rich(1:rich_lvl-1,1,:)),&
                &x=x_axis(1:rich_lvl-1,:),&
                &draw=.false.)
            
            ! output in file
            call draw_GP(plot_title,plot_name,plot_name,&
                &n_sol_requested,1,.false.)
            
            ! user output
            call lvl_ud(-1)
            call writo('Done plotting Eigenvalues')
        end subroutine draw_X_val_rich
    end subroutine term_rich
    
    ! if this Richardson level should be done
    logical function do_rich()
        if (rich_lvl.eq.1) then                                                 ! first level is always done
            do_rich = .true.
        else if (rich_lvl.gt.max_it_rich) then                                  ! maximum level reached
            do_rich = .false.
        else if (rich_conv) then                                                ! not yet converged
            do_rich = .false.
        else
            do_rich = .true.
        end if
    end function do_rich
    
    subroutine start_rich_lvl
        ! Calculate number of normal points for the solution in Richardson loops
        if (rich_lvl.eq.1) then
            n_r_sol = min_n_r_sol
        else
            n_r_sol = 2 * n_r_sol - 1
        end if
        
        ! set use_guess to .false. if user sets no_guess
        if (no_guess) use_guess = .false.
    end subroutine start_rich_lvl
    
    subroutine stop_rich_lvl
        ! Richardson extrapolation
        if (max_it_rich.gt.1) then                                              ! only do this if more than 1 Richardson level
            write(*,*) '!!! FOR RICHARDSON EXTRAPOLATION YOU HAVE TO FIND &
                &CONSTANT EIGENVALUES !!!'
            
            ! user output
            if (rich_lvl.gt.1) call writo('Richardson level '//&
                &trim(i2str(rich_lvl))//' summary')
            call lvl_ud(1)
            
            ! update the x axis of the Eigenvalue plot
            x_axis(rich_lvl,:) = 1.0_dp*n_r_sol
            
            ! decide whether converged or not
            call check_conv()
            
            ! setup possible guess for next Richardson level
            call set_guess()
            
            ! increase level
            rich_lvl = rich_lvl + 1
            
            call writo('')
            call lvl_ud(-1)
        else                                                                    ! if not, Richardson is done
            rich_conv = .true.
        end if
    contains
        ! Decides  whether  the  difference  between the  approximation  of  the
        ! Eigenvalues  in  this Richardson  level  and  the previous  Richardson
        ! level, with the same order of the error, falls below a threshold
        ! [ADD SOURCE]
        subroutine check_conv()
            ! local variables
            real(dp), allocatable :: corr(:)                                    ! correction for order (rich_lvl-1) of error
            integer :: ir                                                       ! counter
            
            if (rich_lvl.gt.1) then                                             ! only do this if in Richardson level higher than 1
                ! set relative correction
                allocate(corr(size(X_val_rich,3)))
                corr = abs(2*(X_val_rich(rich_lvl,rich_lvl-1,:)-&
                    &X_val_rich(rich_lvl-1,rich_lvl-1,:))/&
                    &(X_val_rich(rich_lvl,rich_lvl-1,:)+&
                    &X_val_rich(rich_lvl-1,rich_lvl-1,:)))
                
                ! get maximum and location of maximum for relative correction
                max_rel_error(rich_lvl-1) = maxval(corr)
                loc_max_rel_err(rich_lvl-1,:) = maxloc(corr)
                
                ! user output
                do ir = 1,rich_lvl-1
                    call writo('Richardson level '//trim(i2str(ir))//' -> '//&
                        &trim(i2str(ir+1))//': maximum relative error '//&
                        &trim(r2strt(max_rel_error(ir)))//' for Eigenvalue '//&
                        &trim(i2str(loc_max_rel_err(ir,1))))
                end do
                
                ! check whether tolerance has been reached
                if (max_rel_error(rich_lvl-1).lt.tol_rich) then
                    rich_conv = .true.
                    call writo('tolerance '//trim(r2strt(tol_rich))//&
                        &' reached after '//trim(i2str(rich_lvl))//&
                        &' iterations')
                else
                    call writo('tolerance '//trim(r2strt(tol_rich))//' not yet &
                        &reached')
                end if
            end if
        end subroutine check_conv
        
        ! Decides whether  a guess should be used in  a possible next Richardson
        ! level.
        subroutine set_guess()
            ! decide on guess
            if (rich_lvl.ge.3) then                                             ! start from level 3
                use_guess = max_rel_error(rich_lvl-1).lt.&
                    &max_rel_error(rich_lvl-2)                                  ! set guess for next level to true if max_rel_error is decreasing
            else
                use_guess = .true.                                              ! for first Richardson level, set guess for next level to true
            end if
        end subroutine set_guess
    end subroutine stop_rich_lvl
    
    ! Calculates  the   coefficients  of  the  Eigenvalues   in  the  Richardson
    ! extrapolation. This is done using the recursive formula:
    !   X_val_rich(lvl,ir,:) = X_val_rich(lvl,ir-1,:) +  1/(2^(2 p ir) - 1) * 
    !       (X_val_rich(lvl,ir-1,:) - X_val_rich(lvl-1,ir-1,:)),
    ! where p is norm_disc_prec_sol, the order of the numerical discretization.
    ! as described in [ADD REF]
    integer function calc_rich_ex(X_val) result(ierr)
        use num_vars, only: norm_disc_prec_sol
        
        character(*), parameter :: rout_name = 'calc_rich_ex'
        
        ! input / output
        complex(dp), intent(in) :: X_val(:)                                     ! EV for this Richardson level
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: ir                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(X_val_rich,1).ne.size(X_val_rich,2) .or. &
            &rich_lvl.gt.size(X_val_rich,1)) then
            ierr = 1
            err_msg = 'X_val_rich has to have correct dimensions'
            CHCKERR(err_msg)
        end if
        
        ! do calculations if rich_lvl > 1
        X_val_rich(rich_lvl,1,1:size(X_val)) = X_val                            ! size(X_val) can be less than n_sol_requested
        do ir = 2,rich_lvl
            X_val_rich(rich_lvl,ir,1:size(X_val)) = &
                &X_val_rich(rich_lvl,ir-1,1:size(X_val)) + &
                &1._dp/(2**(2*(ir-1)*norm_disc_prec_sol)-1._dp) * &
                &(X_val_rich(rich_lvl,ir-1,1:size(X_val)) - &
                &X_val_rich(rich_lvl-1,ir-1,1:size(X_val)))
        end do
    end function calc_rich_ex
    
    ! possible extension with Richardson level or nothing if only one level
    elemental character(len=max_str_ln) function rich_info()                    ! full version
        if (max_it_rich.gt.1) then
            rich_info = ' for Richardson level '//trim(i2str(rich_lvl))
        else
            rich_info = ''
        end if
    end function rich_info
    elemental character(len=max_str_ln) function rich_info_short()              ! short version
        if (max_it_rich.gt.1) then
            rich_info_short = '_R_'//trim(i2str(rich_lvl))
        else
            rich_info_short = ''
        end if
    end function rich_info_short
end module rich
