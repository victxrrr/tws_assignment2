program Q3
    use Utils

    implicit none

    character(len=32)     :: buffer
    character(len=256)    :: MSG
    integer               :: N, N0, STATUS = 0, i, m
    integer, parameter    :: n_meas = 10
    logical               :: measurements = .true.
    real(wp)              :: T, beta, beta0, diff, target, table(16,4), start, end
    real(wp), allocatable :: x(:,:), tk(:)

    call utils_set_save(.false.)

    if (command_argument_count() /= 4) then 
        print *, 'Wrong number of arguments'
        stop STATUS
    endif

    call get_command_argument(1, buffer)
    read(buffer,*) N
    N0 = N
    call get_command_argument(2, buffer)
    read(buffer,*) T
    call get_command_argument(3, buffer)
    read(buffer,*) diff
    call get_command_argument(4, buffer)
    read(buffer,*) target

    allocate(x(N+1,5), STAT=STATUS, ERRMSG=MSG)
    call memcheck(STATUS,MSG)
    allocate(tk(N+1), STAT=STATUS, ERRMSG=MSG)
    call memcheck(STATUS,MSG)

    x(1,:) = initialize('input/parameters.in', beta_opt=beta0)

    if (measurements) then 
        table(:,3) = 0._wp
        do i = 1, 16
            diff = 1.0/10.0**i
            do m = 1, n_meas
                beta = beta0
                N = N0
                call cpu_time(start)
                call NewtonMethod(T, N, tk, x, 'Heun', beta)
                call cpu_time(end)
                table(i,3) = table(i,3) + (end - start)/n_meas
            enddo
            table(i,1) = diff
            table(i,2) = beta
            table(i,4) = N
        enddo
        call display(table)
    else 
        call NewtonMethod(T, N, tk, x, 'Heun', beta0)
        print '(a,i3,a,es22.15)', '>>> Optimal beta found in ', N, ' iterations :', beta0
    endif

    deallocate(x, STAT=STATUS, ERRMSG=MSG)
    call memcheck(STATUS,MSG)
    deallocate(tk, STAT=STATUS, ERRMSG=MSG)
    call memcheck(STATUS,MSG)

contains

    subroutine memcheck(STATUS,MSG)
        integer            :: STATUS
        character(len=256) :: MSG
        if (STATUS /= 0) then
            print *, MSG
            stop STATUS
        endif
    end subroutine

    !> Implementation of Newton's method to solve F(beta) - target = 0
    !! @param T simulation horizon
    !! @param N [in] number of discretization steps
    !!          [out] number of required iterations to converge (overwritten)
    !! @param tk array containing the N + 1 equidistant grid points
    !! @param x 2D-array containing the solution approximation at time 0, T/N, 2T/N, ... , T
    !! @param beta [in] the initial guess for beta
    !!             [out] the optimal value of beta computed (overwritten)
    !! @param tol_opt [optional] tolerance on the backward error for the stopping criterion
    !! @param nmax_opt [optional] maximum number of iterations that can be performed
    subroutine NewtonMethod(T, N, tk, x, method, beta, tol_opt, nmax_opt)
        real(wp), intent(in)           :: T
        integer, intent(inout)         :: N
        real(wp), intent(out)          :: tk(:), x(:,:)
        character(*), intent(in)       :: method
        real(wp), intent(inout)        :: beta
        real(wp), intent(in), optional :: tol_opt
        integer, intent(in), optional  :: nmax_opt
        real(wp)                       :: v0, v1, tol, x0(5)
        integer                        :: j, nmax

        if (present(tol_opt)) then
            tol = tol_opt
        else 
            tol = 1e-16
        endif 

        if (present(nmax_opt)) then
            nmax = nmax_opt
        else
            nmax = 1000
        endif

        x0 = x(1,:)
        do j = 1, nmax
            x(1,:) = x0
            v0 = F(T, N, tk, x, method, beta) - target
            x(1,:) = x0
            v1 = F(T, N, tk, x, method, beta + diff) - target
            beta = beta + diff/(1 - v1/v0)
            
            if (abs(diff/(1 - v1/v0)) < tol*abs(beta)) exit 
        enddo
        
        N = j - 1
        if (N == nmax) print *, "ERROR: failed to converge, maximum number of iterations reached."
    end subroutine

    function F(T, N, tk, x, method, beta) result(max)
        real(wp), intent(in)     :: T, beta
        integer, intent(in)      :: N
        real(wp), intent(out)    :: tk(:), x(:,:)
        character(*), intent(in) :: method
        real(wp)                 :: max
        integer                  :: k

        select case(method)
        case('EulerForward')
            call EulerForward(T, N, tk, x, beta_opt=beta)
        case('Heun')
            call Heun(T, N, tk, x, beta_opt=beta)
        case('EulerBackward')
            call EulerBackward(T, N, tk, x, beta_opt=beta)
        case default
            print *, 'ERROR: ', method, ' is not a valid method.'
            stop
        end select

        max = 0._wp
        do k = 1,size(x,1)
            if (x(k,2) + x(k,3) > max) then
                max = x(k,2) + x(k,3)
            endif
        enddo
    end function

    subroutine display(x) 
        real(wp), intent(in) :: x(:,:)
        integer              :: l

        do l=1,size(x,1)
            print '(es12.0, es25.15, f9.2, i5)', x(l,1), x(l,2), x(l,3)*1000, int(x(l,4))
        enddo
    end subroutine
        
end program
    