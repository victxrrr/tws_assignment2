program Q2
    use Utils

    implicit none

    character(len=32)     :: Nin, Tin
    character(len=256)    :: MSG
    integer               :: N, STATUS = 0
    real(wp)              :: T
    real(wp), allocatable :: x(:,:), tk(:)


    if (command_argument_count() /= 2) then 
        print *, 'Wrong number of arguments'
        stop STATUS
    endif

    call get_command_argument(1, Nin)
    call get_command_argument(2, Tin)
    read(Nin,*) N
    read(Tin,*) T

    allocate(x(N+1,5), STAT=STATUS, ERRMSG=MSG)
    call memcheck(STATUS,MSG)

    allocate(tk(N+1), STAT=STATUS, ERRMSG=MSG)
    call memcheck(STATUS,MSG)

    x(1,:) = initialize('input/parameters.in')
    call EulerForward(T, N, tk, x)
    call PrintAndSave(tk, x, 'plot/sol.dat')
    x(1,:) = initialize('input/parameters.in')
    call Heun(T, N, tk, x)
    call PrintAndSave(tk, x, 'plot/sol.dat')
    x(1,:) = initialize('input/parameters.in')
    call EulerBackward(T, N, tk, x)
    call PrintAndSave(tk, x, 'plot/sol.dat')

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
end program