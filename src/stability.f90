program Q4

    implicit none

    integer, parameter ::             &
    sp = selected_real_kind(6, 37),   &
    dp = selected_real_kind(15, 307), &
    wp = dp

    real(wp) :: beta, mu, gamma, alpha, delta, J(5,5), wr(5), wi(5)

    character(len=256)    :: MSG, file
    integer               :: INFO, STATUS = 0
    real(wp), allocatable :: work(:)
    integer               :: lwork, k
    logical               :: unstable = .false.

    if (command_argument_count() /= 1) then 
        print *, 'Wrong number of arguments'
        stop STATUS
    endif

    call get_command_argument(1,file)
    open(7,file=trim(file))
    read(7,*) beta,mu,gamma,alpha,delta
    close(7)

    J = jac_f_eq(beta, mu, gamma, alpha, delta)

    if (wp == sp) then
        call sgeev('N', 'N', 5, J, 5, wr, wi, J, 5, J, 5, wi, -1, INFO)
        call workcheck(INFO)

        lwork = int(wi(1))
        allocate(work(lwork), STAT=STATUS, ERRMSG=MSG)
        call memcheck(STATUS,MSG)

        call sgeev('N', 'N', 5, J, 5, wr, wi, J, 5, J, 5, work, lwork, INFO)
        call errorcheck(INFO)
    else if (wp == dp) then 
        call dgeev('N', 'N', 5, J, 5, wr, wi, J, 5, J, 5, wi, -1, INFO)
        call workcheck(INFO)

        lwork = int(wi(1))
        allocate(work(lwork), STAT=STATUS, ERRMSG=MSG)
        call memcheck(STATUS,MSG)

        call dgeev('N', 'N', 5, J, 5, wr, wi, J, 5, J, 5, work, lwork, INFO)
        call errorcheck(INFO)
    endif

    do k = 1, 5 
        if (wr(k) > 0._wp) unstable = .true.
        print '(a,i2,a,f5.2,a,f5.2,a)', ' Eigenvalue', k, ' :', wr(k), ' +', wi(k), 'j'
    enddo
    if (unstable) then 
        print *, 'The given parameters have led to an outbreak.'
    else
        print *, 'The given parameters have not led to an outbreak.'
    endif

    deallocate(work, STAT=STATUS, ERRMSG=MSG)
    call memcheck(STATUS,MSG)

contains

    function jac_f_eq(beta, mu, gamma, alpha, delta) result(J)
        real(wp), intent(in)     :: beta, mu, gamma, alpha, delta
        real(wp), dimension(5,5) :: J

        J(:,1)   = 0._wp
        J(:,5)   = 0._wp
        J(1,2)   = - beta 
        J(1,3:4) = 0._wp
        J(2,2)   = beta - gamma - delta - alpha
        J(2,3:4) = 0._wp
        J(3,2)   = delta
        J(3,3)   = - (alpha + gamma)
        J(3,4)   = 0._wp
        J(4,2:3) = gamma
        J(4,4)   = - mu
        J(5,2:3) = alpha
        J(5,4)   = 0._wp
    end function

    subroutine memcheck(STATUS,MSG)
        integer, intent(in)            :: STATUS
        character(len=256), intent(in) :: MSG
        if (STATUS /= 0) then
            print *, MSG
            stop STATUS
        endif
    end subroutine

    subroutine workcheck(INFO)
        integer, intent(in) :: INFO
        if (INFO /= 0) then
            print *, 'ERROR: something went wrong when computing the optimal size of the WORK array.'
            stop INFO
        endif
    end subroutine

    subroutine errorcheck(INFO)
        integer, intent(in) :: INFO
    
        if (INFO < 0) then
            print *, 'ERROR: The ', ABS(INFO), '-th argument had an illegal value. (INFO = ', INFO, ')'
        else if (INFO > 0) then
            print *, 'ERROR: the QR algorithm failed to compute all the eigenvalues, and no eigenvectors have been computed;&
            &elements ', ABS(INFO), '+1:N of WR and WI contain eigenvalues which have converged. (INFO = ', INFO, ')'
        end if
    end subroutine
end program

