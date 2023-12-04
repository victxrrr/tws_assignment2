module MySolver 

    implicit none
    private
    public solve
    integer, parameter ::             &
    sp = selected_real_kind(6, 37),   &
    dp = selected_real_kind(15, 307), &
    ep = selected_real_kind(18,4931), &
    qp = selected_real_kind(33,4931)

    interface solve
        module procedure solve_sp
        module procedure solve_dp
    end interface

contains

    !> Solve Ax = b in single precision
    !! @param A matrix of coefficients
    !! @param bx RHS of the cequations, contains x after computation
    subroutine solve_sp(A,bx)
        real(sp), dimension(:,:), intent(inout) :: A
        real(sp), dimension(:), intent(inout)   :: bx
        integer :: info, ipiv(size(A,1))

        call sgesv(size(A,1),1,A,size(A,1),ipiv,bx,size(bx),info)
        call error_msg(info)
    end subroutine

    !> Solve Ax = b in double precision
    !! @param A matrix of coefficients
    !! @param bx RHS of the cequations, contains x after computation
    subroutine solve_dp(A,bx)
        real(dp), dimension(:,:), intent(inout) :: A
        real(dp), dimension(:), intent(inout)   :: bx
        integer :: info, ipiv(size(A,1))

        call dgesv(size(A,1),1,A,size(A,1),ipiv,bx,size(bx),info)
        call error_msg(info)
    end subroutine

    !> Simple return value checker
    subroutine error_msg(INFO)
        integer, intent(in) :: INFO
    
        if (INFO < 0) then
            print *, 'ERROR: The ', ABS(INFO), '-th argument had an illegal value. (INFO = ', INFO, ')'
        else if (INFO > 0) then
            print *, 'ERROR: U(', ABS(INFO), ',', ABS(INFO), ') is exactly zero. The factorization has been completed,&
             & but the factor U is exactly singular, so the solution could not be computed. (INFO = ', INFO, ')'
        end if
    end subroutine    
end module

module Utils
    use MySolver, only:solve

    implicit none   
    save
    private
    public EulerForward, Heun, EulerBackward, initialize, utils_set_verbose, utils_set_save, wp, PrintAndSave

    integer, parameter ::             &
    sp = selected_real_kind(6, 37),   &
    dp = selected_real_kind(15, 307), &
    ep = selected_real_kind(18,4931), &
    qp = selected_real_kind(33,4931), &
    wp = dp

    real(wp) :: beta, mu, gamma, alpha, delta
    logical  :: verbose = .true.
    logical  :: save = .true.

contains

    !> Auxiliary function for initializing the parameters and the state vector of the SIQRD model
    !! @param filename file that contains the initial parameters
    !! @param [optional] beta_opt to store the initial value of beta
    !! @return x0 = [S_0, I_0, Q_0, R_0, D_0]
    function initialize(filename, beta_opt) result(x0)
        character(*), intent(in)        :: filename
        real(wp), intent(out), optional :: beta_opt
        real(wp)                        :: x0(5)
        x0 = 0._wp

        open(7,file=trim(filename))
        read(7,*) beta,mu,gamma,alpha,delta,x0(1),x0(2)
        close(7)

        if (present(beta_opt)) beta_opt = beta
    end function

    !> Auxiliary function for extracting the system parameters from parameters.in  describing the SIQRD model 
    !! and for returning the corresponding evalutation of xk as a function of these parameters
    !! @param xk state vector at time tk
    !! @return eval = f([params], xk)
    function f(xk) result(eval)
        real(wp), intent(in)           :: xk(:)
        real(wp), dimension(size(xk))  :: eval

        eval(1) = -beta*xk(1)*xk(2)/(xk(1) + xk(2) + xk(4)) + mu*xk(4)
        eval(2) = (beta*xk(1)/(xk(1) + xk(2) + xk(4)) - gamma - delta - alpha)*xk(2)
        eval(3) = delta*xk(2) - (gamma + alpha)*xk(3)
        eval(4) = gamma*(xk(2) + xk(3)) - mu*xk(4)
        eval(5) = alpha*(xk(2) + xk(3))

    end function

    !> Auxiliary function for extracting the system paramteres from parameters.in describing the SIQRD model
    !! and for returning the jacobian matrix of f evaluated at xk
    !! @param xk
    !! @return eval = df/dx([params], xk)
    function jac_f(xk) result(eval)
        real(wp), intent(in)                    :: xk(:)
        real(wp), dimension(size(xk), size(xk)) :: eval
        real(wp)                                :: tmp

        tmp = (xk(1) + xk(2) + xk(4))**2

        eval(1,1)   = -beta * (xk(2)*xk(2) + xk(2)*xk(4))/tmp
        eval(2,1)   = -eval(1,1)
        eval(3:,1)  = 0._wp
        eval(1,2)   = -beta * (xk(1)*xk(1) + xk(1)*xk(4))/tmp
        eval(2,2)   = -eval(1,2)- gamma - delta - alpha
        eval(3,2)   = delta
        eval(4,2)   = gamma
        eval(5,2)   = alpha
        eval(1:2,3) = 0._wp
        eval(3,3)   = -(gamma + alpha)
        eval(4,3)   = gamma
        eval(5,3)   = alpha
        eval(1,4)   = beta*xk(1)*xk(2)/tmp + mu
        eval(2,4)   = -(eval(1,4) - mu)
        eval(3,4)   = 0._wp
        eval(4,4)   = -mu
        eval(5,4)   = 0._wp
        eval(:,5)   = 0._wp
   
    end function

    !> Euler's forward method tailored to the SIQRD model
    !! @param T simulation horizon
    !! @param N number of discretization steps
    !! @param tk array containing the N + 1 equidistant grid points
    !! @param x 2D-array containing the solution approximation at time 0, T/N, 2T/N, ... , T
    !! @param beta_opt to specify the value of the beta parameter to be used [optional]
    subroutine EulerForward(T, N, tk, x, beta_opt)
        real(wp), intent(in)           :: T
        integer, intent(in)            :: N
        real(wp), intent(out)          :: tk(:), x(:,:)
        real(wp), intent(in), optional :: beta_opt
        integer                        :: k

        tk(1) = 0._wp
        if (present(beta_opt)) beta = beta_opt
        do k = 1,N
            tk(k+1) = (T/N)*k
            x(k+1,:) = x(k,:) + (T/N)*f(x(k,:))
        enddo
    end subroutine

    !> Heun's method tailored to the SIQRD model
    !! @param T simulation horizon
    !! @param N number of discretization steps
    !! @param tk array containing the N + 1 equidistant grid points
    !! @param x 2D-array containing the solution approximation at time 0, T/N, 2T/N, ... , T
    !! @param beta_opt to specify the value of the beta parameter to be used [optional]
    subroutine Heun(T, N, tk, x, beta_opt)
        real(wp), intent(in)              :: T
        integer, intent(in)               :: N
        real(wp), intent(out)             :: tk(:), x(:,:)
        real(wp), intent(in), optional    :: beta_opt
        integer                           :: k

        tk(1) = 0._wp
        if (present(beta_opt)) beta = beta_opt
        do k = 1,N
            tk(k+1) = (T/N)*k
            x(k+1,:) = x(k,:) + (T/N)*(0.5*f(x(k,:)) + 0.5*f(x(k,:) + (T/N)*f(x(k,:))))
        enddo
    end subroutine

    !> Euler's backward method tailored to the SIQRD model
    !! @param T simulation horizon
    !! @param N number of discretization steps
    !! @param tk array containing the N + 1 equidistant grid points
    !! @param x 2D-array containing the solution approximation at time 0, T/N, 2T/N, ... , T
    !! @param tol_opt [optional] tolerance for the stopping criterion of Newton's method
    !! @param nmax_opt [optional] maximum number of iterations of Newton's method
    !! @param beta_opt [optional] to specify the value of the beta parameter to be used
    subroutine EulerBackward(T, N, tk, x, tol_opt, nmax_opt, beta_opt)
        real(wp), intent(in)                     :: T
        integer, intent(in)                      :: N
        real(wp), intent(out)                    :: tk(:), x(:,:)
        real(wp), intent(in), optional           :: tol_opt
        integer, intent(in), optional            :: nmax_opt
        real(wp), intent(in), optional           :: beta_opt
        real(wp), dimension(size(x,2),size(x,2)) :: A
        real(wp)                                 :: bx(size(x,2)), tol
        integer                                  :: k, i, s, nmax

        if (present(tol_opt)) then 
            tol = tol_opt
        else 
            tol = 1e-10
        endif

        if (present(nmax_opt)) then
            nmax = nmax_opt
        else
            nmax = 100
        endif

        tk(1) = 0._wp
        if (present(beta_opt)) beta = beta_opt
        do k = 1, N
            tk(k+1) = (T/N)*k
            
            ! Euler's method with Backward Error stopping criterion
            x(k+1,:) = x(k,:)
            do s = 1, nmax
                bx = x(k, :) + (T/N)*f(x(k+1,:)) - x(k+1,:)
                if (norm2(bx) < tol*norm2(x(k+1,:))) exit
                A = (T/N)*jac_f(x(k+1,:))
                do i = 1, size(A,1)
                    A(i,i) = A(i,i) - 1._wp
                enddo
    
                call solve(A,bx)
                x(k+1,:) = x(k+1,:) - bx
            enddo
            
            if (s == nmax + 1) then
                print *, "ERROR: Newton's method failed to converge during the Euler Backward method, &
                        &maximum number of iterations reached."
            endif
        enddo
    end subroutine

    !> Auxiliary function for displaying the calculated solution on the flow output and 
    !! saving it in the file filename
    !! @param tk array containing the N + 1 equidistant grid points
    !! @param x 2D-array containing the solution approximation at time 0, T/N, 2T/N, ... , T
    !! @param filename the existing file in which the values are stored
    subroutine PrintAndSave(tk, x, filename)
        real(wp), intent(in)     :: tk(:), x(:,:)
        character(*), intent(in) :: filename
        integer                  :: k

        open(7, file=trim(filename))
        do k = 1, size(x,1)
            if (verbose) write(*, '(f12.2, 5f16.5)') tk(k), x(k, 1), x(k, 2), x(k, 3), x(k, 4), x(k, 5)
            if (save) write(7, '(f12.2, 5f16.5)') tk(k), x(k, 1), x(k, 2), x(k, 3), x(k, 4), x(k, 5)
        enddo
        close(7)
    end subroutine

    subroutine utils_set_verbose(bool)
        logical, intent(in) :: bool
        verbose = bool
    end subroutine

    subroutine utils_set_save(bool)
        logical, intent(in) :: bool
        save = bool
    end subroutine
end module