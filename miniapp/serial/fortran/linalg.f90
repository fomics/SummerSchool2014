! linear algebra subroutines
! Ben Cumming @ CSCS

module linalg

use stats,     only: flops_blas1, iters_cg
use data,      only: discretizationT, options, bndN, bndE, bndS, bndS
use operators, only: diffusion

implicit none

logical    :: cg_initialized=.false.
real (kind=8), allocatable :: r(:), Ap(:), p(:), Fx(:), Fxold(:), v(:), xold(:)

contains

! initialize temporary storage fields used by the cg solver
! I do this here so that the fields are persistent between calls
! to the CG solver. This is useful if we want to avoid malloc/free calls
! on the device for the OpenACC implementation (feel free to suggest a better
! method for doing this)
subroutine cg_init(N)
    ! arguments
    integer,    intent(in)  :: N

    allocate(Ap(N), r(N), p(N), Fx(N), Fxold(N), v(N), xold(N))
    cg_initialized = .true.
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   blas level 1 reductions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computes the inner product of x and y
!       x and y are vectors on length N
real (kind=8) function ss_dot(x, y, N)
    ! arguments
    real (kind=8), intent(in) :: x(N)
    real (kind=8), intent(in) :: y(N)
    integer, intent(in) :: N

    ! local variables
    !real (kind=8) :: ss_dot
    integer :: i

    ! the logic
    ss_dot = 0
    do i = 1, N
        ss_dot = ss_dot + x(i) * y(i)
    enddo

    ! record the number of floating point oporations
    flops_blas1 = flops_blas1 + 2*N

    return
end function

! computes the 2-norm of x
!       x is a vector on length N
real (kind=8) function ss_norm2(x, N)
    ! arguments
    real (kind=8), intent(in) :: x(N)
    integer, intent(in) :: N

    ! local variables
    !real (kind=8) :: ss_dot
    integer :: i

    ! the logic
    ss_norm2 = 0
    do i = 1, N
        ss_norm2 = ss_norm2 + x(i) * x(i)
    enddo
    ss_norm2 = sqrt(ss_norm2)

    flops_blas1 = flops_blas1 + 2*N

    return
end function

! sets entries in a vector to value
!       x is a vector on length N
!       value is th
subroutine ss_fill(x, value, N)
    ! arguments
    real (kind=8), intent(inout) :: x(N)
    real (kind=8), intent(in)    :: value
    integer, intent(in) :: N

    ! local variables
    integer :: i

    ! the logic
    do i = 1, N
        x(i) = value
    enddo
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   blas level 1 vector-vector operations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computes y := alpha*x + y
!       x and y are vectors on length N
!       alpha is a scalar
subroutine ss_axpy(y, alpha, x, N)
    ! arguments
    real (kind=8), intent(in)    :: alpha
    real (kind=8), intent(in)    :: x(N)
    real (kind=8), intent(inout) :: y(N)
    integer, intent(in)          :: N

    ! local variables
    integer :: i

    ! the logic
    do i = 1, N
        y(i) = alpha*x(i) + y(i)
    enddo

    ! update the flops counter
    flops_blas1 = flops_blas1 + 2*N
end subroutine

! computes y = x + alpha*(l-r)
!       y, x, l and r are vectors of length N
!       alpha is a scalar
subroutine ss_add_scaled_diff(y, x, alpha, l, r, N)
    ! arguments
    real (kind=8), intent(in)    :: alpha
    real (kind=8), intent(in)    :: x(N)
    real (kind=8), intent(out)   :: y(N)
    real (kind=8), intent(in)    :: l(N)
    real (kind=8), intent(in)    :: r(N)
    integer, intent(in)          :: N

    ! local variables
    integer :: i

    do i = 1, N
        y(i) = x(i) + alpha * (l(i) - r(i))
    enddo

    ! update the flops counter
    flops_blas1 = flops_blas1 + 3*N
end subroutine

! computes y = alpha*(l-r)
!       y, l and r are vectors of length N
!       alpha is a scalar
subroutine ss_scaled_diff(y, alpha, l, r, N)
    ! arguments
    real (kind=8), intent(in)    :: alpha
    real (kind=8), intent(out)   :: y(N)
    real (kind=8), intent(in)    :: l(N)
    real (kind=8), intent(in)    :: r(N)
    integer, intent(in)          :: N

    ! local variables
    integer :: i

    do i = 1, N
        y(i) = alpha * (l(i) - r(i))
    enddo

    ! update the flops counter
    flops_blas1 = flops_blas1 + 2*N
end subroutine

! computes y := alpha*x
!   alpha is scalar
!   y and x are vectors on length n
subroutine ss_scale(y, alpha, x, N)
    ! arguments
    real (kind=8), intent(in)    :: alpha
    real (kind=8), intent(in)    :: x(N)
    real (kind=8), intent(inout) :: y(N)
    integer, intent(in)          :: N

    ! local variables
    integer :: i

    ! the logic
    do i = 1, N
        y(i) = alpha*x(i)
    enddo

    ! update the flops counter
    flops_blas1 = flops_blas1 + N
end subroutine

! computes linear combination of two vectors y := alpha*x + beta*z
!   alpha and beta are scalar
!   y, x and z are vectors on length n
subroutine ss_lcomb(y, alpha, x, beta, z, N)
    ! arguments
    real (kind=8), intent(inout) :: y(N)
    real (kind=8), intent(in)    :: alpha
    real (kind=8), intent(in)    :: x(N)
    real (kind=8), intent(in)    :: beta
    real (kind=8), intent(in)    :: z(N)
    integer, intent(in)          :: N

    ! local variables
    integer :: i

    ! the logic
    do i = 1, N
        y(i) = alpha*x(i) + beta*z(i)
    enddo

    ! update the flops counter
    flops_blas1 = flops_blas1 + 3*N
end subroutine

! copy one vector into another y := x
!       x and y are vectors of length N
subroutine ss_copy(y, x, N)
    ! arguments
    real (kind=8), intent(in)    :: x(N)
    real (kind=8), intent(inout) :: y(N)
    integer, intent(in)          :: N

    ! local variables
    integer :: i

    ! the logic
    do i = 1, N
        y(i) = x(i)
    enddo
end subroutine

! conjugate gradient solver
! solve the linear system A*x = b for x
! the matrix A is implicit in the objective function for the diffusion equation
!   the value in x constitute the "first guess" at the solution
!   x(N)
!       ON ENTRY contains the initial guess for the solution
!       ON EXIT  contains the solution
subroutine ss_cg(x, b, maxiters, tol, success)
    ! arguments
    real (kind=8), intent(inout)     :: x(options%N)
    real (kind=8), intent(in)        :: b(options%N)
    integer,       intent(in)        :: maxiters
    real (kind=8), intent(in)        :: tol
    logical,       intent(out)       :: success

    ! local parameters
    integer                    :: iter, N, i
    real (kind=8)              :: alpha, rold, rnew, eps, eps_inv
    real (kind=8)              :: one, zero

    ! this is the dimension of the linear system that we are to solve
    N = options%N

    if (.NOT. cg_initialized) then
        call cg_init(N)
    endif

    ! useful constants
    zero = 0.
    one  = 1.

    ! epslion value use for matrix-vector approximation
    eps     = 1.e-8
    eps_inv = 1./eps

    ! allocate memory for temporary storage
    call ss_fill(Fx,    zero, N)
    call ss_fill(Fxold, zero, N)
    call ss_copy(xold, x, N)

    ! matrix vector multiplication is approximated with
    !   A*v = 1/epsilon * ( F(x+epsilon*v) - F(x) )
    !       = 1/epsilon * ( F(x+epsilon*v) - Fxold )
    ! we compute Fxold at startup
    ! we have to keep x so that we can compute the F(x+exps*v)

    call diffusion(x, Fxold)

    ! v = x + epsilon*x
    call ss_scale(v, one+eps, x, N)

    ! Fx = F(v)
    call diffusion(v, Fx)

    ! r = b - A*x
    ! where A*x = (Fx-Fxold)/eps
    call ss_add_scaled_diff(r, b, -eps_inv, Fx, Fxold, N)

    ! p = r
    call ss_copy(p, r, N)

    ! rold = <r,r>
    rold = ss_dot(r, r, N)

    ! check for convergence
    success = .false.
    if( dsqrt(rold)<tol ) then
        success = .true.
        return
    endif

    do iter=1, maxiters
        !Ap = A*p
        call ss_lcomb(v, one, xold, eps, p, N)
        call diffusion(v, Fx)
        call ss_scaled_diff(Ap, eps_inv, Fx, Fxold, N)

        !alpha = rold / p'*Ap
        alpha = rold / ss_dot(p, Ap, N)

        ! x += alpha*p
        call ss_axpy(x, alpha, p, N)

        ! r -= alpha*Ap
        call ss_axpy(r, -alpha, Ap, N)

        ! find new norm
        rnew = ss_dot(r, r, N)

        ! test for convergence
        if( dsqrt(rnew)<tol ) then
            success = .true.
            exit
        endif

        ! p = r + rnew.rold * p
        call ss_lcomb(p, one, r, rnew/rold, p, N)

        rold = rnew
    enddo

    iters_cg = iters_cg + iter

    if (.NOT. success) then
        write(*,*) 'ERROR: CG failed to converge after ', maxiters, ' iterations'
    endif
end subroutine

end module linalg

