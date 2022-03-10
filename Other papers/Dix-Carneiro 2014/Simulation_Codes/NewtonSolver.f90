MODULE Newton
    
Contains    
    
subroutine NewtonSolver(Const1,Const2,Const3,sigma,init,solution)
    
implicit none

integer, parameter :: DOUBLE     = SELECTED_REAL_KIND(p=10)

real(KIND=DOUBLE), intent(in)  :: Const1, Const2, Const3, sigma, init
real(KIND=DOUBLE), intent(out) :: solution

real(KIND=DOUBLE) eps, x, check, f, fcheck, fprime
logical           converged
integer           it


eps = 0.00000001
x = init
it = 0

converged = .false.

do while (.not. converged .and. it <= 500)
    
    f = ((Const1 + Const3*x**sigma)**(1/sigma - 1))*Const3*x**(sigma - 1) - Const2
    fprime = (Const3**2)*(x**(2*sigma - 2))*(1-sigma)*(Const1 + Const3*x**sigma)**(1/sigma - 2) + &
             Const3*(sigma - 1)*(x**(sigma - 2))*(Const1 + Const3*x**sigma)**(1/sigma - 1)
    
    x = x - 0.1*(f / fprime)
    
    fcheck = log( abs(((Const1 + Const3*x**sigma)**(1/sigma - 1))*Const3*x**(sigma - 1)) ) - log( abs(Const2) )
    
    check = abs(fcheck)
    
    converged = (check <= eps)
    
    it = it + 1
    
end do

solution = x

end subroutine

end MODULE Newton
