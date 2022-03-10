MODULE Newton
    
Contains    
    
subroutine NewtonSolver(const1,const2,const3,const4,sigma_prod,init,solution)
    
implicit none

integer, parameter :: DOUBLE     = SELECTED_REAL_KIND(p=10)

real(KIND=DOUBLE), intent(in)  :: const1, const2, const3, const4, sigma_prod, init
real(KIND=DOUBLE), intent(out) :: solution

real(KIND=DOUBLE) eps, x, check, f, fprime
logical           converged
integer           it


eps = 0.000001
x = init
it = 0

converged = .false.

do while (.not. converged .and. it <= 200)
    
    f = (const1*x**(sigma_prod-1)) / (const2 + const3*x**sigma_prod) - const4
    fprime = ( (const1*(sigma_prod-1)*x**(sigma_prod-2))*(const2+const3*x**sigma_prod) - &
               (const1*x**(sigma_prod-1))*const3*sigma_prod*x**(sigma_prod-1) ) / &
               (const2 + const3*x**sigma_prod)**2
    
    x = x - 0.2*(f / fprime)
    
    check = abs((const1*x**(sigma_prod-1)) / (const2 + const3*x**sigma_prod) - const4)
    
    converged = (check <= eps)
    
    it = it + 1
    
end do

solution = x

end subroutine

end MODULE Newton
