Module LinReg_MOD

! *****************************************************
! Computes least square regression coefficients and R2.
! LinReg2 gives the option of feeding the inv(X'X) 
! matrix directly.
! *****************************************************

Contains

subroutine LinReg(Y,X,dim1,dim2,COEF,SST,SSE,info)

! *******************************************************************
! Estimates the linear regression Y(dim1,1) = COEF*X(dim1,dim2) + eps
! Output variables: COEF, SST (total sum of squares)
! SSE (total sum of residuals) and info which is an indicator for 
! whether the cross product matrix is of full rank
! *******************************************************************

implicit none

integer, parameter :: DOUBLE     = SELECTED_REAL_KIND(p=10)

integer          , intent(in)  :: dim1, dim2
real(KIND=DOUBLE), intent(in)  :: Y(:,:), X(:,:)
real(KIND=DOUBLE), intent(out) :: COEF(:,:), SST, SSE
integer          , intent(out) :: info

real(KIND=DOUBLE) XX(dim2,dim2), invXX(dim2,dim2)

COEF  = 0.0

call get_invXX(X,dim1,dim2,XX,invXX,info)
call LinReg2(invXX,X,Y,dim1,dim2,COEF,SSE,SST,1)

end subroutine LinReg




subroutine LinReg2(invXX,X,Y,dim1,dim2,COEF,SSE,SST,flag,eps)

! *******************************************************************
! Estimates the linear regression Y(dim1,1) = COEF*X(dim1,dim2) + eps
! However, we feed this subroutine with inv(X'X)
! Output variables: COEF, SST (total sum of squares)
! SSE (total sum of residuals) and info which is an indicator for 
! whether the cross product matrix is of full rank, eps (residuals)
! are optional
! *******************************************************************

! **************************************************
! Estimates the linear regression Y = COEF*X + eps
! Output variables: COEF, SST (total sum of squares)
! and SSE (total sum of residuals)
! It uses MKL's subroutine dgels
! **************************************************

implicit none

integer, parameter :: DOUBLE     = SELECTED_REAL_KIND(p=10)

integer          , intent(in)  :: dim1, dim2, flag
real(KIND=DOUBLE), intent(in)  :: invXX(:,:), Y(:,:), X(:,:)
real(KIND=DOUBLE), intent(out) :: COEF(:,:), SSE, SST
real(KIND=DOUBLE), intent(out), optional :: eps(:,:)

real(KIND=DOUBLE) XY(dim2,1), resid(dim1,1), resid_2(dim1,1), Y_hat(dim1,1), & 
                  disp(dim1,1), Y_bar

integer one

resid = 0.0

COEF = 0.0

one = 1

call mult_XY(X,Y,dim1,dim2,XY)

call matrix_mult(invXX,XY,dim2,dim2,one,COEF)

if (flag == 1) then

    call matrix_mult(X,COEF,dim1,dim2,one,Y_hat)    
    
    resid = Y - Y_hat
    
    resid_2 = resid **2
    
    SSE = sum(resid_2)
    
    Y_bar = sum(Y) / dim1
    
    disp  = Y - Y_bar
    
    SST = sum(disp**2)
    
    if(present(eps)) then
        eps = resid
    end if

end if

end subroutine LinReg2



subroutine get_invXX(X,m,n,XX,invXX,info)

! *****************
! Computes inv(X'X)
! *****************

implicit none

integer, parameter :: DOUBLE     = SELECTED_REAL_KIND(p=10)

integer            , intent(in) :: m, n
real(KIND=DOUBLE), intent(in) :: X(m,n)
real(KIND=DOUBLE), intent(out) :: XX(n,n), invXX(n,n)


real(KIND=DOUBLE) alpha, beta, U(n,n)
integer i, j, k, info

XX = 0.0

alpha = 1.0
beta  = 0.0

call dgemm('T', 'N', n, n, m, alpha, X, m, X, m, beta, XX, n)



U = XX

call dpotrf('U',n,U,n,info)

invXX = U

call dpotri('U',n,invXX,n,info)

do i = 1, n
    do j = 1, n
        if (i > j) then
            invXX(i,j) = invXX(j,i)
        end if
    end do
end do

end subroutine get_invXX


subroutine mult_XY(X,Y,m,n,XY)

! ************
! Computes X*Y
! ************

implicit none

integer, parameter :: DOUBLE     = SELECTED_REAL_KIND(p=10)

integer            , intent(in) :: m, n
real(KIND = DOUBLE), intent(in) :: X(m,n), Y(m,1)
real(KIND = DOUBLE), intent(out) :: XY(n,1)

real(KIND=DOUBLE) alpha, beta
integer colY
            
alpha = 1.0
beta  = 0.0
colY = 1

call dgemm('T', 'N', n, colY, m, alpha, X, m, Y, m, beta, XY, n)

end subroutine mult_XY


subroutine matrix_mult(A,B,m,n,k,AB)

implicit none

integer, parameter :: DOUBLE     = SELECTED_REAL_KIND(p=10)

integer, intent(in) :: m,n,k
real(KIND=DOUBLE), intent(in) :: A(m,n), B(n,k), AB(m,k)

real(KIND=DOUBLE) alpha, beta
integer colY
            
alpha = 1.0
beta  = 0.0

call dgemm('N', 'N', m, k, n, alpha, A, m, B, n, beta, AB, m)

end subroutine matrix_mult

end module LinReg_MOD

