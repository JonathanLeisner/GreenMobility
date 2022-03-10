
      SUBROUTINE CALFUN (N,X,F)
      
      USE Global_Data
      USE Loss_Function_MOD
      
      IMPLICIT NONE
      
      INTEGER N
      REAL(KIND=DOUBLE) X(N), F
      
      CALL SMM_OBJ_FCN(X,F)      
      
      RETURN
      END 
