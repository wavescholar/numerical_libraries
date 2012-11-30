        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:28 2012
        MODULE DSESRT__genmod
          INTERFACE 
            SUBROUTINE DSESRT(WHICH,APPLY,N,X,NA,A,LDA)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: N
              CHARACTER(LEN=2) :: WHICH
              LOGICAL(KIND=4) :: APPLY
              REAL(KIND=8) :: X(0:N-1)
              INTEGER(KIND=4) :: NA
              REAL(KIND=8) :: A(LDA,0:N-1)
            END SUBROUTINE DSESRT
          END INTERFACE 
        END MODULE DSESRT__genmod
