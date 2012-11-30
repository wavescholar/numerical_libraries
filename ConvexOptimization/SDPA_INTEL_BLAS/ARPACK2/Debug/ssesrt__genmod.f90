        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:22 2012
        MODULE SSESRT__genmod
          INTERFACE 
            SUBROUTINE SSESRT(WHICH,APPLY,N,X,NA,A,LDA)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: N
              CHARACTER(LEN=2) :: WHICH
              LOGICAL(KIND=4) :: APPLY
              REAL(KIND=4) :: X(0:N-1)
              INTEGER(KIND=4) :: NA
              REAL(KIND=4) :: A(LDA,0:N-1)
            END SUBROUTINE SSESRT
          END INTERFACE 
        END MODULE SSESRT__genmod
