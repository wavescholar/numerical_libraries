        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:31 2012
        MODULE ZSORTC__genmod
          INTERFACE 
            SUBROUTINE ZSORTC(WHICH,APPLY,N,X,Y)
              INTEGER(KIND=4) :: N
              CHARACTER(LEN=2) :: WHICH
              LOGICAL(KIND=4) :: APPLY
              COMPLEX(KIND=8) :: X(0:N-1)
              COMPLEX(KIND=8) :: Y(0:N-1)
            END SUBROUTINE ZSORTC
          END INTERFACE 
        END MODULE ZSORTC__genmod
