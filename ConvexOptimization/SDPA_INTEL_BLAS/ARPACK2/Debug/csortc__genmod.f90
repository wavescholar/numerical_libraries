        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:30 2012
        MODULE CSORTC__genmod
          INTERFACE 
            SUBROUTINE CSORTC(WHICH,APPLY,N,X,Y)
              INTEGER(KIND=4) :: N
              CHARACTER(LEN=2) :: WHICH
              LOGICAL(KIND=4) :: APPLY
              COMPLEX(KIND=4) :: X(0:N-1)
              COMPLEX(KIND=4) :: Y(0:N-1)
            END SUBROUTINE CSORTC
          END INTERFACE 
        END MODULE CSORTC__genmod
