        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:21 2012
        MODULE DSORTC__genmod
          INTERFACE 
            SUBROUTINE DSORTC(WHICH,APPLY,N,XREAL,XIMAG,Y)
              INTEGER(KIND=4) :: N
              CHARACTER(LEN=2) :: WHICH
              LOGICAL(KIND=4) :: APPLY
              REAL(KIND=8) :: XREAL(0:N-1)
              REAL(KIND=8) :: XIMAG(0:N-1)
              REAL(KIND=8) :: Y(0:N-1)
            END SUBROUTINE DSORTC
          END INTERFACE 
        END MODULE DSORTC__genmod
