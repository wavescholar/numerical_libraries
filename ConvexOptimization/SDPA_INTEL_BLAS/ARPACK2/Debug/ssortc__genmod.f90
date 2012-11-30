        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:21 2012
        MODULE SSORTC__genmod
          INTERFACE 
            SUBROUTINE SSORTC(WHICH,APPLY,N,XREAL,XIMAG,Y)
              INTEGER(KIND=4) :: N
              CHARACTER(LEN=2) :: WHICH
              LOGICAL(KIND=4) :: APPLY
              REAL(KIND=4) :: XREAL(0:N-1)
              REAL(KIND=4) :: XIMAG(0:N-1)
              REAL(KIND=4) :: Y(0:N-1)
            END SUBROUTINE SSORTC
          END INTERFACE 
        END MODULE SSORTC__genmod
