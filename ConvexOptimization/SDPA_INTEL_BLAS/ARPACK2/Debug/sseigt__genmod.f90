        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:25 2012
        MODULE SSEIGT__genmod
          INTERFACE 
            SUBROUTINE SSEIGT(RNORM,N,H,LDH,EIG,BOUNDS,WORKL,IERR)
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: RNORM
              REAL(KIND=4) :: H(LDH,2)
              REAL(KIND=4) :: EIG(N)
              REAL(KIND=4) :: BOUNDS(N)
              REAL(KIND=4) :: WORKL(3*N)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE SSEIGT
          END INTERFACE 
        END MODULE SSEIGT__genmod
