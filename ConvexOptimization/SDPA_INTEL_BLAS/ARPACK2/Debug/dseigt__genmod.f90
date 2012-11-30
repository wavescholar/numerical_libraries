        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:31 2012
        MODULE DSEIGT__genmod
          INTERFACE 
            SUBROUTINE DSEIGT(RNORM,N,H,LDH,EIG,BOUNDS,WORKL,IERR)
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: RNORM
              REAL(KIND=8) :: H(LDH,2)
              REAL(KIND=8) :: EIG(N)
              REAL(KIND=8) :: BOUNDS(N)
              REAL(KIND=8) :: WORKL(3*N)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE DSEIGT
          END INTERFACE 
        END MODULE DSEIGT__genmod
