        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:31 2012
        MODULE SNEIGH__genmod
          INTERFACE 
            SUBROUTINE SNEIGH(RNORM,N,H,LDH,RITZR,RITZI,BOUNDS,Q,LDQ,   &
     &WORKL,IERR)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: RNORM
              REAL(KIND=4) :: H(LDH,N)
              REAL(KIND=4) :: RITZR(N)
              REAL(KIND=4) :: RITZI(N)
              REAL(KIND=4) :: BOUNDS(N)
              REAL(KIND=4) :: Q(LDQ,N)
              REAL(KIND=4) :: WORKL(N*(N+3))
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE SNEIGH
          END INTERFACE 
        END MODULE SNEIGH__genmod
