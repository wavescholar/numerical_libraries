        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:31 2012
        MODULE DNEIGH__genmod
          INTERFACE 
            SUBROUTINE DNEIGH(RNORM,N,H,LDH,RITZR,RITZI,BOUNDS,Q,LDQ,   &
     &WORKL,IERR)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: RNORM
              REAL(KIND=8) :: H(LDH,N)
              REAL(KIND=8) :: RITZR(N)
              REAL(KIND=8) :: RITZI(N)
              REAL(KIND=8) :: BOUNDS(N)
              REAL(KIND=8) :: Q(LDQ,N)
              REAL(KIND=8) :: WORKL(N*(N+3))
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE DNEIGH
          END INTERFACE 
        END MODULE DNEIGH__genmod
