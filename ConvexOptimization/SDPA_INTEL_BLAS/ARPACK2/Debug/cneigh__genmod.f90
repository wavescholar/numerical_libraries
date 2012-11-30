        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:23 2012
        MODULE CNEIGH__genmod
          INTERFACE 
            SUBROUTINE CNEIGH(RNORM,N,H,LDH,RITZ,BOUNDS,Q,LDQ,WORKL,    &
     &RWORK,IERR)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: RNORM
              COMPLEX(KIND=4) :: H(LDH,N)
              COMPLEX(KIND=4) :: RITZ(N)
              COMPLEX(KIND=4) :: BOUNDS(N)
              COMPLEX(KIND=4) :: Q(LDQ,N)
              COMPLEX(KIND=4) :: WORKL(N*(N+3))
              REAL(KIND=4) :: RWORK(N)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE CNEIGH
          END INTERFACE 
        END MODULE CNEIGH__genmod
