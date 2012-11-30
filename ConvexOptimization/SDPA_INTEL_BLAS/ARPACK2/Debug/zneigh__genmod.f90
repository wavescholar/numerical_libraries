        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:24 2012
        MODULE ZNEIGH__genmod
          INTERFACE 
            SUBROUTINE ZNEIGH(RNORM,N,H,LDH,RITZ,BOUNDS,Q,LDQ,WORKL,    &
     &RWORK,IERR)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: RNORM
              COMPLEX(KIND=8) :: H(LDH,N)
              COMPLEX(KIND=8) :: RITZ(N)
              COMPLEX(KIND=8) :: BOUNDS(N)
              COMPLEX(KIND=8) :: Q(LDQ,N)
              COMPLEX(KIND=8) :: WORKL(N*(N+3))
              REAL(KIND=8) :: RWORK(N)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE ZNEIGH
          END INTERFACE 
        END MODULE ZNEIGH__genmod
