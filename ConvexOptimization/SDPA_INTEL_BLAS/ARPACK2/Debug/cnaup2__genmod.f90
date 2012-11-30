        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:28 2012
        MODULE CNAUP2__genmod
          INTERFACE 
            SUBROUTINE CNAUP2(IDO,BMAT,N,WHICH,NEV,NP,TOL,RESID,MODE,   &
     &IUPD,ISHIFT,MXITER,V,LDV,H,LDH,RITZ,BOUNDS,Q,LDQ,WORKL,IPNTR,WORKD&
     &,RWORK,INFO)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: NEV
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IDO
              CHARACTER(LEN=1) :: BMAT
              CHARACTER(LEN=2) :: WHICH
              REAL(KIND=4) :: TOL
              COMPLEX(KIND=4) :: RESID(N)
              INTEGER(KIND=4) :: MODE
              INTEGER(KIND=4) :: IUPD
              INTEGER(KIND=4) :: ISHIFT
              INTEGER(KIND=4) :: MXITER
              COMPLEX(KIND=4) :: V(LDV,NEV+NP)
              COMPLEX(KIND=4) :: H(LDH,NEV+NP)
              COMPLEX(KIND=4) :: RITZ(NEV+NP)
              COMPLEX(KIND=4) :: BOUNDS(NEV+NP)
              COMPLEX(KIND=4) :: Q(LDQ,NEV+NP)
              COMPLEX(KIND=4) :: WORKL((NEV+NP)*(NEV+NP+3))
              INTEGER(KIND=4) :: IPNTR(13)
              COMPLEX(KIND=4) :: WORKD(3*N)
              REAL(KIND=4) :: RWORK(NEV+NP)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE CNAUP2
          END INTERFACE 
        END MODULE CNAUP2__genmod
