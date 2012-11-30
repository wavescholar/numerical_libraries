        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:31 2012
        MODULE CNAPPS__genmod
          INTERFACE 
            SUBROUTINE CNAPPS(N,KEV,NP,SHIFT,V,LDV,H,LDH,RESID,Q,LDQ,   &
     &WORKL,WORKD)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: KEV
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=4) :: SHIFT(NP)
              COMPLEX(KIND=4) :: V(LDV,KEV+NP)
              COMPLEX(KIND=4) :: H(LDH,KEV+NP)
              COMPLEX(KIND=4) :: RESID(N)
              COMPLEX(KIND=4) :: Q(LDQ,KEV+NP)
              COMPLEX(KIND=4) :: WORKL(KEV+NP)
              COMPLEX(KIND=4) :: WORKD(2*N)
            END SUBROUTINE CNAPPS
          END INTERFACE 
        END MODULE CNAPPS__genmod
