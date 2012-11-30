        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:32 2012
        MODULE ZNAPPS__genmod
          INTERFACE 
            SUBROUTINE ZNAPPS(N,KEV,NP,SHIFT,V,LDV,H,LDH,RESID,Q,LDQ,   &
     &WORKL,WORKD)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: KEV
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: SHIFT(NP)
              COMPLEX(KIND=8) :: V(LDV,KEV+NP)
              COMPLEX(KIND=8) :: H(LDH,KEV+NP)
              COMPLEX(KIND=8) :: RESID(N)
              COMPLEX(KIND=8) :: Q(LDQ,KEV+NP)
              COMPLEX(KIND=8) :: WORKL(KEV+NP)
              COMPLEX(KIND=8) :: WORKD(2*N)
            END SUBROUTINE ZNAPPS
          END INTERFACE 
        END MODULE ZNAPPS__genmod
