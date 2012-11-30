        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:21 2012
        MODULE SSAPPS__genmod
          INTERFACE 
            SUBROUTINE SSAPPS(N,KEV,NP,SHIFT,V,LDV,H,LDH,RESID,Q,LDQ,   &
     &WORKD)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: KEV
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: SHIFT(NP)
              REAL(KIND=4) :: V(LDV,KEV+NP)
              REAL(KIND=4) :: H(LDH,2)
              REAL(KIND=4) :: RESID(N)
              REAL(KIND=4) :: Q(LDQ,KEV+NP)
              REAL(KIND=4) :: WORKD(2*N)
            END SUBROUTINE SSAPPS
          END INTERFACE 
        END MODULE SSAPPS__genmod
