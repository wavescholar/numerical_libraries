        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:26 2012
        MODULE DSAPPS__genmod
          INTERFACE 
            SUBROUTINE DSAPPS(N,KEV,NP,SHIFT,V,LDV,H,LDH,RESID,Q,LDQ,   &
     &WORKD)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: KEV
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: SHIFT(NP)
              REAL(KIND=8) :: V(LDV,KEV+NP)
              REAL(KIND=8) :: H(LDH,2)
              REAL(KIND=8) :: RESID(N)
              REAL(KIND=8) :: Q(LDQ,KEV+NP)
              REAL(KIND=8) :: WORKD(2*N)
            END SUBROUTINE DSAPPS
          END INTERFACE 
        END MODULE DSAPPS__genmod
