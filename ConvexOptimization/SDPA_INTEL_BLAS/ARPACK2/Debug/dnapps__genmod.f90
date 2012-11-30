        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:32 2012
        MODULE DNAPPS__genmod
          INTERFACE 
            SUBROUTINE DNAPPS(N,KEV,NP,SHIFTR,SHIFTI,V,LDV,H,LDH,RESID,Q&
     &,LDQ,WORKL,WORKD)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: KEV
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: SHIFTR(NP)
              REAL(KIND=8) :: SHIFTI(NP)
              REAL(KIND=8) :: V(LDV,KEV+NP)
              REAL(KIND=8) :: H(LDH,KEV+NP)
              REAL(KIND=8) :: RESID(N)
              REAL(KIND=8) :: Q(LDQ,KEV+NP)
              REAL(KIND=8) :: WORKL(KEV+NP)
              REAL(KIND=8) :: WORKD(2*N)
            END SUBROUTINE DNAPPS
          END INTERFACE 
        END MODULE DNAPPS__genmod
