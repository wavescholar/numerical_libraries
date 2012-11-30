        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:29 2012
        MODULE SNAPPS__genmod
          INTERFACE 
            SUBROUTINE SNAPPS(N,KEV,NP,SHIFTR,SHIFTI,V,LDV,H,LDH,RESID,Q&
     &,LDQ,WORKL,WORKD)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: KEV
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: SHIFTR(NP)
              REAL(KIND=4) :: SHIFTI(NP)
              REAL(KIND=4) :: V(LDV,KEV+NP)
              REAL(KIND=4) :: H(LDH,KEV+NP)
              REAL(KIND=4) :: RESID(N)
              REAL(KIND=4) :: Q(LDQ,KEV+NP)
              REAL(KIND=4) :: WORKL(KEV+NP)
              REAL(KIND=4) :: WORKD(2*N)
            END SUBROUTINE SNAPPS
          END INTERFACE 
        END MODULE SNAPPS__genmod
