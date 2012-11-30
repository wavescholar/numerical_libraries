        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:32 2012
        MODULE SNAITR__genmod
          INTERFACE 
            SUBROUTINE SNAITR(IDO,BMAT,N,K,NP,NB,RESID,RNORM,V,LDV,H,LDH&
     &,IPNTR,WORKD,INFO)
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: K
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IDO
              CHARACTER(LEN=1) :: BMAT
              INTEGER(KIND=4) :: NB
              REAL(KIND=4) :: RESID(N)
              REAL(KIND=4) :: RNORM
              REAL(KIND=4) :: V(LDV,K+NP)
              REAL(KIND=4) :: H(LDH,K+NP)
              INTEGER(KIND=4) :: IPNTR(3)
              REAL(KIND=4) :: WORKD(3*N)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE SNAITR
          END INTERFACE 
        END MODULE SNAITR__genmod
