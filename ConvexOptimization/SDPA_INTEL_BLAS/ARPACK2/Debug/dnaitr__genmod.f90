        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:35 2012
        MODULE DNAITR__genmod
          INTERFACE 
            SUBROUTINE DNAITR(IDO,BMAT,N,K,NP,NB,RESID,RNORM,V,LDV,H,LDH&
     &,IPNTR,WORKD,INFO)
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: K
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IDO
              CHARACTER(LEN=1) :: BMAT
              INTEGER(KIND=4) :: NB
              REAL(KIND=8) :: RESID(N)
              REAL(KIND=8) :: RNORM
              REAL(KIND=8) :: V(LDV,K+NP)
              REAL(KIND=8) :: H(LDH,K+NP)
              INTEGER(KIND=4) :: IPNTR(3)
              REAL(KIND=8) :: WORKD(3*N)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DNAITR
          END INTERFACE 
        END MODULE DNAITR__genmod
