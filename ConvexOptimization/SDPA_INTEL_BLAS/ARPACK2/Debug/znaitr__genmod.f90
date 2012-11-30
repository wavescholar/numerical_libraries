        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:20 2012
        MODULE ZNAITR__genmod
          INTERFACE 
            SUBROUTINE ZNAITR(IDO,BMAT,N,K,NP,NB,RESID,RNORM,V,LDV,H,LDH&
     &,IPNTR,WORKD,INFO)
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: K
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IDO
              CHARACTER(LEN=1) :: BMAT
              INTEGER(KIND=4) :: NB
              COMPLEX(KIND=8) :: RESID(N)
              REAL(KIND=8) :: RNORM
              COMPLEX(KIND=8) :: V(LDV,K+NP)
              COMPLEX(KIND=8) :: H(LDH,K+NP)
              INTEGER(KIND=4) :: IPNTR(3)
              COMPLEX(KIND=8) :: WORKD(3*N)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE ZNAITR
          END INTERFACE 
        END MODULE ZNAITR__genmod
