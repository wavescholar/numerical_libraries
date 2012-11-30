        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:34 2012
        MODULE CNAITR__genmod
          INTERFACE 
            SUBROUTINE CNAITR(IDO,BMAT,N,K,NP,NB,RESID,RNORM,V,LDV,H,LDH&
     &,IPNTR,WORKD,INFO)
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: K
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IDO
              CHARACTER(LEN=1) :: BMAT
              INTEGER(KIND=4) :: NB
              COMPLEX(KIND=4) :: RESID(N)
              REAL(KIND=4) :: RNORM
              COMPLEX(KIND=4) :: V(LDV,K+NP)
              COMPLEX(KIND=4) :: H(LDH,K+NP)
              INTEGER(KIND=4) :: IPNTR(3)
              COMPLEX(KIND=4) :: WORKD(3*N)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE CNAITR
          END INTERFACE 
        END MODULE CNAITR__genmod
