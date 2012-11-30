        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:22 2012
        MODULE CGETV0__genmod
          INTERFACE 
            SUBROUTINE CGETV0(IDO,BMAT,ITRY,INITV,N,J,V,LDV,RESID,RNORM,&
     &IPNTR,WORKD,IERR)
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: J
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IDO
              CHARACTER(LEN=1) :: BMAT
              INTEGER(KIND=4) :: ITRY
              LOGICAL(KIND=4) :: INITV
              COMPLEX(KIND=4) :: V(LDV,J)
              COMPLEX(KIND=4) :: RESID(N)
              REAL(KIND=4) :: RNORM
              INTEGER(KIND=4) :: IPNTR(3)
              COMPLEX(KIND=4) :: WORKD(2*N)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE CGETV0
          END INTERFACE 
        END MODULE CGETV0__genmod
