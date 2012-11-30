        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:22 2012
        MODULE SGETV0__genmod
          INTERFACE 
            SUBROUTINE SGETV0(IDO,BMAT,ITRY,INITV,N,J,V,LDV,RESID,RNORM,&
     &IPNTR,WORKD,IERR)
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: J
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IDO
              CHARACTER(LEN=1) :: BMAT
              INTEGER(KIND=4) :: ITRY
              LOGICAL(KIND=4) :: INITV
              REAL(KIND=4) :: V(LDV,J)
              REAL(KIND=4) :: RESID(N)
              REAL(KIND=4) :: RNORM
              INTEGER(KIND=4) :: IPNTR(3)
              REAL(KIND=4) :: WORKD(2*N)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE SGETV0
          END INTERFACE 
        END MODULE SGETV0__genmod
