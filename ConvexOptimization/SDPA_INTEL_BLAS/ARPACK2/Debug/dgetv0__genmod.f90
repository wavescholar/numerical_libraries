        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:25 2012
        MODULE DGETV0__genmod
          INTERFACE 
            SUBROUTINE DGETV0(IDO,BMAT,ITRY,INITV,N,J,V,LDV,RESID,RNORM,&
     &IPNTR,WORKD,IERR)
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: J
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IDO
              CHARACTER(LEN=1) :: BMAT
              INTEGER(KIND=4) :: ITRY
              LOGICAL(KIND=4) :: INITV
              REAL(KIND=8) :: V(LDV,J)
              REAL(KIND=8) :: RESID(N)
              REAL(KIND=8) :: RNORM
              INTEGER(KIND=4) :: IPNTR(3)
              REAL(KIND=8) :: WORKD(2*N)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE DGETV0
          END INTERFACE 
        END MODULE DGETV0__genmod
