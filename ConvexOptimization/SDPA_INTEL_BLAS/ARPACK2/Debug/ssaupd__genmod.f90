        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:30 2012
        MODULE SSAUPD__genmod
          INTERFACE 
            SUBROUTINE SSAUPD(IDO,BMAT,N,WHICH,NEV,TOL,RESID,NCV,V,LDV, &
     &IPARAM,IPNTR,WORKD,WORKL,LWORKL,INFO)
              INTEGER(KIND=4) :: LWORKL
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NCV
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IDO
              CHARACTER(LEN=1) :: BMAT
              CHARACTER(LEN=2) :: WHICH
              INTEGER(KIND=4) :: NEV
              REAL(KIND=4) :: TOL
              REAL(KIND=4) :: RESID(N)
              REAL(KIND=4) :: V(LDV,NCV)
              INTEGER(KIND=4) :: IPARAM(11)
              INTEGER(KIND=4) :: IPNTR(11)
              REAL(KIND=4) :: WORKD(3*N)
              REAL(KIND=4) :: WORKL(LWORKL)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE SSAUPD
          END INTERFACE 
        END MODULE SSAUPD__genmod
