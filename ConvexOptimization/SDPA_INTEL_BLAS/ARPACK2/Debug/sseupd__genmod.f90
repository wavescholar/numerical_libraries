        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:34 2012
        MODULE SSEUPD__genmod
          INTERFACE 
            SUBROUTINE SSEUPD(RVEC,HOWMNY,SELECT,D,Z,LDZ,SIGMA,BMAT,N,  &
     &WHICH,NEV,TOL,RESID,NCV,V,LDV,IPARAM,IPNTR,WORKD,WORKL,LWORKL,INFO&
     &)
              INTEGER(KIND=4) :: LWORKL
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NCV
              INTEGER(KIND=4) :: NEV
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: LDZ
              LOGICAL(KIND=4) :: RVEC
              CHARACTER(LEN=1) :: HOWMNY
              LOGICAL(KIND=4) :: SELECT(NCV)
              REAL(KIND=4) :: D(NEV)
              REAL(KIND=4) :: Z(LDZ,NEV)
              REAL(KIND=4) :: SIGMA
              CHARACTER(LEN=1) :: BMAT
              CHARACTER(LEN=2) :: WHICH
              REAL(KIND=4) :: TOL
              REAL(KIND=4) :: RESID(N)
              REAL(KIND=4) :: V(LDV,NCV)
              INTEGER(KIND=4) :: IPARAM(7)
              INTEGER(KIND=4) :: IPNTR(11)
              REAL(KIND=4) :: WORKD(2*N)
              REAL(KIND=4) :: WORKL(LWORKL)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE SSEUPD
          END INTERFACE 
        END MODULE SSEUPD__genmod
