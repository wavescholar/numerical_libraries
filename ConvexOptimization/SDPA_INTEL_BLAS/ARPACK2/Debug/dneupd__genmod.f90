        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:35 2012
        MODULE DNEUPD__genmod
          INTERFACE 
            SUBROUTINE DNEUPD(RVEC,HOWMNY,SELECT,DR,DI,Z,LDZ,SIGMAR,    &
     &SIGMAI,WORKEV,BMAT,N,WHICH,NEV,TOL,RESID,NCV,V,LDV,IPARAM,IPNTR,  &
     &WORKD,WORKL,LWORKL,INFO)
              INTEGER(KIND=4) :: LWORKL
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NCV
              INTEGER(KIND=4) :: NEV
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: LDZ
              LOGICAL(KIND=4) :: RVEC
              CHARACTER(LEN=1) :: HOWMNY
              LOGICAL(KIND=4) :: SELECT(NCV)
              REAL(KIND=8) :: DR(NEV+1)
              REAL(KIND=8) :: DI(NEV+1)
              REAL(KIND=8) :: Z(LDZ,*)
              REAL(KIND=8) :: SIGMAR
              REAL(KIND=8) :: SIGMAI
              REAL(KIND=8) :: WORKEV(3*NCV)
              CHARACTER(LEN=1) :: BMAT
              CHARACTER(LEN=2) :: WHICH
              REAL(KIND=8) :: TOL
              REAL(KIND=8) :: RESID(N)
              REAL(KIND=8) :: V(LDV,NCV)
              INTEGER(KIND=4) :: IPARAM(11)
              INTEGER(KIND=4) :: IPNTR(14)
              REAL(KIND=8) :: WORKD(3*N)
              REAL(KIND=8) :: WORKL(LWORKL)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DNEUPD
          END INTERFACE 
        END MODULE DNEUPD__genmod
