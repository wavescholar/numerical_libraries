        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:27 2012
        MODULE SNEUPD__genmod
          INTERFACE 
            SUBROUTINE SNEUPD(RVEC,HOWMNY,SELECT,DR,DI,Z,LDZ,SIGMAR,    &
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
              REAL(KIND=4) :: DR(NEV+1)
              REAL(KIND=4) :: DI(NEV+1)
              REAL(KIND=4) :: Z(LDZ,*)
              REAL(KIND=4) :: SIGMAR
              REAL(KIND=4) :: SIGMAI
              REAL(KIND=4) :: WORKEV(3*NCV)
              CHARACTER(LEN=1) :: BMAT
              CHARACTER(LEN=2) :: WHICH
              REAL(KIND=4) :: TOL
              REAL(KIND=4) :: RESID(N)
              REAL(KIND=4) :: V(LDV,NCV)
              INTEGER(KIND=4) :: IPARAM(11)
              INTEGER(KIND=4) :: IPNTR(14)
              REAL(KIND=4) :: WORKD(3*N)
              REAL(KIND=4) :: WORKL(LWORKL)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE SNEUPD
          END INTERFACE 
        END MODULE SNEUPD__genmod
