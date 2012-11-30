        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:33 2012
        MODULE SSAUP2__genmod
          INTERFACE 
            SUBROUTINE SSAUP2(IDO,BMAT,N,WHICH,NEV,NP,TOL,RESID,MODE,   &
     &IUPD,ISHIFT,MXITER,V,LDV,H,LDH,RITZ,BOUNDS,Q,LDQ,WORKL,IPNTR,WORKD&
     &,INFO)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: NEV
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IDO
              CHARACTER(LEN=1) :: BMAT
              CHARACTER(LEN=2) :: WHICH
              REAL(KIND=4) :: TOL
              REAL(KIND=4) :: RESID(N)
              INTEGER(KIND=4) :: MODE
              INTEGER(KIND=4) :: IUPD
              INTEGER(KIND=4) :: ISHIFT
              INTEGER(KIND=4) :: MXITER
              REAL(KIND=4) :: V(LDV,NEV+NP)
              REAL(KIND=4) :: H(LDH,2)
              REAL(KIND=4) :: RITZ(NEV+NP)
              REAL(KIND=4) :: BOUNDS(NEV+NP)
              REAL(KIND=4) :: Q(LDQ,NEV+NP)
              REAL(KIND=4) :: WORKL(3*(NEV+NP))
              INTEGER(KIND=4) :: IPNTR(3)
              REAL(KIND=4) :: WORKD(3*N)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE SSAUP2
          END INTERFACE 
        END MODULE SSAUP2__genmod
