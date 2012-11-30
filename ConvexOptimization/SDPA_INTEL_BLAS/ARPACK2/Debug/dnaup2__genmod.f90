        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:29 2012
        MODULE DNAUP2__genmod
          INTERFACE 
            SUBROUTINE DNAUP2(IDO,BMAT,N,WHICH,NEV,NP,TOL,RESID,MODE,   &
     &IUPD,ISHIFT,MXITER,V,LDV,H,LDH,RITZR,RITZI,BOUNDS,Q,LDQ,WORKL,    &
     &IPNTR,WORKD,INFO)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: NEV
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IDO
              CHARACTER(LEN=1) :: BMAT
              CHARACTER(LEN=2) :: WHICH
              REAL(KIND=8) :: TOL
              REAL(KIND=8) :: RESID(N)
              INTEGER(KIND=4) :: MODE
              INTEGER(KIND=4) :: IUPD
              INTEGER(KIND=4) :: ISHIFT
              INTEGER(KIND=4) :: MXITER
              REAL(KIND=8) :: V(LDV,NEV+NP)
              REAL(KIND=8) :: H(LDH,NEV+NP)
              REAL(KIND=8) :: RITZR(NEV+NP)
              REAL(KIND=8) :: RITZI(NEV+NP)
              REAL(KIND=8) :: BOUNDS(NEV+NP)
              REAL(KIND=8) :: Q(LDQ,NEV+NP)
              REAL(KIND=8) :: WORKL((NEV+NP)*(NEV+NP+3))
              INTEGER(KIND=4) :: IPNTR(13)
              REAL(KIND=8) :: WORKD(3*N)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DNAUP2
          END INTERFACE 
        END MODULE DNAUP2__genmod
