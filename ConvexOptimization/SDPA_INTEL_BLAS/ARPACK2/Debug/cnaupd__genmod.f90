        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:26 2012
        MODULE CNAUPD__genmod
          INTERFACE 
            SUBROUTINE CNAUPD(IDO,BMAT,N,WHICH,NEV,TOL,RESID,NCV,V,LDV, &
     &IPARAM,IPNTR,WORKD,WORKL,LWORKL,RWORK,INFO)
              INTEGER(KIND=4) :: LWORKL
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NCV
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IDO
              CHARACTER(LEN=1) :: BMAT
              CHARACTER(LEN=2) :: WHICH
              INTEGER(KIND=4) :: NEV
              REAL(KIND=4) :: TOL
              COMPLEX(KIND=4) :: RESID(N)
              COMPLEX(KIND=4) :: V(LDV,NCV)
              INTEGER(KIND=4) :: IPARAM(11)
              INTEGER(KIND=4) :: IPNTR(14)
              COMPLEX(KIND=4) :: WORKD(3*N)
              COMPLEX(KIND=4) :: WORKL(LWORKL)
              REAL(KIND=4) :: RWORK(NCV)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE CNAUPD
          END INTERFACE 
        END MODULE CNAUPD__genmod
