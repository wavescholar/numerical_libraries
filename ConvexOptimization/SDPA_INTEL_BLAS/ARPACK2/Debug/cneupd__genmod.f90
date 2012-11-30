        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:21 2012
        MODULE CNEUPD__genmod
          INTERFACE 
            SUBROUTINE CNEUPD(RVEC,HOWMNY,SELECT,D,Z,LDZ,SIGMA,WORKEV,  &
     &BMAT,N,WHICH,NEV,TOL,RESID,NCV,V,LDV,IPARAM,IPNTR,WORKD,WORKL,    &
     &LWORKL,RWORK,INFO)
              INTEGER(KIND=4) :: LWORKL
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: NCV
              INTEGER(KIND=4) :: NEV
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: LDZ
              LOGICAL(KIND=4) :: RVEC
              CHARACTER(LEN=1) :: HOWMNY
              LOGICAL(KIND=4) :: SELECT(NCV)
              COMPLEX(KIND=4) :: D(NEV)
              COMPLEX(KIND=4) :: Z(LDZ,NEV)
              COMPLEX(KIND=4) :: SIGMA
              COMPLEX(KIND=4) :: WORKEV(2*NCV)
              CHARACTER(LEN=1) :: BMAT
              CHARACTER(LEN=2) :: WHICH
              REAL(KIND=4) :: TOL
              COMPLEX(KIND=4) :: RESID(N)
              COMPLEX(KIND=4) :: V(LDV,NCV)
              INTEGER(KIND=4) :: IPARAM(11)
              INTEGER(KIND=4) :: IPNTR(14)
              COMPLEX(KIND=4) :: WORKD(3*N)
              COMPLEX(KIND=4) :: WORKL(LWORKL)
              REAL(KIND=4) :: RWORK(NCV)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE CNEUPD
          END INTERFACE 
        END MODULE CNEUPD__genmod
