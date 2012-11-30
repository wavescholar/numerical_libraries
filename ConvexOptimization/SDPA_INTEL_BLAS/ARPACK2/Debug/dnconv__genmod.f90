        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:34 2012
        MODULE DNCONV__genmod
          INTERFACE 
            SUBROUTINE DNCONV(N,RITZR,RITZI,BOUNDS,TOL,NCONV)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: RITZR(N)
              REAL(KIND=8) :: RITZI(N)
              REAL(KIND=8) :: BOUNDS(N)
              REAL(KIND=8) :: TOL
              INTEGER(KIND=4) :: NCONV
            END SUBROUTINE DNCONV
          END INTERFACE 
        END MODULE DNCONV__genmod
