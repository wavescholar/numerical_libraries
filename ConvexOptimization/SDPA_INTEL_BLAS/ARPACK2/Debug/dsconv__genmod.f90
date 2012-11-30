        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:33 2012
        MODULE DSCONV__genmod
          INTERFACE 
            SUBROUTINE DSCONV(N,RITZ,BOUNDS,TOL,NCONV)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: RITZ(N)
              REAL(KIND=8) :: BOUNDS(N)
              REAL(KIND=8) :: TOL
              INTEGER(KIND=4) :: NCONV
            END SUBROUTINE DSCONV
          END INTERFACE 
        END MODULE DSCONV__genmod
