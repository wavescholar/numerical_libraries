        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:23 2012
        MODULE DSGETS__genmod
          INTERFACE 
            SUBROUTINE DSGETS(ISHIFT,WHICH,KEV,NP,RITZ,BOUNDS,SHIFTS)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: KEV
              INTEGER(KIND=4) :: ISHIFT
              CHARACTER(LEN=2) :: WHICH
              REAL(KIND=8) :: RITZ(KEV+NP)
              REAL(KIND=8) :: BOUNDS(KEV+NP)
              REAL(KIND=8) :: SHIFTS(NP)
            END SUBROUTINE DSGETS
          END INTERFACE 
        END MODULE DSGETS__genmod
