        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:33 2012
        MODULE CNGETS__genmod
          INTERFACE 
            SUBROUTINE CNGETS(ISHIFT,WHICH,KEV,NP,RITZ,BOUNDS)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: KEV
              INTEGER(KIND=4) :: ISHIFT
              CHARACTER(LEN=2) :: WHICH
              COMPLEX(KIND=4) :: RITZ(KEV+NP)
              COMPLEX(KIND=4) :: BOUNDS(KEV+NP)
            END SUBROUTINE CNGETS
          END INTERFACE 
        END MODULE CNGETS__genmod
