        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:34 2012
        MODULE ZNGETS__genmod
          INTERFACE 
            SUBROUTINE ZNGETS(ISHIFT,WHICH,KEV,NP,RITZ,BOUNDS)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: KEV
              INTEGER(KIND=4) :: ISHIFT
              CHARACTER(LEN=2) :: WHICH
              COMPLEX(KIND=8) :: RITZ(KEV+NP)
              COMPLEX(KIND=8) :: BOUNDS(KEV+NP)
            END SUBROUTINE ZNGETS
          END INTERFACE 
        END MODULE ZNGETS__genmod
