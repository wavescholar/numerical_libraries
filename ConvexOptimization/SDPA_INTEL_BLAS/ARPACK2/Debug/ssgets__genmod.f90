        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:31 2012
        MODULE SSGETS__genmod
          INTERFACE 
            SUBROUTINE SSGETS(ISHIFT,WHICH,KEV,NP,RITZ,BOUNDS,SHIFTS)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: KEV
              INTEGER(KIND=4) :: ISHIFT
              CHARACTER(LEN=2) :: WHICH
              REAL(KIND=4) :: RITZ(KEV+NP)
              REAL(KIND=4) :: BOUNDS(KEV+NP)
              REAL(KIND=4) :: SHIFTS(NP)
            END SUBROUTINE SSGETS
          END INTERFACE 
        END MODULE SSGETS__genmod
