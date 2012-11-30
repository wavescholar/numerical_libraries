        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:25 2012
        MODULE SNGETS__genmod
          INTERFACE 
            SUBROUTINE SNGETS(ISHIFT,WHICH,KEV,NP,RITZR,RITZI,BOUNDS,   &
     &SHIFTR,SHIFTI)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: KEV
              INTEGER(KIND=4) :: ISHIFT
              CHARACTER(LEN=2) :: WHICH
              REAL(KIND=4) :: RITZR(KEV+NP)
              REAL(KIND=4) :: RITZI(KEV+NP)
              REAL(KIND=4) :: BOUNDS(KEV+NP)
              REAL(KIND=4) :: SHIFTR(1)
              REAL(KIND=4) :: SHIFTI(1)
            END SUBROUTINE SNGETS
          END INTERFACE 
        END MODULE SNGETS__genmod
