        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:32 2012
        MODULE DNGETS__genmod
          INTERFACE 
            SUBROUTINE DNGETS(ISHIFT,WHICH,KEV,NP,RITZR,RITZI,BOUNDS,   &
     &SHIFTR,SHIFTI)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: KEV
              INTEGER(KIND=4) :: ISHIFT
              CHARACTER(LEN=2) :: WHICH
              REAL(KIND=8) :: RITZR(KEV+NP)
              REAL(KIND=8) :: RITZI(KEV+NP)
              REAL(KIND=8) :: BOUNDS(KEV+NP)
              REAL(KIND=8) :: SHIFTR(1)
              REAL(KIND=8) :: SHIFTI(1)
            END SUBROUTINE DNGETS
          END INTERFACE 
        END MODULE DNGETS__genmod
