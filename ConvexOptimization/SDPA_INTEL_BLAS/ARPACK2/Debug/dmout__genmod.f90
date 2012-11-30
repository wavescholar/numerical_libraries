        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:24 2012
        MODULE DMOUT__genmod
          INTERFACE 
            SUBROUTINE DMOUT(LOUT,M,N,A,LDA,IDIGIT,IFMT)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: IDIGIT
              CHARACTER(*) :: IFMT
            END SUBROUTINE DMOUT
          END INTERFACE 
        END MODULE DMOUT__genmod
