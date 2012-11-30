        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:32 2012
        MODULE SMOUT__genmod
          INTERFACE 
            SUBROUTINE SMOUT(LOUT,M,N,A,LDA,IDIGIT,IFMT)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: A(LDA,*)
              INTEGER(KIND=4) :: IDIGIT
              CHARACTER(*) :: IFMT
            END SUBROUTINE SMOUT
          END INTERFACE 
        END MODULE SMOUT__genmod
