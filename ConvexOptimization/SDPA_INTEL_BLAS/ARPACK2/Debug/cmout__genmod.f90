        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:28 2012
        MODULE CMOUT__genmod
          INTERFACE 
            SUBROUTINE CMOUT(LOUT,M,N,A,LDA,IDIGIT,IFMT)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=4) :: A(LDA,*)
              INTEGER(KIND=4) :: IDIGIT
              CHARACTER(*) :: IFMT
            END SUBROUTINE CMOUT
          END INTERFACE 
        END MODULE CMOUT__genmod
