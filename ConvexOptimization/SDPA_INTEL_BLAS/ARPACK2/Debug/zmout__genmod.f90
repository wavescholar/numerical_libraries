        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:27 2012
        MODULE ZMOUT__genmod
          INTERFACE 
            SUBROUTINE ZMOUT(LOUT,M,N,A,LDA,IDIGIT,IFMT)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: IDIGIT
              CHARACTER(*) :: IFMT
            END SUBROUTINE ZMOUT
          END INTERFACE 
        END MODULE ZMOUT__genmod
