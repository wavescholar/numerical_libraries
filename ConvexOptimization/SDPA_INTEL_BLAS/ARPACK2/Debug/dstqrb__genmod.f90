        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:25 2012
        MODULE DSTQRB__genmod
          INTERFACE 
            SUBROUTINE DSTQRB(N,D,E,Z,WORK,INFO)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: D(N)
              REAL(KIND=8) :: E(N-1)
              REAL(KIND=8) :: Z(N)
              REAL(KIND=8) :: WORK(2*N-2)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DSTQRB
          END INTERFACE 
        END MODULE DSTQRB__genmod
