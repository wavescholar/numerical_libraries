        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:25 2012
        MODULE SSTQRB__genmod
          INTERFACE 
            SUBROUTINE SSTQRB(N,D,E,Z,WORK,INFO)
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: D(N)
              REAL(KIND=4) :: E(N-1)
              REAL(KIND=4) :: Z(N)
              REAL(KIND=4) :: WORK(2*N-2)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE SSTQRB
          END INTERFACE 
        END MODULE SSTQRB__genmod
