        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:22 2012
        MODULE DLAQRB__genmod
          INTERFACE 
            SUBROUTINE DLAQRB(WANTT,N,ILO,IHI,H,LDH,WR,WI,Z,INFO)
              INTEGER(KIND=4) :: LDH
              LOGICAL(KIND=4) :: WANTT
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ILO
              INTEGER(KIND=4) :: IHI
              REAL(KIND=8) :: H(LDH,*)
              REAL(KIND=8) :: WR(*)
              REAL(KIND=8) :: WI(*)
              REAL(KIND=8) :: Z(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DLAQRB
          END INTERFACE 
        END MODULE DLAQRB__genmod
