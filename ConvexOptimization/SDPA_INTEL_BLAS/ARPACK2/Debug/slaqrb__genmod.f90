        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 17 10:02:34 2012
        MODULE SLAQRB__genmod
          INTERFACE 
            SUBROUTINE SLAQRB(WANTT,N,ILO,IHI,H,LDH,WR,WI,Z,INFO)
              INTEGER(KIND=4) :: LDH
              LOGICAL(KIND=4) :: WANTT
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ILO
              INTEGER(KIND=4) :: IHI
              REAL(KIND=4) :: H(LDH,*)
              REAL(KIND=4) :: WR(*)
              REAL(KIND=4) :: WI(*)
              REAL(KIND=4) :: Z(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE SLAQRB
          END INTERFACE 
        END MODULE SLAQRB__genmod
