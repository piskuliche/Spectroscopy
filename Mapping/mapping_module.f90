MODULE mapping
    IMPLICIT NONE

CONTAINS

    ! Here you can define procedures

    SUBROUTINE Read_Empirical_Map

    IMPLICIT NONE

    OPEN(10, file='empirical_map.in')
    READ(10,*) 
    READ(10,*) c0, c1, c2
    READ(10,*)
    READ(10,*) c3, c4, c5
    READ(10,*)
    READ(10,*) b0, b1, b2
    READ(10,*)
    READ(10,*) d0, d1
    READ(10,*)
    READ(10,*) d2, d3
    READ(10,*) 
    READ(10,*) a0, a1
    READ(10,*) 
    READ(10,*) c15
    READ(10,*)
    READ(10,*) T1
        
    END SUBROUTINE Read_Empirical_Map

END MODULE mapping