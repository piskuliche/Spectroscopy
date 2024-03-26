MODULE mapping
    IMPLICIT NONE

CONTAINS

    ! Here you can define procedures

    SUBROUTINE Read_Empirical_Map(filename)

    IMPLICIT NONE

    OPEN(10, file='empirical_map.in')
    READ(10,*) 
    READ(10,*) c0, c1, c2 ! w10 = c0 + c1*E + c2*E^2
    READ(10,*)
    READ(10,*) c3, c4, c5 ! w21 = c3 + c4*E + c5*E^2
    READ(10,*)
    READ(10,*) b0, b1, b2 ! mu' = b0 + b1*E + b2*E^2
    READ(10,*)
    READ(10,*) d0, d1 ! x10 = d0 + d1*E
    READ(10,*)
    READ(10,*) d2, d3 ! x21 = d2 + d3*E
    READ(10,*) 
    READ(10,*) a0, a1 ! alpha = a0 + a1*E
    READ(10,*) 
    READ(10,*) c15 ! alpha_par/alpha_perp, typically 5.6
    READ(10,*)
    READ(10,*) T1 ! Vibrational Relaxation Time

    ! Conver map parameters to atomic units
    c0 = c0/cmiperau; c1 = c1/cmiperau; c2 = c2/cmiperau
    c3 = c3/cmiperau; c4 = c4/cmiperau; c5 = c5/cmiperau

    ! Convert vib. relax time to atomic units
    T1rel = T1rel/fsperau

    ! Convert x01 paramters to atomic units
    d0 = d0/angperau; d1 = d1*cmiperau/angperau
    d2 = d2/angperau; d3 = d3*cmiperau/angperau

        
    END SUBROUTINE Read_Empirical_Map

END MODULE mapping