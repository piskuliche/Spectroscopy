MODULE input_module
    IMPLICIT NONE

CONTAINS

    ! Here you can define procedures

    SUBROUTINE Read_Empirical_Map

    USE map_data
    USE constants


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

    CLOSE(10)
    END SUBROUTINE Read_Empirical_Map

    SUBROUTINE Read_Input
    ! *********************************************************************
    ! This subroutine reads the input file sfg.in and stores the data in the 
    ! data modules.
    ! 
    ! These data modules are time_data, map_data, hist_data, and constants.
    ! The code also does a number of important unit conversions for SFG.
    !
    ! ToDo: 
    !   - Make a generalized module for multiple types of spectroscopy
    ! *********************************************************************
        
        ! Data modules
        USE time_data
        USE map_data
        USE hist_data
        USE constants

        IMPLICIT NONE

        OPEN(10, file='spectra.in', status='old')

        READ(10,*)
        READ(10,*) noh, ndigits ! # OH groups, # of digits in field file extensions
        READ(10,*) 
        READ(10,*) ntimes, dt, ncorr, nskip ! # number of times, timestep, # corelation times, # skip times
        READ(10,*) 
        READ(10,*) w_resol ! freq resolution
        READ(10,*) 
        READ(10,*) flag_hist, nhist
        READ(10,*)
        READ(10,*) w1min, w1max, w3min, w3max ! Ranges for printing the spectra
        READ(10,*)
        READ(10,*) z_c ! center of slab. 
        READ(10,*) 
        READ(10,*) (Tw(j), j = 1, nTw)  ! Waiting time list (fs)

        CLOSE(10)
        ! Grab the Empriical Mapping Parameters
        CALL Read_Empirical_Map

        ! Covert time step to atomic units
        dt = dt/fsperau

        ! Convert T1 by dt
        T1bydt = T1rel/dt

        ! Convert freq. resolution to a.u.
        w_resol = w_resol/cmiperau

        ! Convert histogram parameters to a.u.
        w1min = w1min/cmiperau; w1max = w1max/cmiperau
        w3min = w3min/cmiperau; w3max = w3max/cmiperau

        dw = (w1max - w1min)/real(nhist)
    

    END SUBROUTINE Read_Input

END MODULE input_module