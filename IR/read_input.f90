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

    OPEN(10, file='ir_spectra.in', status='old')

    READ(10,*)
    READ(10,*) noh, ndigits ! # OH groups, # of digits in field file extensions
    READ(10,*) 
    READ(10,*) ntimes, dt, ncorr, nskip ! # number of times, timestep, # corelation times, # skip times
    READ(10,*) 
    READ(10,*) w_resol ! freq resolution
    READ(10,*) 
    READ(10,*) c0, c1, c2 ! frequency map
    READ(10,*) 
    READ(10,*) b0, b1, b2 ! mu' map
    READ(10,*) 
    READ(10,*) d0, d1 ! x01 map
    READ(10,*)
    READ(10,*) a0, a1 ! alpha' map
    READ(10,*) 
    READ(10,*) c15 ! ratio Raman para/perp
    REAd(10,*)
    READ(10,*) T1 ! vib. relaxation time
    READ(10,*) 
    READ(10,*) flag_hist, nhist, wmin, wmax ! Calc freq_dist and spec dense, +params
    READ(10,*)
    READ(10,*) z_c ! center of slab. 

    ! Covert time step to atomic units
    dt = dt/fsperau

    ! Convert freq. resolution to a.u.
    w_resol = w_resol/cmiperau

    ! Convert map parameters to atomic units
    c0 = c0/cmiperau; c1 = c1/cmiperau; c2 = c2/cmiperau

    ! Convert vib. relax time to a.u.
    T1 = T1/fsperau

    ! Convert x01 parameters
    d0 = d0/angperau
    d1 = d1*cmiperau/angperau

    ! Convert histogram parameters to a.u.
    wmin = wmin/cmiperau
    wmax = wmax/cmiperau
    dw = (wmax - wmin)/real(nhist)
    
    CLOSE(10)

END SUBROUTINE Read_Input
