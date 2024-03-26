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
    USE mapping

    IMPLICIT NONE

    OPEN(10, file='ir_spectra.in', status='old')

    READ(10,*)
    READ(10,*) noh, ndigits ! # OH groups, # of digits in field file extensions
    READ(10,*) 
    READ(10,*) ntimes, dt, ncorr, nskip ! # number of times, timestep, # corelation times, # skip times
    READ(10,*) 
    READ(10,*) w_resol ! freq resolution
    READ(10,*) 
    READ(10,*) flag_hist, nhist, wmin, wmax ! Calc freq_dist and spec dense, +params

    CLOSE(10)
    ! Grab the Empriical Mapping Parameters
    CALL Read_Empirical_Map

    ! Covert time step to atomic units
    dt = dt/fsperau

    ! Convert freq. resolution to a.u.
    w_resol = w_resol/cmiperau

    ! Convert histogram parameters to a.u.
    wmin = wmin/cmiperau
    wmax = wmax/cmiperau
    dw = (wmax - wmin)/real(nhist)
    


END SUBROUTINE Read_Input
