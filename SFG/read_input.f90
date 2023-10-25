!   Subroutine to read in the input parameters from spectra.in

  Subroutine read_input

  use time_data
  use map_data
  use hist_data
  use fluc_data
  use constants

  implicit none
 ! include 'constants.f'

  open(10,file='spectra.in',status='old')

  read(10,*)
  read(10,*) noh, ndigits    ! # of OH groups, # of digits in the field file extension
  read(10,*)
  read(10,*) ntimes, dt, ncorr, nskip
  read(10,*)
  read(10,*) w_resol      ! freq resolution (used to pad TCF w/ zeroes) 
  read(10,*)
  read(10,*) c0, c1, c2   ! frequency map
  read(10,*)
  read(10,*) b0, b1, b2   ! mu' map
  read(10,*)
  read(10,*) d0, d1       ! x01 map
  read(10,*)
  read(10,*) a0, a1       ! alpha' map
  read(10,*)
  read(10,*) c15          ! ratio Raman para/perp
  read(10,*) 
  read(10,*) T1           ! vib. relax. time
  read(10,*)
  read(10,*) flag_hist, nhist, wmin, wmax  ! Calculate freq dist & spec dens? + parameters
  read(10,*)
  read(10,*) flag_fluc    ! Calculate the derivative using fluc theory?
  read(10,*)
  read(10,*) z_c          ! Center of the slab


  ! Convert time step to atomic units
  dt = dt/fsperau
  
  ! Convert freq. resol to atomic units
  w_resol = w_resol/cmiperau

  ! Convert map parameters to atomic units
  c0 = c0/cmiperau; c1 = c1/cmiperau; c2 = c2/cmiperau
  
  ! Convert vib. relax. time to atomic units
  T1 = T1/fsperau

  ! Convert x01 parameters
  d0 = d0/angperau; d1 = d1*cmiperau/angperau

  ! Convert histogram parameters to atomic units
  wmin = wmin/cmiperau; wmax = wmax/cmiperau
  dw = (wmax - wmin)/real(nhist)
 
  End Subroutine read_input
