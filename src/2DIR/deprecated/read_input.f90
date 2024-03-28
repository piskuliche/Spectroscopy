!   Subroutine to read in the input parameters from spectra.in

  Subroutine read_input

  use time_data
  use map_data
  use constants

  implicit none
  integer :: j

  open(10,file='2dir.in',status='old')

  read(10,*)
  read(10,*) noh, ndigits    ! # of OH groups, # of digits in the field file extension
  read(10,*)
  read(10,*) ntimes, dt, ncorr, nskip
  read(10,*)
  read(10,*) w_resol      ! freq resolution (used to pad TCF w/ zeroes) 
  read(10,*)
  read(10,*) c0, c1, c2   ! w01 frequency map  (units of cm^-1 for freq, au for field)
  read(10,*)
  read(10,*) c3, c4, c5   ! w12 frequency maps (units of cm^-1 for freq, au for field)
  read(10,*)
  read(10,*) b0, b1, b2   ! mu' map (units of au for field)
  read(10,*)
  read(10,*) d0, d1       ! x01 map (units of cm^-1 for freq, au for x)
  read(10,*)
  read(10,*) d2, d3       ! x12 map (units of cm^-1 for freq, au for x)
  read(10,*)
  read(10,*) T1rel        ! vib. relax. time (in fs)
  read(10,*)
  read(10,*) nTw          ! # of waiting times
  if(nTw.gt.20) then
     write(6,*) ' No more than 20 waiting times allowed'
     stop
  endif
  read(10,*)  
  read(10,*) (Tw(j), j = 1, nTw)        ! waiting time list (in fs)
  read(10,*)
  read(10,*) w1min, w1max, w3min, w3max      ! Ranges in which spectrum will be printed (cm^-1) 
  read(10,*)
  !read(10,*) flag_fluc    ! Calculate the derivative using fluc theory?

  ! Convert time step to atomic units
  dt = dt/fsperau

  ! Convert freq. resol to atomic units
  w_resol = w_resol/cmiperau

  ! Convert map parameters to atomic units
  c0 = c0/cmiperau; c1 = c1/cmiperau; c2 = c2/cmiperau
  c3 = c3/cmiperau; c4 = c4/cmiperau; c5 = c5/cmiperau
  
  ! Convert vib. relax. time to atomic units
  T1rel = T1rel/fsperau
  T1bydt = T1rel/dt

  ! Convert waiting times to atomic units
  Tw = Tw/fsperau

  ! Convert x01 parameters
  d0 = d0/angperau; d1 = d1*cmiperau/angperau; d2 = d2/angperau; d3 = d3*cmiperau/angperau

  ! Convert frequency ranges
  w1min = w1min/cmiperau; w1max = w1max/cmiperau; w3min = w3min/cmiperau; w3max = w3max/cmiperau

  close(10)

  End Subroutine read_input

