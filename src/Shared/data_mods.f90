  module time_data

    implicit none
    integer :: noh, ndigits
    integer :: ntimes, ncorr, nskip
    double precision :: dt, w_resol
    integer :: nTw ! Used for 2D IR Only
    double precision, dimension(20) :: Tw ! Used for 2D IR Only
    double precision :: w1min, w1max, w3min, w3max ! Used for 2D IR Only

  end module time_data

  module map_data

    implicit none
    double precision :: c0, c1, c2   ! freq map parameters for 0-1
    double precision :: c3, c4, c5   ! freq map parameters for 1-2
    double precision :: b0, b1, b2   ! mu' map parameters
    double precision :: d0, d1       ! x01 map parameters
    double precision :: d2, d3       ! x12 map parameters
    double precision :: a0, a1       ! alpha' map parameters
    double precision :: c15          ! ratio Raman para/perp
    double precision :: T1 ! Unconverted T1
    double precision :: T1rel, T1bydt   ! vib. relax. time, vib relax time in units of timesteps
    double precision :: z_c          ! center of slab
    
  end module map_data

  module freq_data

    implicit none

    double precision :: w01_avg, w01_sq_avg
    double precision :: w12_avg, w12_sq_avg

  end module freq_data

  module hist_data

    implicit none
    logical :: flag_hist
    integer :: nhist
    double precision :: dw, wmin, wmax

  end module hist_data

  module constants

    implicit none
    double precision :: fsperau=2.4188843265857d-2
    double precision :: cmiperau=2.1947463d5
    double precision :: angperau=0.5291772083d0
    
  end module constants

