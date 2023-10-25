  module time_data

    implicit none
    integer :: noh, ndigits
    integer :: ntimes, ncorr, nskip
    double precision :: dt, w_resol

  end module time_data

  module map_data

    implicit none
    double precision :: c0, c1, c2   ! freq map parameters
    double precision :: b0, b1, b2   ! mu' map parameters
    double precision :: d0, d1       ! x01 map parameters
    double precision :: a0, a1       ! alpha' map parameters
    double precision :: c15          ! ratio Raman para/perp
    double precision :: T1           ! vib. relax. time               
    double precision :: z_c          ! center of slab

  end module map_data

  module freq_data

    implicit none
    double precision :: w01avg, w01sqavg  ! avg freq and freq^2

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

