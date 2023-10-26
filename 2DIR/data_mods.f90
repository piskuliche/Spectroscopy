  module time_data

    implicit none
    integer :: noh, ndigits
    integer :: ntimes, ncorr, nskip
    double precision :: dt, w_resol
    integer :: nTw
    double precision, dimension(20) :: Tw
    double precision :: w1min, w1max, w3min, w3max

  end module time_data

  module map_data

    implicit none
    double precision :: c0, c1, c2   ! freq map parameters for 0-1
    double precision :: c3, c4, c5   ! freq map parameters for 1-2
    double precision :: b0, b1, b2   ! mu' map parameters
    double precision :: d0, d1       ! x01 map parameters
    double precision :: d2, d3       ! x12 map parameters
    double precision :: T1rel, T1bydt   ! vib. relax. time, vib relax time in units of timesteps

  end module map_data

  module freq_data

    implicit none
    double precision :: w01avg, w12avg, w01sqavg, w12sqavg  ! avg freq and freq^2

  end module freq_data

  module fluc_data

    implicit none
    logical :: flag_fluc
    double precision :: Havg(5), dH2avg(5)     ! Avg energy, Avg energy fluctuation size
    double precision :: dw01avg(5), dw01sqavg(5), dw12avg(5), dw12sqavg(5)    ! fluctuation avg freq and freq^2
    double precision :: d2w01avg(5), d2w01sqavg(5), d2w12avg(5), d2w12sqavg(5)  ! fluctuation^2 avg freq and freq^2

  end module fluc_data


  module constants

    implicit none
    double precision :: fsperau=2.4188843265857d-2
    double precision :: cmiperau=2.1947463d5
    double precision :: angperau=0.5291772083d0
    
  end module constants
