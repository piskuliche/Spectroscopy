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
    double precision :: T1           ! vib. relax. time               

  end module map_data

  module freq_data

    implicit none
    double precision :: w01_avg, w01_sq_avg  ! avg freq and freq^2

  end module freq_data

  module hist_data

    implicit none
    logical :: flag_hist
    integer :: nhist
    double precision :: dw, wmin, wmax

  end module hist_data

  module fluc_data

    implicit none
    logical :: flag_fluc, flag_gg
    double precision :: Havg(5), dH2avg(5)     ! Avg energy, Avg energy fluctuation size
    double precision :: Coul_gg_avg(6), LJ_gg_avg(6)
    double precision :: dwavg(5), dw2avg(5)    ! fluctuation avg freq and freq^2
    double precision :: d2wavg(5), d2w2avg(5)  ! fluctuation^2 avg freq and freq^2

  end module fluc_data


  module constants

    implicit none
    double precision :: fsperau=2.4188843265857d-2
    double precision :: cmiperau=2.1947463d5
    double precision :: angperau=0.5291772083d0
    
  end module constants

