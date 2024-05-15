SUBROUTINE Hist_Calc(w01, mu, w01_dist, spec_dist)
    
    USE time_data
    USE hist_data
    USE constants

    IMPLICIT NONE

    INTEGER :: k, m, iw

    DOUBLE PRECISION :: count, mu_sq
    DOUBLE PRECISION, DIMENSION(ntimes,3) :: mu
    DOUBLE PRECISION, DIMENSION(ntimes) :: w01
    DOUBLE PRECISION, DIMENSION(0:nhist) :: w01_dist, spec_dist

    ! Zero the Arrays 
    w01_dist = 0.0d0; spec_dist = 0.0d0
    count = 0d0
    DO k=1, ntimes
        mu_sq = DOT_PRODUCT(mu(k,:), mu(k,:))
        iw = NINT ( (w01(k) - w1min)/dw )

        IF (iw >= 0 .and. iw <= nhist) THEN
            count = count + 1d0
            w01_dist(iw) = w01_dist(iw) + 1d0
            spec_dist(iw) = spec_dist(iw) + mu_sq
        END IF
    END DO

    ! Normalize the distributions
    w01_dist = w01_dist/count
    spec_dist = spec_dist/count

END SUBROUTINE Hist_Calc

SUBROUTINE Hist_Print(wd_tot, sd_tot)

    USE time_data
    USE hist_data
    USE constants
    USE cli_data

    IMPLICIT NONE

    INTEGER :: k
    DOUBLE PRECISION, DIMENSION(0:nhist) :: wd_tot, sd_tot

    OPEN(23, FILE=trim(tag_output_cli)//'freq_dist.dat')
    OPEN(24, FILE=trim(tag_output_cli)//'spec_dens.dat')

    wd_tot = wd_tot/REAL(noh)
    sd_tot = sd_tot/REAL(noh)

    DO k=0, nhist
        WRITE(23,*) (w1min + REAL(k)*dw)*cmiperau, wd_tot(k)
        WRITE(24,*) (w1min + REAL(k)*dw)*cmiperau, sd_tot(k)
    END DO
    CLOSE(23)
    CLOSE(24)


END SUBROUTINE Hist_Print