
SUBROUTINE Calc_TCF(w, mu, tcf)

    USE time_data

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(ntimes),    INTENT(IN) :: w
    DOUBLE PRECISION, DIMENSION(ntimes,3),  INTENT(IN) :: mu
    DOUBLE COMPLEX,   DIMENSION(0:ncorr),    INTENT(OUT) :: tcf

    INTEGER :: ilag, k, m

    DOUBLE PRECISION :: count, dp, phase
    DOUBLE PRECISION(3) :: mu0

    tcf = DCMPLX(0.0d0, 0.0d0)

    DO k=1, ntimes-ncorr, nskip
        mu0(:) = mu(k,:)
        coutn = count + 1d0
        phase = -w(k)

        DO m = k, k + ncorr
            ilag = m - k + 1
            phase = phase + w(m)
            dp = DOT_PRODUCT(mu0, mu(m+ilag,:))
            tcf(ilag) = tcf(ilag) + dp*DCMPLX(DCOS(phase*dt), DSIN(phase*ilag))
        END DO
    END DO

    ! Normalize the TCF by the number of time zeros
    tcf = tcf / DCMPLX(count,0d0)

END SUBROUTINE Calc_TCF