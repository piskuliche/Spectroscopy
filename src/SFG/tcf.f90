SUBROUTINE Calc_TCF(w, mu, a_ij, z0, tcf)

    use time_data

    IMPLICIT NONE

    ! Arguments
    DOUBLE PRECISION, DIMENSION(ntimes),    INTENT(IN) :: w     ! The frequency
    DOUBLE PRECISION, DIMENSION(ntimes,3),  INTENT(IN) :: mu    ! The dipole moment
    DOUBLE PRECISION, DIMENSION(ntimes),    INTENT(IN) :: a_ij  ! The polarizability
    DOUBLE PRECISION, DIMENSION(ntimes),    INTENT(IN) :: z0    ! The z position
    DOUBLE COMPLEX,   DIMENSION(0:ncorr),   INTENT(OUT) :: tcf  ! the Time Correlation Function

    ! Loop Counters
    INTEGER :: ilag
    INTEGER :: k, m

    ! Local Variables
    DOUBLE PRECISION :: count, dp, phase
    DOUBLE PRECISION :: a_ij0



    tcf = dcmplx(0.0d0, 0.0d0); count = 0d0

    DO k=1, ntimes-ncorr, nskip
        a_ij0 = a_ij(k)*dsign(1d0, z0(k))
        count = count + 1d0
        phase = -w(k)

        DO m=k, k + ncorr
            ilag = m - k
            phase = phase + w(m)
            dp = a_ij0*mu(m,3) !! put in the sign *and only works with a_ijp*
            tcf(ilag) = tcf(ilag) + dp*dcmplx(dcos(phase*dt), dsin(phase*dt))
        ENDDO
    ENDDO

    ! Normalize the tcf by the # time zeros
    tcf = tcf/dcmplx(count,0d0)



END SUBROUTINE Calc_TCF