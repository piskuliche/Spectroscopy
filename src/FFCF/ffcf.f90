SUBROUTINE FFCF_TCF_CALC(w, tcf)

USE time_data

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(ntimes), INTENT(IN) :: w
DOUBLE PRECISION, DIMENSION(0:ncorr), INTENT(OUT) :: tcf

INTEGER :: ilag, k, m

DOUBLE PRECISION :: count

DO k=1, ntimes-ncorr, nskip
    count = count + 1d0
    DO m = k, k + ncorr
        ilag = m - k
        tcf(ilag) = tcf(ilag) + w(k)*w(m)
    ENDDO
ENDDO
tcf = tcf / count

END SUBROUTINE FFCF_TCF_CALC