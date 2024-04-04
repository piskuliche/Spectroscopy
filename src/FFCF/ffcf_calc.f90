PROGRAM FFCF_CALC
! **********************************************************************
! *                                                                    *    
! *  This program calculates the FFCF for a given set of parameters.   *
! *  The FFCF is calculated for a given set of parameters.             *
! *                                                                    *
! **********************************************************************


    use map_data
    use time_data

    IMPLICIT NONE

    INTEGER, PARAMETER :: nperchunk=1000
    INTEGER :: chunk, iper, ioh, nchunks

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: w01
    DOUBLE PRECISION :: ti

    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: ffcf, ffcf_tot

    Call Read_Input

    ALLOCATE(w01(nperchunk,times))
    ALLOCATE(ffcf(0:ncorr)); ALLOCATE(ffcf_tot(0:ncorr))

    ffcf = 0.0d0; ffcf_tot = 0.0d0

    nchunks = CEILING(REAL(noh)/nperchunk)

    DO chunk=1, nchunks
        w01 = 0.0
        DO iper=1, nperchunk
            ioh = (chunk-1)*nperchunk + iper
            IF (ioh > noh) EXIT
            CALL Read_Field(ioh, w01(iper,:))
        ENDDO

        DO ioh=(chunk-1)*nperchunk+1, MIN(chunk*nperchunk, noh)
            iper = ioh - (chunk-1)*nperchunk
            CALL FFCF_TCF_CALC(w01(iper,:))
            ffcf_tot = ffcf_tot + ffcf
        ENDDO
    ENDDO

    ffcf_tot = ffcf_tot/DCMPLX(DFLOAT(noh),0d0)
    w01_avg = w01_avg/DFLOAT(noh*ntimes)

    OPEN(21, file='ffcf.dat')
    DO i=0, ncorr
        ti = float(i)*dt*fsperau
        WRITE(21,*) ti, ffcf_tot(i)
    ENDDO
    CLOSE(21)

    DEALLOCATE(w01)
    DEALLOCATE(ffcf); DEALLOCATE(ffcf_tot)
    
END PROGRAM FFCF_CALC

