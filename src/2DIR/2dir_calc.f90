PROGRAM IR2D_CALC

USE map_data
USE time_data
USE freq_data
USE constants

IMPLICIT NONE

INTEGER, PARAMETER :: nperchunk = 1000

INTEGER :: chunk, iper, ioh, nchunks
INTEGER :: iTw

DOUBLE PRECISION :: tstart, tend, read_time, tcf_time, ta, tb

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: w01, mu01, w12, mu12
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: eOH

DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: tcf_rp, tcf_rp_tot
DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: tcf_np, tcf_np_tot

! *********************************************************************
! I.  Read The Input File, and set the data for the calculations.
!     Also allocate the relevant variables.
! *********************************************************************
CALL read_input

ALLOCATE(w01(nperchunk, ntimes)); ALLOCATE(mu01(nperchunk, ntimes))
ALLOCATE(w12(nperchunk, ntimes)); ALLOCATE(mu12(nperchunk, ntimes))

ALLOCATE(eOH(nperchunk, ntimes, 3))

ALLOCATE(tcf_rp(0:ncorr, 0:ncorr)); ALLOCATE(tcf_rp_tot(0:ncorr, 0:ncorr))
ALLOCATE(tcf_np(0:ncorr, 0:ncorr)); ALLOCATE(tcf_np_tot(0:ncorr, 0:ncorr))


nchunks = ceiling(real(noh)/real(nperchunk))

! *********************************************************************
! 2.  Read The Field, and Calculate the TCFs
! *********************************************************************

DO iTW = 1, nTw
    WRITE(6,*) 'Calculating for Tw = ', Tw(iTW)*fsperau, ' fs'
    CALL flush(6)
    ! Zero the Arrays
    tcf_rp_tot = DCMPLX(0.0d0, 0.0d0); tcf_np_tot = DCMPLX(0.0d0, 0.0d0)
    ! Set the reading time to zero
    read_time = 0.0d0; tcf_time = 0.0d0
    ! Set the average frequencies to zero
    w01_avg = 0.0d0; w01_sq_avg = 0.0d0
    w12_avg = 0.0d0; w12_sq_avg = 0.0d0

    CALL CPU_TIME(tstart)

    DO chunk=1, nchunks
        WRITE(6,*) 'Calculating for chunk = ', chunk, ' of ', nchunks
        ! ***** READ THE FIELD *****
        w01 = 0.0; w12 = 0.0; mu01 = 0.0; mu12 = 0.0; eOH = 0.0
        DO iper=1, nperchunk
            ioh = (chunk-1)*nperchunk + iper
            IF(ioh > noh) EXIT
            CALL read_field(ioh, w01(iper,:), w12(iper,:), mu01(iper,:),  mu12(iper,:), eOH(iper,:,:))
        ENDDO

        ! ***** CALCULATE THE TCFs *****
        DO ioh=(chunk-1)*nperchunk + 1, min(chunk*nperchunk, noh)
            iper = ioh - (chunk-1)*nperchunk
            CALL calc_tcf(ioh, iTw, w01(iper,:), w12(iper,:), mu01(iper,:), mu12(iper,:), eOH(iper,:,:), tcf_rp, tcf_np)

            ! Add the TCFS to the total
            tcf_rp_tot = tcf_rp_tot + tcf_rp
            tcf_np_tot = tcf_np_tot + tcf_np

            CALL CPU_TIME(tb)
            tcf_time = tcf_time + tb - ta
        ENDDO

        ! ***** Calculate Averages *****
        tcf_rp_tot = tcf_rp_tot/DCMPLX(DBLE(noh),0d0)
        tcf_np_tot = tcf_np_tot/DCMPLX(DBLE(noh),0d0)

        w01_avg = w01_avg/DBLE(noh*ntimes); w01_sq_avg = w01_sq_avg/DBLE(noh*ntimes)
        w12_avg = w12_avg/DBLE(noh*ntimes); w12_sq_avg = w12_sq_avg/DBLE(noh*ntimes)

        CALL CPU_TIME(ta)

        ! ***** Calculate the Spectra *****
        CALL spec_calc(iTw, tcf_rp_tot, tcf_np_tot)

        CALL CPU_TIME(tb)

        CALL CPU_TIME(tend)
        write(6,'(A,F12.2,A)') ' Waiting time = ',Tw(iTw)*fsperau,' fs'
        write(6,'(A,F12.2,A,F10.2,A)') ' Read  cpu time = ',read_time,' s, ',read_time/60d0,' min'
        write(6,'(A,F12.2,A,F10.2,A)') ' TCF   cpu time = ',tcf_time,' s, ',tcf_time/60d0,' min'
        write(6,'(A,F12.2,A,F10.2,A)') ' FFT   cpu time = ',tb-ta,' s, ',(tb-ta)/60d0,' min'
        write(6,'(A,F12.2,A,F10.2,A)') ' Total cpu time = ',tend-tstart,' s, ',(tend-tstart)/60d0,' min'
    ENDDO
ENDDO

DEALLOCATE(w01); DEALLOCATE(mu01)
DEALLOCATE(w12); DEALLOCATE(mu12)
DEALLOCATE(eOH)

DEALLOCATE(tcf_rp); DEALLOCATE(tcf_rp_tot)
DEALLOCATE(tcf_np); DEALLOCATE(tcf_np_tot)

END PROGRAM IR2D_CALC