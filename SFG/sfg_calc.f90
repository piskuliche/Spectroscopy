PROGRAM SFG_CALC
! *********************************************************************
! This program calculates the SFG spectra using the empirical mapping
! approach developed by Skinner and co-workers. 
! 
! 
! *********************************************************************

    USE map_data
    USE time_data
    USE freq_data
    USE hist_data

    IMPLICIT NONE

    INTEGER, PARAMETER :: nperchunk = 1000
    INTEGER :: ioh, chunk, nchunks, iper

    DOUBLE PRECISION :: ta, tb, tstart, tend, read_time, tcf_time
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: w01, z0, a_ss, a_sp, a_pp
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: mu, eOH
    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: tcf, tcf_tot


! *********************************************************************
! I.  Read The Input File, and set the data for the calculations.
!     Also allocate the relevant variables.
! *********************************************************************
    CALL Read_Input
    WRITE(*,*) nperchunk, nchunks

    ALLOCATE(w01(nperchunk,ntimes)); ALLOCATE(mu(nperchunk,ntimes,3))
    ALLOCATE(z0(nperchunk,ntimes))
    ALLOCATE(a_ss(nperchunk, ntimes))
    ALLOCATE(a_sp(nperchunk, ntimes))
    ALLOCATE(a_pp(nperchunk, ntimes))

    ALLOCATE(tcf(0:ncorr)); ALLOCATE(tcf_tot(0:ncorr))

    tcf_tot = dcmplx(0.0d0, 0.0d0); read_time=0.0d0; tcf_time = 0.0d0

! *********************************************************************
!   II. Calculate TCFS
! *********************************************************************
    DO chunk=1, nchunks

        w01 = 0.0; mu = 0.0; eOH = 0.0
        DO iper=1, nperchunk
            WRITE(*,*) iper
            ioh = (chunk-1)*nperchunk + iper
            IF (ioh > noh) EXIT
            CALL Read_Field(ioh, w01(iper,:), mu(iper,:,:), eOH(iper,:,:), &
                         a_ss(iper,:), a_sp(iper,:), a_pp(iper,:), z0(iper,:))
            WRITE(*,*) w01(iper,1) 
        END DO 

        DO ioh=1, noh
            ! Some sort of histogramming?
            CALL Calc_TCF(w01(iper,:), mu(iper,:,:), a_ss(iper,:), z0(iper,:), tcf(:))
            tcf_tot = tcf_tot + tcf
        ENDDO ! io

    ENDDO ! chunk

    ! normalize by the number of oh groups
    tcf_tot = tcf_tot/dcmplx(dfloat(noh),0d0)

! *********************************************************************
! III. Calculate the SFG Spectra
! *********************************************************************

    CALL Spec_Calc(tcf_tot)

! *********************************************************************
! IV. Cleanup Calculation
! *********************************************************************

    DEALLOCATE(w01)
    DEALLOCATE(mu)
    DEALLOCATE(z0)
    DEALLOCATE(a_ss); DEALLOCATE(a_sp); DEALLOCATE(a_pp)

END PROGRAM SFG_CALC