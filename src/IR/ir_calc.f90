PROGRAM IR_CALC
! *********************************************************************
! This program calculates the infrared spectra using the empirical mapping
! approach developed by Skinner and co-workers. 
! 
! The equation for the IR spectra is given by:
! 
! 
! Copyright, Zeke Piskulich, 2024. All rights reserved.
! *********************************************************************

    USE map_data
    USE time_data
    USE freq_data
    USE hist_data
    USE input_module
    USE cli_data
    USE CLI
    IMPLICIT NONE

    INTEGER, PARAMETER :: nperchunk = 1000
    INTEGER :: chunk, iper, ioh, nchunks

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: w01
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: mu01, eOH
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: w01_dist, spec_dist, wd_tot, sd_tot

    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: tcf, tcf_tot

    CALL Read_CLI_Arguments

! *********************************************************************
! I.  Read The Input File, and set the data for the calculations.
!     Also allocate the relevant variables.
! *********************************************************************

    CALL Read_Input

    CALL Apply_CLI_Args

    ALLOCATE(w01(nperchunk,ntimes))
    ALLOCATE(mu01(nperchunk,ntimes,3))
    ALLOCATE(eOH(nperchunk,ntimes,3))
    ALLOCATE(z0(nperchunk,ntimes))
    ALLOCATE(efield(ntimes))

    ALLOCATE(tcf(0:ncorr)); ALLOCATE(tcf_tot(0:ncorr))

    ALLOCATE(w01_dist(0:nhist)); ALLOCATE(spec_dist(0:nhist))
    ALLOCATE(wd_tot(0:nhist)); ALLOCATE(sd_tot(0:nhist))

    w01_dist = 0.0d0; spec_dist = 0.0d0; wd_tot = 0.0d0; sd_tot = 0.0d0

    tcf_tot = dcmplx(0.0d0, 0.0d0)

    nchunks = ceiling(real(noh)/real(nperchunk))


! *********************************************************************
!   II. Calculate TCFS
! *********************************************************************

    DO chunk=1, nchunks
        w01 = 0.0; w12 = 0.0; mu01 = 0.0; mu12 = 0.0; eOH = 0.0
        DO iper=1, nperchunk
            ioh = (chunk-1)*nperchunk + iper
            if (ioh > noh) EXIT
            CALL Read_Field_File(ioh, efield(:), eOH(iper,:,:), z0(iper,:))
            CALL Get_w01(efield(:),  w01(iper,:))
            CALL Get_01_Dipole(efield(:), w01(iper,:),  eOH(iper,:,:), mu01(iper,:,:))
        END DO

        DO ioh=(chunk-1)*nperchunk+1, min(chunk*nperchunk, noh)
                iper = ioh - (chunk-1)*nperchunk
                ! Calculate the Spectral Density and Histogram
                CALL Hist_Calc(w01(iper,:), mu01(iper,:,:), w01_dist(:), spec_dist(:))
                wd_tot = wd_tot + w01_dist
                sd_tot = sd_tot + spec_dist
                ! Calculate the TCFS
                CALL Calc_TCF(w01(iper,:), mu01(iper,:,:), tcf(:))
                tcf_tot = tcf_tot + tcf
        ENDDO !
    END DO

    ! Normalize by the number of oh groups
    tcf_tot = tcf_tot/dcmplx(dfloat(noh),0d0)
    w01_avg = w01_avg/dfloat(noh*ntimes)
    w01_sq_avg = w01_sq_avg/dfloat(noh*ntimes)

! *********************************************************************
! III. Calculate the IR Spectra
! *********************************************************************
    CALL Hist_Print(wd_tot, sd_tot)

    CALL Spec_Calc(tcf_tot)

! *********************************************************************
! IV. Cleanup Calculation
! *********************************************************************

    DEALLOCATE(w01)
    DEALLOCATE(mu01)

END PROGRAM IR_CALC