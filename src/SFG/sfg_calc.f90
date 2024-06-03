PROGRAM SFG_CALC
! *********************************************************************
! This program calculates the SFG spectra using the empirical mapping
! approach developed by Skinner and co-workers. 
! 
! 
! Copyright, Zeke Piskulich, 2023.
! *********************************************************************

    USE map_data
    USE time_data
    USE freq_data
    USE hist_data
    USE input_module
    USE cli_data
    USE CLI
    USE output_module

    IMPLICIT NONE

    INTEGER, PARAMETER :: nperchunk = 1000
    INTEGER :: ioh, chunk, nchunks, iper
    INTEGER :: i

    DOUBLE PRECISION :: ta, tb, tstart, tend, read_time, tcf_time
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: w01, z0, a_ss, a_sp, a_pp
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: mu01, eOH
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: efield
    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: tcf, tcf_tot
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: multi_a_ss
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: spec_dist, sd_tot

    CALL Read_CLI_Arguments

! *********************************************************************
! I.  Read The Input File, and set the data for the calculations.
!     Also allocate the relevant variables.
! *********************************************************************
    CALL Read_Input

    CALL Apply_CLI_Args

    ALLOCATE(w01(nperchunk,ntimes)); ALLOCATE(mu01(nperchunk,ntimes,3))
    ALLOCATE(eOH(nperchunk,ntimes,3))
    ALLOCATE(z0(nperchunk,ntimes))
    ALLOCATE(efield(ntimes))
    
    ALLOCATE(a_ss(nperchunk, ntimes))
    ALLOCATE(a_sp(nperchunk, ntimes))
    ALLOCATE(a_pp(nperchunk, ntimes))
    ALLOCATE(multi_a_ss(ntimes,3))
    ALLOCATE(spec_dist(0:nhist)); ALLOCATE(sd_tot(0:nhist))
    

    ALLOCATE(tcf(0:ncorr)); ALLOCATE(tcf_tot(0:ncorr))

    multi_a_ss = 0.0d0
    tcf_tot = dcmplx(0.0d0, 0.0d0); read_time=0.0d0; tcf_time = 0.0d0
    spec_dist = 0.0d0; sd_tot = 0.0d0

    nchunks = ceiling(real(noh)/real(nperchunk))

    w01_avg = 0.0; w01_sq_avg = 0.0; tcf_tot = DCMPLX(0.0d0, 0.0d0)
! *********************************************************************
!   II. Calculate TCFS
! *********************************************************************
    DO chunk=1, nchunks

        w01 = 0.0; mu01 = 0.0; eOH = 0.0; efield = 0.0
        DO iper=1, nperchunk
            ioh = (chunk-1)*nperchunk + iper
            IF (ioh > noh) EXIT
            CALL Read_Field_File(ioh, efield(:), eOH(iper,:,:), z0(iper,:))
            CALL Get_w01(efield(:),  w01(iper,:))
            CALL Get_01_Dipole(efield(:), w01(iper,:),  eOH(iper,:,:), mu01(iper,:,:))
            CALL Get_Transtion_Pol_Polarizations(efield(:), w01(iper,:), eOH(iper,:,:),&
                                            &    a_ss(iper,:), a_sp(iper,:), a_pp(iper,:))
        END DO 

        DO ioh=(chunk-1)*nperchunk+1, min(chunk*nperchunk, noh)
            iper = ioh - (chunk-1)*nperchunk
            ! Calculate the Spectral Density
            multi_a_ss = 0.0d0
            DO i=1, 3
                multi_a_ss(:,i) = a_ss(iper,:)
            END DO
            CALL Spec_Dist_1D(w01(iper,:), mu01(iper,:,:), multi_a_ss(iper,:), spec_dist(:))
            sd_tot = sd_tot + spec_dist
            ! Some sort of histogramming?
            CALL Calc_TCF(w01(iper,:), mu01(iper,:,:), a_ss(iper,:), z0(iper,:), tcf(:))
            tcf_tot = tcf_tot + tcf

        ENDDO ! io

    ENDDO ! chunk

    ! normalize by the number of oh groups
    tcf_tot = tcf_tot/dcmplx(dfloat(noh),0d0)
    w01_avg = w01_avg/dfloat(noh*ntimes)
    w01_sq_avg = w01_sq_avg/dfloat(noh*ntimes)

! *********************************************************************
! III. Calculate the SFG Spectra
! *********************************************************************
    CALL Spectral_Print_1D(sd_tot, 'sfg_ss_spec_dist.dat')
    CALL Spec_Calc(tcf_tot)

! *********************************************************************
! IV. Cleanup Calculation
! *********************************************************************

    DEALLOCATE(w01)
    DEALLOCATE(mu01)
    DEALLOCATE(z0)
    DEALLOCATE(a_ss); DEALLOCATE(a_sp); DEALLOCATE(a_pp)
    DEALLOCATE(multi_a_ss)
    DEALLOCATE(spec_dist); DEALLOCATE(sd_tot)

END PROGRAM SFG_CALC