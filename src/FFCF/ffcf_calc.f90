PROGRAM FFCF_CALC
! *********************************************************************
! This program calculates the FFCF using the empirical mapping
! approach developed by Skinner and co-workers. 
! 
! 
! Copyright, Zeke Piskulich, 2024.
!
! Adapted with permission from the original code by:
! Dr. Ashley Borkowski, Dr. Hasini Senanayake, and Prof. Ward Thompson
! *********************************************************************


    USE map_data
    USE time_data
    USE constants
    USE freq_data
    USE input_module
    USE cli_data
    USE CLI
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: nperchunk=1000

    ! Variables for Loop Control
    INTEGER :: chunk, iper, ioh, nchunks
    INTEGER :: i

    ! Variables for Data Storage
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)     :: efield
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: z0
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: eOH

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: w01
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ffcf, ffcf_tot
    DOUBLE PRECISION :: save_w01_avg

    ! Other Variables
    DOUBLE PRECISION :: ti
    LOGICAL :: cli_freq

    ! ///////////////////////////////////////////////////////////////////
    
    cli_freq = .False.

    CALL Read_CLI_Arguments

    WRITE(*,*) 'Calculating FFCF for the given parameters'
    ! Read Input Parameters
    CALL Read_Input

    ! Apply the CLI Arg Overrides
    CALL Apply_CLI_Args

    ! Check if frequency was read from the command line
    IF (avfreq_cli /= -1) THEN
        cli_freq = .TRUE.
        save_w01_avg = avfreq_cli
    ELSE
        IF (flag_z_range) THEN
            WRITE(*,*) 'Must provide frequency when using z-range'
            STOP
        ENDIF
    ENDIF


    ! Allocate Memory
    ALLOCATE(w01(nperchunk,ntimes))
    ALLOCATE(efield(ntimes)); ALLOCATE(eOH(nperchunk,ntimes,3))
    ALLOCATE(z0(nperchunk,ntimes))
    ALLOCATE(ffcf(0:ncorr)); ALLOCATE(ffcf_tot(0:ncorr))

    ! Initialize Data Storage
    ffcf = 0.0d0; ffcf_tot = 0.0d0

    nchunks = CEILING(REAL(noh)/nperchunk)
    

    IF (.NOT. cli_freq) THEN
        ! Loop over all the chunks
        DO chunk=1, nchunks
            w01 = 0.0
            DO iper=1, nperchunk
                ioh = (chunk-1)*nperchunk + iper
                IF (ioh > noh) EXIT
                CALL Read_Field_File(ioh, efield(:), eOH(iper,:,:), z0(iper,:))
                CALL Get_w01(efield(:),  w01(iper,:))
            ENDDO
        ENDDO
        save_w01_avg = w01_avg/DFLOAT(noh*ntimes)
    ENDIF

    WRITE(*,*) 'Average w01 = ', save_w01_avg*cmiperau
    OPEN(21, file=trim(tag_output_cli)//'w01_avg.dat')
    WRITE(21,*) save_w01_avg*cmiperau
    CLOSE(21)

    w01_avg = 0.0

    DO chunk=1, nchunks
        w01 = 0.0
        DO iper=1, nperchunk
            ioh = (chunk-1)*nperchunk + iper
            IF (ioh > noh) EXIT
            CALL Read_Field_File(ioh, efield(:), eOH(iper,:,:), z0(iper,:))
            CALL Get_w01(efield(:),  w01(iper,:))
        ENDDO

        DO ioh=(chunk-1)*nperchunk+1, MIN(chunk*nperchunk, noh)
            iper = ioh - (chunk-1)*nperchunk
            CALL FFCF_TCF_CALC(w01(iper,:)-save_w01_avg, ffcf(:))
            ffcf_tot = ffcf_tot + ffcf(:)
        ENDDO
    ENDDO

    ffcf_tot = ffcf_tot/DFLOAT(noh)

    ! Write the FFCF to a file
    OPEN(21, file=trim(tag_output_cli)//'ffcf.dat')
    DO i=0, ncorr
        ti = float(i)*dt*fsperau
        WRITE(21,*) ti, ffcf_tot(i)
    ENDDO
    ! Close the ffcf file
    CLOSE(21)

    DEALLOCATE(w01)
    DEALLOCATE(ffcf); DEALLOCATE(ffcf_tot)
    
END PROGRAM FFCF_CALC

