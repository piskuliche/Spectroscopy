PROGRAM FFCF_CALC
! **********************************************************************
! *                                                                    *    
! *  This program calculates the FFCF for a given set of parameters.   *
! *  The FFCF is calculated for a given set of parameters.             *
! *                                                                    *
! *  Copyright (C) 2024, Zeke Piskulich, All Rights Reserved           *
! **********************************************************************


    USE map_data
    USE time_data
    USE constants
    USE freq_data
    USE input_module
    USE CLI
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: nperchunk=1000

    ! Variables for Loop Control
    INTEGER :: chunk, iper, ioh, nchunks
    INTEGER :: i

    ! Variables for Data Storage
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: w01
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ffcf, ffcf_tot
    DOUBLE PRECISION :: save_w01_avg

    ! Other Variables
    DOUBLE PRECISION :: ti
    LOGICAL :: file_exists

    ! ///////////////////////////////////////////////////////////////////
    
    file_exists = .False.

    CALL Read_CLI_Arguments

    WRITE(*,*) 'Calculating FFCF for the given parameters'
    ! Read Input Parameters
    CALL Read_Input

    INQUIRE(FILE='ffcf.in', EXIST=file_exists)
    IF (file_exists) THEN
        OPEN(9, file='ffcf.in', status='old', action='read')
        READ(9,*) 
        READ(9,*) ncorr
        READ(9,*)
        READ(9,*) save_w01_avg
        CLOSE(9)
    ENDIF

    ! Allocate Memory
    ALLOCATE(w01(nperchunk,ntimes))
    ALLOCATE(ffcf(0:ncorr)); ALLOCATE(ffcf_tot(0:ncorr))

    ! Initialize Data Storage
    ffcf = 0.0d0; ffcf_tot = 0.0d0

    nchunks = CEILING(REAL(noh)/nperchunk)
    

    IF (.NOT. file_exists) THEN
        ! Loop over all the chunks
        DO chunk=1, nchunks
            w01 = 0.0
            DO iper=1, nperchunk
                ioh = (chunk-1)*nperchunk + iper
                IF (ioh > noh) EXIT
                CALL Read_Field(ioh, w01(iper,:))
            ENDDO
        ENDDO
        save_w01_avg = w01_avg/DFLOAT(noh*ntimes)
    ENDIF

    write(*,*) 'Average w01 = ', save_w01_avg
    w01_avg = 0.0

    DO chunk=1, nchunks
        w01 = 0.0
        DO iper=1, nperchunk
            ioh = (chunk-1)*nperchunk + iper
            IF (ioh > noh) EXIT
            CALL Read_Field(ioh, w01(iper,:))
        ENDDO

        DO ioh=(chunk-1)*nperchunk+1, MIN(chunk*nperchunk, noh)
            iper = ioh - (chunk-1)*nperchunk
            CALL FFCF_TCF_CALC(w01(iper,:)-save_w01_avg, ffcf(:))
            ffcf_tot = ffcf_tot + ffcf(:)
        ENDDO
    ENDDO

    ffcf_tot = ffcf_tot/DFLOAT(noh)

    ! Write the FFCF to a file
    OPEN(21, file='ffcf.dat')
    DO i=0, ncorr
        ti = float(i)*dt*fsperau
        WRITE(21,*) ti, ffcf_tot(i)
    ENDDO
    ! Close the ffcf file
    CLOSE(21)

    DEALLOCATE(w01)
    DEALLOCATE(ffcf); DEALLOCATE(ffcf_tot)
    
END PROGRAM FFCF_CALC

