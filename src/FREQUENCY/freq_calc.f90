PROGRAM FFCF_CALC
! *********************************************************************
! This program calculates the average frequency using the empirical mapping
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

    ! Other Variables
    DOUBLE PRECISION    :: save_w01_avg, w01_sum
    INTEGER             :: count_w01

    ! ///////////////////////////////////////////////////////////////////


    CALL Read_CLI_Arguments

    WRITE(*,*) 'Calculating Frequency'
    ! Read Input Parameters
    CALL Read_Input

    ! Apply the CLI Arg Overrides
    CALL Apply_CLI_Args

    ! Check if frequency was read from the command line

    ! Allocate Memory
    ALLOCATE(w01(nperchunk,ntimes))
    ALLOCATE(efield(ntimes)); ALLOCATE(eOH(nperchunk,ntimes,3))
    ALLOCATE(z0(nperchunk,ntimes))

    nchunks = CEILING(REAL(noh)/nperchunk)
    w01_sum = 0.0
    count_w01 = 0

    ! Loop over all the chunks
    DO chunk=1, nchunks
        w01 = 0.0
        DO iper=1, nperchunk
            ioh = (chunk-1)*nperchunk + iper
            IF (ioh > noh) EXIT
            CALL Read_Field_File(ioh, efield(:), eOH(iper,:,:), z0(iper,:))
            CALL Get_w01(efield(:),  w01(iper,:))

            IF (flag_z_range) THEN
                DO i=1, ntimes
                    IF (ABS(z0(iper,i)) > zmin_cli .AND. ABS(z0(iper,i)) < zmax_cli) THEN
                        w01_sum = w01_sum + w01(iper,i)
                        count_w01 = count_w01 + 1
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    save_w01_avg = w01_sum / REAL(count_w01) * cmiperau
    write(*,*) 'Average w01 = ', save_w01_avg
    

    ! Write the frequency to a file
    OPEN(21, file=trim(tag_output_cli)//'_w01.dat')
        WRITE(21,*) save_w01_avg
    CLOSE(21)

    DEALLOCATE(w01)
    DEALLOCATE(efield); DEALLOCATE(eOH); DEALLOCATE(z0)
    
END PROGRAM FFCF_CALC

