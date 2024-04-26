MODULE CLI
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE Read_CLI_Arguments
        USE cli_data
        IMPLICIT NONE

        INTEGER :: num_args
        CHARACTER(len=100) :: arg

        INTEGER :: i

        ! Default values
        ncorr_cli = -1
        nskip_cli = -1
        input_file_cli = "spectra.in"
        map_file_cli = "empirical_map.in"
        tag_output_cli = ""
        avfreq_cli = -1
        flag_z_range = .FALSE.
        zmin_cli = -1
        zmax_cli = -1

        num_args = COMMAND_ARGUMENT_COUNT()

        DO i = 1, num_args
            CALL GET_COMMAND_ARGUMENT(i, arg)
            SELECT CASE (arg)
                CASE ("-h")
                    CALL Print_Help
                CASE ("-ncorr")
                    CALL GET_COMMAND_ARGUMENT(i+1, arg)
                    READ(arg, *) ncorr_cli 
                CASE ("-nskip")
                    CALL GET_COMMAND_ARGUMENT(i+1, arg)
                    READ(arg, *) nskip_cli
                CASE ("-in")
                    CALL GET_COMMAND_ARGUMENT(i+1, arg)
                    READ(arg, *) input_file_cli
                CASE ("-map")
                    CALL GET_COMMAND_ARGUMENT(i+1, arg)
                    READ(arg, *) map_file_cli
                CASE ("-tag")
                    CALL GET_COMMAND_ARGUMENT(i+1, arg)
                    READ(arg, *) tag_output_cli
                CASE ("-avfreq")
                    CALL GET_COMMAND_ARGUMENT(i+1, arg)
                    READ(arg, *) avfreq_cli
                CASE ("-zmin")
                    CALL GET_COMMAND_ARGUMENT(i+1, arg)
                    READ(arg, *) zmin_cli
                    flag_z_range = .TRUE.
                CASE ("-zmax")
                    CALL GET_COMMAND_ARGUMENT(i+1, arg)
                    READ(arg, *) zmax_cli
                    flag_z_range = .TRUE.
            END SELECT
        END DO
        WRITE(*,*) "-----------------------"
        WRITE(*,*) "Command Line Arguments:"
        WRITE(*,*) "-----------------------"
        WRITE(*,*) "ncorr_cli = ", ncorr_cli
        WRITE(*,*) "nskip_cli = ", nskip_cli
        WRITE(*,*) "input_file_cli = ", TRIM(ADJUSTL(input_file_cli))
        WRITE(*,*) "map_file_cli = ", TRIM(ADJUSTL(map_file_cli))
        WRITE(*,*) "tag_output_cli = ", TRIM(ADJUSTL(tag_output_cli))
        WRITE(*,*) "avfreq_cli = ", avfreq_cli
        WRITE(*,*) "flag_z_range = ", flag_z_range
        WRITE(*,*) "zmin_cli = ", zmin_cli
        WRITE(*,*) "zmax_cli = ", zmax_cli
        WRITE(*,*) "-----------------------"

        IF (avfreq_cli == -1) THEN
            IF (flag_z_range) THEN
                WRITE(*,*) 'Must provide frequency when using z-range'
                STOP
            ENDIF
        ENDIF

    END SUBROUTINE

    SUBROUTINE Print_Help
    ! Prints help message
        WRITE(*,*) "Usage: ./program [options]"
        WRITE(*,*) "Options:"
        WRITE(*,*) "  -h, --help: Print this help message"
        WRITE(*,*) "  -ncorr: Correlation Length"
        WRITE(*,*) "  -nskip: Number of configs between origins"
        WRITE(*,*) "  -in: Input file"
        WRITE(*,*) "  -map: Map file"
        WRITE(*,*) "  -tag: Tag output"
        WRITE(*,*) "  -avfreq: Average frequency [cm^-1]"
        WRITE(*,*) "  -zmin: Minimum redshift"
        WRITE(*,*) "  -zmax: Maximum redshift"
        STOP
    END SUBROUTINE

    SUBROUTINE Apply_CLI_Args
    ! Apply command line arguments to the global variables
    ! This subroutine should be called after Read_CLI_Arguments
    ! to apply the values to the global variables
    ! To override input file variables, use AFTER Read_Input.
        USE cli_data
        USE constants
        USE time_data
        IMPLICIT NONE

        ! Correlation time
        IF (ncorr_cli /= -1) THEN
            ncorr = ncorr_cli
        ENDIF

        ! Time Origin Skips
        IF (nskip_cli /= -1) THEN
            nskip = nskip_cli
        ENDIF

        ! Average Frequency
        IF (avfreq_cli /= -1) THEN
            avfreq_cli = avfreq_cli/cmiperau
        ENDIF
        
    END SUBROUTINE

END MODULE CLI