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

        WRITE(*,*) "ncorr_cli = ", ncorr_cli
        WRITE(*,*) "nskip_cli = ", nskip_cli
        WRITE(*,*) "input_file_cli = ", TRIM(input_file_cli)
        WRITE(*,*) "map_file_cli = ", TRIM(map_file_cli)
        WRITE(*,*) "tag_output_cli = ", TRIM(tag_output_cli)
        WRITE(*,*) "avfreq_cli = ", avfreq_cli
        WRITE(*,*) "flag_z_range = ", flag_z_range
        WRITE(*,*) "zmin_cli = ", zmin_cli
        WRITE(*,*) "zmax_cli = ", zmax_cli

    END SUBROUTINE

    SUBROUTINE Print_Help
        WRITE(*,*) "Usage: ./program [options]"
        WRITE(*,*) "Options:"
        WRITE(*,*) "  -h, --help: Print this help message"
        WRITE(*,*) "  -ncorr: Correlation Length"
        WRITE(*,*) "  -nskip: Number of configs between origins"
        WRITE(*,*) "  -in: Input file"
        WRITE(*,*) "  -map: Map file"
        WRITE(*,*) "  -tag: Tag output"
        STOP
    END SUBROUTINE

END MODULE CLI