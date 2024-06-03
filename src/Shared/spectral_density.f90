MODULE output_module
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE Spec_Dist_1D(w01, weight1, weight2,  spectral_density)
    ! Calcualtes the spectral density using the dot product of two terms.
    ! 
    ! For the spectral density of different types of spectroscopy, different weights are used.
    ! List:
    !  1. IR: weight1 = dipole_moment, weight2 = dipole_moment
    !  2. Raman: weight1 = polarizability, weight2 = polarizability
    !  3. SFG: weight1 = dipole_moment, weight2 = polarizability
    !

        USE time_data
        USE hist_data
        USE constants

        IMPLICIT NONE

        INTEGER :: k, m, iw

        DOUBLE PRECISION :: count, dot_spec
        DOUBLE PRECISION, DIMENSION(ntimes,3) :: weight1, weight2
        DOUBLE PRECISION, DIMENSION(ntimes) :: w01
        DOUBLE PRECISION, DIMENSION(0:nhist) :: spectral_density

        ! Zero the Arrays 
        spectral_density = 0.0d0
        count = 0d0
        DO k=1, ntimes

            dot_spec = DOT_PRODUCT(weight1(k,:), weight2(k,:))
            iw = NINT ( (w01(k) - w1min)/dw )

            IF (iw >= 0 .and. iw <= nhist) THEN
                count = count + 1d0
                spectral_density(iw) = spectral_density(iw) + dot_spec
            END IF
        END DO

        ! Normalize the Spectral Density
        spec_density = spec_density/count

    END SUBROUTINE Spec_Dist_1D

    SUBROUTINE Freq_Dist_1D(w01, w01_dist)
    ! Calculates the frequency distribution using values of w01
    ! 
    ! Args:
    !   w01: Array of w01 values
    !   w01_dist: Array to store the frequency distribution
    !
    ! Returns:
    !   None
    ! 
        USE time_data
        USE hist_data
        USE constants

        IMPLICIT NONE

        INTEGER :: k, iw
        DOUBLE PRECISION :: count
        DOUBLE PRECISION, DIMENSION(ntimes) :: w01
        DOUBLE PRECISION, DIMENSION(0:nhist) :: w01_dist

        w01_dist = 0.0d0
        count = 0d0
        DO k=1, ntimes
            iw = NINT((w01(k) - w1min)/dw )

            IF (iw >= 0 .and. iw <= nhist) THEN
                count = count + 1d0
                w01_dist(iw) = w01_dist(iw) + 1d0
            END IF
        ENDDO
        w01_dist = w01_dist/count

    END SUBROUTINE Freq_Dist_1D

    SUBROUTINE Spectral_Print_1D(spectral_feature, output_file_name)
    ! Prints 1D spectra to a file, described by output_file_name.
    ! 
    ! Args:
    !   spectral: Array of spectral data
    !   output_file_name: Name of the output file
    !
    ! Returns:
    !   None
    !

        USE time_data
        USE hist_data
        USE constants
        USE cli_data

        IMPLICIT NONE

        INTEGER :: k
        DOUBLE PRECISION, DIMENSION(0:nhist) :: spectral_feature
        CHARACTER(LEN=100) :: output_file_name

        OPEN(24, FILE=trim(tag_output_cli)//trim(output_file_name))

        spectral_feature = spectral_feature/REAL(noh)

        DO k=0, nhist
            WRITE(24,*) (w1min + REAL(k)*dw)*cmiperau, spectral_feature(k)
        END DO

        CLOSE(24)


    END SUBROUTINE Spectral_Print_1D
END MODULE output_module