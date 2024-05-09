! Currently not implemented or used.

SUBROUTINE Read_Field_File(ioh, efield, eOH, z0)
! This subroutine reads the field files
! 
! This program requires the HDF5 library for reading.
! 
! Copyright, Zeke Piskulich, 2017
! 
    USE time_data
    USE map_data
    USE freq_data
    USE HDF5
    USE ieee_arithmetic, ONLY: ieee_is_finite

    IMPLICIT NONE

    INTEGER,                                INTENT(IN) :: ioh
    DOUBLE PRECISION, DIMENSION(ntimes),    INTENT(OUT) :: efield
    DOUBLE PRECISION, DIMENSION(ntimes,3),  INTENT(OUT) :: eOH
    DOUBLE PRECISION, DIMENSION(ntimes),    INTENT(OUT) :: z0

    INTEGER :: k
    REAL, DIMENSION(ntimes) :: etmp, z0_tmp
    REAL, DIMENSION(ntimes, 3) :: eoh_tmp

    ! HDF5 Variables
    CHARACTER(len=8), PARAMETER :: filename = 'field.h5'
    INTEGER(HID_T) :: file_id, dataset_id
    INTEGER :: ERROR_FLAG

    INTEGER(HSIZE_T), DIMENSION(1) :: dot_dims 
    INTEGER(HSIZE_T), DIMENSION(2) :: eoh_dims 

    CHARACTER(len=20) :: dataset_name

    dot_dims = (/ntimes/)
    eoh_dims = (/3, ntimes/)


!I. Read the field file *******************************************************
    ! Open the hdf5 library
    CALL h5open_f(ERROR_FLAG)
    ! Open the hdf5 file
    CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ERROR_FLAG)
    ! Read the value of the dot product ***
    ! Set dataset name
    WRITE(dataset_name, '(A, I0)') 'dot_', ioh

    ! Open the existing dataset
    CALL h5dopen_f(file_id, dataset_name, dataset_id, ERROR_FLAG)
    ! Read the dataset
    CALL h5dread_f(dataset_id, H5T_NATIVE_REAL, etmp, dot_dims, ERROR_FLAG)
    ! Close the dataset
    CALL h5dclose_f(dataset_id, ERROR_FLAG)

    ! Read the value of the unit vector eoh_tmp ***
    ! Set dataset name
    WRITE(dataset_name, '(A, I0)') 'eoh_', ioh

    ! Open the existing dataset
    CALL h5dopen_f(file_id, dataset_name, dataset_id, ERROR_FLAG)
    ! Read the dataset
    CALL h5dread_f(dataset_id, H5T_NATIVE_REAL, eoh_tmp, eoh_dims, ERROR_FLAG)
    ! Close the dataset
    CALL h5dclose_f(dataset_id, ERROR_FLAG)


    ! Set dataset name
    WRITE(dataset_name, '(A, I0)') 'z0_', ioh

    ! Open the existing dataset
    CALL h5dopen_f(file_id, dataset_name, dataset_id, ERROR_FLAG)
    ! Read the dataset
    CALL h5dread_f(dataset_id, H5T_NATIVE_REAL, z0_tmp, dot_dims, ERROR_FLAG)
    ! Close the dataset
    CALL h5dclose_f(dataset_id, ERROR_FLAG)

    ! Close the file
    CALL h5fclose_f(file_id, ERROR_FLAG)
    ! Close the library
    CALL h5close_f(ERROR_FLAG)

! II. Calculate the field parameters ******************************************
    DO k=1,ntimes
        z0(k) = z0_tmp(k) - z_c
        eOH(k,:) = eoh_tmp(k,:)
        efield(k) = etmp(k)
    ENDDO
    WRITE(*,*) eoh_tmp(1,1), eoh_tmp(1,2), eoh_tmp(1,3)
    WRITE(*,*) SHAPE(eoh_tmp)
    WRITE(*,*) eoh_tmp(2,1), eoh_tmp(2,2), eoh_tmp(2,3)
    STOP
END SUBROUTINE Read_Field_File

SUBROUTINE Get_w01(efield, w01)
    USE freq_data
    USE time_data
    USE map_data
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(IN) :: efield
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(OUT) :: w01

    INTEGER :: k
    DO k=1, ntimes
        w01(k) = c0 + c1*efield(k) + c2*efield(k)**2
        w01_avg = w01_avg + w01(k)
        w01_sq_avg = w01_sq_avg + w01(k)**2
    ENDDO
END SUBROUTINE Get_w01

SUBROUTINE Get_w12(efield, w12)
    USE freq_data
    USE time_data
    USE map_data
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(IN) :: efield
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(OUT) :: w12

    INTEGER :: k

    DO k=1, ntimes
        ! ** 21 Transition Frequencies ***
        w12(k) = c3 + c4*efield(k) + c5*efield(k)**2     
        w12_avg = w12_avg + w12(k)
        w12_sq_avg = w12_sq_avg + w12(k)**2
    ENDDO
END SUBROUTINE Get_w12

SUBROUTINE Get_mu01_Prime(efield, w01, mu01prime)
    USE freq_data
    USE time_data
    USE map_data
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(IN) :: efield, w01
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(OUT) :: mu01prime

    ! Loop Variables
    INTEGER :: k

    ! Temporary Variables
    DOUBLE PRECISION :: x01tmp
    DOUBLE PRECISION :: muprime

    DO k=1, ntimes
        ! *** Position Matrix Elements ***
        x01tmp = d0 + d1*w01(k)

        ! *** Transition Dipole Matrix Elements ***
        muprime = b0 + b1*efield(k) + b2*efield(k)**2
        mu01prime(k) = muprime*x01tmp
    ENDDO
END SUBROUTINE Get_mu01_Prime

SUBROUTINE Get_mu12_Prime(efield, w12, mu12prime)
    USE freq_data
    USE time_data
    USE map_data
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(IN) :: efield, w12
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(OUT) :: mu12prime

    ! Loop Variables
    INTEGER :: k

    ! Temporary Variables
    DOUBLE PRECISION :: x12tmp
    DOUBLE PRECISION :: muprime

    DO k=1, ntimes
        ! *** Position Matrix Elements ***
        x12tmp = d2 + d3*w12(k)

        ! *** Transition Dipole Matrix Elements ***
        muprime = b0 + b1*efield(k) + b2*efield(k)**2
        mu12prime(k) = muprime*x12tmp
    ENDDO
END SUBROUTINE Get_mu12_Prime


SUBROUTINE Get_01_Dipole(efield, w01, eOH, mu01)
    USE freq_data
    USE time_data
    USE map_data
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(ntimes),    INTENT(IN) :: efield, w01
    DOUBLE PRECISION, DIMENSION(ntimes,3),  INTENT(IN) :: eOH
    DOUBLE PRECISION, DIMENSION(ntimes, 3), INTENT(OUT) :: mu01

    DOUBLE PRECISION, DIMENSION(ntimes) :: mu01prime

    ! Loop Variables
    INTEGER :: k

    mu01 = 0.0
    CALL Get_mu01_Prime(efield, w01, mu01prime)
    DO k=1, ntimes
        mu01(k,:) = mu01prime(k)*eOH(k,:)
    ENDDO
END SUBROUTINE Get_01_Dipole

SUBROUTINE Get_12_Dipole(efield, w12, eOH, mu12)
    USE freq_data
    USE time_data
    USE map_data
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(ntimes),    INTENT(IN) :: efield,  w12
    DOUBLE PRECISION, DIMENSION(ntimes,3),  INTENT(IN) :: eOH
    DOUBLE PRECISION, DIMENSION(ntimes, 3), INTENT(OUT) ::  mu12

    DOUBLE PRECISION, DIMENSION(ntimes) :: mu12prime

    ! Loop Variables
    INTEGER :: k

    mu12 = 0.0
    CALL Get_mu12_Prime(efield, w12, mu12prime)
    DO k=1, ntimes
        mu12(k,:) = mu12prime(k)*eOH(k,:)
    ENDDO
END SUBROUTINE Get_12_Dipole

SUBROUTINE Get_Transition_Pol_Para_Perp(efield, w01, eOH, a_para, a_perp)
    USE freq_data
    USE time_data
    USE map_data
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(ntimes),    INTENT(IN) :: efield
    DOUBLE PRECISION, DIMENSION(ntimes),    INTENT(IN) :: w01
    DOUBLE PRECISION, DIMENSION(ntimes,3),  INTENT(IN) :: eOH
    DOUBLE PRECISION, DIMENSION(ntimes,3),  INTENT(OUT) :: a_para, a_perp
    ! Loop Variables
    INTEGER :: k
    ! Temporary Variables
    DOUBLE PRECISION :: x01tmp
    DOUBLE PRECISION :: alpha, apara, aperp

    DO k=1, ntimes
        x01tmp = d0 + d1*w01(k)

        alpha = x01tmp*(a0 + a1*efield(k))

        apara = 3d0/(1d0 + 2d0/c15) * alpha
        aperp = 3d0/(2d0 + c15) * alpha

        a_para(k,:) = aperp+(apara-aperp)*eOH(k,:)**2
        a_perp(k,1) = (aperp-apara)*eOH(k,1)*eOH(k,2)
        a_perp(k,2) = (aperp-apara)*eOH(k,2)*eOH(k,3)
        a_perp(k,3) = (aperp-apara)*eOH(k,3)*eOH(k,1)
    ENDDO
END SUBROUTINE Get_Transition_Pol_Para_Perp

SUBROUTINE Get_Transtion_Pol_Polarizations(efield, w01, eOH, a_ss, a_sp, a_pp)
    USE map_data
    USE freq_data
    USE time_data
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(ntimes),    INTENT(IN) :: efield, w01
    DOUBLE PRECISION, DIMENSION(ntimes,3),  INTENT(IN) :: eOH
    DOUBLE PRECISION, DIMENSION(ntimes),    INTENT(OUT):: a_ss, a_sp, a_pp

    INTEGER :: k

    DOUBLE PRECISION :: x01tmp
    DOUBLE PRECISION :: axx, ayy, azz, axy, ayz, azx
    DOUBLE PRECISION, DIMENSION(ntimes,3) :: a_para, a_perp

    a_ss = 0.0; a_sp = 0.0; a_pp = 0.0
    DO k=1, ntimes
        x01tmp = d0 + d1*w01(k)

        CALL Get_Transition_Pol_Para_Perp(efield, w01, eOH, a_para, a_perp)
        axx = a_para(k,1); ayy = a_para(k,2); azz = a_para(k,3)
        axy = a_perp(k,1); ayz = a_perp(k,2); azx = a_perp(k,3)

        ! Final polarizabilities
        a_ss(k) = x01tmp*(axx + 2.0d0*axy + ayy)
        a_sp(k) = x01tmp*(azx + ayz)
        a_pp(k) = x01tmp*(azz)
    ENDDO 
END SUBROUTINE Get_Transtion_Pol_Polarizations