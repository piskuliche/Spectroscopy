SUBROUTINE read_field(ioh, w01, w12, mu01, mu12, eOH)
    USE time_data
    USE map_data
    USE freq_data
    USE HDF5
    USE ieee_arithmetic, ONLY: ieee_is_finite

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ioh
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(OUT) :: w01, w12
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(OUT) :: mu01, mu12
    DOUBLE PRECISION, DIMENSION(ntimes,3), INTENT(OUT) :: eoh

    INTEGER :: j, k
    REAL, DIMENSION(ntimes) :: etmp
    REAL, DIMENSION(ntimes, 3) :: eoh_tmp
    DOUBLE PRECISION ::  muprime, x01tmp, x12tmp

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

    ! Close the file
    CALL h5fclose_f(file_id, ERROR_FLAG)
    ! Close the library
    CALL h5close_f(ERROR_FLAG)

! II. Calculate the field parameters ******************************************
    DO k=1, ntimes

      
        w01(k) = c0 + c1*etmp(k) + c2*etmp(k)**2
        IF (.NOT. ieee_is_finite(w01(k))) THEN
         WRITE(6,*) ioh, k
         STOP 'w01 is not finite'
        ENDIF
        w01_avg = w01_avg + w01(k)
        w01_sq_avg = w01_sq_avg + w01(k)**2

        w12(k) = c3 + c4*etmp(k) + c5*etmp(k)**2     
        w12_avg = w12_avg + w12(k)
        w12_sq_avg = w12_sq_avg + w12(k)**2

        ! convert unit vector to transition dipole moment
        muprime = b0 + b1*etmp(k) + b2*etmp(k)**2

        x01tmp = d0 + d1*w01(k)     
        mu01(k) = muprime*x01tmp

        x12tmp = d2 + d3*w12(k)     
        mu12(k) = muprime*x12tmp

        eoh(k,:) = eoh_tmp(k,:)
    END DO


END SUBROUTINE read_field



