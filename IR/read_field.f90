SUBROUTINE read_field(ioh, w01, mu, eOH)
! This subroutine reads the field files and calculates various calculated quantities
!   including the ss, sp, and pp components of the polarizability. 
! 
! This program requires the HDF5 library for reading.
! 
! Copyright, Zeke Piskulich, 2024
! 
    USE time_data
    USE map_data
    USE freq_data
    USE HDF5
    USE ieee_arithmetic, ONLY: ieee_is_finite

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ioh
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(OUT) :: w01
    DOUBLE PRECISION, DIMENSION(ntimes,3), INTENT(OUT) :: eoh, mu

    INTEGER :: j, k
    REAL, DIMENSION(ntimes) :: etmp, z0_tmp
    REAL, DIMENSION(ntimes, 3) :: eoh_tmp
    DOUBLE PRECISION ::  muprime, x01tmp, mutmp, xtmp

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

    DO k=1,ntimes
        ! Calculate distance from center of slab
        eOH(k,:) = eoh_tmp(k,:)
        
        ! Covert to freq
        w01(k) = c0 + c1*etmp(k) + c2*etmp(k)**2
        w01_avg = w01_avg + w01(k)
        w01_sq_avg = w01_sq_avg + w01(k)**2

        ! Convert unit vector to transition dipole
        xtmp = d0 + d1*w01(k)
        mutmp = b0 + b1*etmp(k) + b2*etmp(k)**2

        mu(k,:) = eOH(k,:)*mutmp*xtmp

    ENDDO


END SUBROUTINE read_field
