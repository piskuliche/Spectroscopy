SUBROUTINE read_field(ioh, w01, mu, eOH, a_ss, a_sp, a_pp, z0)
    USE time_data
    USE map_data
    USE freq_data
    USE HDF5
    USE ieee_arithmetic, ONLY: ieee_is_finite

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ioh
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(OUT) :: w01
    DOUBLE PRECISION, DIMENSION(ntimes,3), INTENT(OUT) :: eoh, mu
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(OUT) :: a_ss, a_sp, a_pp, z0

    INTEGER :: j, k
    REAL, DIMENSION(ntimes) :: etmp, z0_tmp
    REAL, DIMENSION(ntimes, 3) :: eoh_tmp
    DOUBLE PRECISION ::  muprime, x01tmp, mutmp, xtmp
    DOUBLE PRECISION :: alpha, apara, aperp
    DOUBLE PRECISION :: axx, ayy, azz, axy, ayz, azx

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
        ! Calculate distance from center of slab
        z0(k) = z0_tmp(k) - z_c
        
        ! Covert to freq
        w01(k) = c0 + c1*etmp(k) + c2*etmp(k)**2
        w01avg = w01avg + w01(k)
        w01sqavg = w01sqavg + w01(k)**2
        ! DO I NEED w2avg?

        ! Convert unit vector to transition dipole
        xtmp = d0+d1*w01(k)
        mutmp = b0 + b1*etmp(k) + b2*etmp(k)**2

        mu(k,:) = eOH(k,:)*mutmp*xtmp

        alpha = xtmp*(a0 + a1*etmp(k))

        ! This next section should be enclosed in a function

        apara = 3d0/(1d0 + 2d0/c15) * alpha
        aperp = 3d0/(2d0 + c15) * alpha

        axx = aperp + (apara - aperp)*eOH(k,1)*eOH(k,1)
        ayy = aperp + (apara - aperp)*eOH(k,2)*eOH(k,2)
        azz = aperp + (apara - aperp)*eOH(k,3)*eOH(k,3)
        axy = (aperp - apara)*eOH(k,1)*eOH(k,2)
        ayz = (aperp - apara)*eOH(k,2)*eOH(k,3)
        azx = (aperp - apara)*eOH(k,3)*eOH(k,1)

        ! Final polarizabilities
        a_ss(k) = xtmp*(axx + 2.0d0*axy + ayy)
        a_sp(k) = xtmp*(azx + ayz)
        a_pp(k) = xtmp*(azz)

    ENDDO


END SUBROUTINE read_field
