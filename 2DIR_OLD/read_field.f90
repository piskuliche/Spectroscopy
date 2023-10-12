SUBROUTINE read_field(ioh, w01, w12, mu01, mu12, eOH)
    USE time_data
    USE map_data
    USE freq_data
    USE HDF5

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ioh
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(OUT) :: w01, w12
    DOUBLE PRECISION, DIMENSION(ntimes), INTENT(OUT) :: mu01, mu12
    DOUBLE PRECISION, DIMENSION(ntimes,3), INTENT(OUT) :: eoh

    INTEGER :: j, k
    REAL, DIMENSION(2000) :: etmp
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
      WRITE(*,*) w01(k)
        w01(k) = c0 + c1*etmp(k) + c2*etmp(k)**2
        w01avg = w01avg + w01(k)
        w01sqavg = w01sqavg + w01(k)**2

        w12(k) = c3 + c4*etmp(k) + c5*etmp(k)**2     
        w12avg = w12avg + w12(k)
        w12sqavg = w12sqavg + w12(k)**2

        ! convert unit vector to transition dipole moment
        muprime = b0 + b1*etmp(k) + b2*etmp(k)**2

        x01tmp = d0 + d1*w01(k)     
        mu01(k) = muprime*x01tmp

        x12tmp = d2 + d3*w12(k)     
        mu12(k) = muprime*x12tmp

        eoh(k,:) = eoh_tmp(k,:)
    END DO
    WRITE(6,*) "**************"
    WRITE(6,*) 'w01avg tmp= ', w01avg


END SUBROUTINE read_field


!   Subroutine to read in the field parameters from field.xxx
!     and calculate the derivatives of the averages

  Subroutine read_field_fluc(ioh, w01, w12, mu01, mu12, eOH, dH)

  use time_data
  use map_data
  use freq_data
  use fluc_data

  implicit none
  integer :: ioh, j, k, n
  character*4 :: ext
  double precision :: etmp, muprime, x01tmp, x12tmp
  double precision, dimension(ntimes) :: w01, w12
  double precision, dimension(ntimes) :: mu01, mu12
  double precision, dimension(ntimes,3) :: eOH
  double precision, dimension(ntimes,4) :: dH



  ! Open the field field file for OH group iOH
  if(ndigits.eq.4) then
     write(ext,'(I0.4)') ioh 
  elseif(ndigits.eq.3) then
     write(ext,'(I0.3)') ioh
  else
     write(6,*) ' Invalid ndigits, only 3 or 4 are allowed'
     stop
  endif

  ! Open the field file
  open(ioh+100,file='field_files/field.'//ext//'',status='old')

  do k = 1 , ntimes

     ! read in data from file (here mu = e_OH)
     read(ioh+100,*) etmp, (eOH(k,j), j=1,3)
     
     ! convert etmp to frequency
     w01(k) = c0 + c1*etmp + c2*etmp**2     
     w01avg = w01avg + w01(k)
     w01sqavg = w01sqavg + w01(k)**2

     w12(k) = c3 + c4*etmp + c5*etmp**2     
     w12avg = w12avg + w12(k)
     w12sqavg = w12sqavg + w12(k)**2


     ! calculate the fluctuation weighted averages
     do n = 1, 5
        dw01avg(n) = dw01avg(n) + dH(k,n)*w01(k)
        dw01sqavg(n) = dw01sqavg(n) + dH(k,n)*w01(k)**2
        d2w01avg(n) = d2w01avg(n) + (dH(k,n)**2 - dH2avg(n))*w01(k)
        d2w01sqavg(n) = d2w01sqavg(n) + (dH(k,n)**2 - dH2avg(n))*w01(k)**2

        dw12avg(n) = dw12avg(n) + dH(k,n)*w12(k)
        dw12sqavg(n) = dw12sqavg(n) + dH(k,n)*w12(k)**2
        d2w12avg(n) = d2w12avg(n) + (dH(k,n)**2 - dH2avg(n))*w12(k)
        d2w12sqavg(n) = d2w12sqavg(n) + (dH(k,n)**2 - dH2avg(n))*w12(k)**2
     enddo

     ! convert unit vector to transition dipole moment
     muprime = b0 + b1*etmp + b2*etmp**2

     x01tmp = d0 + d1*w01(k)     
     mu01(k) = muprime*x01tmp

     x12tmp = d2 + d3*w12(k)
     mu12(k) = muprime*x12tmp

  enddo
  rewind(ioh+100)
  close(ioh+100)

  End subroutine read_field_fluc
