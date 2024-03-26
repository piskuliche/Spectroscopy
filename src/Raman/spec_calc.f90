SUBROUTINE Spec_Calc(vv_tcf_tot, vh_tcf_tot)
! ******************************************************************
! This subroutine calculates the Raman spectrum from the total TCF.  
!                                                                  
!  INPUTS:                                                          
!   tcf_tot - The total TCF.    
!                                    
!  OUTPUTS:                                                         
!   None.                                                          
! *******************************************************************  

    USE time_data
    USE hist_data
    USE constants
    USE map_data
    USE freq_data

    IMPLICIT NONE
    INCLUDE '../Shared/fftw3.f'

    INTEGER :: i, p, nt
    INTEGER*8 :: plan
    DOUBLE PRECISION :: w_spec, ti, pi

    DOUBLE COMPLEX, DIMENSION(0:ncorr) :: vv_tcf_tot, vh_tcf_tot
    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: input, vv_spec, vh_spec

    pi = 4.0d0*datan(1.0d0)

    ! I. Multiply the toal TCF by the
    ! ...A. The Vibrational Relaxation Factor
    DO i = 0, ncorr
        ti = float(i)*dt
        vv_tcf_tot(i) = vv_tcf_tot(i)*DCMPLX(DEXP(-0.5d0*ti/T1),0d0)
        vh_tcf_tot(i) = vh_tcf_tot(i)*DCMPLX(DEXP(-0.5d0*ti/T1),0d0)
    END DO
    WRITE(6,*) '<w> = ',w01_avg*cmiperau

    ! ...B. The average frequency exp{-i*<w>*t}
    DO i = 0, ncorr
        ti = float(i)*dt
        vv_tcf_tot(i) = vv_tcf_tot(i)*DCMPLX(DCOS(ti*w01_avg), -DSIN(ti*w01_avg))
        vh_tcf_tot(i) = vh_tcf_tot(i)*DCMPLX(DCOS(ti*w01_avg), -DSIN(ti*w01_avg))
    END DO 

    ! II. Write out the total TCF
    OPEN(21, file='vv_tcf.dat')
    OPEN(22, file="vh_tcf.dat")
    DO i = 0, ncorr
        ti = float(i)*dt*fsperau
        WRITE(21,*) ti, REAL(vv_tcf_tot(i)), AIMAG(vv_tcf_tot(i))
        WRITE(21,*) ti, REAL(vh_tcf_tot(i)), AIMAG(vh_tcf_tot(i))
    END DO
    CLOSE(21)
    CLOSE(22)

    ! III. Initialize Fourier Transform
    ! ...A. Calculate the time domain grid for the freq resolution
    p = INT(dlog(2.0*pi/(dt*w_resol))/dlog(2.0d0))
    nt = 2**p
    WRITE(6,*) ' Time Grid power, 2^ = ',p
    WRITE(6,*) ' Time Grid Size      = ',nt

    ALLOCATE(input(nt)); ALLOCATE(vv_spec(nt)); ALLOCATE(vh_spec(nt))
    
    ! PART 1 - VV TCF
    CALL dfftw_plan_dft_1d(plan, nt, input, vv_spec, FFTW_FORWARD, FFTW_ESTIMATE)

    ! ... Pad TCF with zeros to get the frequency resolution.
    input = dcmplx(0.0d0,0.0d0)
    DO i = 1, min(nt, ncorr)
        input(i) = vv_tcf_tot(i-1)
    END DO

    CALL dfftw_execute(plan)

    ! ... Write out the VV Raman Spectrum
    ! Raman Spectra is split in the middle because the average frequency
    ! has been subtracted from the phase.
    OPEN(23, file='vv_spectrum.dat')
    DO i=nt/2 + 1, nt
        w_spec = (2d0*pi/dt)*dfloat(i-1-nt)/dfloat(nt) + w01_avg
        WRITE(23,*) w_spec*cmiperau, DREAL(vv_spec(i))
    ENDDO
    DO i=1, nt/2
        w_spec = (2d0*pi/dt)*dfloat(i-1)/dfloat(nt) + w01_avg
        WRITE(23,*) w_spec*cmiperau, DREAL(vv_spec(i))
    ENDDO
    CLOSE(23)

    ! ... Deallocate the VV Raman Spectrum
    CALL dfftw_destroy_plan(plan)

    ! PART 2 - VV TCF
    CALL dfftw_plan_dft_1d(plan, nt, input, vh_spec, FFTW_FORWARD, FFTW_ESTIMATE)

    ! ... Pad TCF with zeros to get the frequency resolution.
    input = dcmplx(0.0d0,0.0d0)
    DO i = 1, min(nt, ncorr)
        input(i) = vh_tcf_tot(i-1)
    END DO

    CALL dfftw_execute(plan)

    ! ... Write out the VH Raman Spectrum
    ! Raman Spectra is split in the middle because the average frequency
    ! has been subtracted from the phase.
    OPEN(24, file='vh_spectrum.dat')
    DO i=nt/2 + 1, nt
        w_spec = (2d0*pi/dt)*dfloat(i-1-nt)/dfloat(nt) + w01_avg
        WRITE(24,*) w_spec*cmiperau, DREAL(vh_spec(i))
    ENDDO
    DO i=1, nt/2
        w_spec = (2d0*pi/dt)*dfloat(i-1)/dfloat(nt) + w01_avg
        WRITE(24,*) w_spec*cmiperau, DREAL(vh_spec(i))
    ENDDO
    CLOSE(24)

    ! ... Deallocate the VV Raman Spectrum
    CALL dfftw_destroy_plan(plan)

    ! IV. Depolarized Spectrum
    OPEN(25, file='depol_spectrum.dat')
    DO i=nt/2 + 1, nt
        w_spec = (2d0*pi/dt)*dfloat(i-1-nt)/dfloat(nt) + w01_avg
        WRITE(25,*) w_spec*cmiperau, DREAL(10d0*vh_spec(i))
    ENDDO
    DO i=1, nt/2
        w_spec = (2d0*pi/dt)*dfloat(i-1)/dfloat(nt) + w01_avg
        WRITE(25,*) w_spec*cmiperau, DREAL(10d0*vh_spec(i))
    ENDDO
    CLOSE(25)

    ! V. Isotropic Spectrum
    OPEN(26, file='iso_spectrum.dat')
    DO i=nt/2 + 1, nt
        w_spec = (2d0*pi/dt)*dfloat(i-1-nt)/dfloat(nt) + w01_avg
        WRITE(26,*) w_spec*cmiperau, DREAL(vv_spec(i) - 4d0/3d0*vh_spec(i))
    ENDDO
    DO i=1, nt/2
        w_spec = (2d0*pi/dt)*dfloat(i-1)/dfloat(nt) + w01_avg
        WRITE(26,*) w_spec*cmiperau, DREAL(vv_spec(i) - 4d0/3d0*vh_spec(i))
    ENDDO
    CLOSE(26)


    DEALLOCATE(input); DEALLOCATE(vv_spec); DEALLOCATE(vh_spec)

END SUBROUTINE Spec_Calc
