SUBROUTINE Spec_Calc(tcf_tot)
! ******************************************************************
! This subroutine calculates the SFG spectrum from the total TCF.  
!                                                                  
!  INPUTS:                                                          
!   tcf_tot - The total TCF.    
!                                    
!  OUTPUTS:                                                         
!   None.                                                          
! *******************************************************************  
    USE time_data
    USE map_data
    USE hist_data
    USE constants

    IMPLICIT NONE
    INCLUDE 'fftw3.f'

    DOUBLE COMPLEX, DIMENSION(0:ncorr), INTENT(INOUT) :: tcf_tot

    INTEGER :: i, p, nt
    INTEGER*8 :: plan

    DOUBLE PRECISION :: w_spec, ti, pi
    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: input, ir_spectra

    pi = 4.0d0*datan(1.0d0)

    ! Write out the average frequency
    WRITE(6,*) '<w> = ', w01_avg*cmiperau

    ! I. Multiply the toal TCF by the vibrational relaxation factor.
    DO i = 0, ncorr
        ti = float(i)*dt
        tcf_tot(i) = tcf_tot(i)*dcmplx(dexp(-0.5d0*ti/T1),0d0)
    END DO

    ! II. Write out the total TCF
    OPEN(21, file='ir_tcf.dat')
    DO i = 0, ncorr
        ti = float(i)*dt*fsperau
        WRITE(21,*) ti, REAL(tcf_tot(i)), AIMAG(tcf_tot(i))
    END DO
    CLOSE(21)

    ! III. Initialize Fourier Transform
    ! ...A. Calculate the time domain grid for the freq resolution
    p = INT(dlog(2.0*pi/(dt*w_resol))/dlog(2.0d0))
    nt = 2**p
    WRITE(6,*) ' Time Grid power, 2^ = ',p
    WRITE(6,*) ' Time Grid Size      = ',nt

    ALLOCATE(input(nt)); ALLOCATE(ir_spectra(nt))

    CALL dfftw_plan_dft_1d(plan, nt, input, ir_spectra, FFTW_FORWARD, FFTW_ESTIMATE)

    ! ...B. Pad TCF with zeros to get the frequency resolution.
    input = dcmplx(0.0d0,0.0d0)
    DO i = 0, ncorr
        input(i) = tcf_tot(i)
    END DO

    CALL dfftw_execute(plan)

    ! IV. Write out the IR spectra
    ! IR Spectra is split in the middle because the average frequency
    ! has been subtracted from the phase.
    OPEN(22, file='ir_spectra.dat')
    DO i=nt/2 + 1, nt
        w_spec = (2d0*pi/dt)*dfloat(i-1)/dfloat(nt)
        WRITE(22,*) w_spec*cmiperau, DREAL(ir_spectra(i))
    ENDDO
    DO i=1, nt/2
        w_spec = (2d0*pi/dt)*dfloat(i-1)/dfloat(nt)
        WRITE(22,*) w_spec*cmiperau, DREAL(ir_spectra(i))
    ENDDO
    CLOSE(22)

    ! V. Deallocate memory
    CALL dfftw_destroy_plan(plan)
    DEALLOCATE(input); DEALLOCATE(ir_spectra)


END SUBROUTINE Spec_Calc