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
    USE freq_data
    USE hist_data
    USE constants

    IMPLICIT NONE
    INCLUDE '../Shared/fftw3.f'

    DOUBLE COMPLEX, DIMENSION(0:ncorr), INTENT(INOUT) :: tcf_tot

    INTEGER :: i, p, nt
    INTEGER*8 :: plan

    DOUBLE PRECISION :: w_spec, ti, pi
    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: input, sfg_spec

    pi = 4.0d0*datan(1.0d0)

    ! Write out the average frequency
    write(6,*) ' <w> = ', w01_avg*cmiperau

    ! I. Multiply the total TCF by the vib relaxation factor
    DO i=0, ncorr
        ti = float(i)*dt
        tcf_tot(i) = tcf_tot(i)*DCMPLX(DEXP(-0.5d0*ti/T1),0d0)
    ENDDO 

    ! II. Write out the total tcf
    OPEN(21, file='sfg_tcf.dat')
    DO i=0, ncorr
        WRITE(21,*) float(i)*dt*fsperau, REAL(tcf_tot(i)), AIMAG(tcf_tot(i))
    ENDDO
    CLOSE(21)

    ! III. Initialize Fourier Transform
    ! ...A: Calculate time domain grid for freq resolution
    p = int(dlog(2.0*pi/(dt*w_resol))/dlog(2.0d0))
    nt = 2**p
    WRITE(6,*) ' Time Grid power, 2^ = ',p
    WRITE(6,*) ' Time Grid Size      = ',nt

    ALLOCATE(input(nt)); ALLOCATE(sfg_spec(nt))

    CALL dfftw_plan_dft_1d(plan, nt, input, sfg_spec, FFTW_FORWARD, FFTW_ESTIMATE)

    ! ...B: Pad TCF with zeroes to get frequency resoltuion
    input = DCMPLX(0.0d0,0.0d0)
    DO i=1, min(nt, ncorr)
        input(i) = tcf_tot(i-1)
    ENDDO
    
    CALL dfftw_execute(plan)
    
    OPEN(22,file='sfg_spectrum.dat')
    ! IV. Write out the spectrum
    DO i=1, nt
        w_spec = (2d0*pi/dt)*dfloat(i-1)/dfloat(nt)
        IF (w_spec .ge. w1min .and. w_spec .le. w1max) THEN
            WRITE(22,*) w_spec*cmiperau, DREAL(sfg_spec(i)), -DIMAG(sfg_spec(i))
        ENDIF
    ENDDO

    CLOSE(22)

    ! V. Cleanup the calculations

    CALL dfftw_destroy_plan(plan)
    DEALLOCATE(input); DEALLOCATE(sfg_spec)

END SUBROUTINE Spec_Calc