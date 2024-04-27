!   Subroutine to Fourier Transform the total TCF to calculate the IR spectrum

SUBROUTINE spec_calc(iTw, tcf_rp_tot, tcf_np_tot)

   USE time_data
   USE map_data
   USE freq_data
   USE constants
   USE cli_data
   
   IMPLICIT NONE
   INCLUDE '../Shared/fftw3.f'
   INTEGER :: i, j, ineg, jneg, k, k_rp, p, nt1, nt3, iTw
   INTEGER*8 :: plan_rp, plan_np
   DOUBLE PRECISION :: w_spec1, w_spec3, ti, pi, fact1, fact3
   CHARACTER*5 :: cTw
   
   DOUBLE COMPLEX, DIMENSION(0:ncorr,0:ncorr) :: tcf_rp_tot, tcf_np_tot
   DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: input_rp, input_np, spec_rp, spec_np
   

   pi = 4.0d0*datan(1.0d0)
   
   ! Write out the average frequencies
   WRITE(6,*) ' <w01>   = ',w01_avg*cmiperau
   WRITE(6,*) ' <w12>   = ',w12_avg*cmiperau
   WRITE(6,*) ' <w01^2> = ',w01_sq_avg*cmiperau**2
   WRITE(6,*) ' <w12^2> = ',w12_sq_avg*cmiperau**2
   
   ! Initialize parameters for Fourier Transform
   
   ! Calculate the time domain grid needed to get the desired frequency resolution
   p = int(dlog(2.0d0*pi/(dt*w_resol))/dlog(2.0d0))
   nt1 = 2**p
   nt3 = 2**p
   WRITE(6,*) ' Time Grid power, 2^ = ',p
   WRITE(6,*) ' Time Grid Size t1   = ',nt1
   WRITE(6,*) ' Time Grid Size t3   = ',nt3
   
   ALLOCATE(input_rp(nt1,nt3)); ALLOCATE(spec_rp(nt1,nt3))
   ALLOCATE(input_np(nt1,nt3)); ALLOCATE(spec_np(nt1,nt3))
   
   CALL dfftw_plan_dft_2d(plan_rp, nt1, nt3, input_rp, spec_rp, FFTW_BACKWARD, FFTW_ESTIMATE)
   CALL dfftw_plan_dft_2d(plan_np, nt1, nt3, input_np, spec_np, FFTW_BACKWARD, FFTW_ESTIMATE)
   
   ! Write waiting time to a CHARACTER string
   WRITE(cTw,'(I0.5)') NINT(Tw(iTw)*fsperau)
   
   ! Pad the calculated TCF with zeroes to get the desired freq resolution before FTing
   input_rp = dcmplx(0.0d0,0.0d0); input_np = dcmplx(0.0d0,0.0d0)
   DO i = 0, MIN(nt1/2 - 1,ncorr - 1)
      DO j = 0, MIN(nt3/2 - 1,ncorr - 1)
         ! Fill the first quadrant (t1, t3 > 0) as given
         input_rp(i+1,j+1) = tcf_rp_tot(i,j)
         input_np(i+1,j+1) = tcf_np_tot(i,j)
         
         ! Symmetrize the other three quadrants
         ineg = nt1 - i
         jneg = nt3 - j
         input_rp(i+1, jneg) = tcf_rp_tot(i,j) 
         input_np(i+1, jneg) = tcf_np_tot(i,j) 
         input_rp(ineg, j+1) = tcf_rp_tot(i,j) 
         input_np(ineg, j+1) = tcf_np_tot(i,j) 
         input_rp(ineg, jneg) = tcf_rp_tot(i,j) 
         input_np(ineg, jneg) = tcf_np_tot(i,j) 
      ENDDO
   ENDDO
   
   CALL dfftw_execute(plan_rp)
   CALL dfftw_execute(plan_np)
      
   ! Write out the calculated spectrum
   OPEN(22,file=trim(tag_output_cli)//'2dir_imag_'//trim(cTw)//'.dat')
   OPEN(23,file=trim(tag_output_cli)//'2dir_real_'//trim(cTw)//'.dat')
      
!!! N/A    ! When the average frequency is subtracted from the phase the spectrum must be split in the middle
!!! N/A    ! First print out the frequencies below the average
   fact1 = (2d0*pi/dt)/DBLE(nt1); fact3 = (2d0*pi/dt)/DBLE(nt3)
   DO k = 1, nt1/2 
      w_spec1 = fact1*DBLE(k-1)
      k_rp = MODULO(nt1 + 1 - k, nt1 + 1) + FLOOR(DBLE(nt1 + 1 - k)/DBLE(nt1+1))
      DO j = 1, nt3/2 
         w_spec3 = fact3*DBLE(j-1)
         IF(w_spec1.ge.w1min.and.w_spec1.le.w1max.and.w_spec3.ge.w3min.and.w_spec3.le.w3max) THEN
            WRITE(22,'(2F20.10,2F22.10)') w_spec1*cmiperau, w_spec3*cmiperau, dimag(spec_rp(k_rp,j) + spec_np(k,j))
            WRITE(23,'(2F20.10,2F22.10)') w_spec1*cmiperau, w_spec3*cmiperau, dreal(spec_rp(k_rp,j) + spec_np(k,j))
         ENDIF
      ENDDO
   ENDDO
   
   close(22); close(23)
   
   CALL dfftw_destroy_plan(plan_rp)
   CALL dfftw_destroy_plan(plan_np)
   
   DEALLOCATE(input_rp); DEALLOCATE(input_np)
   
End Subroutine spec_calc

