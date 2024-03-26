!   Subroutine to Fourier Transform the total TCF to calculate the IR spectrum

   Subroutine spec_calc(iTw, tcf_rp_tot, tcf_np_tot)
  
      use time_data
      use map_data
      use freq_data
      use constants
      
      implicit none
      include '../Shared/fftw3.f'
      integer :: i, j, ineg, jneg, k, k_rp, p, nt1, nt3, iTw
      integer*8 :: plan_rp, plan_np
      double precision :: w_spec1, w_spec3, ti, pi, fact1, fact3
      character*5 :: cTw
      
      double complex, dimension(0:ncorr,0:ncorr) :: tcf_rp_tot, tcf_np_tot
      double complex, allocatable, dimension(:,:) :: input_rp, input_np, spec_rp, spec_np
      

      pi = 4.0d0*datan(1.0d0)
      
      ! Write out the average frequencies
      write(6,*) ' <w01>   = ',w01_avg*cmiperau
      write(6,*) ' <w12>   = ',w12_avg*cmiperau
      write(6,*) ' <w01^2> = ',w01_sq_avg*cmiperau**2
      write(6,*) ' <w12^2> = ',w12_sq_avg*cmiperau**2
      
      ! Initialize parameters for Fourier Transform
      
      ! Calculate the time domain grid needed to get the desired frequency resolution
      p = int(dlog(2.0d0*pi/(dt*w_resol))/dlog(2.0d0))
      nt1 = 2**p
      nt3 = 2**p
      write(6,*) ' Time Grid power, 2^ = ',p
      write(6,*) ' Time Grid Size t1   = ',nt1
      write(6,*) ' Time Grid Size t3   = ',nt3
      
      allocate(input_rp(nt1,nt3)); allocate(spec_rp(nt1,nt3))
      allocate(input_np(nt1,nt3)); allocate(spec_np(nt1,nt3))
      
      call dfftw_plan_dft_2d(plan_rp, nt1, nt3, input_rp, spec_rp, FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_plan_dft_2d(plan_np, nt1, nt3, input_np, spec_np, FFTW_BACKWARD, FFTW_ESTIMATE)
      
      ! Write waiting time to a character string
      write(cTw,'(I0.5)') nint(Tw(iTw)*fsperau)
      
      ! Pad the calculated TCF with zeroes to get the desired freq resolution before FTing
      input_rp = dcmplx(0.0d0,0.0d0); input_np = dcmplx(0.0d0,0.0d0)
      do i = 0, min(nt1/2 - 1,ncorr - 1)
         do j = 0, min(nt3/2 - 1,ncorr - 1)
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
         enddo
      enddo
      
      call dfftw_execute(plan_rp)
      call dfftw_execute(plan_np)
         
      ! Write out the calculated spectrum
      open(22,file='2dir_imag_'//trim(cTw)//'.dat')
      open(23,file='2dir_real_'//trim(cTw)//'.dat')
      !         open(24,file='2dtcf_imag_'//trim(cTw)//'.dat')
      !         open(25,file='2dtcf_real_'//trim(cTw)//'.dat')

         ! Write out the TCF

!         do k = 0, ncorr
!            do j = 0, ncorr
!               write(24,'(3F16.8)') dble(k)*dt*fsperau, dble(j)*dt*fsperau, dimag(tcf_tot(k,j,iTw))
!               write(25,'(3F16.8)') dble(k)*dt*fsperau, dble(j)*dt*fsperau, dreal(tcf_tot(k,j,iTw))
!            enddo
!         enddo
         
!!! N/A    ! When the average frequency is subtracted from the phase the spectrum must be split in the middle
!!! N/A    ! First print out the frequencies below the average
      fact1 = (2d0*pi/dt)/dble(nt1); fact3 = (2d0*pi/dt)/dble(nt3)
      do k = 1, nt1/2 
         w_spec1 = fact1*dble(k-1)
         k_rp = modulo(nt1 + 1 - k, nt1 + 1) + floor(dble(nt1 + 1 - k)/dble(nt1+1))
         do j = 1, nt3/2 
            w_spec3 = fact3*dble(j-1)
            if(w_spec1.ge.w1min.and.w_spec1.le.w1max.and.w_spec3.ge.w3min.and.w_spec3.le.w3max) then
               write(22,'(2F16.6,2F18.6)') w_spec1*cmiperau, w_spec3*cmiperau, dimag(spec_rp(k_rp,j) + spec_np(k,j))
               write(23,'(2F16.6,2F18.6)') w_spec1*cmiperau, w_spec3*cmiperau, dreal(spec_rp(k_rp,j) + spec_np(k,j))
            endif
            
         enddo
      enddo
      
      close(22); close(23)
      
      call dfftw_destroy_plan(plan_rp)
      call dfftw_destroy_plan(plan_np)
      deallocate(input_rp,input_np)
      
    End Subroutine spec_calc



!   Subroutine to Fourier Transform the total TCF to calculate the IR spectrum
!     and the energy fluctuated versions

   Subroutine spec_calc_fluc(iTw, tcf_rp_tot, tcf_np_tot, tcfH_rp_tot, tcfH_np_tot)
  
      use time_data
      use map_data
      use freq_data
      use constants
      
      implicit none
      include 'fftw3.f'
      integer :: i, j, ineg, jneg, k, k_rp, p, nt1, nt3, iTw, iH
      integer*8 :: plan_rp, plan_np
      double precision :: w_spec1, w_spec3, ti, pi, fact1, fact3
      character*5 :: cTw
      character*2, dimension(5) :: Hlab
      
      double complex, dimension(0:ncorr,0:ncorr) :: tcf_rp_tot, tcf_np_tot
      double complex, dimension(0:ncorr,0:ncorr,5) :: tcfH_rp_tot, tcfH_np_tot
      double complex, allocatable, dimension(:,:) :: input_rp, input_np, spec_rp, spec_np
      

      pi = 4.0d0*datan(1.0d0)
      Hlab(1) = 'H'; Hlab(2) = 'KE'; Hlab(3) = 'LJ'; Hlab(4) = 'Co'; Hlab(5) = 'V'
      
      ! Write out the average frequencies
      write(6,*) ' <w01>   = ',w01_avg*cmiperau
      write(6,*) ' <w12>   = ',w12_avg*cmiperau
      write(6,*) ' <w01^2> = ',w01_sq_avg*cmiperau**2
      write(6,*) ' <w12^2> = ',w12_sq_avg*cmiperau**2
      
      ! Initialize parameters for Fourier Transform
      
      ! Calculate the time domain grid needed to get the desired frequency resolution
      p = int(dlog(2.0d0*pi/(dt*w_resol))/dlog(2.0d0))
      nt1 = 2**p
      nt3 = 2**p
      write(6,*) ' Time Grid power, 2^ = ',p
      write(6,*) ' Time Grid Size t1   = ',nt1
      write(6,*) ' Time Grid Size t3   = ',nt3
      
      allocate(input_rp(nt1,nt3)); allocate(spec_rp(nt1,nt3))
      allocate(input_np(nt1,nt3)); allocate(spec_np(nt1,nt3))

      
      call dfftw_plan_dft_2d(plan_rp, nt1, nt3, input_rp, spec_rp, FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_plan_dft_2d(plan_np, nt1, nt3, input_np, spec_np, FFTW_BACKWARD, FFTW_ESTIMATE)
      
      ! Write the waiting time to a character string  
      write(cTw,'(I0.5)') nint(Tw(iTw)*fsperau)
         
      ! Pad the calculated TCF with zeroes to get the desired freq resolution before FTing
      input_rp = dcmplx(0.0d0,0.0d0); input_np = dcmplx(0.0d0,0.0d0)
      do i = 0, min(nt1/2 - 1,ncorr - 1)
         do j = 0, min(nt3/2 - 1,ncorr - 1)
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
         enddo
      enddo
      
      call dfftw_execute(plan_rp)
      call dfftw_execute(plan_np)
         
      ! Write out the calculated spectrum
      open(23,file='2dir_real_'//trim(cTw)//'.dat')
      open(24,file='2dir_imag_'//trim(cTw)//'.dat')

      !         open(24,file='2dtcf_imag_'//trim(cTw)//'.dat')
      !         open(25,file='2dtcf_real_'//trim(cTw)//'.dat')

         ! Write out the TCF

!         do k = 0, ncorr
!            do j = 0, ncorr
!               write(24,'(3F16.8)') dble(k)*dt*fsperau, dble(j)*dt*fsperau, dimag(tcf_tot(k,j,iTw))
!               write(25,'(3F16.8)') dble(k)*dt*fsperau, dble(j)*dt*fsperau, dreal(tcf_tot(k,j,iTw))
!            enddo
!         enddo
         
      fact1 = (2d0*pi/dt)/dble(nt1); fact3 = (2d0*pi/dt)/dble(nt3)
      do k = 1, nt1/2 
         w_spec1 = fact1*dble(k-1)
         k_rp = modulo(nt1 + 1 - k, nt1 + 1) + floor(dble(nt1 + 1 - k)/dble(nt1+1))
         do j = 1, nt3/2 
            w_spec3 = fact3*dble(j-1)
            if(w_spec1.ge.w1min.and.w_spec1.le.w1max.and.w_spec3.ge.w3min.and.w_spec3.le.w3max) then
               write(23,'(4F16.8)') w_spec1*cmiperau, w_spec3*cmiperau, dreal(spec_rp(k_rp,j) + spec_np(k,j))
               write(24,'(4F16.8)') w_spec1*cmiperau, w_spec3*cmiperau, dimag(spec_rp(k_rp,j) + spec_np(k,j))
            endif
            
         enddo
      enddo
      
      close(23); close(24) 

      do iH = 1, 4
         ! Pad the calculated TCF with zeroes to get the desired freq resolution before FTing
         input_rp = dcmplx(0.0d0,0.0d0); input_np = dcmplx(0.0d0,0.0d0)
         do i = 0, min(nt1/2 - 1,ncorr - 1)
            do j = 0, min(nt3/2 - 1,ncorr - 1)
               ! Fill the first quadrant (t1, t3 > 0) as given
               input_rp(i+1,j+1) = tcfH_rp_tot(i,j,iH)
               input_np(i+1,j+1) = tcfH_np_tot(i,j,iH)
                  
               ! Symmetrize the other three quadrants
               ineg = nt1 - i
               jneg = nt3 - j
               input_rp(i+1, jneg) = tcfH_rp_tot(i,j,iH) 
               input_np(i+1, jneg) = tcfH_np_tot(i,j,iH) 
               input_rp(ineg, j+1) = tcfH_rp_tot(i,j,iH) 
               input_np(ineg, j+1) = tcfH_np_tot(i,j,iH) 
               input_rp(ineg, jneg) = tcfH_rp_tot(i,j,iH) 
               input_np(ineg, jneg) = tcfH_np_tot(i,j,iH) 
            enddo
         enddo
            
         call dfftw_execute(plan_rp)
         call dfftw_execute(plan_np)
         
         ! Write out the calculated spectrum
         open(23,file='2dir'//trim(Hlab(iH))//'_real_'//trim(cTw)//'.dat')
         open(24,file='2dir'//trim(Hlab(iH))//'_imag_'//trim(cTw)//'.dat')

         fact1 = (2d0*pi/dt)/dble(nt1); fact3 = (2d0*pi/dt)/dble(nt3)
         do k = 1, nt1/2 
            w_spec1 = fact1*dble(k-1)
            k_rp = modulo(nt1 + 1 - k, nt1 + 1) + floor(dble(nt1 + 1 - k)/dble(nt1+1))
            do j = 1, nt3/2 
               w_spec3 = fact3*dble(j-1)
               if(w_spec1.ge.w1min.and.w_spec1.le.w1max.and.w_spec3.ge.w3min.and.w_spec3.le.w3max) then
                  write(23,'(4F16.8)') w_spec1*cmiperau, w_spec3*cmiperau, dreal(spec_rp(k_rp,j) + spec_np(k,j))
                  write(24,'(4F16.8)') w_spec1*cmiperau, w_spec3*cmiperau, dimag(spec_rp(k_rp,j) + spec_np(k,j))
               endif
            enddo
         enddo
         
      enddo  ! Loop over energy fluctuation components
         
      close(23); close(24) 
      
      call dfftw_destroy_plan(plan_rp)
      call dfftw_destroy_plan(plan_np)
      deallocate(input_rp,input_np)
      
    End Subroutine spec_calc_fluc



