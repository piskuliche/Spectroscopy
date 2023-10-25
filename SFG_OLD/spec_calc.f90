!   Subroutine to Fourier Transform the total TCF to calculate the SFG spectrum

  Subroutine spec_calc(tcf_tot)

  use time_data
  use map_data
  use freq_data
  use hist_data
  use constants

  implicit none
  include 'fftw3.f'
  integer :: i, p, nt
  integer*8 :: plan
  double precision :: w_spec, ti, pi
  double complex, dimension(0:ncorr) :: tcf_tot
  double complex, allocatable, dimension(:) :: input, sfg_spec

  pi = 4.0d0*datan(1.0d0)

  ! Write out the average frequency
  write(6,*) ' <w>  = ',wavg*cmiperau

  ! Multiply the total TCF by the vibrational relaxation factor
  do i = 0, ncorr
     ti = float(i)*dt
     tcf_tot(i) = tcf_tot(i)*dcmplx(dexp(-0.5d0*ti/T1),0d0)
  enddo

!  ! Multiply the total TCF by exp{-i*<w>*t}
!  do i = 0, ncorr
!     ti = float(i)*dt
!     tcf_tot(i) = tcf_tot(i)*dcmplx(dcos(ti*wavg),-dsin(ti*wavg))
!  enddo

  ! Write out the total TCF
  open(21,file='sfg_tcf.dat')
  do i = 0, ncorr
     write(21,*) float(i)*dt*fsperau, real(tcf_tot(i)), aimag(tcf_tot(i))
  enddo
  close(21)


  ! Initialize parameters for Fourier Transform

  ! Calculate the time domain grid needed to get the desired frequency resolution
  p = int(dlog(2.0*pi/(dt*w_resol))/dlog(2.0))
  nt = 2**p
  write(6,*) ' Time Grid power, 2^ = ',p
  write(6,*) ' Time Grid Size      = ',nt

  allocate(input(nt)); allocate(sfg_spec(nt))

  call dfftw_plan_dft_1d(plan, nt, input, sfg_spec, FFTW_FORWARD, FFTW_ESTIMATE)
  
  ! Pad the calculated TCF with zeroes to get the desired freq resolution before FTing
  input = dcmplx(0.0d0,0.0d0)
  do i = 1, min(nt,ncorr)
     input(i) = tcf_tot(i-1)
  enddo
  
  call dfftw_execute(plan)

  ! Write out the calculated spectrum
  open(22,file='sfg_spectrum.dat')

  ! When the average frequency is subtracted from the phase the spectrum must be split in the middle
  ! First print out the frequencies below the average
!  do i = nt/2 + 1, nt
!     w_spec = (2d0*pi/dt)*dfloat(i-1-nt)/dfloat(nt) !+ wavg
!     if(w_spec.ge.wmin.and.w_spec.le.wmax) then
!        write(22,*) w_spec*cmiperau, dreal(sfg_spec(i)), -dimag(sfg_spec(i)) ! dreal = imag part and dimag = - real part
!     endif
!  enddo
  ! Then print out the frequencies above the average
!  do i = 1, nt/2 
  do i = 1, nt
     w_spec = (2d0*pi/dt)*dfloat(i-1)/dfloat(nt) !+ wavg
     if(w_spec.ge.wmin.and.w_spec.le.wmax) then
        write(22,*) w_spec*cmiperau, dreal(sfg_spec(i)), -dimag(sfg_spec(i)) ! dreal = imag part and dimag = - real part
     endif
  enddo

  call dfftw_destroy_plan(plan)
  deallocate(input)

   End Subroutine spec_calc



!   Subroutine to Fourier Transform the total TCF to calculate the SFG spectrum

  Subroutine spec_calc_fluc(tcf_tot, tcfH_tot)

  use time_data
  use map_data
  use freq_data
  use hist_data
  use fluc_data
  use constants

  implicit none
  include 'fftw3.f'
  integer :: i, k, p, nt
  integer*8 :: plan
  double precision :: w_spec, ti, pi
  double complex, dimension(0:ncorr) :: tcf_tot
  double complex, dimension(0:ncorr,8) :: tcfH_tot
  double complex, allocatable, dimension(:) :: input, sfg_spec, sfg_deriv

  pi = 4.0d0*datan(1.0d0)

  ! Write out the average frequency, freq squared and the fluct.-weighted values
  write(6,*) ' <w>  = ',wavg*cmiperau,' cm^-1'
  write(6,*) ' <w2> = ',w2avg*cmiperau**2,' (cm^-1)^2'
  write(6,*) ' <dH*w>     = ',dwavg(1)*cmiperau,' cm^-1*kcal/mol '
  write(6,*) ' <dKE*w>    = ',dwavg(2)*cmiperau,' cm^-1*kcal/mol '
  write(6,*) ' <dLJ*w>    = ',dwavg(3)*cmiperau,' cm^-1*kcal/mol '
  write(6,*) ' <dCoul*w>  = ',dwavg(4)*cmiperau,' cm^-1*kcal/mol '
  write(6,*) ' <dV*w>     = ',dwavg(5)*cmiperau,' cm^-1*kcal/mol '
!  write(6,*) ' <dH_top*w> = ',dwavg(6)*cmiperau,' cm^-1*kcal/mol '
!  write(6,*) ' <dH_bot*w> = ',dwavg(7)*cmiperau,' cm^-1*kcal/mol '
!  write(6,*) ' <dH_t-b*w> = ',dwavg(8)*cmiperau,' cm^-1*kcal/mol '
  write(6,*) ' <dH*w2>    = ',dw2avg(1)*cmiperau**2,' (cm^-1)^2*kcal/mol '
  write(6,*) ' <dKE*w2>   = ',dw2avg(2)*cmiperau**2,' (cm^-1)^2*kcal/mol '
  write(6,*) ' <dLJ*w2>   = ',dw2avg(3)*cmiperau**2,' (cm^-1)^2*kcal/mol '
  write(6,*) ' <dCoul*w2> = ',dw2avg(4)*cmiperau**2,' (cm^-1)^2*kcal/mol '
  write(6,*) ' <dV*w2>    = ',dw2avg(5)*cmiperau**2,' (cm^-1)^2*kcal/mol '
  write(6,*) ' <d2H*w>     = ',d2wavg(1)*cmiperau,' cm^-1*(kcal/mol)^2 '
  write(6,*) ' <d2KE*w>    = ',d2wavg(2)*cmiperau,' cm^-1*(kcal/mol)^2 '
  write(6,*) ' <d2LJ*w>    = ',d2wavg(3)*cmiperau,' cm^-1*(kcal/mol)^2 '
  write(6,*) ' <d2Coul*w>  = ',d2wavg(4)*cmiperau,' cm^-1*(kcal/mol)^2 '
  write(6,*) ' <d2V*w>     = ',d2wavg(5)*cmiperau,' cm^-1*(kcal/mol)^2 '
  write(6,*) ' <d2H*w2>    = ',d2w2avg(1)*cmiperau**2,' (cm^-1)^2*(kcal/mol)^2 '
  write(6,*) ' <d2KE*w2>   = ',d2w2avg(2)*cmiperau**2,' (cm^-1)^2*(kcal/mol)^2 '
  write(6,*) ' <d2LJ*w2>   = ',d2w2avg(3)*cmiperau**2,' (cm^-1)^2*(kcal/mol)^2 '
  write(6,*) ' <d2Coul*w2> = ',d2w2avg(4)*cmiperau**2,' (cm^-1)^2*(kcal/mol)^2 '
  write(6,*) ' <d2V*w2>    = ',d2w2avg(5)*cmiperau**2,' (cm^-1)^2*(kcal/mol)^2 '

  ! Multiply the total TCF by the vibrational relaxation factor
  do i = 0, ncorr
     ti = float(i)*dt
     tcf_tot(i) = tcf_tot(i)*dcmplx(dexp(-0.5d0*ti/T1),0d0)
     tcfH_tot(i,:) = tcfH_tot(i,:)*dcmplx(dexp(-0.5d0*ti/T1),0d0)
  enddo

!  ! Multiply the total TCF by exp{-i*<w>*t}
!  do i = 0, ncorr
!     ti = float(i)*dt
!     tcf_tot(i) = tcf_tot(i)*dcmplx(dcos(ti*wavg),-dsin(ti*wavg))
!     tcfH_tot(i,:) = tcfH_tot(i,:)*dcmplx(dcos(ti*wavg),-dsin(ti*wavg))
!  enddo

  ! Write out the total TCF
  open(21,file='sfg_tcf.dat')
  open(25,file='dH_sfg_tcf.dat')
  open(26,file='dKE_sfg_tcf.dat')
  open(27,file='dLJ_sfg_tcf.dat')
  open(28,file='dCoul_sfg_tcf.dat')
!  open(29,file='dHtop_sfg_tcf.dat')
!  open(30,file='dHbot_sfg_tcf.dat')
!  open(31,file='dHtb_sfg_tcf.dat')
  do i = 0, ncorr
     write(21,*) float(i)*dt*fsperau, real(tcf_tot(i)), aimag(tcf_tot(i))
     write(25,*) float(i)*dt*fsperau, real(tcfH_tot(i,1)), aimag(tcfH_tot(i,1))
     write(26,*) float(i)*dt*fsperau, real(tcfH_tot(i,2)), aimag(tcfH_tot(i,2))
     write(27,*) float(i)*dt*fsperau, real(tcfH_tot(i,3)), aimag(tcfH_tot(i,3))
     write(28,*) float(i)*dt*fsperau, real(tcfH_tot(i,4)), aimag(tcfH_tot(i,4))
!     write(29,*) float(i)*dt*fsperau, real(tcfH_tot(i,6)), aimag(tcfH_tot(i,6))
!     write(30,*) float(i)*dt*fsperau, real(tcfH_tot(i,7)), aimag(tcfH_tot(i,7))
!     write(31,*) float(i)*dt*fsperau, real(tcfH_tot(i,8)), aimag(tcfH_tot(i,8))
  enddo
  close(21); close(25); close(26); close(27); close(28)!; close(29); close(30); close(31)

  ! Initialize parameters for Fourier Transform

  ! Calculate the time domain grid needed to get the desired frequency resolution
  p = int(dlog(2.0*pi/(dt*w_resol))/dlog(2.0))
  nt = 2**p
  write(6,*) ' Time Grid power, 2^ = ',p
  write(6,*) ' Time Grid Size      = ',nt

  allocate(input(nt)); allocate(sfg_spec(nt)); ; allocate(sfg_deriv(nt))

!!!!!!!! Calculate the FT to get the Spectrum 

  call dfftw_plan_dft_1d(plan, nt, input, sfg_spec, FFTW_FORWARD, FFTW_ESTIMATE)
  
  ! Pad the calculated TCF with zeroes to get the desired freq resolution before FTing
  input = dcmplx(0.0d0,0.0d0)
  do i = 1, min(nt,ncorr)
     input(i) = tcf_tot(i-1)
  enddo
  
  call dfftw_execute(plan)

  ! Write out the calculated spectrum
  open(22,file='sfg_spectrum.dat')

  ! When the average frequency is subtracted from the phase the spectrum must be split in the middle
  ! First print out the frequencies below the average
!  do i = nt/2 + 1, nt
!     w_spec = (2d0*pi/dt)*dfloat(i-1-nt)/dfloat(nt) + wavg
!!     w_spec = (2d0*pi/dt)*dfloat(nt + 1 - i)/dfloat(nt) !+ wavg
!     if(w_spec.ge.wmin.and.w_spec.le.wmax) then
!        write(22,*) w_spec*cmiperau, dreal(sfg_spec(i)), -dimag(sfg_spec(i))
!     endif
!  enddo
  ! Then print out the frequencies above the average
!  do i = 1, nt/2 
  do i = 1, nt 
     w_spec = (2d0*pi/dt)*dfloat(i-1)/dfloat(nt) !+ wavg
     if(w_spec.ge.wmin.and.w_spec.le.wmax) then
        write(22,*) w_spec*cmiperau, dreal(sfg_spec(i)), -dimag(sfg_spec(i))
     endif
  enddo
  close(22)

  call dfftw_destroy_plan(plan)

!!!!!!!! End SFG Spectrum Calc


!!!!!!!! Calculate the FT to get the Derivative of the Spectrum

  open(51,file='dH_sfg_spectrum.dat')
  open(52,file='dKE_sfg_spectrum.dat')
  open(53,file='dLJ_sfg_spectrum.dat')
  open(54,file='dCoul_sfg_spectrum.dat')
  open(55,file='dV_sfg_spectrum.dat')
!  open(56,file='dHtop_sfg_spectrum.dat')
!  open(57,file='dHbot_sfg_spectrum.dat')
!  open(58,file='dHt-b_sfg_spectrum.dat')
!  do k = 1, 8
  call dfftw_plan_dft_1d(plan, nt, input, sfg_deriv, FFTW_FORWARD, FFTW_ESTIMATE)
  do k = 1, 5

     ! Pad the calculated TCF with zeroes to get the desired freq resolution before FTing
     input = dcmplx(0.0d0,0.0d0)
     do i = 1, min(nt,ncorr)
        input(i) = tcfH_tot(i-1,k)
     enddo
  
     call dfftw_execute(plan)

     ! Write out the calculated spectrum derivative

     ! When the average frequency is subtracted from the phase the spectrum must be split in the middle
     ! First print out the frequencies below the average
!     do i = nt/2 + 1, nt
!        w_spec = (2d0*pi/dt)*dfloat(i-1-nt)/dfloat(nt) + wavg
!!        w_spec = (2d0*pi/dt)*dfloat(nt + 1 - i)/dfloat(nt) !+ wavg
!        if(w_spec.ge.wmin.and.w_spec.le.wmax) then
!           write(50+k,*) w_spec*cmiperau, dreal(sfg_deriv(i)), -dimag(sfg_deriv(i))
!        endif
!     enddo
     ! Then print out the frequencies above the average
!     do i = 1, nt/2
     do i = 1, nt
        w_spec = (2d0*pi/dt)*dfloat(i-1)/dfloat(nt) !+ wavg
        if(w_spec.ge.wmin.and.w_spec.le.wmax) then
           write(50+k,*) w_spec*cmiperau, dreal(sfg_deriv(i)), -dimag(sfg_deriv(i))
        endif
     enddo
     close(50+k)
  enddo

  call dfftw_destroy_plan(plan)

!!!!!!!! End SFG Derivative Calc

  deallocate(input)

   End Subroutine spec_calc_fluc

