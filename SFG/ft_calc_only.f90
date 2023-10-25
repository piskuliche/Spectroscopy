!   Program to Fourier Transform the total TCF & derivatives to calculate the SFG spectrum

Program spec_calc_fluc

  !use constants

  implicit none
  include 'fftw3.f'
  integer :: i, k, p, nt, ncorr
  integer*8 :: plan
  double precision :: dt, wmin, wmax, w_resol, T1
  double precision :: w_spec, ti, pi, tmp, tmpr, tmpi

  ! Constants
  double precision :: fsperau=2.4188843265857d-2
  double precision :: cmiperau=2.1947463d5

  double complex, allocatable, dimension(:) :: tcf_tot
  double complex, allocatable, dimension(:,:) :: tcfH_tot
  double complex, allocatable, dimension(:) :: input, sfg_spec, sfg_deriv

  pi = 4.0d0*datan(1.0d0)

  ! Read in input
  open(11,file='ft_only.in',status='old')
  read(11,*)
  read(11,*) ncorr, dt
  read(11,*) 
  read(11,*) T1, wmin, wmax, w_resol
  close(11)

  ! Convert units
  dt = dt/fsperau
  T1 = T1/fsperau
  wmin = wmin/cmiperau; wmax = wmax/cmiperau
  w_resol = w_resol/cmiperau

  ! Allocate arrays
  allocate(tcf_tot(0:ncorr)); allocate(tcfH_tot(0:ncorr,4))

  ! Zero the arrays
  tcf_tot = dcmplx(0d0,0d0); tcfH_tot = dcmplx(0d0,0d0) 

  ! First read in the SFG TCF and the derivative TCFs
  open(21,file='sfg_tcf.dat',status='old')
  !open(21,file='real_sfg_tcf.dat',status='old')
  !open(22,file='imag_sfg_tcf.dat',status='old')
  open(23,file='dH_sfg_tcf.dat',status='old')
  !open(23,file='real_dH_sfg_tcf.dat',status='old')
  !open(24,file='imag_dH_sfg_tcf.dat',status='old')
  open(25,file='dKE_sfg_tcf.dat',status='old')
  open(27,file='dLJ_sfg_tcf.dat',status='old')
  open(29,file='dCoul_sfg_tcf.dat',status='old')
  do i = 0, ncorr
     read(21,*) tmp, tmpr, tmpi
     !read(22,*) tmp, tmpi
     tcf_tot(i) = dcmplx(tmpi, tmpr)  !! Not sure why the parts need to be switched to get the right answer!!!!

     ! Read in total derivative of TCF
     read(23,*) tmp, tmpr, tmpi
     !read(24,*) tmp, tmpi
     tcfH_tot(i,1) = dcmplx(tmpi, tmpr)

     ! Read in KE derivative of TCF
     read(25,*) tmp, tmpr, tmpi
     tcfH_tot(i,2) = dcmplx(tmpi, tmpr)
     ! Read in LJ derivative of TCF
     read(27,*) tmp, tmpr, tmpi
     tcfH_tot(i,3) = dcmplx(tmpi, tmpr)
     ! Read in Coul derivative of TCF
     read(29,*) tmp, tmpr, tmpi
     tcfH_tot(i,4) = dcmplx(tmpi, tmpr)
  enddo
  !close(21); close(22); close(23); close(24)
  close(21); close(23); close(25); close(27); close(29)

  ! Multiply the total TCF by the vibrational relaxation factor
  do i = 0, ncorr
     ti = float(i)*dt
     tcf_tot(i) = tcf_tot(i)*dcmplx(dexp(-0.5d0*ti/T1),0d0)
     tcfH_tot(i,:) = tcfH_tot(i,:)*dcmplx(dexp(-0.5d0*ti/T1),0d0)
  enddo


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
  open(22,file='ft_sfg_spectrum.dat')

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

  open(51,file='ft_dH_sfg_spectrum.dat')
  open(52,file='ft_dKE_sfg_spectrum.dat')
  open(53,file='ft_dLJ_sfg_spectrum.dat')
  open(54,file='ft_dCoul_sfg_spectrum.dat')
!  open(55,file='dV_sfg_spectrum.dat')
!  open(56,file='dHtop_sfg_spectrum.dat')
!  open(57,file='dHbot_sfg_spectrum.dat')
!  open(58,file='dHt-b_sfg_spectrum.dat')
!  do k = 1, 8
  call dfftw_plan_dft_1d(plan, nt, input, sfg_deriv, FFTW_FORWARD, FFTW_ESTIMATE)
  do k = 1, 4

     ! Pad the calculated TCF with zeroes to get the desired freq resolution before FTing
     input = dcmplx(0.0d0,0.0d0)
     do i = 1, min(nt,ncorr)
        input(i) = tcfH_tot(i-1,k)
     enddo
  
     call dfftw_execute(plan)

     ! Write out the calculated spectrum derivative

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

End Program spec_calc_fluc

