!    Program to calculate the vibrational spectra

  Program spectra

  use map_data
  use time_data
  use freq_data
  use hist_data
  use fluc_data

  implicit none
  integer :: k,ioh
  double precision :: apara,aperp,tstart, tend, read_time, tcf_time, ta, tb
  double precision, allocatable, dimension(:) :: w,a_ss,a_sp,a_pp,zO
  double precision, allocatable, dimension(:,:) :: e, mu, dH
  double precision, allocatable, dimension(:) :: wdist, specdist, wd_tot, sd_tot
  double precision, allocatable, dimension(:,:) :: wdistH, specdistH, wdH_tot, sdH_tot
  double precision, allocatable, dimension(:,:) :: wdistH2, specdistH2, wdH2_tot, sdH2_tot
  double complex, allocatable, dimension(:) :: tcf, tcf_tot
  double complex, allocatable, dimension(:,:) :: tcfH, tcfH_tot

  ! Read in the parameters for the calculation from spectra.in
  call read_input

!  write(6,*) "The input has been read"

  ! Allocate all arrays that are used regardless of flags
  allocate(w(ntimes)); allocate(mu(ntimes,3)) 
  allocate(e(ntimes,3)); allocate(zO(ntimes))
  allocate(a_ss(ntimes)); allocate(a_sp(ntimes)); allocate(a_pp(ntimes))

  allocate(tcf(0:ncorr)); allocate(tcf_tot(0:ncorr))

  ! Allocate arrays use for frequency distribution/spectral density
!  if(flag_hist) then
     allocate(wdist(0:nhist)); allocate(specdist(0:nhist)); allocate(wd_tot(0:nhist)); allocate(sd_tot(0:nhist))
     wd_tot = 0.0d0; sd_tot = 0.0d0
!     if(flag_fluc) then
        allocate(wdistH(0:nhist,8)); allocate(specdistH(0:nhist,8))
        allocate(wdistH2(0:nhist,8)); allocate(specdistH2(0:nhist,8))
        allocate(wdH_tot(0:nhist,8)); allocate(sdH_tot(0:nhist,8))
        allocate(wdH2_tot(0:nhist,8)); allocate(sdH2_tot(0:nhist,8))
        wdH_tot = 0.0d0; sdH_tot = 0.0d0; wdH2_tot = 0.0d0; sdH2_tot = 0.0d0
!     endif
!  endif

  ! Allocate arrays used for flag_fluc and read in the energies
!  if(flag_fluc) then
       allocate(tcfH(0:ncorr,8)); allocate(tcfH_tot(0:ncorr,8)); allocate(dH(ntimes,8))
       tcfH_tot = dcmplx(0.0d0,0.0d0)
       ! Read in the energies
       if(flag_fluc) then
          call read_eners(dH)
       endif
       dwavg = 0.0d0; dw2avg = 0.0d0; d2wavg = 0.0d0; d2w2avg = 0.0d0
!  endif

  tcf_tot = dcmplx(0.0d0,0.0d0); read_time = 0.0d0; tcf_time = 0.0d0
  wavg = 0.0d0; w2avg = 0.0d0
!  write(6,*) ' Made it to loop over OH groups'

  ! Loop over the OH groups
  call cpu_time(tstart)


  ! set the number of OpenMP threads
  call omp_set_num_threads(56)

!$omp parallel do private(ioh,ta,tb,w,mu,e,a_ss,a_sp,a_pp,zO) &
!$omp private(wdist,specdist,wdistH,specdistH,wdistH2,specdistH2) &
!$omp private(tcf,tcfH) &
!$omp reduction(+:read_time,tcf_time,wd_tot,sd_tot,wdH_tot,sdH_tot,wdH2_tot,sdH2_tot,tcf_tot,tcfH_tot)
  do ioh = 1, noh
     ! Read in the field file for OH number ioh
     call cpu_time(ta)

     if(flag_fluc) then
        call read_field_fluc(ioh, w, mu, e, a_ss, a_sp, a_pp, zO, dH)
     else
        call read_field(ioh, w, mu, e, a_ss, a_sp, a_pp, zO)
     endif
!     write(6,*) ' Finished field read, ioh = ',ioh
!     call flush(6)

     call cpu_time(tb)
     read_time = read_time + tb - ta

     ! Calculate the frequency distribution and spectral density for the OH number ioh & add to total
     if(flag_hist) then
       if(flag_fluc) then
           call hist_calc_fluc(w, mu, a_ss, zO, wdist, specdist, dH, wdistH, specdistH, wdistH2, specdistH2)
           wd_tot = wd_tot + wdist
           sd_tot = sd_tot + specdist
           wdH_tot = wdH_tot + wdistH
           sdH_tot = sdH_tot + specdistH
           wdH2_tot = wdH2_tot + wdistH2
           sdH2_tot = sdH2_tot + specdistH2
        else
           call hist_calc(w, mu, a_ss, zO, wdist, specdist)
           wd_tot = wd_tot + wdist
           sd_tot = sd_tot + specdist
        endif
     endif

!     write(6,*) ' Finished Hist calc, ioh = ',ioh
!     call flush(6)

     ! Calculate the SFG spectrum TCF for the OH number ioh 
     call cpu_time(ta)
     if(flag_fluc) then
        call calc_tcf_fluc(w, mu, a_ss, zO, dH, tcf, tcfH) ! a_sp and a_pp?

        ! Add TCF to the total TCF
        tcf_tot = tcf_tot + tcf
        tcfH_tot = tcfH_tot + tcfH
     else
       call calc_tcf(w, mu, a_ss, zO, tcf)

       ! Add TCF to the total TCF
       tcf_tot = tcf_tot + tcf
     endif

!     write(6,*) ' Finished TCF calc, ioh = ',ioh
!     call flush(6)

     call cpu_time(tb)
     tcf_time = tcf_time + tb - ta

  enddo  ! end loop over OH groups
!$omp end parallel do

  tcf_tot = tcf_tot/dcmplx(dfloat(noh),0d0)
  wavg = wavg/dfloat(noh*ntimes); w2avg = w2avg/dfloat(noh*ntimes)
  dwavg = dwavg/dfloat(noh*ntimes); dw2avg = dw2avg/dfloat(noh*ntimes)
  d2wavg = d2wavg/dfloat(noh*ntimes); d2w2avg = d2w2avg/dfloat(noh*ntimes)

  if(flag_fluc) tcfH_tot = tcfH_tot/dcmplx(dfloat(noh),0d0)

!  write(6,*) ' OH loop'
!  call flush(6)

  ! Printout the distributions if calculated
  if(flag_hist) then
     if(flag_fluc) then
        call hist_print_fluc(wd_tot, sd_tot, wdH_tot, sdH_tot, wdH2_tot, sdH2_tot)
     else
        call hist_print(wd_tot, sd_tot)
     endif
  endif

!  write(6,*) ' Finished hist print'
!  call flush(6)

  ! Calculate the SFG spectrum from the total TCF
  call cpu_time(ta)

  if(flag_fluc) then
     call spec_calc_fluc(tcf_tot, tcfH_tot)
  else
     call spec_calc(tcf_tot)
  endif

!  write(6,*) ' Finished spec_calc'

  call cpu_time(tb)

  call cpu_time(tend)

  write(6,'(A,F10.4,A,F10.2,A)') ' Read  cpu time = ',read_time,' s, ',read_time/60d0,' min'
  write(6,'(A,F10.4,A,F10.2,A)') ' TCF   cpu time = ',tcf_time,' s, ',tcf_time/60d0,' min'
  write(6,'(A,F10.4,A,F10.2,A)') ' FFT   cpu time = ',tb-ta,' s, ',(tb-ta)/60d0,' min'
  write(6,'(A,F10.4,A,F10.2,A)') ' Total cpu time = ',tend-tstart,' s, ',(tend-tstart)/60d0,' min'

  deallocate(w); deallocate(mu); deallocate(e)
  deallocate(a_ss); deallocate(a_sp); deallocate(a_pp); deallocate(zO)

  End Program spectra
