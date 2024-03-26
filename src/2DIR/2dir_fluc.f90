!    Program to calculate the 2D-IR vibrational spectra
!
  Program spectra

  use map_data
  use time_data
  use freq_data
  use fluc_data
  use constants

  implicit none
  integer, parameter :: nperchunk = 1000
  integer :: ioh, iTw, chunk, nchunks, iper
  double precision :: tstart, tend, read_time, tcf_time, ta, tb
  double precision, allocatable, dimension(:,:) :: w01, w12, mu01, mu12
  double precision, allocatable, dimension(:, :,:) :: eOH
  double precision, allocatable, dimension(:,:) :: dH
  double complex, allocatable, dimension(:,:) :: tcf_rp, tcf_np, tcf_rp_tot, tcf_np_tot
  double complex, allocatable, dimension(:,:,:) :: tcfH_rp, tcfH_np, tcfH_rp_tot, tcfH_np_tot

  ! Read in the parameters for the calculation from spectra.in
  call read_input
  write(6,*) ' Read in input parameters'

  ! Allocate arrays for calculation of the spectra
  allocate(w01(nperchunk,ntimes)); allocate(mu01(nperchunk,ntimes)); allocate(w12(nperchunk,ntimes))
   allocate(mu12(nperchunk,ntimes)); allocate(eOH(nperchunk,ntimes,3))
  allocate(tcf_rp(0:ncorr,0:ncorr)); allocate(tcf_rp_tot(0:ncorr,0:ncorr))
  allocate(tcf_np(0:ncorr,0:ncorr)); allocate(tcf_np_tot(0:ncorr,0:ncorr))

  ! Allocated arrays for calculation of the spectra derivative
  if(flag_fluc) then
     allocate(dH(ntimes,5))
     allocate(tcfH_rp(0:ncorr,0:ncorr,5)); allocate(tcfH_rp_tot(0:ncorr,0:ncorr,5))
     allocate(tcfH_np(0:ncorr,0:ncorr,5)); allocate(tcfH_np_tot(0:ncorr,0:ncorr,5))
     call read_eners(dH)
  !     dw01avg = 0.0d0; dw12avg = 0.0d0; dw01sqavg = 0.0d0; dw12sqavg = 0.0d0
  !     d2w01avg = 0.0d0; d2w01sqavg = 0.0d0; d2w12avg = 0.0d0; d2w12sqavg = 0.0d0
  endif

  write(6,*) ' Allocated arrays'

  ! set the number of OpenMP threads
  !call omp_set_num_threads(20)

  nchunks = ceiling(real(noh)/real(nperchunk))


  ! Loop over the waiting times
  do iTw = 1, nTw

     write(6,'(A,F12.2,A)') ' Waiting time = ',Tw(iTw)*fsperau,' fs'
     call flush(6)

     ! Zero TCFs etc.
     tcf_rp_tot = dcmplx(0.0d0,0.0d0); tcf_np_tot = dcmplx(0.0d0,0.0d0) 
     if(flag_fluc) then
        tcfH_rp_tot = dcmplx(0.0d0,0.0d0); tcfH_np_tot = dcmplx(0.0d0,0.0d0)
     endif
     read_time = 0.0d0; tcf_time = 0.0d0
     w01_avg = 0.0d0; w01_sq_avg = 0.0d0; w12_avg = 0.0d0; w12_sq_avg = 0.0d0     

     ! Loop over the OH groups
     call cpu_time(tstart)

   

   DO chunk=1, nchunks
      w01 = 0.0; w12 = 0.0; mu01 = 0.0; mu12 = 0.0; eOH = 0.0
      DO iper=1, nperchunk
        ioh = (chunk-1)*nperchunk + iper
        if(ioh > noh) exit
        call read_field(ioh, w01(iper,:), w12(iper,:), mu01(iper,:), mu12(iper,:), eOH(iper,:,:))
      END DO 
      WRITE(*,*) "t", w01(1,1), w01_avg

!$omp parallel do private(ioh,ta,tb,tcf_rp,tcf_np) &
!$omp shared(w01,w12,mu01,mu12,eOH) &
!$omp reduction(+:read_time,tcf_time,tcf_rp_tot,tcf_np_tot)
     do ioh=(chunk-1)*nperchunk+1, min(chunk*nperchunk, noh)
      iper = ioh - (chunk-1)*nperchunk
        
        if(flag_fluc) then
           call calc_tcf_fluc(ioh, iTw, w01(iper,:), w12(iper,:), mu01(iper,:), mu12(iper,:), eOH(iper,:,:), &
                  & dH, tcf_rp, tcf_np, tcfH_rp, tcfH_np)
           
           ! Add TCF to the total TCF
           tcf_rp_tot = tcf_rp_tot + tcf_rp
           tcf_np_tot = tcf_np_tot + tcf_np
           tcfH_rp_tot = tcfH_rp_tot + tcfH_rp
           tcfH_np_tot = tcfH_np_tot + tcfH_np
        else
           call calc_tcf(ioh, iTw, w01(iper,:), w12(iper,:), mu01(iper,:), mu12(iper,:), eOH(iper,:,:), tcf_rp, tcf_np)
           
           ! Add TCF to the total TCF
           tcf_rp_tot = tcf_rp_tot + tcf_rp
           tcf_np_tot = tcf_np_tot + tcf_np
        endif
        
        call cpu_time(tb)
        tcf_time = tcf_time + tb - ta     
        
        write(6,'(A,I4,A,F10.2,A)') ' ioh = ',ioh,' time = ',tb-ta,' s'
        call flush(6)
        
     enddo  ! end loop over OH groups
!$omp end parallel do
   END DO ! end loop over chunks

     write(6,*) ' Finished TCF calc'
     call flush(6)

     tcf_rp_tot = tcf_rp_tot/dcmplx(dble(noh),0d0); tcf_np_tot = tcf_np_tot/dcmplx(dble(noh),0d0)
     w01_avg = w01_avg/dble(noh*ntimes); w01_sq_avg = w01_sq_avg/dble(noh*ntimes)
     w12_avg = w12_avg/dble(noh*ntimes); w12_sq_avg = w12_sq_avg/dble(noh*ntimes)
     dw01avg = dw01avg/dble(noh*ntimes); dw01sqavg = dw01sqavg/dble(noh*ntimes)
     dw12avg = dw12avg/dble(noh*ntimes); dw12sqavg = dw12sqavg/dble(noh*ntimes)
     d2w01avg = d2w01avg/dble(noh*ntimes); d2w01sqavg = d2w01sqavg/dble(noh*ntimes)
     d2w12avg = d2w12avg/dble(noh*ntimes); d2w12sqavg = d2w12sqavg/dble(noh*ntimes)
     if(flag_fluc) then
        tcfH_rp_tot = tcfH_rp_tot/dcmplx(dble(noh),0d0)
        tcfH_np_tot = tcfH_np_tot/dcmplx(dble(noh),0d0)
     endif

     ! Calculate (and output) the 2D-IR spectra from the total TCF
     call cpu_time(ta)
     if(flag_fluc) then
        call spec_calc_fluc(iTw, tcf_rp_tot, tcf_np_tot, tcfH_rp_tot, tcfH_np_tot)
     else
        call spec_calc(iTw, tcf_rp_tot,tcf_np_tot)
     endif
     call cpu_time(tb)

     call cpu_time(tend)
     write(6,'(A,F12.2,A)') ' Waiting time = ',Tw(iTw)*fsperau,' fs'
     write(6,'(A,F12.2,A,F10.2,A)') ' Read  cpu time = ',read_time,' s, ',read_time/60d0,' min'
     write(6,'(A,F12.2,A,F10.2,A)') ' TCF   cpu time = ',tcf_time,' s, ',tcf_time/60d0,' min'
     write(6,'(A,F12.2,A,F10.2,A)') ' FFT   cpu time = ',tb-ta,' s, ',(tb-ta)/60d0,' min'
     write(6,'(A,F12.2,A,F10.2,A)') ' Total cpu time = ',tend-tstart,' s, ',(tend-tstart)/60d0,' min'

  enddo  ! End loop over waiting times

  deallocate(w01); deallocate(mu01); deallocate(w12); deallocate(mu12); deallocate(eOH)
  deallocate(tcf_rp); deallocate(tcf_np); deallocate(tcf_rp_tot); deallocate(tcf_np_tot)

  End Program spectra


