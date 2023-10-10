! Subroutine to calculate the IR spectral time correlation function

   subroutine calc_tcf(ioh, iTw, w01, w12, mu01, mu12, eOH, tcf_rp, tcf_np)

     use time_data
     use map_data
     
     implicit none
     integer :: ioh
     double precision, dimension(ntimes) :: mu01, mu12
     double precision, dimension(ntimes) :: w01, w12
     double precision, dimension(ntimes,3) :: eOH
     double precision, dimension(0:ntimes) :: phase01, phase12
     integer :: it0, it1, it2, it3, t1, t2, t3  ! "i" indicates absolute timestep in trajectory, no "i" means relative timestep
     integer :: k, iTw
     integer, dimension(1) :: t2max
     double precision :: count, mu01_t0, orient, relax1, relax2
     double precision :: dip1, dip3, dip4, dip6, dip7, dip8
     double precision :: np1, np3, np4, np6, np7, np8
     double complex, dimension(0:ncorr,0:ncorr) :: tcf_rp, tcf_np
     
     ! First calculate the running integral of the 01 and 12 frequencies to use in the phases

     phase01 = 0d0; phase12 = 0d0
     do k = 1, ntimes
        phase01(k) = phase01(k-1) + w01(k)
        phase12(k) = phase12(k-1) + w12(k)
     enddo

     t2max(1) = nint(maxval(Tw)/dt)
     t2 = nint(Tw(iTw)/dt)    ! calculate the # of timesteps in the waiting time

     tcf_rp = dcmplx(0.0d0,0.0d0); tcf_np = dcmplx(0.0d0,0.0d0); count = 0d0
     do it0 = 1, ntimes - 2*ncorr - t2max(1), nskip  ! Loop over time zeroes indexed by integer it0
        mu01_t0 = mu01(it0)
        count = count + 1d0      

        do it1 = it0, it0 + ncorr  ! Loop over first time delay, measured relative to it0 by t1
           t1 = it1 - it0

           it2 = it1 + t2             ! This is the point in the trajectory corresponding to Tw (or t2) past it1 
           
           do it3 = it2, it2 + ncorr  ! Loop over third time delay, measured relative to it2 by t3
              t3 = it3 - it2 
                 
              ! Assume unpolarized -- average of XXXX, YYYY and ZZZZ
              orient = (eOH(it0,1)*eOH(it1,1)*eOH(it2,1)*eOH(it3,1) + eOH(it0,2)*eOH(it1,2)*eOH(it2,2)*eOH(it3,2) &
                   + eOH(it0,3)*eOH(it1,3)*eOH(it2,3)*eOH(it3,3) )/3d0 
              
              relax1 = dexp(-0.5d0*dble(t3 + 2*t2 + t1)/T1bydt)
              relax2 = dexp(-0.5d0*dble(3*t3 + 2*t2 + t1)/T1bydt)
              !relax1 = dexp(-dble(t2)/T1bydt)
              !relax2 = dexp(-dble(t2)/T1bydt)
              
              dip1 = mu01_t0*mu01(it1)*mu01(it2)*mu01(it3)*relax1
              dip3 = -mu01_t0*mu01(it1)*mu12(it2)*mu12(it3)*relax2
              dip4 = dip1
              dip6 = dip3
              !dip7 = mu01_t0*mu12(it1)*mu12(it2)*mu01(it3)*relax1
              !dip8 = -mu01_t0*mu12(it1)*mu01(it2)*mu12(it3)*relax2
              
              np1 = dt*(  phase01(it1) - phase01(it0) - phase01(it3) + phase01(it2) )
              np3 = dt*(  phase01(it1) - phase01(it0) - phase12(it3) + phase12(it2) )
              np4 = dt*( -phase01(it1) + phase01(it0) - phase01(it3) + phase01(it2) )
              np6 = dt*( -phase01(it1) + phase01(it0) - phase12(it3) + phase12(it2) )
              !np7 = dt*( -phase01(it3) + phase01(it0) - phase12(it2) + phase12(it1) )
              !np8 = dt*( -phase01(it2) + phase01(it0) - phase12(it3) + phase12(it1) )
              
              tcf_rp(t1, t3) = tcf_rp(t1, t3) + ( dcmplx(2d0*dip1*dcos(np1),2d0*dip1*dsin(np1)) &
                   + dcmplx(dip3*dcos(np3),dip3*dsin(np3)) )*dcmplx(orient,0d0)
              
              tcf_np(t1, t3) = tcf_np(t1, t3) + ( dcmplx(2d0*dip4*dcos(np4),2d0*dip4*dsin(np4)) &
                   + dcmplx(dip6*dcos(np6),dip6*dsin(np6)) )*dcmplx(orient,0d0) 
              !+ dcmplx(dip7*dcos(np7),dip7*dsin(np7)) &
              !+ dcmplx(dip8*dcos(np8),dip8*dsin(np8)) )*dcmplx(orient,0d0)

           enddo  ! End loop over it3 index
        enddo   ! End loop over it1 index

     enddo  ! End loop over time zeroes, it0
     
     ! Normalize the TCF by the number of time zeroes
     tcf_rp = tcf_rp/dcmplx(count,0d0); tcf_np = tcf_np/dcmplx(count,0d0)

   end subroutine calc_tcf



! Subroutine to calculate the IR spectral time correlation function and
!   correlations with the energy fluctuations

   subroutine calc_tcf_fluc(ioh, iTw, w01, w12, mu01, mu12, eOH, dH, tcf_rp, tcf_np, tcfH_rp, tcfH_np)

     use time_data
     use map_data
     
     implicit none
     integer :: ioh
     double precision, dimension(ntimes) :: mu01, mu12
     double precision, dimension(ntimes) :: w01, w12
     double precision, dimension(ntimes,3) :: eOH
     double precision, dimension(ntimes,5) :: dH
     double precision, dimension(0:ntimes) :: phase01, phase12
     integer :: it0, it1, it2, it3, t1, t2, t3  ! "i" indicates absolute timestep in trajectory, no "i" means relative timestep
     integer :: j, k, iTw
     integer, dimension(1) :: t2max
     double precision :: count, mu01_t0, orient, relax1, relax2
     double precision :: dip1, dip3, dip4, dip6, dip7, dip8
     double precision :: np1, np3, np4, np6, np7, np8
     double complex, dimension(0:ncorr,0:ncorr) :: tcf_rp, tcf_np
     double complex, dimension(0:ncorr,0:ncorr,5) :: tcfH_rp, tcfH_np
     
     ! First calculate the running integral of the 01 and 12 frequencies to use in the phases

     phase01 = 0d0; phase12 = 0d0
     do k = 1, ntimes
        phase01(k) = phase01(k-1) + w01(k)
        phase12(k) = phase12(k-1) + w12(k)
     enddo

     t2max(1) = nint(maxval(Tw)/dt)
     t2 = nint(Tw(iTw)/dt)    ! calculate the # of timesteps in the waiting time
           
     tcf_rp = dcmplx(0.0d0,0.0d0); tcf_np = dcmplx(0.0d0,0.0d0); count = 0d0
     tcfH_rp = dcmplx(0.0d0,0.0d0); tcfH_np = dcmplx(0.0d0,0.0d0)
     do it0 = 1, ntimes - 2*ncorr - t2max(1), nskip  ! Loop over time zeroes indexed by integer it0
        mu01_t0 = mu01(it0)
        count = count + 1d0      

        do it1 = it0, it0 + ncorr  ! Loop over first time delay, measured relative to it0 by t1
           t1 = it1 - it0

           it2 = it1 + t2             ! This is the point in the trajectory corresponding to Tw (or t2) past it1 
           
           do it3 = it2, it2 + ncorr  ! Loop over third time delay, measured relative to it2 by t3
              t3 = it3 - it2 
                 
              ! Assume unpolarized -- average of XXXX, YYYY and ZZZZ
              orient = (eOH(it0,1)*eOH(it1,1)*eOH(it2,1)*eOH(it3,1) + eOH(it0,2)*eOH(it1,2)*eOH(it2,2)*eOH(it3,2) &
                   + eOH(it0,3)*eOH(it1,3)*eOH(it2,3)*eOH(it3,3) )/3d0 
              
              relax1 = dexp(-0.5d0*dble(t3 + 2*t2 + t1)/T1bydt)
              relax2 = dexp(-0.5d0*dble(3*t3 + 2*t2 + t1)/T1bydt)
              !relax1 = dexp(-dble(t2)/T1bydt)
              !relax2 = dexp(-dble(t2)/T1bydt)
              
              dip1 = mu01_t0*mu01(it1)*mu01(it2)*mu01(it3)*relax1
              dip3 = -mu01_t0*mu01(it1)*mu12(it2)*mu12(it3)*relax2
              dip4 = dip1
              dip6 = dip3
              !dip7 = mu01_t0*mu12(it1)*mu12(it2)*mu01(it3)*relax1
              !dip8 = -mu01_t0*mu12(it1)*mu01(it2)*mu12(it3)*relax2
              
              np1 = dt*(  phase01(it1) - phase01(it0) - phase01(it3) + phase01(it2) )
              np3 = dt*(  phase01(it1) - phase01(it0) - phase12(it3) + phase12(it2) )
              np4 = dt*( -phase01(it1) + phase01(it0) - phase01(it3) + phase01(it2) )
              np6 = dt*( -phase01(it1) + phase01(it0) - phase12(it3) + phase12(it2) )
              !np7 = dt*( -phase01(it3) + phase01(it0) - phase12(it2) + phase12(it1) )
              !np8 = dt*( -phase01(it2) + phase01(it0) - phase12(it3) + phase12(it1) )

              tcf_rp(t1, t3) = tcf_rp(t1, t3) +  ( dcmplx(2d0*dip1*dcos(np1),2d0*dip1*dsin(np1)) &
                   + dcmplx(dip3*dcos(np3),dip3*dsin(np3)) )*dcmplx(orient,0d0)

              tcf_np(t1, t3) = tcf_np(t1, t3) + ( dcmplx(2d0*dip4*dcos(np4),2d0*dip4*dsin(np4)) &
                   + dcmplx(dip6*dcos(np6),dip6*dsin(np6)) )*dcmplx(orient,0d0)
              !+ dcmplx(dip7*dcos(np7),dip7*dsin(np7)) &
              !+ dcmplx(dip8*dcos(np8),dip8*dsin(np8)) )*dcmplx(orient,0d0)

              do j = 1, 5
                 tcfH_rp(t1, t3, j) = tcfH_rp(t1, t3, j) +  ( dcmplx(2d0*dip1*dcos(np1),2d0*dip1*dsin(np1)) &
                      + dcmplx(dip3*dcos(np3),dip3*dsin(np3)) )*dcmplx(orient*dH(it0,j),0d0)                 
                 tcfH_np(t1, t3, j) = tcfH_np(t1, t3, j) + ( dcmplx(2d0*dip4*dcos(np4),2d0*dip4*dsin(np4)) &
                      + dcmplx(dip6*dcos(np6),dip6*dsin(np6)) )*dcmplx(orient*dH(it0,j),0d0)
                 !+ dcmplx(dip7*dcos(np7),dip7*dsin(np7)) &
                 !+ dcmplx(dip8*dcos(np8),dip8*dsin(np8)) )*dcmplx(orient*dH(it0,j),0d0)
              enddo
              
           enddo  ! End loop over it3 index
        enddo   ! End loop over it1 index

     enddo  ! End loop over time zeroes, it0
     
     ! Normalize the TCF by the number of time zeroes
     tcf_rp = tcf_rp/dcmplx(count,0d0); tcf_np = tcf_np/dcmplx(count,0d0)
     tcfH_rp = tcfH_rp/dcmplx(count,0d0); tcfH_np = tcfH_np/dcmplx(count,0d0)

   end subroutine calc_tcf_fluc



