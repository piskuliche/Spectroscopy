! Subroutine to calculate the IR spectral time correlation function

SUBROUTINE calc_tcf(ioh, iTw, w01, w12, mu01, mu12, eOH, tcf_rp, tcf_np)

   USE time_data
   USE map_data
   
   IMPLICIT NONE
   INTEGER :: ioh
   DOUBLE PRECISION, DIMENSION(ntimes) :: mu01, mu12
   DOUBLE PRECISION, DIMENSION(ntimes) :: w01, w12
   DOUBLE PRECISION, DIMENSION(ntimes,3) :: eOH
   DOUBLE PRECISION, DIMENSION(0:ntimes) :: phase01, phase12
   INTEGER :: it0, it1, it2, it3, rel_t1, rel_t2, rel_t3  ! "i" indicates absolute timestep in trajectory, no "i" means relative timestep
   INTEGER :: k, iTw
   INTEGER, DIMENSION(1) :: t2max
   DOUBLE PRECISION :: count, mu01_t0, orient, relax1, relax2
   DOUBLE PRECISION :: dip1, dip3, dip4, dip6, dip7, dip8
   DOUBLE PRECISION :: np1, np3, np4, np6, np7, np8
   DOUBLE complex, DIMENSION(0:ncorr,0:ncorr) :: tcf_rp, tcf_np
   
   ! First calculate the running integral of the 01 and 12 frequencies to USE in the phases

   phase01 = 0d0; phase12 = 0d0
   DO k = 1, ntimes
      phase01(k) = phase01(k-1) + w01(k)
      phase12(k) = phase12(k-1) + w12(k)
   ENDDO

   t2max(1) = NINT(MAXVAL(Tw)/dt)
   rel_t2 = NINT(Tw(iTw)/dt)    ! calculate the # of timesteps in the waiting time

   tcf_rp = dcmplx(0.0d0,0.0d0); tcf_np = dcmplx(0.0d0,0.0d0); count = 0d0
   DO it0 = 1, ntimes - 2*ncorr - t2max(1), nskip  ! Loop over time zeroes indexed by INTEGER it0
      mu01_t0 = mu01(it0)
      count = count + 1d0      

      DO it1 = it0, it0 + ncorr  ! Loop over first time delay, measured relative to it0 by rel_t1
         rel_t1 = it1 - it0

         it2 = it1 + rel_t2             ! This is the point in the trajectory corresponding to Tw (or rel_t2) past it1 
         
         DO it3 = it2, it2 + ncorr  ! Loop over third time delay, measured relative to it2 by rel_t3
            rel_t3 = it3 - it2 
               
            ! Assume unpolarized -- average of XXXX, YYYY and ZZZZ
            orient = (eOH(it0,1)*eOH(it1,1)*eOH(it2,1)*eOH(it3,1) + eOH(it0,2)*eOH(it1,2)*eOH(it2,2)*eOH(it3,2) &
                  + eOH(it0,3)*eOH(it1,3)*eOH(it2,3)*eOH(it3,3) )/3d0 
            
            relax1 = dexp(-0.5d0*DBLE(rel_t3 + 2*rel_t2 + rel_t1)/T1bydt)
            relax2 = dexp(-0.5d0*DBLE(3*rel_t3 + 2*rel_t2 + rel_t1)/T1bydt)
            !relax1 = dexp(-DBLE(rel_t2)/T1bydt)
            !relax2 = dexp(-DBLE(rel_t2)/T1bydt)
            
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
            
            tcf_rp(rel_t1, rel_t3) = tcf_rp(rel_t1, rel_t3) + ( dcmplx(2d0*dip1*dcos(np1),2d0*dip1*dsin(np1)) &
                  + dcmplx(dip3*dcos(np3),dip3*dsin(np3)) )*dcmplx(orient,0d0)
            
            tcf_np(rel_t1, rel_t3) = tcf_np(rel_t1, rel_t3) + ( dcmplx(2d0*dip4*dcos(np4),2d0*dip4*dsin(np4)) &
                  + dcmplx(dip6*dcos(np6),dip6*dsin(np6)) )*dcmplx(orient,0d0) 
            !+ dcmplx(dip7*dcos(np7),dip7*dsin(np7)) &
            !+ dcmplx(dip8*dcos(np8),dip8*dsin(np8)) )*dcmplx(orient,0d0)

         ENDDO  ! End loop over it3 index
      ENDDO   ! End loop over it1 index

   ENDDO  ! End loop over time zeroes, it0
   
   ! Normalize the TCF by the number of time zeroes
   tcf_rp = tcf_rp/dcmplx(count,0d0); tcf_np = tcf_np/dcmplx(count,0d0)

end subroutine calc_tcf

