! Subroutine to calculate the SFG spectral time correlation function

   subroutine calc_tcf(w, mu, a_ij, zO, tcf)

     use time_data
     
     implicit none
     
     double precision, dimension(ntimes,3) :: mu
     double precision, dimension(ntimes) :: w, a_ij, zO
     integer :: ilag
     integer :: k,m
     double precision :: count, dp, phase
     double precision :: a_ij0
     double complex, dimension(0:ncorr) :: tcf
     
     tcf = dcmplx(0.0d0,0.0d0); count = 0d0

     do k = 1, ntimes-ncorr, nskip
        a_ij0 = a_ij(k)*dsign(1d0,zO(k))
        count = count + 1d0      
        phase = -w(k)
       
        do m = k, k + ncorr
           ilag = m - k
           phase = phase + w(m)
           dp = a_ij0*mu(m,3) !! put in the sign! *and only works with a_ijp*        
           tcf(ilag) = tcf(ilag) + dp*dcmplx(dcos(phase*dt),dsin(phase*dt))        
        enddo
     enddo
     
     ! Normalize the TCF by the number of time zeroes
     tcf = tcf/dcmplx(count,0d0)     

   end subroutine calc_tcf



! Subroutine to calculate the SFG spectral time correlation function

   subroutine calc_tcf_fluc(w, mu, a_ij, zO, dH, tcf, tcfH)

     use time_data
     
     implicit none
     
     double precision, dimension(ntimes,3) :: mu
     double precision, dimension(ntimes) :: w, a_ij, zO
     double precision, dimension(ntimes,8) :: dH
     integer :: ilag
     integer :: k, m, n
     double precision :: count, dp, phase, sz
     double precision, dimension(3) :: mu0
     double precision :: a_ij0
     double complex, dimension(0:ncorr) :: tcf
     double complex, dimension(0:ncorr,8) :: tcfH
     
     tcf = dcmplx(0.0d0,0.0d0); tcfH = dcmplx(0.0d0,0.0d0); count = 0d0

     do k = 1, ntimes-ncorr, nskip
        sz = dsign(1d0,zO(k))
        a_ij0 = a_ij(k)*sz
        count = count + 1d0      
        phase = -w(k)
       
        do m = k, k + ncorr
           ilag = m - k
           phase = phase + w(m)
           dp = a_ij0*mu(m,3)          
           tcf(ilag) = tcf(ilag) + dp*dcmplx(dcos(phase*dt),dsin(phase*dt))

           ! Calculate the derivative TCF for each energy category
           do n = 1, 5 
              tcfH(ilag,n) = tcfH(ilag,n) + dH(k,n)*dp*dcmplx(dcos(phase*dt),dsin(phase*dt))
           enddo
!           tcfH(ilag,6) = tcfH(ilag,6) + 0.5d0*( dH(k,6)*(sz+1d0) + dH(k,7)*(1d0-sz) )*dp*dcmplx(dcos(phase*dt),dsin(phase*dt))
!           tcfH(ilag,7) = tcfH(ilag,7) + 0.5d0*( dH(k,7)*(sz+1d0) + dH(k,6)*(1d0-sz) )*dp*dcmplx(dcos(phase*dt),dsin(phase*dt))
!           tcfH(ilag,8) = tcfH(ilag,8) + dH(k,8)*dp*dcmplx(dcos(phase*dt),dsin(phase*dt))
        enddo
     enddo
     
     ! Normalize the TCF by the number of time zeroes
     tcf = tcf/dcmplx(count,0d0)
     tcfH = tcfH/dcmplx(count,0d0)

   end subroutine calc_tcf_fluc

