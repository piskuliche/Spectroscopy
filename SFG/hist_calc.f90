! Subroutine to calculate the frequency distribution & spectral density
!    for a single OH group

   subroutine hist_calc(w, mu, a_ss, zO, wdist, specdist)

     use time_data
     use hist_data
     use constants
     
     implicit none
     
     integer :: k, m, iw

     double precision :: count, a_mu
     double precision, dimension(ntimes,3) :: mu
     double precision, dimension(ntimes) :: w, a_ss, zO
     double precision, dimension(0:nhist) :: wdist, specdist
     
     ! Zero the arrays
     wdist = 0.0d0; specdist = 0.0d0; count = 0.0d0

     do k = 1, ntimes

        a_mu = a_ss(k)*mu(k,3)  !  this is just alpha_ss m_z for SSP polarization
        iw = nint( (w(k) - wmin)/dw )

        if(iw.ge.0.and.iw.le.nhist) then
           count = count + 1d0
           wdist(iw) = wdist(iw) + 1d0
           specdist(iw) = specdist(iw) + a_mu
        endif

     enddo
     
     ! Normalize the distributions by the number of times
     wdist = wdist/count; specdist = specdist/count

   end subroutine hist_calc


! Subroutine to calculate the frequency distribution & spectral density
!    for a single OH group INCLUDING the derivative and its components 
!    from fluctuation theory

   subroutine hist_calc_fluc(w, mu, a_ss, zO, wdist, specdist, dH, wdistH, specdistH, wdistH2, specdistH2)

     use time_data
     use hist_data
     use fluc_data
     use constants
     
     implicit none
     
     integer :: k, m, n, iw

     double precision :: count, a_mu
     double precision, dimension(ntimes,3) :: mu
     double precision, dimension(ntimes) :: w, a_ss, zO
     double precision, dimension(0:nhist) :: wdist, specdist
     double precision, dimension(ntimes,8) :: dH
     double precision, dimension(0:nhist,8) :: wdistH, specdistH
     double precision, dimension(0:nhist,8) :: wdistH2, specdistH2
     
     ! Zero the arrays
     wdist = 0.0d0; specdist = 0.0d0; wdistH = 0.0d0; specdistH = 0.0d0; wdistH2 = 0.0d0; specdistH2 = 0.0d0 
     count = 0.0d0

     do k = 1, ntimes

        a_mu = a_ss(k)*mu(k,3)  !  this is just alpha_ss m_z for SSP polarization
        iw = nint( (w(k) - wmin)/dw )

        if(iw.ge.0.and.iw.le.nhist) then
           count = count + 1d0
           wdist(iw) = wdist(iw) + 1d0
           specdist(iw) = specdist(iw) + a_mu
           
!           do n = 1, 8
           do n = 1, 5
!              wdistH(iw,n) = wdistH(iw,n) + (dH(k,n) + Havg(n))
              wdistH(iw,n) = wdistH(iw,n) + dH(k,n)
              specdistH(iw,n) = specdistH(iw,n) + a_mu*dH(k,n)
!              do m = 1, 8
!                 wdistH2(iw,n,m) = wdistH2(iw,n,m) + (dH(k,n) + Havg(n))*(dH(k,m) + Havg(m))
              wdistH2(iw,n) = wdistH2(iw,n) + dH(k,1)*dH(k,n)
!                 specdistH2(iw,n,m) = specdistH2(iw,n,m) + a_mu*(dH(k,n) + Havg(n))*(dH(k,m) + Havg(m))
              specdistH2(iw,n) = specdistH2(iw,n) + a_mu*dH(k,1)*dH(k,n)
!              enddo
           enddo

        endif

     enddo
     
     ! Normalize the distributions by the number of times
     wdist = wdist/count; specdist = specdist/count
     wdistH = wdistH/count; specdistH = specdistH/count
     wdistH2 = wdistH2/count; specdistH2 = specdistH2/count

   end subroutine hist_calc_fluc


! Subroutine to calculate the frequency distribution & spectral density
!    for a single OH group

   subroutine hist_print(wd_tot, sd_tot)

     use time_data
     use hist_data
     use constants

     implicit none
     integer :: k
     double precision, dimension(0:nhist) :: wd_tot, sd_tot

     open(23,file='freq_dist.dat')
     open(24,file='spec_dens.dat')

     wd_tot = wd_tot/real(noh)
     sd_tot = sd_tot/real(noh)

     do k = 0, nhist
        write(23,*) (wmin + real(k)*dw)*cmiperau, wd_tot(k)
        write(24,*) (wmin + real(k)*dw)*cmiperau, sd_tot(k)
     enddo
     close(23)
     close(24)

   end subroutine hist_print


! Subroutine to calculate the frequency distribution & spectral density
!    for a single OH group including the derivatives from fluctuation theory

   subroutine hist_print_fluc(wd_tot, sd_tot, wdH_tot, sdH_tot, wdH2_tot, sdH2_tot)

     use time_data
     use hist_data
     use constants

     implicit none
     integer :: k, n
     double precision, dimension(0:nhist) :: wd_tot, sd_tot
     double precision, dimension(0:nhist,8) :: wdH_tot, sdH_tot
     double precision, dimension(0:nhist,8) :: wdH2_tot, sdH2_tot

     open(23,file='freq_dist.dat')
     open(24,file='spec_dens.dat')
     open(25,file='dH_freq_dist.dat')
     open(26,file='dH_spec_dens.dat')
     open(27,file='dH2_freq_dist.dat')
     open(28,file='dH2_spec_dens.dat')

     wd_tot = wd_tot/real(noh); sd_tot = sd_tot/real(noh)
     wdH_tot = wdH_tot/real(noh); sdH_tot = sdH_tot/real(noh)
     wdH2_tot = wdH2_tot/real(noh); sdH2_tot = sdH2_tot/real(noh)

     do k = 0, nhist
        write(23,*) (wmin + real(k)*dw)*cmiperau, wd_tot(k)
        write(24,*) (wmin + real(k)*dw)*cmiperau, sd_tot(k)
!        write(25,'(15F12.5)') (wmin + real(k)*dw)*cmiperau, (wdH_tot(k,n), n=1,8)
!        write(26,'(15F12.5)') (wmin + real(k)*dw)*cmiperau, (sdH_tot(k,n), n=1,8)
!        write(27,'(15E17.6)') (wmin + real(k)*dw)*cmiperau, (wdH2_tot(k,n), n=1,8) 
!        write(28,'(15E17.6)') (wmin + real(k)*dw)*cmiperau, (sdH2_tot(k,n), n=1,8) 
        write(25,'(15F12.5)') (wmin + real(k)*dw)*cmiperau, (wdH_tot(k,n), n=1,5)
        write(26,'(15F12.5)') (wmin + real(k)*dw)*cmiperau, (sdH_tot(k,n), n=1,5)
        write(27,'(15E17.6)') (wmin + real(k)*dw)*cmiperau, (wdH2_tot(k,n), n=1,5) 
        write(28,'(15E17.6)') (wmin + real(k)*dw)*cmiperau, (sdH2_tot(k,n), n=1,5) 
     enddo
     close(23); close(24); close(25); close(26); close(27); close(28)

   end subroutine hist_print_fluc
