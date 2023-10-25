!   Subroutine to read in the field parameters from field.xxx

  Subroutine read_field(ioh,w,mu,e,a_ss,a_sp,a_pp,zO)

  use time_data
  use map_data
  use freq_data

  implicit none
  integer :: ioh, j, k 
  character*4 :: ext
  double precision :: etmp, mutmp, xtmp,alpha
  double precision :: axx,ayy,azz,axy,ayz,azx
  double precision, dimension(ntimes) :: w,a_ss,a_sp,a_pp,zO
  double precision, dimension(ntimes,3) :: mu, e
  double precision :: apara,aperp


  ! Open the field field file for OH group iOH
  if(ndigits.eq.4) then
     write(ext,'(I0.4)') ioh
  elseif(ndigits.eq.3) then
     write(ext,'(I0.3)') ioh 
  else
     write(6,*) ' Invalid ndigits, only 3 or 4 are allowed'
     stop
  endif

  ! Open the field file
  open(ioh+100,file='field_files/field.'//ext//'',status='old')

  do k = 1 , ntimes

     ! read in data from file (here e = e_OH)
     read(ioh+100,*) etmp, (e(k,j), j=1,3), zO(k)
     zO(k) = zO(k) - z_c
     
     ! convert etmp to frequency
     w(k) = c0 + c1*etmp + c2*etmp**2      
     wavg = wavg + w(k)                                       
     w2avg = w2avg + w(k)**2
   
     ! convert unit vector to transition dipole moment
     xtmp = d0 + d1*w(k)                
     mutmp = b0 + b1*etmp + b2*etmp**2
     
     mu(k,:) = e(k,:)*mutmp*xtmp

     alpha = xtmp*(a0 + a1*etmp) 

     apara = 3d0/(1d0 + 2d0/c15)*alpha
     aperp = 3d0/(2d0 + c15)*alpha
     axx = aperp + (apara - aperp)*e(k,1)*e(k,1)
     ayy = aperp + (apara - aperp)*e(k,2)*e(k,2)
     azz = aperp + (apara - aperp)*e(k,3)*e(k,3)
     axy = (aperp - apara)*e(k,1)*e(k,2)
     ayz = (aperp - apara)*e(k,2)*e(k,3)
     azx = (aperp - apara)*e(k,3)*e(k,1)
     a_ss(k) = xtmp*(axx + 2.0d0*axy + ayy)
     a_sp(k) = xtmp*(azx + ayz)
     a_pp(k) = xtmp*(azz)
  enddo

  close(ioh+100)

  End subroutine read_field


!   Subroutine to read in the field parameters from field.xxx
!     and calculate the derivatives of the averages

  Subroutine read_field_fluc(ioh,w,mu,e,a_ss,a_sp,a_pp,zO,dH)

  use time_data
  use map_data
  use freq_data
  use fluc_data

  implicit none
  integer :: ioh, j, k, n
  character*4 :: ext
  double precision :: apara,aperp
  double precision :: axx,ayy,azz,axy,ayz,azx
  double precision :: etmp, mutmp, xtmp,alpha
  double precision, dimension(ntimes) :: w,a_ss,a_sp,a_pp,zO
  double precision, dimension(ntimes,3) :: mu, e
  double precision, dimension(ntimes,4) :: dH

  ! Open the field field file for OH group iOH
  if(ndigits.eq.4) then
     write(ext,'(I0.4)') ioh
  elseif(ndigits.eq.3) then
     write(ext,'(I0.3)') ioh
  else
     write(6,*) ' Invalid ndigits, only 3 or 4 are allowed'
     stop
  endif

  ! Open the field file
  open(ioh+100,file='field_files/field.'//ext//'',status='old')

  do k = 1 , ntimes

     ! read in data from file (here e = e_OH)
     read(ioh+100,*) etmp, (e(k,j), j=1,3), zO(k)
     zO(k) = zO(k) - z_c
     
     ! convert etmp to frequency
     w(k) = c0 + c1*etmp + c2*etmp**2       
     wavg = wavg + w(k)                 
     w2avg = w2avg + w(k)**2

     ! calculate the fluctuation weighted averages
     do n = 1, 8
        dwavg(n) = dwavg(n) + dH(k,n)*w(k)
        dw2avg(n) = dw2avg(n) + dH(k,n)*w(k)**2
        d2wavg(n) = d2wavg(n) + (dH(k,n)**2 - dH2avg(n))*w(k)
        d2w2avg(n) = d2w2avg(n) + (dH(k,n)**2 - dH2avg(n))*w(k)**2
     enddo

     ! convert unit vector to transition dipole moment
     xtmp = d0 + d1*w(k)               
     mutmp = b0 + b1*etmp + b2*etmp**2
     
     mu(k,:) = e(k,:)*mutmp*xtmp
     
     alpha = xtmp*(a0 + a1*etmp)
     apara = 3d0/(1d0 + 2d0/c15)*alpha
     aperp = 3d0/(2d0 + c15)*alpha
     axx = aperp + (apara - aperp)*e(k,1)*e(k,1)
     ayy = aperp + (apara - aperp)*e(k,2)*e(k,2)
     azz = aperp + (apara - aperp)*e(k,3)*e(k,3)
     axy = (aperp - apara)*e(k,1)*e(k,2)
     ayz = (aperp - apara)*e(k,2)*e(k,3)
     azx = (aperp - apara)*e(k,3)*e(k,1)
     a_ss(k) = xtmp*(axx + 2.0d0*axy + ayy)
     a_sp(k) = xtmp*(azx + ayz)
     a_pp(k) = xtmp*(azz)
  enddo
  close(ioh+100)

  End subroutine read_field_fluc
