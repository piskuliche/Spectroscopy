!   Subroutine to read in the energies at each step to do fluc theory calculations

  Subroutine read_eners(dH)

  use time_data
  use fluc_data

  implicit none
  integer :: i, k
  double precision :: rspc, ewald, count
  double precision, dimension(ntimes,8) :: dH

  ! Open the field file
  open(30,file='ener_tot.dat',status='old')
  open(41,file='ener_avg.dat')

  call flush(6)
  Havg = 0.0d0; count = 0d0
  do i = 1 , ntimes - ncorr, nskip
     count = count + 1d0

     ! read in the total energies and calculate the average
     ! Here, dH(:,6), dH(:,7), dH(:,8) are top, bottom, and top-bottom energies
     !read(30,*) dH(i,1), dH(i,2), dH(i,3), rspc, ewald, dH(i,6), dH(i,7), dH(i,8)
     read(30,*) dH(i,1), dH(i,2), dH(i,3), rspc, ewald
     dH(i,4) = rspc + ewald
     dH(i,5) = dH(i,3) + dH(i,4)

     Havg(:) = Havg(:) + dH(i,:)

  enddo
  Havg = Havg/count
  write(6,*) ' <H>      = ',Havg(1),' (kcal/mol)'
  write(6,*) ' <KE>     = ',Havg(2),' (kcal/mol)'
  write(6,*) ' <LJ>     = ',Havg(3),' (kcal/mol)'
  write(6,*) ' <Coul>   = ',Havg(4),' (kcal/mol)'
  write(6,*) ' <V>      = ',Havg(5),' (kcal/mol)'
!  write(6,*) ' <H_top>  = ',Havg(6),' (kcal/mol)'
!  write(6,*) ' <H_bot>  = ',Havg(7),' (kcal/mol)'
!  write(6,*) ' <H_t-b>  = ',Havg(8),' (kcal/mol)'

  write(41,'(10F20.7)') Havg(1), Havg(2), Havg(3), Havg(4), Havg(5)
  close(41)

  ! Subtract the averages to determine the fluctuations
!  do k = 1, 8
!     dH(:,k) = dH(:,k) - Havg(k)
!  enddo

!  do k = 1, 8
  do k = 1, 5
     dH2avg(k) = 0d0
     do i = 1, ntimes
        dH2avg(k) = dH2avg(k) + dH(i,k)**2
     enddo
  enddo
  dH2avg = dH2avg/dfloat(ntimes)
  write(6,*) ' <H^2>     = ',dH2avg(1)+Havg(1)**2,' (kcal/mol)^2'
  write(6,*) ' <KE^2>    = ',dH2avg(2)+Havg(2)**2,' (kcal/mol)^2'
  write(6,*) ' <LJ^2>    = ',dH2avg(3)+Havg(3)**2,' (kcal/mol)^2'
  write(6,*) ' <Coul^2>  = ',dH2avg(4)+Havg(4)**2,' (kcal/mol)^2'
  write(6,*) ' <V^2>     = ',dH2avg(5)+Havg(5)**2,' (kcal/mol)^2'

  write(6,*) ' <dH^2>    = ',dH2avg(1),' (kcal/mol)^2'
  write(6,*) ' <dKE^2>   = ',dH2avg(2),' (kcal/mol)^2'
  write(6,*) ' <dLJ^2>   = ',dH2avg(3),' (kcal/mol)^2'
  write(6,*) ' <dCoul^2> = ',dH2avg(4),' (kcal/mol)^2'
  write(6,*) ' <dV^2>    = ',dH2avg(5),' (kcal/mol)^2'

  close(30)

  End subroutine read_eners
