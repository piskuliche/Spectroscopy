!   Subroutine to read in the energies at each step to do fluc theory calculations

  Subroutine read_eners(dH)

  use time_data
  use fluc_data

  implicit none
  integer :: i, k
  double precision :: rspc, ewald
  double precision, dimension(ntimes,5) :: dH


  ! Open the field file
  open(30,file='ener_tot.dat',status='old')

  Havg = 0.0d0
  do i = 1 , ntimes

     ! read in the total energies and calculate the average
     read(30,*) dH(i,1), dH(i,2), dH(i,3), rspc, ewald
     dH(i,4) = rspc + ewald
     dH(i,5) = dH(i,3) + dH(i,4)

     Havg(:) = Havg(:) + dH(i,:)

  enddo
  Havg = Havg/dble(ntimes)
  write(6,*) ' <H>      = ',Havg(1),' (kcal/mol)'
  write(6,*) ' <KE>     = ',Havg(2),' (kcal/mol)'
  write(6,*) ' <LJ>     = ',Havg(3),' (kcal/mol)'
  write(6,*) ' <Coul>   = ',Havg(4),' (kcal/mol)'
  write(6,*) ' <V>      = ',Havg(5),' (kcal/mol)'

!!! For simulations with NVE trajectories, we need to compute the average later
  ! Subtract the averages to determine the fluctuations
!  do k = 1, 5
!     dH(:,k) = dH(:,k) - Havg(k)
!  enddo
!!! 

  do k = 1, 5
     dH2avg(k) = 0d0
     do i = 1, ntimes
        dH2avg(k) = dH2avg(k) + (dH(i,k)- Havg(k))**2
     enddo
  enddo
  dH2avg = dH2avg/dble(ntimes)
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
