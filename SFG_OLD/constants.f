!       Fondamental Physical Constants

!       Latest update: Jan 03 2001 from NIST website

!    Energy
      double precision  jperau ! = Hartree
      parameter         (jperau=4.35974381d-18)
      double precision  evsi
      parameter         (evsi=1.602176462d-19)
      double precision  jpercal
      parameter         (jpercal=4.1868d0)
!    Planck
      double precision  hbar
      parameter         (hbar=1.054571596d-34)
!    Length
      double precision  mperau ! = Bohr
      parameter         (mperau=0.5291772083d-10)
!    Pressure
      double precision  Paperatm
      parameter         (Paperatm=101325d0)
!    Boltzmann
      double precision  kbsi
      parameter         (kbsi=1.3806503d-23)
!    Mass
      double precision  kgperau
      parameter         (kgperau=9.10938188d-31)
!    Electric Dipole Moment
!     source = www.unc.edu/~rowlett/units/dictD.html
      double precision  siperdebye 
      parameter         (siperdebye=3.33564d-30)
!    Misc
      double precision  avogadro
      parameter         (avogadro=6.02214199d23)
      double complex    aye
      parameter         (aye=(0d0,1d0))
      double precision  lightvel
      parameter         (lightvel=299792458d0)

      double precision  pi
      parameter         (pi=3.14159265359d0)

!    Resulting conversion factors

!    Energy
      double precision  evperau
      parameter         (evperau=jperau/evsi)
      double precision  kjperau
      parameter         (kjperau=jperau*1d-3*avogadro)
      double precision   kjperkcal,kcalperev,kcalperau
      parameter         (kjperkcal=jpercal)
      parameter         (kcalperev=evsi/jpercal*avogadro*1d-3)
      parameter         (kcalperau=kjperau/jpercal)
!    Length
      !double precision  angperau,cmperau,dmperau
      double precision  cmperau,dmperau
      !parameter         (angperau=mperau*1d10)
      parameter         (cmperau=mperau*1d2,dmperau=mperau*1d1)
!    Pressure
      double precision  atmperau
      parameter         (atmperau=jperau/Paperatm/mperau**3)
!    Boltzmann
      double precision  kb
      parameter         (kb=kbsi/jperau)
!    Electric Dipole Moment
      double precision  dipolesiperau
      parameter         (dipolesiperau=evsi*mperau)
      double precision  debyeperau
      parameter         (debyeperau=evsi*mperau/siperdebye)
!    Mass
      double precision  gmolperau
      parameter         (gmolperau=kgperau*1d3*avogadro)
!    Time
      !double precision  sperau,fsperau,psperau
      double precision  sperau,psperau
      parameter         (sperau=hbar/jperau)
      !parameter         (fsperau=sperau*1d15)
      parameter         (psperau=sperau*1d12)
!    Wavenumbers (cm-1) per a.u.
      double precision  wnperau
      parameter         (wnperau=1d-2/(sperau*2d0*pi*lightvel))

!    Angle
      double precision  degperrad
      parameter         (degperrad = 57.2957795131d0)

!    Force
      double precision   forcedlperau
      parameter          (forcedlperau = 2.01553132145d-6)


      integer nit2max
      parameter (nit2max=10)
