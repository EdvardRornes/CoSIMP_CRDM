      program DMCR_sigtot
*******************************************************************************
*** This program calculates the sigtot for a DM+p->DM+DM+p process          ***
***                                                                         ***
*** Author: Edvard Rørnes, 2025-01-07                                       ***
*******************************************************************************
      implicit none
      
      include 'dsidtag.h'    ! contains information about which 
                             ! particle module is used.
c   include 'dsddcom.h'    ! for access to Deff, rho0
      include 'dsmpconst.h' ! For m_p

ccc functions
      real*8 dsddDMCRsigtot
      
ccc output files
      character*80 outfile

      real*8 tstart,tfinish       ! aux
      integer ierr,iwarn
      real*8 sigtot(5)
      real*8 sigsi, sigvan
      real*8 mdm(5), mN, Tdm, Tdmmin, Tdmmax
      integer npoints, i, j

      call CPU_TIME(tstart)
      call dsinit

c... USER INPUT:
c**************************************
      outfile = 'data/DMCR_sigtot.dat'
      
c... fiducial scattering and annihilation cross sections;
c... result will not depend on these
      sigsi  = 1.d-22 ! cm^2
      sigvan = 3.d-26  ! cm^3 s^-1
c... define DM masses to be tabulated (in GeV)
      mdm(1) = 1.d-3
      mdm(2) = 1.d-2
      mdm(3) = 1.d-1
      mdm(4) = 1.d0
      mdm(5) = 1.d1
      mN = 122.05
      mN = m_p
c... range of kinetic energies to tabulate [GeV]
      Tdmmin = mdm(1)+3.d0*mdm(1)**2/(2.d0*mN)!+1.d-5
      Tdmmax = 1.d8
      Tdmmin = ((2.d0)*mdm(1)+mN)**2
      Tdmmin = (3.d0*mdm(1)+mN)**2
      npoints = 500

      open (unit=10,file=outfile,status='unknown')
      write(10,*) '# DM masses to be tabulated [GeV]'
      write(10,*) mdm
      write(10,*) 'Tkin [GeV]  | DM sigtot/Cross section [cm^-4 s^-1 GeV^-1] '//
     &            'for DM masses in above order '
      write(10,*)

      do i=1,npoints
        Tdm = Tdmmin*(Tdmmax/Tdmmin)**((i-1.)/(npoints-1.))
        do j=1,5
      !     if (Tdm.lt.mdm(j)+3.d0*mdm(j)**2/(2.d0*mN)) then
      !       sigtot(j) = 0.d0
      !       cycle
      !     endif
          call dsgivemodel_generic_wimp(mdm(j), sigvan, 5, sigsi*1.d36)
          call dsmodelsetup(ierr, iwarn)
          sigtot(j)=dsddDMCRsigtot(Tdm) ! NB: rholocal and Deff are set in dsdd_init!
          if (sigtot(j).ne.sigtot(j)) then
            write(*,*) 'WARNING: dsddDMCRsigtot returned NaN!'
            write(*,*) 'Setting result of the integral to zero.'
            sigtot(j) = 0.d0
          endif
        enddo 
        write(10,*) Tdm, sigtot
      enddo

      close(10)

c... Final output
      write(*,*)
      write(*,*) 'Outputs written to ', outfile
      write(*,*)

      call CPU_TIME(tfinish)
      write (*,*)
      write (*,*) '----------------------------------------------------------'
      write (*,*) 'The DarkSUSY DMCR_sigtot program has finished successfully.'
      write (*,'(A,F5.3,A)') ' Total time needed: ', tfinish - tstart, ' seconds'
      write (*,*) 'Particle module that was used: ', moduletag
      write (*,*) '----------------------------------------------------------'
      write (*,*)
      end program DMCR_sigtot