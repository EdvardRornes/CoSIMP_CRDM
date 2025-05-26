      program DMCR_target_recoil
*******************************************************************************
*** This program calculates the recoil energy of a particle N from an       ***
*** incident DM particle given the DM flux given by dsddDMCRflux for a 2->3 ***
*** process which has been upscattered by a CR.                             ***
***                                                                         ***
*** Author: Edvard Rørnes, 2025-01-07                                       ***
*******************************************************************************
      implicit none
      include 'dsidtag.h'    ! contains information about which 
                             ! particle module is used.
      include 'dsddcom.h'    ! for access to Deff, rho0
      include 'dsnuclides.h' ! also includes dsmpconst.h

ccc functions
      real*8 dsddDMCRdgammadt
      
ccc output files
      character*80 outfile

      real*8 tstart, tfinish, tcurrent    ! aux
      integer ierr, iwarn
      real*8 recoil(5)
      real*8 sigsi, sigvan
      real*8 mdm(5)
      real*8 mN, TN, TNmin, TNmax
      integer npoints, i, j, percent_complete

      call CPU_TIME(tstart)
      call dsinit

c... USER INPUT:
c**************************************
      outfile = 'data/DMCR_target_recoil_ff.dat'
      
c... Fiducial scattering and annihilation cross sections;
c... The results in DDCR_limits will not depend on these
      sigsi  = 1.d-24  ! cm^2 GeV^2
      sigvan = 3.d-26  ! cm^3 s^-1
c... Xenon mass
      mN = 122.05
      ! mN = m_p
c... Define DM masses to be calculated (in GeV)
      mdm(1) = 1.d-3
      mdm(2) = 1.d-2
      mdm(3) = 1.d-1
      mdm(4) = 1.d0
      mdm(5) = 1.d1
c... Range of kinetic energies to tabulate [GeV]
      TNmin = 1.d-6
      TNmax = 1.d1
      npoints = 500
c... This is the max range the program allows without encountering NaNs
c TNmin = 1.d-26
c TNmax = 1.d10

c... Write target recoil to second file
      open(unit=14, file=outfile, status='unknown')
      write(14,*) '# DM masses to be tabulated [GeV]'
      write(14,*) mdm
      write(14,*) 'Tkin [GeV]  | Target recoil [s^-1 GeV^-1] '//
     &            'for DM masses in above order '
      write(14,*)

      do i=1,npoints
        TN = TNmin*(TNmax/TNmin)**((i-1.)/(npoints-1.))
        do j = 1, 5
c... Factor of 1.d36 converts cm^2 to pb
          call dsgivemodel_generic_wimp(mdm(j), sigvan, 5, sigsi*1.d36)
          call dsmodelsetup(ierr, iwarn)
          recoil(j) = dsddDMCRdgammadt(TN, mN, 0.d0)
        enddo
      !   write(*,*) TN, recoil
        write(14,*) TN, recoil
c... Calculate and print progress
        if (mod(i, max(1, npoints/100)).eq.0.or.i.eq.npoints) then
          call CPU_TIME(tcurrent)
          percent_complete = int(100.d0*i/npoints)
          if (tcurrent-tstart.lt.1.d1) then
            write(*, '(I6, A, F6.4, A)') percent_complete, 
     &      '%  time: ', tcurrent - tstart, ' seconds'
          elseif (tcurrent-tstart.gt.1.d1.and.tcurrent-tstart.lt.1.d2) then
            write(*, '(I6, A, F6.3, A)') percent_complete, 
     &      '%  time: ', tcurrent - tstart, ' seconds'
          elseif (tcurrent-tstart.gt.1.d2.and.tcurrent-tstart.lt.1.d3) then
            write(*, '(I6, A, F6.2, A)') percent_complete, 
     &      '%  time: ', tcurrent - tstart, ' seconds'
          else
            write(*, '(I6, A, F6.1, A)') percent_complete, 
     &      '%  time: ', tcurrent - tstart, ' seconds'
          endif
        endif
      enddo
      close(14)

c... Final output
      write(*,*)
      write(*,*) 'Outputs written to ', outfile

      call CPU_TIME(tfinish)
      write (*,*)
      write (*,*) '----------------------------------------------------'//
     &            '--------------'
      write (*,*) 'The DarkSUSY DMCR_target_recoil program has finished'//
     &            ' successfully.'
      write (*,'(A,F7.3,A)') ' Total time needed: ', tfinish-tstart, ' sec'
      write (*,*) 'Particle module that was used: ', moduletag
      write (*,*) '----------------------------------------------------'//
     &            '--------------'
      write (*,*)
      end program DMCR_target_recoil