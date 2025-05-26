      program DMCR_EnergyLoss
*******************************************************************************
*** This program calculates the EnergyLoss of DM particles that results     ***
*** from cosmic ray upscattering for a 2->3 process where the incident CR   ***
*** is an oxygen nucleus.                                                   ***
***                                                                         ***
*** Author: Edvard Rørnes, 2025-01-07                                       ***
*******************************************************************************
      use omp_lib  ! Import OpenMP functions
      implicit none
      
      include 'dsidtag.h'    ! contains information about which particle module is used.
      include 'dsddcom.h'    ! for access to Deff, rho0
      include 'dsnuclides.h' ! also includes dsmpconst.h

ccc functions
      real*8 dsddDMCREnergyLoss
      
ccc output files
      character*200 outfile

      real*8 tstart, tfinish, tcurrent
      integer ierr,iwarn
      real*8 EnergyLoss(5)
      real*8 sigsi, sigvan, sigtot(5)
      real*8 mdm(5), s, Tdm, Tdmmin, Tdmmax, mN, x
      integer npoints, i, j, percent_complete

      tstart = OMP_GET_WTIME()
      call dsinit

c... USER INPUT:
c**************************************
      outfile = 'data/DMCR_EnergyLoss'

      if (CRDM_form_factor) then
        outfile = trim(outfile) // '_FF'
      endif
      if (CRDM_cs) then
        outfile = trim(outfile) // '_cs'
      endif
      outfile = trim(outfile) // '_O.dat'
      
c... fiducial scattering and annihilation cross sections;
c... The result in DDCR_limits will not depend on these
      sigsi  = 1.d-10 ! cm^2
      sigvan = 3.d-26  ! cm^3 s^-1
c... define DM masses to be tabulated [GeV]
      mdm(1) = 1.d-3
      mdm(2) = 1.d-2
      mdm(3) = 1.d-1
      mdm(4) = 1.d0
      mdm(5) = 1.d1
c... range of kinetic energies to tabulate [GeV]
      mN = mNCRau(1)*atomicmassunit ! convert masses to GeV

c... This should probably be done automatically in the future if used further...
c... Is only done to sample more accurately close to threshold energies.

c... Sample low energy for mdm(1)
      Tdmmin = 3.d0*mdm(1)**2/(2.d0*mN)+mdm(1)+1.d-6
      ! Tdmmax = 3.d0*mdm(2)**2/(2.d0*mN)+mdm(2)
c... Sample low energy for mdm(2)
      ! Tdmmin = 3.d0*mdm(2)**2/(2.d0*mN)+mdm(2)+1.d-5
      ! Tdmmax = 3.d0*mdm(3)**2/(2.d0*mN)+mdm(3)
c... Sample low energy for mdm(3)
      ! Tdmmin = 3.d0*mdm(3)**2/(2.d0*mN)+mdm(3)+1.d-4
      ! Tdmmax = 3.d0*mdm(4)**2/(2.d0*mN)+mdm(4)
c... Sample low energy for mdm(4)
      ! Tdmmin = 3.d0*mdm(4)**2/(2.d0*mN)+mdm(4)+1.d-3
      ! Tdmmax = 3.d0*mdm(5)**2/(2.d0*mN)+mdm(5)
c... Sample low energy for mdm(5)
      ! Tdmmin = 3.d0*mdm(5)**2/(2.d0*mN)+mdm(5)+1.d-2
      ! Tdmmax = 3.d0*1.d2**2/(2.d0*mN)+1.d2
c... Finish the range
      ! Tdmmin = 3.d0*1.d2**2/(2.d0*mN)+1.d2+1.d-1
      ! Tdmmax = 1.d4
      ! Tdmmin = 1.d4
      ! Tdmmin = 1.d4

      Tdmmax = 1.d5
      npoints = 20
c... To set EnergyLoss = NaN when s<smin
      x = 0

      write(*,*) 'Calculating average energy loss and writing results to file ', outfile
      if (CRDM_form_factor) then
        write(*,*)
        write(*,*) 'With form factors this can take a while...'
      endif
      write(*,*)
      open (unit=12,file=outfile,status='unknown')
      write(12,*) '# DM masses [GeV]'
      write(12,*) mdm
      write(12,*) 'Incoming Tdm [GeV] | Energy Loss [GeV] '//
     &            'for DM masses in above order '
      write(12,*)

      do i=1,npoints
        Tdm = Tdmmin*(Tdmmax/Tdmmin)**((i-1.)/(npoints-1.))
        do j=1,5
          EnergyLoss(j) = 0.d0
          if (Tdm.le.(3.d0*mdm(j)**2/(2.d0*mN)+mdm(j))) then
            EnergyLoss(j) = 0/x
            if (CRDM_form_factor) then
              write(*,*) 'Too low energy for mdm =', mdm(j), 'Setting EnergyLoss = NaN'
            endif
            cycle
          endif
          s = (mN+mdm(j))**2 + 2.d0*mN*Tdm
          if (CRDM_form_factor) then
            write(*,*) 'CM energy =', sqrt(s), 'CF Energy =', (Tdm+mdm(j)+mN), 'min CM energy =', (mN + 2.d0*mdm(j))
            write(*,*) 'Incoming Tdm =', Tdm, 'mdm =', mdm(j)
          endif

          call dsgivemodel_generic_wimp(mdm(j), sigvan, 5, sigsi*1.d36)
          call dsmodelsetup(ierr, iwarn)
c... Calculate the energy loss times cross section
          CRDM_EnergyLoss = .true.
          EnergyLoss(j) = dsddDMCREnergyLoss(s, mN)
          if (CRDM_form_factor) then
            write(*,*) 'ω*σ =', EnergyLoss(j)
          endif

c... Recalculate the cross section without the energy loss
          if (EnergyLoss(j).ne.0.d0) then
            CRDM_EnergyLoss = .false.
            sigtot(j) = dsddDMCREnergyLoss(s, mN)
c... Divide by the cross section to get the energy loss
            EnergyLoss(j) = EnergyLoss(j)/sigtot(j)
            if (CRDM_form_factor) then
            write(*,*) 'ω   =', EnergyLoss(j), 'ω/Tdm =', EnergyLoss(j)/Tdm
          endif
          endif
        enddo
      write(12,*) Tdm, EnergyLoss, sigtot
c... Calculate and print progress
        if (mod(i, max(1, npoints/100)).eq.0.or.i.eq.npoints) then
          tcurrent = OMP_GET_WTIME()
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
          elseif (tcurrent-tstart.gt.1.d3.and.tcurrent-tstart.lt.1.d4) then
            write(*, '(I6, A, F6.1, A)') percent_complete,
     &      '%  time: ', tcurrent - tstart, ' seconds'
          else
            write(*, '(I6, A, F9.1, A)') percent_complete,
     &      '%  time: ', tcurrent - tstart, ' seconds'
          endif
        endif
      enddo

      close(12)

c... Final output
      write(*,*)
      write(*,*) 'Outputs written to ', outfile
      write(*,*)

      tfinish = OMP_GET_WTIME()
      write (*,*)
      write (*,*) '---------------------------------------------------------'
      write (*,*) 'The DarkSUSY DMCR_EnergyLoss program has finished successfully.'
      write (*,'(A,F10.3,A)') ' Total time needed: ', tfinish - tstart, ' sec'
      write (*,*) 'Particle module that was used: ', moduletag
      write (*,*) '---------------------------------------------------------'
      write (*,*)
      end program DMCR_EnergyLoss