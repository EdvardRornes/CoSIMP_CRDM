      program DMCR_EnergyLoss
*******************************************************************************
*** This program calculates the average EnergyLoss of DM particles in a     ***  
*** 2->2 elastic scattering process with mN as the target nucleon.          *** 
***                                                                         ***
*** Author: Edvard Rørnes, 2025-01-07                                       ***
*******************************************************************************
      use omp_lib  ! Import OpenMP functions
      implicit none
      
      include 'dsidtag.h'    ! contains information about which particle module is used.
      include 'dsddcom.h'    ! for access to Deff, rho0
      include 'dsnuclides.h' ! also includes dsmpconst.h

ccc functions
      real*8 dsddDMCREnergyLoss2t2
      
ccc output files
      character*80 outfile, ini_name

      real*8 tstart,tfinish,tcurrent                    ! aux
      integer ierr,iwarn
      real*8 EnergyLoss(5)
      real*8 sigsi, sigvan, sigtot
      real*8 mdm(5), s, Tdm, Tdmmin, Tdmmax, mcr, x
      integer npoints, i, j, percent_complete

      tstart = OMP_GET_WTIME()
      call dsinit

c... USER INPUT:
c**************************************
      ini_name = 'data/DMCR_EnergyLoss'
      outfile = trim(ini_name)  ! Remove trailing spaces before assigning

      if (CRDM_form_factor) then
        outfile = trim(outfile) // '_FF'
      endif
      if (CRDM_cs) then
        outfile = trim(outfile) // '_cs'
      endif
      outfile = trim(outfile) // '_O_2t2.dat'
      
c... fiducial scattering and annihilation cross sections;
c... The result in DDCR_limits will not depend on these
      sigsi  = 1.d0 ! cm^2
      sigvan = 3.d-26  ! cm^3 s^-1
c... define DM masses to be tabulated [GeV]
      mdm(1) = 1.d-3
      mdm(2) = 1.d-2
      mdm(3) = 1.d-1
      mdm(4) = 1.d0
      mdm(5) = 1.d1
c... range of kinetic energies to tabulate [GeV]
      mcr = mNCRau(4)*atomicmassunit ! convert masses to GeV
      Tdmmin = 1.d-5
      Tdmmax = 1.d5
      npoints = 200
c... Just to set EnergyLoss = NaN when s<smin
      x = 0

      write(*,*) 'Calculating flux and writing results to file ', outfile
      if (CRDM_form_factor) then
        write(*,*)
        write(*,*) 'With form factors this can take a while...'
      endif
      write(*,*)
      open (unit=10,file=outfile,status='unknown')
      write(10,*) '# DM masses [GeV]'
      write(10,*) mdm
      write(10,*) 'Tkin [GeV]  | DM flux [cm^-2 s^-1 GeV^-1] '//
     &            'for DM masses in above order '
      write(10,*)

      do i=1,npoints
        Tdm = Tdmmin*(Tdmmax/Tdmmin)**((i-1.)/(npoints-1.))
        do j=1,5
          EnergyLoss(j) = 0.d0
          s = (mcr+mdm(j))**2 + 2.d0*mcr*Tdm
          ! write(*,*) 's =', s, 'CF Energy squared =', (Tdm+mdm(j)+mcr)**2
          ! write(*,*) 'Incoming Tdm =', Tdm, 'mdm = ', mdm(j)

          call dsgivemodel_generic_wimp(mdm(j), sigvan, 5, sigsi*1.d36)
          call dsmodelsetup(ierr, iwarn)
c... Calculate the energy loss times cross section
          CRDM_EnergyLoss = .true.
          EnergyLoss(j) = dsddDMCREnergyLoss2t2(s, mcr)
          ! write(*,*) 'ω*σ =', EnergyLoss(j)

c... Recalculate the cross section without the energy loss
          if (EnergyLoss(j).ne.0.d0) then
            CRDM_EnergyLoss = .false.
            sigtot = dsddDMCREnergyLoss2t2(s, mcr)
c... Divide by the cross section to get the energy loss
            EnergyLoss(j) = EnergyLoss(j)/sigtot
            ! write(*,*) 'ω   =', EnergyLoss(j)
          endif

          if (EnergyLoss(j).ne.EnergyLoss(j)) then
            write(*,*) 'WARNING: dsddDMCREnergyLoss returned NaN!'
            write(*,*) 'Setting EnergyLoss to zero.'
            EnergyLoss(j) = 0.d0
          endif
          if (EnergyLoss(j).lt.0.d0) then
            write(*,*) 'WARNING: dsddDMCREnergyLoss returned a negative value!'
            write(*,*) 'Setting EnergyLoss to zero.'
            EnergyLoss(j) = 0.d0
          endif
          if (EnergyLoss(j).gt.(1.d300)) then
            write(*,*) 'WARNING: dsddDMCREnergyLoss returned Infinity!'
            write(*,*) 'Setting EnergyLoss to zero.'
            EnergyLoss(j) = 0.d0
          endif
        enddo
      write(10,*) Tdm, EnergyLoss
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

      close(10)
      
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