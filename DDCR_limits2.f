      program DDCR_limits2
*******************************************************************************
*** This program calculates the upper and lower limits for the Co-SIMP CRDM ***
*** scenario when attenuation is neglected.                                 ***
*******************************************************************************
      use omp_lib  ! Import OpenMP functions
      implicit none
      include 'dsddcom.h'   ! needed for direct access to variables stored in
                            ! direct detection-related common blocks
                            ! Doesnt seem to do anything currently?
      include 'dsidtag.h'   ! contains information about which particle module is used.
      include 'dsnuclides.h' ! also includes dsmpconst.h

      real*8 tstart, tcurrent, tfinish, count
      real*8 mdm, sigstart, mmin, mmax, sigtot, target_sig, sigsi, sigvan, s, m_O
      character*80 outfile     ! output file
      character*20 addon, experiment
      logical upper
      integer ierr, iwarn, i, npoints, percent_complete

ccc functions
      real*8 dsddDMCRcountrate, dsddDMCRsigtot, dsddDMCREnergyLoss

      call dsinit
      tstart = OMP_GET_WTIME()


********************
***  USER INPUT  ***
********************
      upper = .false. ! If upper=.true. we calculate an estimated upper limit due
                     ! to attenuation. Otherwise we calculate the lower limit.
      unitarity = .true.
c... choose which experiment (and determine output filename correspondingly)
      experiment = 'MiniBoone' ! see src/dd/dsddDMCRcountrate.f for options
      ! experiment = 'Xenon1t'   ! see src/dd/dsddDMCRcountrate.f for options
      ! experiment = 'Borexino'   ! see src/dd/dsddDMCRcountrate.f for options
      ! experiment = 'Darwin'   ! see src/dd/dsddDMCRcountrate.f for options
      experiment = 'LunarExample'
      experiment = 'BGO'

c... limit type
      sigstart = 1.0d-26 ! start value to search for limit; this should not be too
                         ! far off the true limit for the bisection to work
      

c...model settings
      addon = '_CoSIMP'
*************************
***  END USER INPUT   ***
*************************

      write (*,*) '-------------------------------------------------------'
      write (*,*) 
c... determine output filename based on user input above and settings in dsddcrdm_init.f
      outfile = 'data/DMCR'
      if (upper) then
        outfile = trim(outfile) // '_upper'
      else
        outfile = trim(outfile) // '_lower'
      endif
      outfile = trim(outfile) // '_limits'
      if (CRDM_cs) then
        outfile = trim(outfile) // '_cs'
      endif
      if (CRDM_form_factor) then
        outfile = trim(outfile) // '_FF'
      endif

c... Set desired max cross section to approximate attenuation
      if (experiment.eq.'MiniBoone') then 
        target_sig = 1.d-26
        outfile = trim(outfile) // '_MiniBoone'
      endif
      if (experiment.eq.'Xenon1t') then
        target_sig = 1.d-28
        outfile = trim(outfile) // '_Xenon1t'
      endif
      if (experiment.eq.'Borexino') then 
        target_sig = 1.d-26
        outfile = trim(outfile) // '_Borexino'
      endif
      if (experiment.eq.'Darwin') then
        target_sig = 1.d-28
        outfile = trim(outfile) // '_Darwin'
      endif

      outfile = trim(outfile) // trim(addon) // '.dat'        

      write(*,*) 'Calculating upper limits and writing results to file ', outfile
      write(*,*)
      open (unit=15,file=outfile,status='unknown')

c... scan setup
      mmin     = 1.d-4  ! minimal and maximal DM masses (GeV), mmin<1d-6 causes numerical errors
      mmax     = 1.d1   ! With form factors the program is very slow for mdm > 10 GeV
      npoints = 100

      do i=1,npoints
        mdm = mmin*(mmax/mmin)**((i-1.)/(npoints-1.))

        call dsgivemodel_generic_wimp(mdm, sigvan, 5, sigstart*1.d36)
        call dsmodelsetup(ierr, iwarn)
        if (upper) then
          s = (3.d0*mdm+m_p)**2
          s = (2.d0*mdm+2.d0*m_p)**2
c... For the upper limit we put a threshold on the cross section at s = smin + mdm
          sigtot = dsddDMCREnergyLoss(s, m_p)

          if (sigtot.ne.sigtot .or. sigtot.gt.1.d300) then
            write(*,*) 'WARNING: dsddDMCREnergyLoss returned NaN!'
            write(*,*) 'Setting result of the integral to zero.'
            sigsi = 0.d0
          elseif (sigtot.le.0.d0) then
            write(*,*) 'WARNING: dsddDMCREnergyLoss <= 0!'
            write(*,*) 'Setting result of the integral to zero.'
            sigsi = 0.d0
          else
c... If no errors, acquire normalization which gives sigtot=target_sig
            sigsi = target_sig*sigstart/sigtot
          endif

c... To verify that it was done correctly
          call dsgivemodel_generic_wimp(mdm, sigvan, 5, sigsi*1.d36)
          call dsmodelsetup(ierr, iwarn)
          sigtot = dsddDMCREnergyLoss(s, m_p)
          write(15,*) mdm, sigsi, sigtot

c... otherwise, we run the full program
        else
          s = (3.d0*mdm+m_p)**2
          ! s = (2.d0*mdm+2.d0*m_p)**2
          count = dsddDMCRcountrate(experiment)
          sigsi = sigstart/sqrt(count)  ! This gives count=1 when excluding attenuation
          if (count.eq.0.d0) cycle

c... new setup for new sigsi
          call dsgivemodel_generic_wimp(mdm, sigvan, 5, sigsi*1.d36)
          call dsmodelsetup(ierr, iwarn)
c... calculate the lower limit for the cross section at sqrt(s) = sqrt(smin) + mdm
          sigtot = dsddDMCRsigtot(s)
          ! count = dsddDMCRcountrate(experiment)
          write(*,*) '2', mdm, sigsi, sigtot!, count
          write(15,*) mdm, sigsi, sigtot
        endif

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

      close(15)

      tfinish = OMP_GET_WTIME()

      write (*,*)
      write (*,*) '-----------------------------------------------------------'
      write (*,*) 'The DarkSUSY DDCR_limits program has finished successfully.'
      write(*,*)  'Results written to file ', outfile
      write (*,'(A, F9.1)') ' Total time needed (in seconds): ',tfinish-tstart
      write (*,*) 'Particle module that was used:  ', moduletag
      write (*,*) '-----------------------------------------------------------'
      write (*,*)
      end

      
