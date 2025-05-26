      program DMCR_flux
*******************************************************************************
*** This program calculates the flux of DM particles that results from      ***
*** cosmic ray upscattering for a 2->3 process where the incident CR is a   ***
*** proton.                                                                 ***
***                                                                         ***
*** Author: Edvard Rørnes, 2025-01-07                                       ***
*******************************************************************************
      use omp_lib  ! Import OpenMP functions
      implicit none
      
      include 'dsidtag.h'    ! contains information about which particle module is used.
      include 'dsddcom.h'    ! for access to Deff, rho0
      include 'dsnuclides.h' ! also includes dsmpconst.h

ccc functions
      real*8 dsddDMCRflux, dsddDMCRflux_Tab, dsddDMCRflux_2loop
      
ccc output files
      character*80 outfile, ini_name

      real*8 tstart,tfinish,tcurrent                    ! aux
      integer ierr,iwarn
      real*8 flux(5)
      real*8 sigsi, sigvan
      real*8 mdm(5), Tdm, Tdmmin, Tdmmax, mcr, mdmmax, mdmmin
      integer npoints, i, j, percent_complete

      tstart = OMP_GET_WTIME()
      call dsinit

c... USER INPUT:
c**************************************
      ini_name = 'data/DMCR_flux'
      outfile = trim(ini_name)  ! Remove trailing spaces before assigning

      if (CRDM_form_factor) then
        outfile = trim(outfile) // '_FF'
      endif
      if (CRDM_cs) then
        outfile = trim(outfile) // '_cs'
      endif
      outfile = trim(outfile) // '.dat'
      
c... fiducial scattering and annihilation cross sections;
c... The result in DDCR_limits will not depend on these
      sigsi  = 1.d-20 ! cm^2
      sigvan = 3.d-26  ! cm^3 s^-1
c... define DM masses to be tabulated [GeV]
      mdm(1) = 1.d-3
      mdm(2) = 1.d-2
      mdm(3) = 1.d-1
      mdm(4) = 1.d0
      mdm(5) = 1.d1
c... range of kinetic energies to tabulate [GeV]
c... This is to tabulate the flux. We never need energies below this for Co-SIMP
c... interactions in the detector.
      ! mcr = mNCRau(3)*atomicmassunit ! convert masses to GeV
      ! mcr = mNCRau(4)*atomicmassunit ! convert masses to GeV
      mdmmin = 1.d-4
      mdmmax = 1.d1
      ! Tdmmin = 3.d0*mdm(1)**2/(2.d0*mcr)+mdm(1)+1.d-6
      Tdmmin = 1.d-6
      Tdmmax = 1.d6
      npoints = 100

      ! do j=1,100
      !   mdm(j) = mdmmin*(mdmmax/mdmmin)**((j-1.d0)/(npoints-1.d0))
      ! enddo

      write(*,*) 'Calculating flux and writing results to file ', outfile
      if (CRDM_form_factor) then
        write(*,*)
        write(*,*) 'With form factors this can take a while...'
      endif
      write(*,*)
      open (unit=13,file=outfile,status='unknown')
      write(13,*) '# DM masses [GeV]'
      write(13,*) mdm
      write(13,*) 'Tkin [GeV]  | DM flux [cm^-2 s^-1 GeV^-1] '//
     &            'for DM masses in above order '
      write(13,*)

      do i=1,npoints
        Tdm = Tdmmin*(Tdmmax/Tdmmin)**((i-1.d0)/(npoints-1.d0))
        do j=1,5
          call dsgivemodel_generic_wimp(mdm(j), sigvan, 5, sigsi*1.d36)
          call dsmodelsetup(ierr, iwarn)
          ! if (CRDM_form_factor) then
          !   write(*,*) 'Calculating flux for mdm =', mdm(j), 'with Tdm =', Tdm
          ! endif
          flux(j) = dsddDMCRflux(Tdm, rholocal*Deff)
          ! flux(j) = dsddDMCRflux_Tab(Tdm, rholocal*Deff)
          ! flux(j) = dsddDMCRflux_2loop(Tdm, rholocal*Deff)

          if (flux(j).ne.flux(j)) then
            write(*,*) 'WARNING: dsddDMCRflux returned NaN!'
            write(*,*) 'Setting flux to zero.'
            flux(j) = 0.d0
          endif
          if (flux(j).lt.0.d0) then
            write(*,*) 'WARNING: dsddDMCRflux returned a negative value!'
            write(*,*) 'Setting flux to zero.'
            flux(j) = 0.d0
          endif
          if (flux(j).gt.(1.d300)) then
            write(*,*) 'WARNING: dsddDMCRflux returned Infinity!'
            write(*,*) 'Setting flux to zero.'
            flux(j) = 0.d0
          endif
          ! if (CRDM_form_factor) then
          !   write(*,*) 'flux for mdm', mdm(j), ' =', flux(j), 'Tdm =', Tdm
          ! endif
        enddo
      ! write(*,*) Tdm, flux
      write(13,*) Tdm, flux
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

      close(13)

c... Final output
      write(*,*)
      write(*,*) 'Outputs written to ', outfile
      write(*,*) 'CRtypes used are'
      do i=1,4
        write(*,*) 'CRtype', i, ' = ', CR_inc(i)
      enddo
      write(*,*)

      tfinish = OMP_GET_WTIME()
      write (*,*)
      write (*,*) '---------------------------------------------------------'
      write (*,*) 'The DarkSUSY DMCR_flux program has finished successfully.'
      write (*,'(A,F10.3,A)') ' Total time needed: ', tfinish - tstart, ' sec'
      write (*,*) 'Particle module that was used: ', moduletag
      write (*,*) '---------------------------------------------------------'
      write (*,*)
      end program DMCR_flux