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

ccc functions
      real*8 dscrISRflux
      
ccc output files
      character*80 outfile, ini_name

      real*8 tstart,tfinish,tcurrent
      real*8 flux(4)
      real*8 Tcr, Tcrmin, Tcrmax
      integer npoints, i, j, percent_complete

      tstart = OMP_GET_WTIME()
      call dsinit

c... USER INPUT:
c**************************************
      ini_name = 'data/DMCR_CRflux.dat'
      outfile = trim(ini_name)
      
c... fiducial scattering and annihilation cross sections;
c... The result in DDCR_limits will not depend on these
c... define DM masses to be tabulated [GeV]
c... range of kinetic energies to tabulate [GeV]
      Tcrmin = 1.d-6
      Tcrmax = 1.d6
      npoints = 1000

      write(*,*) 'Calculating flux and writing results to file ', outfile
      if (CRDM_form_factor) then
        write(*,*)
        write(*,*) 'With form factors this can take a while...'
      endif
      write(*,*)
      open (unit=10,file=outfile,status='unknown')
      write(10,*) 'Tkin [GeV]  | CR flux [cm^-2 s^-1 GeV^-1] '//
     &            'for CR in order p, He, C, O '
      write(10,*)

      do i=1,npoints
        Tcr = Tcrmin*(Tcrmax/Tcrmin)**((i-1.)/(npoints-1.))
        do j=1,4
          flux(j) = dscrISRflux(Tcr, j)
        enddo
      !$OMP END PARALLEL DO
      ! write(*,*) Tcr, flux
      write(10,*) Tcr, flux
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
      write (*,*) '-----------------------------------------------------------'
      write (*,*) 'The DarkSUSY DMCR_CRflux program has finished successfully.'
      write (*,'(A,F10.3,A)') ' Total time needed: ', tfinish - tstart, ' sec'
      write (*,*) 'Particle module that was used: ', moduletag
      write (*,*) '-----------------------------------------------------------'
      write (*,*)
      end program DMCR_flux