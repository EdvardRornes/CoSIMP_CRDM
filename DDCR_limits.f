      program DDCR_limits
c
c     This program is an example DarkSUSY main program to calculate
c     various limits on the DM-nucleon scattering cross section that
c     result from the cosmic-ray induced high-velocity flux of DM
c     particles.
c
c     It also demonstrates how to use the same main program for
c     different, user-provided scattering cross sections. If compiled
c     with 'make DDCR_limits', this program directly links to the
c     genericWIMP module and thus a constant scattering cross section
c     [producing the same limits as in 1810.10543].
c     When compiled with 'make DDCR_limits_q2', it uses instead the
c     relativistic scattering cross sections provided in
c     user_replaceables/dsddsigmarel.f (The examples in that file
c     implement the specific examples for non-trivial q^2-dependences
c     discussed in 2209.03360. Note that a q^2-dependent scattering
c     rate requires a re-interpretation of published direct detection
c     limits [see discussion in 2209.03360]; the limit returned by this
c     program is the *unrescaled* cross section in the highly
c     non-relativistic limit.
c
c     Author: Torsten Bringmann, 2018-11-23
c     mod TB 2020-10-30: added q2 example
c     mod HK, TB 2022: made mass scan adaptive
c     mod ER 2025-04-29: Edited to use 2->3 logic with 2-loop
c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none
      include 'dsddcom.h'   ! needed for direct access to variables stored in
                            ! direct detection-related common blocks
      include 'dsidtag.h'   ! contains information about which particle module is used.
      include 'dsmpconst.h' ! for m_p

      real*8 tstart, tcurrent, tfinish    ! aux
      integer ierr, iwarn, i
      character*80 outfile     ! output file
      real*8 mdm, sigstart     ! user input
      logical sigmamax
      character*20 lightmedtag, experiment
c... only needed for light mediator/finite DM size
      real*8 mmed, lam
      common /mymed/ mmed, lam
      save /mymed/
c... (further) input parameters to determine scan behaviour
      real*8 mmin, mmax, minminc, stepsize, eps, siglarge, sigtot !, siglargev
      integer npoints
      logical contsearch
c... needed for the actual scan over cross sections
      real*8 count, countlo, counthi, siglo, sighi, counttmp, sigsi, sigvan,slopesign
c...  needed for adaptive mass scan
      real*8 sigsav, mfailed, minc, mincref,sigincred, atten_lim
      integer iback
      logical overshoot, undershoot, newmass, approx_sigmamax

ccc functions
      real*8 dsddDMCRcountrate, dsddDMCRsigtot, dsddDMCRsigtot_2loop, dsddDMCRinvert_cpsi
      

      call CPU_TIME(tstart)
      call dsinit


********************
***  USER INPUT  ***
********************
c... choose which experiment (and determine output filename correspondingly)
      experiment = 'Borexino' ! see src/dd/dsddDMCRcountrate.f for options
      experiment = 'MiniBoone' ! see src/dd/dsddDMCRcountrate.f for options
      ! experiment = 'Xenon1t' ! see src/dd/dsddDMCRcountrate.f for options
      ! experiment = 'Darwin' ! see src/dd/dsddDMCRcountrate.f for options

c... limit type
      sigmamax = .false. ! set to .true. to obtain maximal sigma that can be
                         ! probed (due to soil absorption)
      ! sigmamax = .true.
      sigstart = 1.0d-33 ! start value to search for limit; this should not be too
                         ! far off the true limit for the bisection to work
      if (experiment.eq.'MiniBoone') sigstart=5.3d-28
c... We do not need a different sigstart since mdm=mdmmin has the same upper
c... and lower limit. This is different from previously.
      if (sigmamax) then 
        sigstart = 3.8817d-27
      endif
      CRDM_2loop = .true.
      approx_sigmamax = .false. ! Finds sigsi which gives the 2loop sigma = atten_lim
      atten_lim = 1.d-28 ! attenuation limit for the experiment

c...model settings
c      lightmedtag = ''
      lightmedtag = 'CoSIMP'! to guarantee unique file names, indicating the
                                    ! choice for the mediator mass
                                    ! lightmedtag='' is recommended when compiling
                                    ! with 'make DDCR_limits'
      mmed = 1.d-1     ! fixed DM mass [GeV]
*************************
***  END USER INPUT   ***
*************************


      write (*,*) '-------------------------------------------------------'
      write (*,*) 

c... determine output filename based on user input above
      if (sigmamax) then
        if (CRDM_cs) then
          outfile='data/DMCR_limitsMAX_cs_'
        else
          outfile='data/DMCR_limitsMAX_'
        endif
      else
        if (CRDM_cs) then
          outfile='data/DMCR_limits_cs_'
        else
          outfile='data/DMCR_limits_'
        endif
      endif           
      outfile=outfile(1:index(outfile,' ')-1)//
     &           experiment(1:index(experiment,' ')-1) 
      outfile=outfile(1:index(outfile,' ')-1)//'_'//
     &lightmedtag(1:index(lightmedtag,' ')-1)//'_DIS.dat' 
          
           
      open(unit=21,file=outfile)
      write(*,*) 'Determining limits and writing results to file ', outfile
      write(*,*) 'Please wait...'
      write(*,*)

c... scan setup
      if (CRDM_2loop) then
        if (experiment.eq.'Borexino') mmin = 5.736d-2   ! Minimum for Borexino (quenching makes it higher!)
        if (experiment.eq.'Borexino') mmin = 5.739d-2   ! Minimum for Borexino (quenching makes it higher!)
        if (experiment.eq.'Xenon1t') mmin = 0.009079d0   ! Minimum for Xenon1t
        if (experiment.eq.'Xenon1t') mmin = 0.00911d0   ! Xenon1t both = .false.
        if (experiment.eq.'Darwin') mmin = 0.0041789d0 ! Darwin with both = .true.
        if (experiment.eq.'Darwin') mmin = 0.0041791d0 ! Darwin with both = .false.
        if (experiment.eq.'MiniBoone') mmin = 0.08532d0  ! Minimum for Miniboone
        if (experiment.eq.'MiniBoone') mmin = 0.08539d0  ! Minimum for Miniboone
      else
        mmin = 1.d-4
      endif
      ! mmin = 1.d0
      
      mmax     = 5.d0   ! With form factors the program is very slow for mdm > 10 GeV
      npoints = 40     ! (if curves are not too steep; will be increased otherwise)
      minminc = 1.02d0  ! Do not decrease relative mass stepsize below this
      sigincred = 1.3d0 ! Reduce stepsize if increase in sigma is larger than this
      siglarge = 1d-10  ! ignore adaptive settings above this mass
                        ! (e.g. because it is anyway outside the plotting range)
      contsearch = .false. ! continue with next mass step when no solution is found
                           !  (otherwise reduce mass again)
      eps      = 1.0d-3 ! required precision for limiting cross section
c      n_att1   = 120    ! override default [40] setting in dsddcrdm_init;
c                        ! increases runtime, but may be needed for q2-dependent amplitudes
      sigvan = 3.d-26   ! cm^3 s^-1; results do not depend on this choice
      sigsav = sigstart*1d2
      sigtot = 0.d0
      mfailed  = mmax


      count = 0d0
      mincref = 1.0000001*(mmax/mmin)**((1.)/(1.*npoints-1.))
      minc = mincref
      mdm  = mmin
      i = 1

c... loop over DM masses
50      sighi = 1d-10
        siglo = 1d-40
        overshoot = .false.
        undershoot = .false.
        newmass = .true.
        iback = 0
        sigsi = 1d-1*sigsav ! this should make sure that we are far enough away
                            ! from exponential transition
        stepsize = 5.d0     ! initial stepsize for sigsi


c... The code below tries to find the scattering cross section that exactly
c... corresponds to the reported experimental limit (for a given DM mass).
c... It does so by i) setting up a model with dsgivemodel_generic_wimp
c...               ii) calculating the DMCR recoil rate with dsddDMCRcountrate
c...               iii) narrowing in on scattering rate that corresponds
c...                    to published limit
c...               iv) adaptively setting next mass, then starting from i)

c... loop of cross section
100     call dsgivemodel_generic_wimp(mdm, sigvan, 5, sigsi*1.0d36)
c... For Borexino, we want for better comparison in any case only SD couplings
c... (-> no CRDM component from scttering on He)
        if (experiment.eq.'Borexino'.or.experiment.eq.'Borexino_SD') then
          call dsgivemodel_generic_wimp_opt('spin',0.5d0) ! SD couplings require
                                                          ! DM particle with spin!
          call dsgivemodel_generic_wimp_opt('sigsip',0.0d36)
          call dsgivemodel_generic_wimp_opt('sigsin',0.0d36)
          call dsgivemodel_generic_wimp_opt('sigsdp',sigsi*1.0d36)
        endif  
        call dsmodelsetup(ierr,iwarn)
        counttmp = count
        if (CRDM_2loop.and.approx_sigmamax) then
          sigtot = dsddDMCRsigtot_2loop(m_p, (1.5d0*mdm+m_p)**2)
          count = sigtot/atten_lim
        else
          count = dsddDMCRcountrate(experiment)
        endif
        write(*,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') 
     &       'mdm =', mdm, '   sigsi =', sigsi, '   count =', count, 
     &       '    siglow:', siglo, '   sighigh:', sighi, '   msq:', dsddDMCRinvert_cpsi()
        if (count.gt.1.0d0) then
          overshoot = .true. 
          sighi = sigsi
          counthi = count
        endif
        slopesign = (count-counttmp)/log(stepsize)! are we on a rising or a falling slope?
        if (newmass) then
          slopesign = 0d0
          newmass = .false.
        endif
        if (count.lt.1.0d0.and.(.not.undershoot).and.((sigmamax.and.slopesign.lt.0.0d0)
     &      .or.((.not.sigmamax).and.slopesign.gt.0.0d0))) then
          undershoot = .true. 
          siglo = sigsi
          countlo = count   
        endif
        if (sigsi.lt.1.d-50.or.sigsi.gt.1.d0) goto 250
c... change stepsize adaptively if accuracy goal is not yet met
        if ((abs(sighi-siglo)/siglo).gt.eps) then
        
          if (.not.(overshoot.and.undershoot)) then ! dont (yet) enter enter bisection
            if (sigmamax) then
              if ((stepsize.gt.1.d0.and.count.lt.1.0d0.and.count.lt.counttmp).or.
     &            (stepsize.lt.1.d0.and.count.gt.1.0d0).or.
     &            (stepsize.lt.1.d0.and.count.lt.1.0d0.and.count.lt.counttmp).or.
     &            (stepsize.gt.1.d0.and.count.lt.1.d-6))
     &             stepsize=1.d0/stepsize**0.61 ! change search direction and decrease stepsize
            else
              if ((stepsize.lt.1.d0.and.count.lt.1.0d0.and.count.lt.counttmp).or.
     &            (stepsize.gt.1.d0.and.count.gt.1.0d0).or.
     &            (stepsize.gt.1.d0.and.count.lt.1.0d0.and.count.lt.counttmp)) 
     &             stepsize=1.d0/stepsize**0.71
            endif
            sigsi = sigsi*stepsize
            if (abs(1.0-stepsize).lt.eps) then  ! no solution found...
              counttmp = count
              goto 200
            endif
          else ! bisect
            if (count.gt.1.0d0) then
              sighi = sigsi
              counthi = count
            else
              siglo = sigsi
              countlo = count         
            endif
            sigsi = sqrt(siglo*sighi)
            if (mdm.gt.mmin.and.(abs(sighi-siglo)/siglo).lt.0.9*(sigincred-1d0).and.
     &          sighi.lt.siglarge.and.minc.gt.minminc.and.(mfailed.ge.mmax).and.
     &          (sigsi/sigsav.gt.sigincred.or.sigsav/sigsi.gt.sigincred)) then
               mdm = mdm/minc
               minc = minc**0.7
               mdm = mdm*minc
               sigsav = sigsi
               goto 50
            endif
          endif
          goto 100 ! re-calculate countrate with new value of sigsi
        endif
        if (mdm.eq.mmin) then
          sigsav=sigsi
        endif
                
        if (sigsi.lt.siglarge.and.minc.gt.minminc.and.(mfailed.ge.mmax).and.
     &     (sigsi/sigsav.gt.sigincred.or.sigsav/sigsi.gt.sigincred)) then
           mdm = mdm/minc
           minc = minc**0.7
           mdm = mdm*minc
           sigsav = sigsi
           goto 50
        endif
      
 200    if (overshoot.and.undershoot) then
c... when closing in on the maximal mass, we first need to double-check that we
c..  are really on the correct branch, and that the bisection has converged
          if (mfailed.lt.mmax) then
            if (abs(1.0-stepsize).lt.eps) goto 300
          endif
c          sigtot = dsddDMCRsigtot((1.d3)**2)   ! sig at sqrt(s)= 1 TeV
          if (CRDM_2loop) then
            sigtot = dsddDMCRsigtot_2loop(m_p, ((1.d0+1.d-10)*mdm+m_p)**2)
          else
            sigtot = dsddDMCRsigtot((3.d0*mdm+m_p)**2)   ! sig at Tdm = Tdmmin + mdm
          endif

          call CPU_TIME(tcurrent)
          if (CRDM_2loop.and.sigtot.gt.0.d0) then
            write(21,*) mdm, dsddDMCRinvert_cpsi(), sigtot
          else
            write(21,*) mdm, sigsi, sigtot
          endif
          if (mdm.gt.mmax) goto 300 ! This region is already excluded and not tabulated
          if (mdm.gt.8.d0) goto 300 ! Slightly stricter than the above
c write(*,*)
          write(*,'(I4,A,ES12.4,A,ES12.4,A,ES12.4,A,F8.1,A)') i, 
     &     '    mdm:', mdm, '    sigsi:', sigsi, '    sigtot:', sigtot, 
     &     '    Time elapsed:', tcurrent-tstart, ' seconds'
c write(*,*)
          if (mdm.gt.mmin.and.(mfailed.ge.mmax).and.
     6        (sigsi/sigsav.lt.0.97*sigincred.or.sigsi/sigsav.gt.1.03/sigincred)) then
            iback=iback+1
            if (iback.ge.1.and.minc.lt.mincref) minc = minc**1.3
            if (iback.ge.3) then
              minc = mincref
              iback = 0
            endif
          endif
          !sigsav = sigsi
          if (mfailed.ge.mmax.and.sigsi.gt.siglarge) then ! outside the plotting range...
            minc = mincref**3
          endif
        else
          i = i-1
          mfailed = mdm
          write(*,*) 'Did not manage to find sigsi for DM mass =', mdm
          if (contsearch) then
            write(*,*) '[continue with next mass]'
            sigsi = sigstart
          else
            mdm=mdm/minc
          endif
          sigsav = sigsi
        endif   

 250    if (minc.lt.minminc**0.03.and.mfailed.lt.mmax) goto 300
        if ((.not.contsearch).and.mfailed.lt.mmax) minc=minc**0.5
        mdm = mdm*minc
        
        i=i+1
        if (mdm/minc.lt.mmax) goto 50

 300  close(21)

      call CPU_TIME(tfinish)
      write (*,*)
      write (*,*) '-----------------------------------------------------------'
      write (*,*) 'The DarkSUSY DDCR_limits program has finished successfully.'
      write(*,*)  'Results written to file ', outfile
      write (*,*) 'Total time needed (in seconds): ',tfinish-tstart
      write (*,*) 'Particle module that was used:  ', moduletag
      write (*,*) '-----------------------------------------------------------'
      write (*,*)
      end

      
