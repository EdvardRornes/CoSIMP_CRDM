*******************************************************************************
*** Function dsddDMCRsigCR provides the differential scattering cross     ***
*** section of cosmic rays on DM for a 2->3 process, d\sigma_CR / dTkin.    ***
*** Currently this program only uses CR as the incoming CR.                 ***
***                                                                         ***
***  Input:                                                                 ***
***    Tdm     - kinetic energy of DM particle in CF [GeV]                  ***
***              (recoil energy after scattering)                           ***
***    Tcr     - initial kinetic energy of CR particle in CF [GeV]          ***
***                                                                         ***
***  Output:  d\sigma_CR / dTkin, units: cm**2/ GeV                         ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2020-10-30                                                         ***
*** mod 2021-12 (H.Kolesova): added more elements + form factor suppresion  ***
*** mod 2025-01 (E.Rørnes): adapted to 2->3 Co-SIMP CRDM                    ***
*******************************************************************************
      real*8 function dsddDMCRsigCRff(Tdm, Tcr, CRtype)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h
      include 'dsddcom.h'

      real*8 Tdm, Tcr, Tcrmin, Tdmmax, mcr, mdm, s, t13, s_msq, msq, rescale
      real*8 s45min, s45max, prefactor, intres, lamb_s_mdm_mcr, scom, s45min_aux, s45max_aux
      integer CRtype, CRtypecom

      real*8 C, sip, sin, sdp, sdn, sigCR(NCRelements), sigij(27,27)  ! Normalization

      integer dsidnumber ! to make sure cross section is only called once per model
      integer idold, ierr, i
      data idold/-123456789/
      save idold, sigCR

c... For the alternative integration routine
      integer limit, neval, ier, last
      parameter (limit = 10) ! max number of sample points from the integrand
                             ! 10 works great for the most part. 20 is very slow
                             ! 5 seems to be not so accurate.
      real*8 epsabs, epsrel, abserr
      real*8 alist(limit), blist(limit), rlist(limit), elist(limit)
      integer iord(limit)

      common/dsddDMCRsigCRff_block/ mdm, mcr, scom, t13, lamb_s_mdm_mcr, rescale, CRtypecom

c... functions
      external integrand_sigCRs45
      real*8 dsmwimp, dsddDMCRinvert_cpsi

      mdm = dsmwimp()
      dsddDMCRsigCRff = 0.d0
      
      if (mdm.le.0) return

      if (idold.ne.dsidnumber()) then ! new model -> calc Q2=0d0 cross sections *once*
        call dsddsigmanucleon(0.d0, 0.d0, sip, sin, sdp, sdn, ierr)
        do i=1,NCRelements
          sigCR(i) = 0.d0
          ierr = 0
          call dsddsigma(0.d0, 0.d0, ANCR(i), ZNCR(i), sigij, ierr)
          if (ierr.eq.0.and.sigij(1, 1).gt.0.d0) then
            sigCR(i) = sigij(1, 1)! SI (neglect SD)
            if (i.eq.1) sigCR(1) = sigCR(1) + sigij(4, 4) ! add SD for protons
          else ! approximate coherently enhanced DM-nucleus cross section
            mcr = mNaU(i)*atomicmassunit ! convert masses to GeV
c... Mass factor removed below as this is already included in the kinematics.
            sigCR(i) = (sip+sin)/2.d0*ANCR(i)**2   ! *(mcr*(m_p+mdm)/m_p/(mcr+mdm))**2
          endif
        enddo
        idold=dsidnumber()
      endif

      if (CRtype.ge.1.and.CRtype.le.NCRelements) then
        mcr = mNCRau(CRtype)*atomicmassunit
        rescale = 1.d5*1.d1**CRtype ! NB!! If CRtype indices are ever changed, 
                                    !      this no longer makes sense!
      else
        write(*,*) 'warning in dsddDMCRsigCR:'
        write(*,*) 'unimplemented option CRtype = ',CRtype
        write(*,*) 'Setting cross section to zero'
        return
      endif  

c... Minimum CR kinetic energy for the process to occur
      Tcrmin = 3.d0/2.d0*mdm + mcr

      if (Tcr.le.Tcrmin) then
        write(*,*) 'Tcr =', Tcr, '< Tcr_min =', Tcrmin
        write(*,*) 'Returning dsddDMCRsigCR = 0'
        return
      endif
c... determine kinematics
      s = (mcr+mdm)**2 + 2.d0*mdm*Tcr  ! CM energy squared
      scom = s
      CRtypecom = CRtype
      t13 = -2.d0*mdm*Tdm              ! (p1-p3)^2

c... lambda(s,mdm^2,mcr^2)
      lamb_s_mdm_mcr = s**2 - 2.d0*s*(mdm**2+mcr**2) + (mdm**2-mcr**2)**2
      if (lamb_s_mdm_mcr.le.0.d0) return

      C = 0.d0  ! This is just a dynamical normalization.
                ! Its units will vary depending on choice for msq.
      ierr = 0

      if (CRDM_2loop) then
        C = dsddDMCRinvert_cpsi()  ! This function provides 2->3 constant C based 
                                   ! on the sigmaNR the program is currently using
      else
        C = (sip+sin)/2.d0
        if (C.eq.0.d0) C = sdp
      endif

c... Integration limits
      s45min_aux = (mcr + mdm)**2
      s45max_aux = (s+t13+mcr**2)/2.d0 + (s-mcr**2)*(t13-mdm**2)/(2.d0*mdm**2)
     &       + sqrt(lamb_s_mdm_mcr*(t13**2 - 4.d0*t13*mdm**2))/(2.d0*mdm**2)

      s45min = s45min_aux !+ 1.d-8*(s45max_aux-s45min_aux)
      s45max = s45max_aux !- 1.d-8*(s45max_aux-s45min_aux)
      
c... Maximum Tdm given Tcr
      Tdmmax = - mdm + (s+mdm**2-mcr**2)*(s-s45min+mdm**2)/(4.d0*s*mdm)
     &       + sqrt((s**2 - 2.d0*s*(mdm**2+s45min) + (mdm**2-s45min)**2)
     &       * lamb_s_mdm_mcr)/(4.d0*s*mdm)

      if (Tdm.gt.Tdmmax) then 
        write(*,*) 'Tdm > Tdm_max, Tdm =', Tdm, 'Tdm_max =', Tdmmax
        write(*,*) 'Returning dsddDMCRsigCRff = 0'
        return
      endif

      if (s45max.lt.s45min) then
        write(*,*) 'WARNING: s45max < s45min! Returning dsddDMCRsigCRff = 0'
        write(*,'(A,ES12.2,2(A,ES12.2))') 'mdm', mdm, '  Tdm', Tdm, '  Tcr', Tcr
        return
      endif

c... 1/2=symmetry factor, 16*atan(1)=4*pi, 4*F=2*sqrt(lambda(s,mcr,mdm))
      prefactor = 1.d0/2.d0*8.d0*s*mdm/(16.d0*atan(1.d0))**4
     &          /(2.d0*lamb_s_mdm_mcr**(3.d0/2.d0))

      if (CRDM_cs) then
c... Here we use the proton mass instead of the nucleus mass since the dynamics
c... relate to the proton. Coherent enhancement takes care of the rest.
        s_msq = (m_p+mdm)**2 + 2.d0*mdm*Tcr
        C = C/s_msq
      endif
      msq = C
      

c... The below was an attempt to speed up the code. However it seems like the
c... common blocks are more expensive even if it requires more computation to
c... not do this.
    !   I = mcr**2-3.d0*mdm**2-s+t13
    !   J = s*(-mcr**2-2.d0*mdm**2+t13)+mcr**2*mdm**2+mdm**4-mdm**2*t13+s**2
    !   K = mcr**2+mdm**2-s
    !   L = -mcr**4+s*(-mcr**2-3.d0*mdm**2+t13)-mdm**4-mdm**2*t13+2.d0*s**2
    !   M = mcr**4*mdm**2-3.d0*mcr**2*mdm**4+s*(mcr**4+mcr**2*mdm**2-2.d0*mcr**2*
    ! &     t13-5.d0*mdm**4+2.d0*mdm**2*t13)+2.d0*mdm**6+mdm**4*t13+s**2*(4.d0*
    !  &    mdm**2-t13)-s**3
    !   N = mcr**2+2.d0*mdm**2+s
    !   O = s*(mcr**2-mdm**2)-mcr**2*mdm**2+mdm**4-s45**2
    !   P = s*(mcr**2-mdm**2)-mcr**2*mdm**2+mdm**4
    !   Q = -mcr**2*s**2+s*(-mcr**4+2.d0*mcr**2*mdm**2+mdm**4)-2.d0*mdm**6+mcr**2
    !  &    *mdm**4
    !   R = 2.d0*mcr**2*mdm**2-mcr**2*t13+mdm**2*t13+s*t13
    !   T = s*(mcr**2*t13+2.d0*mdm**2*t13-t13**2)+mcr**2*mdm**2*t13-mdm**4*t13
    !  &    -mdm**2*s45**2-s**2*t13

      ! call dgadap(s45min, s45max, integrand_sigCRs45, CRDM_acc/5.d1, intres)
      epsabs = 0.d0
      epsrel = CRDM_acc/5.d1
      call dqagse(integrand_sigCRs45, s45min, s45max, epsabs, epsrel, limit, 
     &       intres, abserr, neval, ier, alist, blist, rlist, elist, iord, last)

      if (intres.ne.intres) then
        write(*,*) 'WARNING: dsddDMCRsigCRff integral returned NaN!'//
     &             ' Returning dsddDMCRsigCRff = 0'
        return
      elseif (intres.lt.0.d0) then
        write(*,*) 'WARNING: dsddDMCRsigCRff integral < 0!'//
     &             ' Returning dsddDMCRsigCRff = 0'
        return
      endif

      dsddDMCRsigCRff = prefactor*intres
     &                * msq
     &                / rescale

      return
      end


c Define a helper function for the integrand
      real*8 function integrand_sigCRs45(s45)
      implicit none
      real*8 mdm, mcr, s45, s, scom, t13, C13sq, lamb_s_s45_mdm, lamb_s_mdm_mcr, rescale, S13
      real*8 dsddDMCRsigCRs35int
      integer CRtypecom
      common/dsddDMCRsigCRff_block/ mdm, mcr, scom, t13, lamb_s_mdm_mcr, rescale, CRtypecom

      ! s45 = exp(s45)
      integrand_sigCRs45 = 0.d0
      s = scom
      lamb_s_s45_mdm = s**2-2.d0*s*(mdm**2+s45)+(mdm**2-s45)**2
      if (lamb_s_s45_mdm.le.0.d0) return

      C13sq = (-mcr**2*(mdm**2+s-s45)+mdm**4-mdm**2*(2.d0*s+s45)
     &    +s*(s-s45+2.d0*t13))**2/(lamb_s_mdm_mcr*lamb_s_s45_mdm)

c... C13sq ~ 1 corresponds to backwards scattering in the CF and is heavily CRflux suppressed
c... Thus it basically doesnt matter what we set S13 to in this case.
      if (C13sq.ge.1.d0) then
        S13 = 1.d-5 ! this is roughly the error of sqrt(C13sq)
        ! return
      else
        S13 = sqrt(1-C13sq)
      endif

      integrand_sigCRs45 = dsddDMCRsigCRs35int(s, t13, s45, CRtypecom)
     &                   / sqrt(lamb_s_s45_mdm)
     &                   / S13
     &                   * rescale

      return
      end