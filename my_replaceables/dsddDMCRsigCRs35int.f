*******************************************************************************
*** Performs the s35 integral for the CRDM cross section                    ***
***                                                                         ***
*** Input:                                                                  ***
***      s      - CM energy squared [GeV^2]                                 ***
***      t13    - (p1-p3)^2 [GeV^2]                                         ***
***      s45    - (p4+p5)^2 [GeV^2]                                         ***
***      CRtype - CR type {1,2,3,4}                                         ***
***                                                                         ***
*** Output:                                                                 ***
***      dsddDMCRsigCRs35int - dsig/(dt13 ds45)/prefactor [cm^2]            ***
***                                                                         ***
*** Note:  The units are off since some kinematical prefactors are omitted  ***
***        and instead included in sigCRff.                                 ***
*** Note2: This is only valid in the CM frame with the momentum labeling    ***
***        DM(p1)CR(p2)->DM(p3)DM(p4)CR(p5)                                 ***
***                                                                         ***
*** author: edvard-rornes@hotmail.com                                       ***
*** date 2025-05-14                                                         ***
*******************************************************************************
      real*8 function dsddDMCRsigCRs35int(s, t13, s45, CRtype)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h
      include 'dsddcom.h'

      real*8 mcr, mdm, s, t13, s45, lamb_s_s45_mdm, lamb_s45_mdm_mcr, rescale
      real*8 s35min, s35max, s35min_aux, s35max_aux, intres, scom, t13com, s45com
      real*8 mdm_2, mdm_4, mcr_2, mcr_4
      integer CRtype, CRtypecom
      real*8 A, B, D, E, F, G, H

c... For the alternative integration routine. 
c.   (set CRDM_high_acc = .true. in dsddcrdm_init.f to use this)
      integer limit, neval, ier, last
      parameter (limit = 10) ! max number of sample points from the integrand
                             ! 10 works great for the most part. 20 is very slow
                             ! 5 seems to be not so accurate.
      real*8 epsabs, epsrel, abserr
      real*8 alist(limit), blist(limit), rlist(limit), elist(limit)
      integer iord(limit)

      common/dsddDMCRsigCRs35_block/ mdm, mcr, mdm_2, mdm_4, mcr_2, scom, t13com, s45com,
     &                               lamb_s_s45_mdm, rescale, CRtypecom
      common/s35_vars/ D, E, F, G, H

c... functions
      external integrand_sigCRs35
      real*8 dsmwimp

      mdm = dsmwimp()
      dsddDMCRsigCRs35int = 0.d0
      
      if (mdm.le.0.d0) return
      if (CRtype.ge.1.and.CRtype.le.NCRelements) then
        mcr = mNCRau(CRtype)*atomicmassunit
        rescale = 1.d1**CRtype ! NB: If CRtype indices are ever changed, 
                               !     this MUST be updated!!!
      else
        write(*,*) 'warning in dsddDMCRsigCRs35int:'
        write(*,*) 'unimplemented option CRtype = ',CRtype
        write(*,*) 'Setting cross section to zero'
        return
      endif

c... I believe this is the only place where it is ideal to do this.
      mdm_2 = mdm**2
      mdm_4 = mdm_2**2
      mcr_2 = mcr**2
      mcr_4 = mcr_2**2

c... Pass to common block
      scom = s
      t13com = t13
      s45com = s45
      CRtypecom = CRtype
      lamb_s_s45_mdm = s**2-2.d0*s*(mdm_2+s45)+(mdm_2-s45)**2
      lamb_s45_mdm_mcr = s45**2-2.d0*s45*(mdm_2+mcr_2)+(mdm_2-mcr_2)**2

c... Integration limits (NB: These are not the same A and B as in t14int function!)
      A = mcr_2*(-mdm_2+s+s45)+mdm_4-mdm_2*(s-2.d0*s45)+s45*(s-s45)
      B = sqrt(lamb_s45_mdm_mcr*lamb_s_s45_mdm)

      if ((B.ne.B) .or. (B.eq.0.d0)) then
      !   write(*,*) 'WARNING: dsddDMCRsigCRs35int: B = ', B
        return
      endif

      s35min_aux = (A-B)/(2.d0*s45)
      s35max_aux = (A+B)/(2.d0*s45)

c... The t14 integral suffers close to the integration limits. Thus we shave 
c... off a tiny bit to avoid NaNs. The result is largely independent of this.
      s35min = s35min_aux! + 1.d-8*(s35max_aux-s35min_aux)
      s35max = s35max_aux! - 1.d-8*(s35max_aux-s35min_aux)

      if (s35min.lt.(mdm+mcr)**2) then
        write(*,*) 'Warning: s35min < (mdm+mcr)^2! Setting s35min = (mdm+mcr)^2'
        write(*,*) 's35min = ', s35min, '(mdm+mcr)^2 = ', (mdm+mcr)**2
        s35min = (mdm+mcr)**2
      endif
      if (s35min.ge.s35max) then
        write(*,*) 'Warning: s35min >= s35max!'
        write(*,*) 's35min = ', s35min, 's35max = ', s35max
        return
      endif


c... We pass these to the integrand like this to avoid redundant calculations
c... in the innermost integrand. This speeds up the function by ~ 5-10 times.      
      D = mcr_2*(mdm_2-s+s45)+mdm_4+mdm_2*(-2.d0*s-3.d0*s45-t13)+t13*(s+
     &    s45)+s*(s-s45)
      E = mcr_4*(mdm_2+s-s45)+mcr_2*(-3.d0*mdm_4+mdm_2*s-s*(s45+2.d0*t13)
     &    +s45**2)+2.d0*mdm**6+mdm_4*(-5.d0*s-s45+t13)+mdm_2*(4.d0*s**2-3.d0*
     &    s*s45+2.d0*s*t13+s45**2-s45*t13)-s*(s-s45)*(s-s45+t13)
      F = mcr_2*(-mdm_2+s+s45)+mdm_4-mdm_2*(s-2.d0*s45)+s45*(s-s45)
      G = -mcr_4*s+mcr_2*(mdm_4+mdm_2*(2.d0*s-s45)+s*(s45-s))-2.d0*mdm**6+
     &    mdm_4*(s+s45)-mdm_2*s*s45
      H = sqrt(t13*(mcr_2*(mdm_2+s-s45)+s45*(mdm_2+s)-(mdm_2-s)**2)-mdm_2
     &         *(mcr_2-s45)**2-s*t13**2)
c... if H is 0 then t14 integral will yield 0.
      if (H.eq.0.d0 .or. H.ne.H) return
      
      if (CRDM_high_acc) then
        epsabs = 0.d0
        epsrel = CRDM_acc/5.d1
        call dqagse(integrand_sigCRs35, s35min, s35max, epsabs, epsrel, limit, 
     &      intres, abserr, neval, ier, alist, blist, rlist, elist, iord, last)
      else
        call dgadap(s35min, s35max, integrand_sigCRs35, CRDM_acc/5.d0, intres)
      endif


      if (intres.ne.intres) then
        write(*,*) 'WARNING: dsddDMCRsigCRs35int integral returned NaN!'//
     &             ' Returning dsddDMCRsigCRs35int = 0'
        return
      elseif (intres.lt.0.d0) then
        write(*,*) 'WARNING: dsddDMCRsigCRs35int integral < 0!'//
     &             ' Returning dsddDMCRsigCRs35int = 0'
        write(*,*) s35min, s35max, intres, s45, s, s35max-s35min
        return
      endif

      dsddDMCRsigCRs35int = intres
     &                    / rescale

      return
      end

c Define a helper function for the integrand
      real*8 function integrand_sigCRs35(s35)
      implicit none
      real*8 mdm, mcr, scom, t13com, s45com, s35, s, s45, C34sq, S34
      real*8 mdm_2, mdm_4, mcr_2
      real*8 lamb_s_s45_mdm, lamb_s_s35_mdm, rescale
      real*8 dsddDMCRsigCRt14int
      integer CRtypecom
      common/dsddDMCRsigCRs35_block/ mdm, mcr, mdm_2, mdm_4, mcr_2, scom, t13com, 
     &                               s45com, lamb_s_s45_mdm, rescale, CRtypecom
      
      integrand_sigCRs35 = 0.d0

      s45 = s45com
      s = scom
      lamb_s_s35_mdm = s**2-2.d0*s*(mdm_2+s35)+(mdm_2-s35)**2

      if (lamb_s_s35_mdm.le.0.d0) return

      C34sq = (s*(-2.d0*mcr_2-s+s35)+mdm_4+mdm_2*(2.d0*s-s35-s45)
     &      +s45*(s+s35))**2/(lamb_s_s35_mdm*lamb_s_s45_mdm)

      ! if (C34sq.lt.0.d0.or.C34sq.ge.1.d0) return
      if (C34sq.ge.1.d0) then
        S34 = 1.d-5 ! this is roughly the error of sqrt(C34sq)
        ! return
      else
        S34 = sqrt(1.d0-C34sq)
      endif

      integrand_sigCRs35 = dsddDMCRsigCRt14int(s, t13com, s45, s35, CRtypecom)
     &                   / sqrt(lamb_s_s35_mdm)
     &                   / S34
     &                   * rescale
     
      return
      end
