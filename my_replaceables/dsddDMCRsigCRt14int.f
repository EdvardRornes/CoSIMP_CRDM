*******************************************************************************
*** Performs the t14 integral for the CRDM cross section                    ***
***                                                                         ***
*** Input:                                                                  ***
***      s      - CM energy squared [GeV^2]                                 ***
***      t13    - (p1-p3)^2 [GeV^2]                                         ***
***      s45    - (p4+p5)^2 [GeV^2]                                         ***
***      s35    - (p3+p5)^2 [GeV^2]                                         ***
***      CRtype - CR type {1,2,3,4}                                         ***
***                                                                         ***
*** Output:                                                                 ***
***      dsddDMCRsigCRs14int - dsig/(dt13 ds45 ds35)/prefactor [cm^2]       ***
***                                                                         ***
*** Note:  The units are off since some kinematical prefactors are omitted  ***
***        and instead included in sigCRff.                                 ***
*** Note2: This is only valid in the CM frame with the momentum labeling    ***
***        DM(p1)CR(p2)->DM(p3)DM(p4)CR(p5)                                 ***
***                                                                         ***
*** author: edvard-rornes@hotmail.com                                       ***
*** date 2025-05-14                                                         ***
*******************************************************************************
      real*8 function dsddDMCRsigCRt14int(s, t13, s45, s35, CRtype)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h
      include 'dsddcom.h'

      real*8 mcr, mdm, s, t13, A, B, s45, s35, lamb_s_s45_mdm
      real*8 t14min, t14max, t14min_aux, t14max_aux, intres, scom, t13com, s45com, s35com
      real*8 D, E, F, G, H
      integer CRtype, CRtypecom

c... For the alternative integration routine
c.   (set CRDM_high_acc = .true. in dsddcrdm_init.f to use this)
      integer limit, neval, ier, last
      parameter (limit = 10) ! max number of sample points from the integrand
                             ! 10 works great for the most part. 20 is very slow
                             ! 5 seems to be not so accurate.
      real*8 epsabs, epsrel, abserr
      real*8 alist(limit), blist(limit), rlist(limit), elist(limit)
      integer iord(limit)

      common/dsddDMCRsigCRt14_block/ mdm, mcr, scom, t13com, s45com, s35com, 
     &                               A, B, lamb_s_s45_mdm, CRtypecom
      common/s35_vars/ D, E, F, G, H

c... functions
      external integrand_sigCRt14
      real*8 dsmwimp

      mdm = dsmwimp()
      dsddDMCRsigCRt14int = 0.d0
      
      if (mdm.le.0.d0) return

      if (CRtype.ge.1.and.CRtype.le.NCRelements) then
        mcr = mNCRau(CRtype)*atomicmassunit
      else
        write(*,*) 'warning in dsddDMCRsigCRt14int:'
        write(*,*) 'unimplemented option CRtype = ',CRtype
        write(*,*) 'Setting cross section to zero'
        return
      endif

c... Pass to integrand
      CRtypecom = CRtype
      scom = s
      t13com = t13
      s45com = s45
      s35com = s35
      lamb_s_s45_mdm = mdm**4-2.d0*mdm**2*(s+s45)+(s-s45)**2
      if (lamb_s_s45_mdm.le.0.d0) then
        write(*,*) 'WARNING: dsddDMCRsigCRt14int: lamb_s_s45_mdm = ', lamb_s_s45_mdm
        return
      endif

c... Components of the integration limits are used in the integrand as well. 
c... Thus, for numerical purposes, it is much faster to calculate them like this.
    !   A = mcr**4*(mdm**2+s-s45)+mcr**2*(-3.d0*mdm**4+mdm**2*(s+s35)-s*(s35+s45
    !  &  +2.d0*t13)+s45*(s35+s45))+2.d0*mdm**6+mdm**4*(-5.d0*s+s35-s45+t13)
    !  &  +mdm**2*(4.d0*s**2+s*(-2.d0*s35-3.d0*s45+2.d0*t13)-t13*(s35+s45)+s45*
    !  &  (s45-3.d0*s35))+t13*(-s**2+s*(s35+s45)+s35*s45)-s*(s-s45)*(s-s35-s45)
      
    !   C = -mcr**4*s+mcr**2*(mdm**4+mdm**2*(2.d0*s-s35-s45)-s**2+s
    !  &  *(s35+s45)+s35*s45)-2.d0*mdm**6+mdm**4*(s+s35+s45)-mdm**2*(s*(s35+s45)
    !  &  -2.d0*s35*s45)-s35*s45*(-s+s35+s45)
    
    !   B = 2.d0*sqrt(C)*sqrt(t13*(mcr**2*(mdm**2+s-s45)+s45*(mdm**2+s)-(mdm**2
    !  &    -s)**2)-mdm**2*(mcr**2-s45)**2-s*t13**2)

c... We can however make this even faster by computing as much as possible
c... in the outermost integrand. (see dsddDMCRsigCRs35int.f)
c...  (NB: These are not the same A and B as in s35int function!)
      A = s35*D + E
      B = 2.d0*sqrt(G+s35*F-s45*s35**2)*H

      if (B.le.0.d0 .or. B.ne.B) return

c... Integration limits
      t14min_aux = (A-B)/lamb_s_s45_mdm
      t14max_aux = (A+B)/lamb_s_s45_mdm

c... Numerical issues occur at t14min, thus we add a small offset based on the
c... difference between t14max and t14min. This heavily speeds up this function
c... with only a marginal affect on the result.
      t14min = t14min_aux + 1.d-4*(t14max_aux-t14min_aux)
      t14max = t14max_aux - 1.d-4*(t14max_aux-t14min_aux)

      if (t14max.gt.0.d0) then
    !     write(*,*) 'Warning: Unphysical t14max: ', t14max, 
    !  &             '(it must be <= 0 when m1=m4) Setting t14max = 0.'
        t14max = 0.d0
      endif
      if (t14min.ne.t14min .or. t14max.ne.t14max) then
        write(*,*) 'WARNING: dsddDMCRsigCRt14int: t14min or t14max = NaN!'
        write(*,*) B
        return
      endif

      if (CRDM_high_acc) then
        epsabs = 0.d0
        epsrel = CRDM_acc/5.d1
        call dqagse(integrand_sigCRt14, t14min, t14max, epsabs, epsrel, limit, 
     &      intres, abserr, neval, ier, alist, blist, rlist, elist, iord, last)
      else
        call dgadap(t14min, t14max, integrand_sigCRt14, CRDM_acc/5.d0, intres)
      endif

      if (intres.ne.intres) then
        write(*,*) 'WARNING: dsddDMCRsigCRt14int integral returned NaN!'//
     &             ' Returning dsddDMCRsigCRt14int = 0'
        write(*,*) 'Tdm =', -t13/mdm/2.d0, 'mdm =', mdm, 'mcr =', mcr
        return
      elseif (intres.lt.0.d0) then
    !     write(*,*) 'WARNING: dsddDMCRsigCRt14int integral < 0!'//
    !  &             ' Returning dsddDMCRsigCRt14int = 0'
        return
      endif

      dsddDMCRsigCRt14int = intres
      
c... This is just a weird normalization which helps with numerical accuracy
c... This is an alternative to the more accurate integration routine which is
c... faster with similar accuracy. (found after thesis was submitted)
      ! if (-t13.gt.1.d2) dsddDMCRsigCRt14int = dsddDMCRsigCRt14int/(t13com**2)
        

      return
      end

c Define a helper function for the integrand
      real*8 function integrand_sigCRt14(t14)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h
      include 'dsddcom.h'
      real*8 t14, mcr, mdm, Q2, ff, lambda2, lamb_s_s45_mdm
      real*8 scom, t13com, s45com, s35com, S14, C14sq, A, B, omega
      integer CRtype, CRtypecom, ierr
      common/dsddDMCRsigCRt14_block/ mdm, mcr, scom, t13com, s45com, s35com, 
     &                               A, B, lamb_s_s45_mdm, CRtypecom

      integrand_sigCRt14 = 0.d0

      C14sq = ((lamb_s_s45_mdm*t14-A)/B)**2

      if (C14sq.ge.1.d0) then
        S14 = 1.d-5 ! this is roughly the error of sqrt(C14sq)
        ! return
      else
        S14 = sqrt(1.d0-C14sq)
      endif
      ! S14 = sqrt(1.d0-C14sq)
      ff = 1.d0

      integrand_sigCRt14 = 1.d0/S14
 
      if (CRDM_form_factor) then
        Q2 = mdm**2-mcr**2-scom+s35com+s45com-t13com-t14
        ! if (Q2.gt.1d1) then
        !   integrand_sigCRt14 = 0.d0
        !   return
        ! endif
        CRtype = CRtypecom
        if (CRtype.eq.1 .and. (Q2.gt.0.d0)) then !protons -- Dipole suppression
          lambda2 = 0.71d0 ! GeV^2 (H. Kolesova: 0.77 GeV -> lambda2=0.71 GeV^2)
          ! if (Q2/lambda2.gt.1.d0) return  ! Test excluding large momentum transfers
          ff = 1.d0/(1.d0 + Q2/lambda2)**4
        ! elseif (CRtype.eq.2) then  ! Test dipole suppression for He
        !   lambda2 = 0.41d0**2
        !   ff = 1.d0/(1.d0 + Q2/lambda2)**4
        elseif (CRtype.ge.2.and.CRtype.le.NCRelements .and. (Q2.gt.0.d0)) then
        ! default: cascade of FF (first FB, SOG if error)
          call dsddffsi(sqrt(Q2),ANCR(CRtype),ZNCR(CRtype),ff,ierr)
        elseif (CRtype.ge.1.and.CRtype.le.NCRelements .and. (Q2.le.0.d0)) then
          ff = 1.d0
        else
          write(*,*) 'WARNING: CRtype = ', CRtype, ' is not implemented!'
          return
        endif
        integrand_sigCRt14 = integrand_sigCRt14*ff
      endif

c... This is for calculating the minimal average energy loss of outgoing dm
      if (CRDM_EnergyLoss) then
        omega = min(s45com-t13com,s35com-t14)/(2.d0*mcr)-mcr/2.d0
        integrand_sigCRt14 = integrand_sigCRt14*omega
      endif
c... This is just a weird normalization which helps with numerical accuracy
      ! if (-t13com.gt.1.d2) integrand_sigCRt14 = integrand_sigCRt14*(t13com**2)
      
      return
      end
