*******************************************************************************
*** Function dsddDMCREnergyLoss computes the average energy loss for an     ***
*** elastic scattering process chi N -> chi N. It integrates over the       ***
*** outgoing dark matter energy T_chi using the                             ***
*** full 2->3 differential cross section with form factor suppression.      ***
***                                                                         ***
*** Input:                                                                  ***
***      s   - CM energy squared [GeV^2]                                    ***
***      mN  - target nucleon mass [GeV]                                    ***
***                                                                         ***
*** Output:                                                                 ***
***      result    - <EnergyLoss> [GeV]                                     ***
***                                                                         ***
*** author: edvard-rornes@hotmail.com                                       ***
*** date 2025-05-14                                                         ***
*******************************************************************************

      real*8 function dsddDMCREnergyLoss2t2(s, mN)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h
      include 'dsddcom.h'
      real*8 mN, mdm, tmin, tmax, intres, s, Term1, Term2, mcrcom
      common/dsddDMCREnergyLoss2t2Block/ mcrcom ! Pass to integrand_sigtot

c... functions
      external integrand_EnergyLoss2t2
      real*8 dsmwimp

      dsddDMCREnergyLoss2t2 = 0.d0

      mdm = dsmwimp()

      if (mdm.lt.1.d-50) return

      mcrcom = mN

      Term1 = 2*mdm**2 - (s+mdm**2-mN**2)**2/(2.d0*s)
      Term2 = (s**2-2.d0*s*(mN**2+mdm**2)+(mN**2-mdm**2)**2)/(2.d0*s)
      tmin = Term1 - Term2
      tmax = Term1 + Term2

      call dgadap(tmin, tmax, integrand_EnergyLoss2t2, CRDM_acc/5.d0, intres)

      if (intres.ne.intres) then
        write(*,*) 'WARNING: 2t2 energy loss integral returned NaN!'//
     &             ' Returning dsddDMCREnergyLoss2t2 = 0'
        return
      endif
      if (intres.lt.0.d0) then
        write(*,*) 'WARNING: 2t2 energy loss integral < 0!'//
     &             ' Returning dsddDMCREnergyLoss2t2 = 0'
        return
      endif
      
      dsddDMCREnergyLoss2t2 = intres

      return
      end


c Define a helper function for the integrand
      real*8 function integrand_EnergyLoss2t2(t)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h
      include 'dsddcom.h'
      real*8 t, mcrcom, ff, Q2, lambda2
      integer CRtype, ierr
      common/dsddDMCREnergyLoss2t2Block/ mcrcom

      integrand_EnergyLoss2t2 = 1.d0
      if (CRDM_EnergyLoss) then
        integrand_EnergyLoss2t2 = integrand_EnergyLoss2t2*(-t/(2.d0*mcrcom))
      endif

      if (CRDM_form_factor) then
        Q2 = -t
        CRtype = 4
        if (CRtype.eq.1) then !protons -- Dipole suppression
          lambda2 = 0.71d0 ! GeV^2 (H. Kolesova: 0.77 GeV -> lambda2=0.71 GeV^2)
          ! if (Q2/lambda2.gt.1.d0) return  ! Test excluding large momentum transfers
          ff = 1.d0/(1.d0 + Q2/lambda2)**4
        elseif (CRtype.ge.2.and.CRtype.le.NCRelements) then
        ! default: cascade of FF (first FB, SOG if error)
          call dsddffsi(sqrt(Q2),ANCR(CRtype),ZNCR(CRtype),ff,ierr)
        endif
        integrand_EnergyLoss2t2 = integrand_EnergyLoss2t2*ff
      endif

      return
      end

