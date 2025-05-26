*******************************************************************************
*** Function dsddDMCREnergyLoss computes the average energy loss times      ***
*** cross section for the Co-SIMP scattering process chi N -> chi chi N.    ***
*** It integrates over the outgoing dark matter energy T_chi using the      ***
*** full 2->3 differential cross section with form factor suppression.      ***
***                                                                         ***
*** Note: if CRDM_EnergyLoss is set to .false., this function will return   ***
***       the integrated cross section.                                     ***
***                                                                         ***
*** Input:                                                                  ***
***      s   - CM energy squared [GeV^2]                                    ***
***      mN  - target nucleon mass [GeV]                                    ***
***                                                                         ***
*** Output:                                                                 ***
***      result  - <EnergyLoss>*sigma [cm^2 GeV]                            ***
***                (if CRDM_EnergyLoss=.true.)                              ***
***              - sigma [cm^2]                                             ***
***                (if CRDM_EnergyLoss=.false.)                             ***
***                                                                         ***
*** author: edvard-rornes@hotmail.com                                       ***
*** date 2025-05-14                                                         ***
*******************************************************************************

      real*8 function dsddDMCREnergyLoss(s, mN)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h
      include 'dsddcom.h'
      real*8 mN, mdm, TN, Tdmmin, Tdmmax, intres, s, Term1, Term2
      common/dsddDMCREnergyLossBlock/ TN ! Pass to integrand

c... functions
      external integrand_EnergyLoss
      real*8 dsmwimp

      dsddDMCREnergyLoss = 0.d0

      mdm = dsmwimp()

      if (mdm.lt.1.d-50) return
      
      TN = (s-(mN+mdm)**2)/(2.d0*mdm)
      
      Term1 = (s+mdm**2-mN**2)*(s+mdm**2-(mdm+mN)**2)/(4.d0*s*mdm)-mdm
      Term2 = sqrt((s**2 - 2.d0*s*(mN**2 + mdm**2)
     &      + (mN**2 - mdm**2)**2)*(s**2 - 2.d0*s*(mdm**2 + (mdm+mN)**2)
     &      + (mdm**2 - (mdm+mN)**2)**2))/(4.d0*s*mdm)

      Tdmmin = Term1 - Term2
      Tdmmax = Term1 + Term2
      
      if (Tdmmin.lt.0.d0) Tdmmin = 1.d-5

c... Note that this Tdm which is being integrated over is a Tdm out, unlike
c... the one in DDCR_EnergyLoss which is Tdm in.
      call dgadap(Tdmmin, Tdmmax, integrand_EnergyLoss, CRDM_acc/5.d0, intres)

      if (intres.ne.intres) then
        write(*,*) 'WARNING: Energy loss integral returned NaN!'//
     &             ' Returning Energy loss = 0'
        write(*,*) 's=',s,'mdm=',mdm,'mN=',mN,'Tdmmin=',Tdmmin,'Tdmmax=',Tdmmax
        return
      endif
      if (intres.lt.0.d0) then
        write(*,*) 'WARNING: Energy loss integral < 0!'//
     &             ' Returning Energy loss = 0'
        return
      endif
      
      dsddDMCREnergyLoss = intres
     &                   * 1.d-20   ! Undo normalization in integrand

      return
      end


c Define a helper function for the integrand
      real*8 function integrand_EnergyLoss(Tdm)
      implicit none
      real*8 dsddDMCRsigCRff
      real*8 TN, Tdm
      common/dsddDMCREnergyLossBlock/ TN

      integrand_EnergyLoss = dsddDMCRsigCRff(Tdm, TN, 1)
     &                     * 1.d20  ! For better numerical result


      return
      end

