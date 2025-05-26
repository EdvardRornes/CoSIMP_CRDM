*******************************************************************************
*** Function dsddDMCRsigtarget_2loop provides the differential scattering cross   ***
*** section of a (potentially relativistic) DM on a target (detector)       ***
*** nucleus, d\sigma_N / dTr, where Tr is the nuclear recoil energy.        ***
***                                                                         ***
***  Input:                                                                 ***
***    TN    - kinetic energy of nucleus in lab frame [GeV]                 ***
***            (recoil energy after scattering)                             ***
***    Tdm   - initial kinetic energy of DM particle in lab frame [GeV]     ***
***                                                                         ***
***  'Hidden' input (common block flags):                                   ***
***    targetoption - determines which of the implemented concrete          ***
***                   target options should be used. Typically set in       ***
***                   dsddDMCRcountrate                                     ***
***                                                                         ***
***  Output:  d\sigma_N / dTr, units: cm^2 GeV^-1                           ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2019-03-10                                                         ***
*** mod  2019-10-30 (added light mediator examples)                         ***
*** mod  2022-05-05 removed soil scattering, added call to general sigma    ***
*** mod  2025-05-26 (E.Rørnes): minor 2-loop changes                        ***
*******************************************************************************
      real*8 function dsddDMCRsigtarget_2loop(TN, mN, Tdm)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h
      include 'dsddcom.h'

      real*8 TN, mN, mdm, s, s_msq, C, sigv32, cut_off_sq, msq_2loop
      real*8 prefactor, t, msq, Tdm, lamb_s_mdm_mN, TNmax, tmp, c_psi_sq

      
      real*8 c_psi, sip, sin, sdp, sdn, Q2, lambda2, ff, c_psicom, sigma0

      integer ierr, dsidnumber ! to make sure cross section is only called once per model
      integer idold
      data idold/-123456789/
      save idold, sip, sin, sdp, sdn

      common/attencom/ c_psicom

c... functions
      real*8 dsmwimp, dsddsigmarel, dsddTrmax

      mdm = dsmwimp()

      dsddDMCRsigtarget_2loop = 0.d0
      c_psi = 0.d0  ! This is just a dynamical normalization.
                ! Its units will vary depending on choice for msq.

c... determine kinematics
      s  = (mN + mdm)**2 + 2.d0*mN*Tdm   ! cms energy squared
      t = -2.d0*mN*TN                  ! (p_1-p_3)^2 with m1 = m3 = mN
      Q2 = -t                          ! (timelike) momentum transfer
c... lambda(s,mdm^2,mN^2)
      lamb_s_mdm_mN = s**2 - 2.d0*s*(mdm**2+mN**2) + (mdm**2-mN**2)**2

      TNmax = dsddTRmax(Tdm, mdm, mN)
      if (TN.gt.TNmax) return

      if (idold.ne.dsidnumber()) then ! new model -> calc Q2=0d0 cross sections *once*
        call dsddsigmanucleon(0.0d0, 0.0d0, sip, sin, sdp, sdn, ierr)
      endif

c... Now implement energy-independent, isotropic scattering
      C = 0.0d0
      ierr = 0
      if ((targetoption/10).eq.1) then     ! 'xenon1t' or 'Darwin'
        C = (sip+sin)/2.   ! NB: This is not the actual Xe cross section,
                           ! but allows to compare to reported limits 
                           ! *per nucleon*, including form factors
        ff = 1.d0
      else  ! 'Borexino' or 'MiniBoone'
        C = sip+sdp        ! proton target
        lambda2 = 0.71d0 ! GeV^2 (H. Kolesova: 0.77 GeV -> lambda2=0.71 GeV^2)
        ff = 1.d0/(1.d0 + Q2/lambda2)**4 ! Dipole form factor suppression
      endif

      sigma0 = C

      tmp = sigma0/TNmax  ! Highly non-rel limit

      if ((targetoption/10).eq.1) then
c... We have not performed the 2->2 loop diagram for scalars, thus we use
c... spin 1/2 for both cases. 
        dsddDMCRsigtarget_2loop = tmp*dsddsigmarel(Q2,s,mN,1)
      else
        dsddDMCRsigtarget_2loop = tmp*dsddsigmarel(Q2,s,mN,1)
      !   write(*,*) dsddsigmarel(Q2,s,mN,1)
      endif

      if (CRDM_form_factor) then
        dsddDMCRsigtarget_2loop = dsddDMCRsigtarget_2loop*ff
      endif

      return
      end

