*******************************************************************************
*** Function dsddDMCRsigtarget provides the differential scattering cross   ***
*** section of a (potentially relativistic) DM on a target (detector)       ***
*** nucleus, d\sigma_N / dTr, where Tr is the nuclear recoil energy.        ***
***                                                                         ***
***  Input:                                                                 ***
***    TN    - kinetic energy of nucleus in lab frame [GeV]                 ***
***            (recoil energy after scattering)                             ***
***    mN    - target nucleon mass [GeV]                                     ***
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
*** mod 2025-01 (E.Rørnes): adapted to 2->3 Co-SIMP CRDM                    ***
*******************************************************************************
      real*8 function dsddDMCRsigtarget(TN, mN, Tdm)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h
      include 'dsddcom.h'

      real*8 TN, TNmax, mN, mdm, s, s_msq, t13, msq, Tdm, Tdmmin, lamb_s_mdm_mN
      real*8 prefactor, s45min, s45max, intres

      
      real*8 C, sip, sin, sdp, sdn, Q2, lambda2

      integer ierr, dsidnumber ! to make sure cross section is only called once per model
      integer idold
      data idold/-123456789/
      save idold, sip, sin, sdp, sdn

      common/dsddDMCRsigtarget_mass/ mdm ! Pass to integrand_sigtarget

c... functions
      external integrand_sigtarget
      real*8 dsmwimp

      mdm = dsmwimp()

      dsddDMCRsigtarget = 0.d0
      C = 0.d0  ! This is just a dynamical normalization.
                ! Its units will vary depending on choice for msq.

c... determine kinematics
      s  = (mN + mdm)**2 + 2.d0*mN*Tdm   ! cms energy squared
      t13 = -2.d0*mN*TN                  ! (p_1-p_3)^2 with m1 = m3 = mN
      Q2 = -t13                          ! (timelike) momentum transfer
c... lambda(s,mdm^2,mN^2)
      lamb_s_mdm_mN = s**2 - 2.d0*s*(mdm**2+mN**2) + (mdm**2-mN**2)**2

c... Minimum DM kinetic energy in the lab frame to produce another DM particle.
      Tdmmin = 3.d0*mdm**2/(2.d0*mN) + mdm

      if (Tdm.le.Tdmmin) then
        write(*,*) 'Tdm <= Tdm_min, Tdm =', Tdm, 'Tdm_min =', Tdmmin
        write(*,*) 'Returning dsddDMCRsigtarget = 0'
        return
      endif

      if (idold.ne.dsidnumber()) then ! new model -> calc Q2=0d0 cross sections *once*
        call dsddsigmanucleon(0.0d0, 0.0d0, sip, sin, sdp, sdn, ierr)
      endif

      C = (sip+sin)/2.d0
      if (CRDM_form_factor.and.((targetoption/10).ne.1)) then
        lambda2 = 0.71d0 ! GeV^2 (H. Kolesova: 0.77 GeV -> lambda2=0.71 GeV^2)
        C = C/(1d0 + Q2/lambda2)**4
      endif

c... Maximum recoil energy possible given for any given Tdm
      TNmax = (s + mN**2 - mdm**2)*(s + mN**2 - 4.d0*mdm**2)
     &      / (4.d0*s*mN) + sqrt(lamb_s_mdm_mN*(s**2-2.d0*s*(mN**2+4.d0*mdm**2)
     &      + (mN**2 - 4.d0*mdm**2)**2))/(4.d0*s*mN) - mN

      if (TN.gt.TNmax) then
        write(*,*) 'TN > TN_max, TN =', TN, 'TN_max =', TNmax
        write(*,*) 'Returning dsddDMCRsigtarget = 0'
        return
      endif

c... Integration limits
      s45min = 4.d0*mdm**2
      s45max = (s + t13 + mdm**2)/2.d0 - (s-mdm**2)*(mN**2-t13)/(2.d0*mN**2)
     &      + sqrt(lamb_s_mdm_mN*(t13**2 - 4.d0*t13*mN**2))/(2.d0*mN**2)

      if (s45max.lt.s45min) then
        write(*,*) 'WARNING: s45max < s45min! Returning dsddDMCRsigtarget = 0'
        write(*,*) 's45min =', s45min, 's45max =', s45max
        write(*,*) 'TN =', TN, 'Tdm =', Tdm
        return
      endif

c... 16*atan(1)=4*pi, 4F=2*sqrt(lamb_s_mdm_mN), symmetry factor included.
c... Note no sqrt on lamb since there is another sqrt(lamb) factor in the prefactor.
      prefactor = 1.d0/2.d0*mN/(16.d0*atan(1.d0))**3/(2.d0*lamb_s_mdm_mN)

      if (CRDM_cs) then
c... Here we use the proton mass instead of the nucleus mass since the dynamics
c... relate to the proton. Coherent enhancement takes care of the rest.
        s_msq = (m_p + mdm)**2 + 2.d0*m_p*Tdm
        C = C/s_msq
      endif
      msq = C

c... This is the analytical result of the integral however as s45max >>> mdm^2 
c... this overflows. Thus we use numerical integration instead.
      if (4*mdm**2/s45max.lt.1.d-7) then
c... First order expansion
        intres = sqrt(s45max*(s45max - 4.d0*mdm**2))- 2*mdm**2*(log(s45max)-2.d0*log(mdm))
      endif
      intres = sqrt(s45max*(s45max - 4.d0*mdm**2))
     &       - 4.d0*mdm**2*atanh(sqrt(1.d0 - 4.d0*mdm**2/s45max))

      ! call dgadap(s45min, s45max, integrand_sigtarget, CRDM_acc/5.d0, intres)

      if (intres.ne.intres) then
        write(*,*) 'WARNING: dsddDMCRsigtarget integral result returned NaN!'//
     &             'Returning dsddDMCRsigtarget = 0'
        return
      endif
      if (intres.lt.0.d0) then
        write(*,*) 'WARNING: dsddDMCRsigtarget integral result < 0!'//
     &             'Returning dsddDMCRsigtarget = 0'
        return
      endif

      dsddDMCRsigtarget = prefactor*intres
     &                  * msq
      
      return
      end

c Define a helper function for the integrand
      real*8 function integrand_sigtarget(s45)
      implicit none
      real*8 s45, mdm
      common/dsddDMCRsigtarget_mass/ mdm

c integrand_sigtarget = sqrt(s45**2 - 4.d0*s45*mdm**2) / s45
c... This is slightly faster with no noticeable numerical errors
      integrand_sigtarget = sqrt(1.d0 - 4.d0*mdm**2/s45)
      

      return
      end

