*******************************************************************************
*** Provides the total cross section for DM-nucleon elastic scattering      ***
***                                                                         ***
*** Input:                                                                  ***
***    mN    - target nucleon mass [GeV]                                    ***
***    s     - CM energy squared [GeV^2]                                    ***
***                                                                         ***
*** Output:                                                                 ***
***    dsddDMCRsigtot_2loop - total cross section [cm^2]                    ***
***                                                                         ***
*** author: edvard-rornes@hotmail.com                                       ***
*** date 2025-05-23                                                         ***
*******************************************************************************
      real*8 function dsddDMCRsigtot_2loop(mN, s)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h
      include 'dsddcom.h'

      real*8 mN, mdm, s, C, sigv32, cut_off_sq, msq, lamb_s_mdm_mN
      real*8 prefactor, A, B, tmin, tmax, intres, mNcom

      
      real*8 c_psi, sip, sin, sdp, sdn

      integer ierr, dsidnumber ! to make sure cross section is only called once per model
      integer idold
      data idold/-123456789/
      save idold, sip, sin, sdp, sdn

      common/dsddDMCR_sigtot_2loop_block/ mNcom

c... functions
      real*8 integrand_sigtot_2loop
      external integrand_sigtot_2loop
      real*8 dsmwimp

      mdm = dsmwimp()

      dsddDMCRsigtot_2loop = 0.d0
      c_psi = 0.d0  ! This is just a dynamical normalization.
                    ! Its units will vary depending on choice for msq.

      mNcom = mN
c... lambda(s,mdm^2,mN^2)
      lamb_s_mdm_mN = s**2 - 2.d0*s*(mdm**2+mN**2) + (mdm**2-mN**2)**2

      if (idold.ne.dsidnumber()) then ! new model -> calc Q2=0d0 cross sections *once*
        call dsddsigmanucleon(0.0d0, 0.0d0, sip, sin, sdp, sdn, ierr)
      endif

      C = (sip+sin)/2.d0
      if (C.eq.0.d0) C = sdp ! 2->2 normalization
      if (CRDM_cs) then
        C = C/s
      endif
      msq = C

      sigv32 = 1.d0/2.d0/(8.d0*mdm**2*m_p)*sqrt(s**2-2.d0*s*(m_p**2+mdm**2)
     &       + (m_p**2-mdm**2)**2)/(8.d0*(4.d0*atan(1.d0))*s)*msq
      cut_off_sq = sqrt(sqrt(3.d0)/(16.d0*atan(1.d0)*mdm*sigv32))
      c_psi = m_p/((16.d0*atan(1.d0))**4*cut_off_sq)*(1.d0-mdm**2/cut_off_sq)
     &      * log((cut_off_sq+m_p**2)/(4.d0*mdm**2))

      c_psi = 2*sqrt(pi*((sip+sin)/2.d0+sdp))*(mdm+m_p)/m_p
     
      prefactor = c_psi**2/(16.d0*(4.d0*atan(1.d0))*lamb_s_mdm_mN)


      A = 2.d0*mN**2-(s+mN**2-mdm**2)**2/(2.d0*s)
      B = lamb_s_mdm_mN/(2.d0*s)
      tmin = A - B
      tmax = A + B

      call dgadap(tmin, tmax, integrand_sigtot_2loop, CRDM_acc/5.d0, intres)

      ! write(*,*) msq_2loop, C

      dsddDMCRsigtot_2loop = prefactor
     &                     * intres
      
      return
      end

c Define a helper function for the integrand
      real*8 function integrand_sigtot_2loop(t)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h
      include 'dsddcom.h'
      real*8 t, mNcom, Q2, ff, lambda2
      common/dsddDMCR_sigtot_2loop_block/ mNcom

      integrand_sigtot_2loop = 4.d0*mNcom**2-t
 
      if (CRDM_form_factor) then
        Q2 = -t
        lambda2 = 0.71d0 ! GeV^2 (H. Kolesova: 0.77 GeV -> lambda2=0.71 GeV^2)
        ff = 1.d0/(1.d0 + Q2/lambda2)**4
        integrand_sigtot_2loop = integrand_sigtot_2loop*ff
      endif
      
      return
      end