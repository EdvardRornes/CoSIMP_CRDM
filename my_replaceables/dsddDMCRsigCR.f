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
*** mod 2025-01 (E.Rørnes): adapted to 2->3 case and removed form factors   ***
*******************************************************************************
      real*8 function dsddDMCRsigCR(Tdm, Tcr, CRtype)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h
      include 'dsddcom.h'

      real*8 Tdm, Tcr, Tcrmin, Tdmmax, mcr, mdm, s, t13, s_msq, msq
      real*8 s45min, s45max, prefactor, intres, lamb_s_mdm_mcr
      integer CRtype

      real*8 C, sip, sin, sdp, sdn, sigCR(NCRelements), sigij(27,27)

      integer dsidnumber ! to make sure cross section is only called once per model
      integer idold, ierr, i
      data idold/-123456789/
      save idold, sigCR

      common/dsddDMCRsigCR_mass/ mdm, mcr  ! Pass to integrand_sigCR

c... functions
      external integrand_sigCR, dsddDMCRinvert_cpsi
      real*8 dsmwimp, dsddDMCRinvert_cpsi

      mdm = dsmwimp()
      dsddDMCRsigCR = 0.d0
      
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
      else
        write(*,*) 'warning in dsddDMCRsigCR:'
        write(*,*) 'unimplemented option CRtype = ',CRtype
        write(*,*) 'Setting cross section to zero'
        return
      endif  

c... Minimum CR kinetic energy for the process to occur
      Tcrmin = 3.d0/2.d0*mdm + mcr

      if (Tcr.le.Tcrmin) then
        write(*,*) 'Tcr < Tcr_min, Tcr =', Tcr, 'Tcr_min =', Tcrmin
        write(*,*) 'Returning dsddDMCRsigCR = 0'
        return
      endif

c... determine kinematics
      s  = (mcr+mdm)**2 + 2.d0*mdm*Tcr  ! CM energy squared
      t13 = -2.d0*mdm*Tdm               ! (p1-p3)^2

c... lambda(s,mdm^2,mcr^2)
      lamb_s_mdm_mcr = s**2 - 2.d0*s*(mdm**2+mcr**2) + (mdm**2-mcr**2)**2

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
      write(*,*) 'hi'
c... Maximum Tdm given Tcr
      Tdmmax = - mdm + (s+mdm**2-mcr**2)*(s-mcr**2-2.d0*mdm*mcr) 
     &       / (4.d0*s*mdm) + sqrt(lamb_s_mdm_mcr*(s**2 - 2.d0*s*(2.d0*mdm**2
     &       + mcr**2 + 2.d0*mcr*mdm)+(mcr**2+2.d0*mcr*mdm)**2))/(4.d0*s*mdm)

      if (Tdm.gt.Tdmmax) then 
        write(*,*) 'Tdm > Tdm_max, Tdm =', Tdm, 'Tdm_max =', Tdmmax
        write(*,*) 'Returning dsddDMCRsigCR = 0'
        return
      endif

c... Integration limits
      s45min = (mcr + mdm)**2
      s45max = (s+t13+mcr**2)/2.d0 - (s-mcr**2)*(mdm**2-t13)/(2.d0*mdm**2)
     &       + sqrt(lamb_s_mdm_mcr*(t13**2 - 4.d0*t13*mdm**2))/(2.d0*mdm**2)

c... Not sure if this is simply from numerical errors or if it is a real issue
      if (s45max.lt.s45min) then
        write(*,*) 'WARNING: s45max < s45min! Returning dsddDMCRsigCR = 0'
        write(*,'(A,ES12.2,2(A,ES12.2))') 'mdm', mdm, '  Tdm', Tdm, '  Tcr', Tcr
        return
      endif

c... 1/2=symmetry factor, 16*atan(1)=4*pi, 4*F=2*sqrt(lambda(s,mcr,mdm))
c... Note no sqrt on lamb since there is another sqrt factor in the prefactor.
      prefactor = 1.d0/2.d0*mdm/(16.d0*atan(1.d0))**3/(2.d0*lamb_s_mdm_mcr)

      if (CRDM_cs) then
c... Here we use the proton mass instead of the nucleus mass since the dynamics
c... relate to the proton. Coherent enhancement takes care of the rest.
        s_msq = (m_p+mdm)**2 + 2.d0*mdm*Tcr
        C = C/s_msq
      endif
      msq = C

c... This is the analytical result of the integral. 
c... However, the code is faster and more stable when numerically integrating.
c   intres = sqrt(s45max**2 - 2.d0*s45max*(mdm**2+mcr**2)+(mdm**2-mcr**2)**2)
c  &       - abs(mdm**2-mcr**2)*atanh(((mdm**2-mcr**2)**2-(mdm**2+mcr**2)*s45max)
c  &       / sqrt((mdm**2-mcr**2)**2*s45max**2-2.d0*(mdm**2+mcr**2)
c  &       * (mdm**2-mcr**2)**2*s45max +(mdm**2-mcr**2)**4))

      call dgadap(s45min, s45max, integrand_sigCR, CRDM_acc/5.d0, intres)

      if (intres.ne.intres) then
        write(*,*) 'WARNING: dsddDMCRsigCR integral returned NaN!'//
     &             ' Returning dsddDMCRsigCR = 0'
        return
      elseif (intres.lt.0.d0) then
        write(*,*) 'WARNING: dsddDMCRsigCR integral  integral < 0!'//
     &             ' Returning dsddDMCRsigCR = 0'
        return
      endif

      dsddDMCRsigCR = prefactor*intres
     &              * msq

      return
      end

c Define a helper function for the integrand
      real*8 function integrand_sigCR(s45)
      implicit none
      real*8 s45, mcr, mdm
      common/dsddDMCRsigCR_mass/ mdm, mcr

      integrand_sigCR = sqrt(s45**2 - 2.d0*s45*(mcr**2+mdm**2)
     &                     + (mcr**2-mdm**2)**2)/s45
     
      return
      end