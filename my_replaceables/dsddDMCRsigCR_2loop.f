*******************************************************************************
*** Function dsddDMCRsigCR provides the differential scattering cross       ***
*** section of cosmic rays on DM, d\sigma_CR / dTkin.                       ***
***                                                                         ***
***  Input:                                                                 ***
***    Tkin    - kinetic energy of DM particle [GeV]                        ***
***              (recoil energy after scattering)                           ***
***    Tcr     - initial kinetic energy of CR particle [GeV]                ***
***    CRtype  - CR type (1=proton, 2=helium, 3=carbon, 4=oxygen)           ***
***                                                                         ***
***  Output:  d\sigma_CR / dTkin, units: cm**2/ GeV                         ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2020-10-30                                                         ***
*** mod 2021-12 (H.Kolesova): added more elements + form factor suppresion  ***
*** mod 2025-05-26 (E.Rørnes): edited for Co-SIMP CRDM 2-loop               ***
*******************************************************************************
      real*8 function dsddDMCRsigCR_2loop(Tkin,Tcr,CRtype)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h
      include 'dsio.h'

      real*8 Tkin, Tcr
      integer CRtype
      
      real*8 Q2, s, Tkinmax, mcr, mdm, tmp  ! general kinematic quantities
      real*8 sigCR(NCRelements), sigma0, sigij(27,27), ff, lambda2 ! cross sections
      real*8 sip, sin, sdp, sdn
      integer i, ierr, j, k

      integer dsidnumber ! to make sure cross section is only called once per model
      integer idold
      data idold/-123456789/
      save idold, sigCR

c... functions
      real*8 dsmwimp, dsddTrmax, dsddsigmarel

      mdm = dsmwimp()

      if (idold.ne.dsidnumber()) then ! new model -> calc Q2=0d0 cross sections *once*
        call dsddsigmanucleon(0.d0,0.d0,sip,sin,sdp,sdn,ierr)
        do i=1,NCRelements
          sigCR(i) = 0.d0
          ierr = 0
          call dsddsigma(0.0d0,0.0d0,ANCR(i),ZNCR(i),sigij,ierr)
          if (ierr.eq.0.and.sigij(1,1).gt.0d0) then
            sigCR(i) = sigij(1,1)! SI (neglect SD)
            if (i.eq.1) sigCR(1) = sigCR(1) + sigij(4,4) + sdp ! add SD for protons
          else ! approximate coherently enhanced DM-nucleus cross section
            mcr = mNaU(i)*atomicmassunit ! convert masses to GeV
            sigCR(i) = ((sip+sin)/2.d0 +sdp)*ancr(i)**2*(mcr*(m_p+mdm)/m_p/(mcr+mdm))**2
          endif
        enddo
        idold=dsidnumber()
      endif


      dsddDMCRsigCR_2loop=0.0d0
      
c... determine kinematics
      if (CRtype.ge.1.and.CRtype.le.NCRelements) then
        mcr = mNCRau(CRtype)*atomicmassunit
      else
        if (prtlevel.gt.1) then
          write(*,*) 'warning in dsddDMCRsigCRdTkin:'
          write(*,*) 'unimplemented option CRtype = ',CRtype
          write(*,*) 'Setting cross section to zero'
        endif  
        return
      endif  
      Q2 = 2*mdm*Tkin
      s  = (mcr+mdm)**2 + 2*mdm*Tcr      ! cms energy squared
      Tkinmax = dsddTrmax(Tcr, mcr, mdm) ! maximal recoil energy of DM particle
      if (Tkin.ge.Tkinmax) return


c... Now implement energy-independent, isotropic scattering
      sigma0 = 0.d0
      ierr = 0
      sigma0 = sigCR(CRtype)
      if (CRtype.eq.1) then !protons -- Dipole suppression
        lambda2 = 0.71d0 ! GeV^2 (H. Kolesova: 0.77 GeV -> lambda2=0.71 GeV^2)
        ff = 1.d0/(1.d0 + Q2/lambda2)**4 ! Dipole suppression of cross section
      elseif (CRtype.ge.2.and.CRtype.le.NCRelements) then
        ! default: cascade of FF (first FB, SOG if error)
        call dsddffsi(sqrt(Q2),ANCR(CRtype),ZNCR(CRtype),ff,ierr)
c          call dsddffsog(sqrt(Q2),ANCR(CRtype),ZNCR(CRtype),ff,ierr) !'SOG' form factor
c          call dsddffgauss(sqrt(Q2),ANCR(CRtype),ZNCR(CRtype),ff,ierr)
c          call dsddffh41(sqrt(Q2),ANCR(CRtype),ZNCR(CRtype),ff,ierr) !Helm form factor
       endif
       if (ierr.eq.0) sigma0 = sigma0*ff

      tmp = sigma0/Tkinmax ! highly non-relativistic result

c... finally add Q2- and s-dependence (for full relativistic result)
      dsddDMCRsigCR_2loop = tmp*dsddsigmarel(Q2,s,mcr,1)

      return
      end

