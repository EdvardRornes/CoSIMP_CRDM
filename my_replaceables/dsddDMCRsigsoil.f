*******************************************************************************
*** Function dsddDMCRsigsoil provides the differential scattering cross     ***
*** section of a (potentially relativistic) DM on a target nucleus, taking  ***
*** into account both elastic and inelastic energy losses.                  ***
***                                                                         ***
***  Input:                                                                 ***
***    Tdm     - initial kinetic energy of DM particle [GeV]                ***
***    omega   - energy loss of DM particle [GeV]                           ***
***              (=recoil energy of nucleus only for elastic scattering)    ***
***    iel     - index for the nucleus that the DM particles are scattering ***
***              on (see list in dsnuclides.h)                              ***
***                                                                         ***
***  Output:  d\sigma_N / domega, units: cm**2/ GeV                         ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2022-05-05                                                         ***
*******************************************************************************
      real*8 function dsddDMCRsigsoil(Tdm,omega,iel)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h
      include 'dsio.h'
      include 'dsddcom.h'
      
      real*8 omega,Tdm
      integer iel, i
      
      real*8 Q2, Q2ref(4),Trmax, mdm, s, mN, mNp, tmp  ! general kinematic quantities
      real*8 Trmaxnuc, snuc, res_NRnuc                 ! added for description of nucleon scattering
      real*8 sigN(Nelements), sigij(27,27)             ! cross section
      real*8 sip, sin, sdp, sdn, ff, Ifac
      real*8 res_NR,res_el,res_inel, omegaref
      integer ierr,p

      integer dsidnumber     ! to make sure cross section is only called once per model
      integer idold
      data idold/-123456789/
      save idold, sip, sin, sdp, sdn, sigN

c... functions
      real*8 dsmwimp, dsddTrmax, dsddsigmarel, dsddnu_inel
      
      dsddDMCRsigsoil=0.0d0
      res_el=0d0
      res_inel=0d0
      mdm=dsmwimp()
      if (omega.ge.Tdm) return
      
c... determine kinematics
      if (iel.gt.0.and.iel.le.Nelements) then
        mN = mNaU(iel)*atomicmassunit ! convert masses to GeV
      else
        if (prtlevel.gt.1) then
          write(*,*) 'WARNING in dsddDMCRsigsoil:'
          write(*,*) 'unimplemented overburden eleement iel = ',iel
          write(*,*) 'Setting cross section to zero.'
        endif  
        return
      endif  
      Q2 = 2*mN*omega
      s  = (mN+mdm)**2 + 2*mN*Tdm  ! cms energy squared

c... calc q->0 scattering cross only once per model point
      if (idold.ne.dsidnumber()) then
        call dsddsigmanucleon(0.0d0,0.0d0,sip,sin,sdp,sdn,ierr)
        do i=1,Nelements
          sigN(i) = 0d0
          ierr = 0
          call dsddsigma(0.0d0,0.0d0,AN(i),ZN(i),sigij,ierr)
          if (ierr.eq.0.and.sigij(1,1).gt.0d0) then
            sigN(i) = sigij(1,1) ! SI (neglect SD in soil attenuation)
          else ! approximate DM-nucleus cross section
               !  (note that we do not have access to gp and gn,
               !   hence we cannot do better than this here...)
            mNp = mNaU(i)*atomicmassunit ! convert masses to GeV
            sigN(i) = ((sip+sin)/2.+sdp)*an(i)**2*(mNp*(m_p+mdm)/m_p/(mNp+mdm))**2
          endif
        enddo
        idold=dsidnumber()
      endif


c... Coherent elastic scattering
      ierr = 0
      Trmax = dsddTrmax(Tdm, mdm, mN) ! maximal recoil energy of nucleus
                                      ! NB: *must* be kept as function can be called
                                      !     with larger value of Q2!
      res_NR = sigN(iel)/Trmax ! highly non-relativistic result

      if (omega.gt.Trmax) goto 100 ! energy loss too large for elastic scattering

      call dsddffsi(sqrt(Q2),an(iel),zn(iel),ff,ierr) ! cascade of FF (first FB, SOG if error)
      if (ierr.ne.0) then
        if (prtlevel.gt.1) then
          write(*,*) 'WARNING in dsddDMCRsigsoil:'
          write(*,*) 'problem with form factor for eleement iel = ',iel
        endif
      endif
c      ff=1d0 ! DEBUG
      res_el = res_NR*ff*dsddsigmarel(Q2,s,mN,s2N(iel)) ! full relativistic result,
                                                        ! including form factor

c... Inelastic scattering
100   if (CRDM_inelastic.and.omega.gt.0.001) then
        res_inel = 0d0
        snuc = ((m_p+m_n)/2d0+mdm)**2 + (m_p+m_n)*Tdm ! cms energy squared in case of nucleon scattering
        Trmaxnuc = dsddTrmax(Tdm, mdm, (m_p+m_n)/2d0) ! maximal recoil energy of nucleon
        res_NRnuc = (sip+sin)/(2.*Trmaxnuc)+sdp/Trmaxnuc
        Q2ref(1) = (m_p+m_n)*omega              ! [guess for] peak contribution from QE
        Q2ref(2) = (m_p+m_n)*(omega-0.29)       ! [guess for] peak contribution from delta
        Q2ref(3) = (m_p+m_n)*(omega-0.40)       ! [guess for] peak contribution from other res
        Q2ref(4) = 0.3*(m_p+m_n)*(omega-1.0)    ! [guess for] peak contribution from DIS
        do p=1,4
          if (Q2ref(p).gt.(m_p+m_n)*Trmaxnuc) Q2ref(p) = 0.999*(m_p+m_n)*Trmaxnuc
          if (Q2ref(p).lt.1d-10) Q2ref(p)=1d-10
        enddo
        if (tdm.lt.CRDM_inel_Eref) then
          do p=1,4
            Ifac = dsddnu_inel(Tdm,omega,an(iel),zn(iel),p,ierr)
            tmp=res_NRnuc*Ifac*dsddsigmarel(Q2ref(p),snuc,(m_p+m_n)/2d0,1)
            if (prtlevel.gt.1.and.tmp.lt.0d0)
     &        write(*,*) 'ERROR in dsddDMCRsigsoil: negative inelastic cross section!'
            if (tmp.gt.0d0) res_inel = res_inel + tmp
          enddo
        else ! no simulations -- rescale to CRDM_inel_Eref (<=largest tabulated energy): only DIS...
          omegaref = omega*(CRDM_inel_Eref/tdm)**0.25
          Ifac = dsddnu_inel(CRDM_inel_Eref,omegaref,an(iel),zn(iel),0,ierr)
          res_inel = res_NRnuc*Ifac*dsddsigmarel(Q2ref(4),snuc,(m_p+m_n)/2d0,1)
          if (res_inel.lt.0d0) then
            res_inel = 0d0
            if (prtlevel.gt.1) write(*,*) 'ERROR in dsddDMCRsigsoil: negative inelastic cross section!'
          endif
        endif
      endif

      dsddDMCRsigsoil = res_el + res_inel

      return
      end

