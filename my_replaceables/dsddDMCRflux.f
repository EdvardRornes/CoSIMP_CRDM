*******************************************************************************
*** Function dsddDMCRflux provides the local differential flux of DM        ***
*** particles that results from cosmic rays impinging on DM in the          ***
*** diffusive halo, see Bringmann & Pospelov (2018).                        ***
***                                                                         ***
***  Input:                                                                 ***
***    Tkin    - kinetic energy of DM particle [GeV]                        ***
***    rhoDeff - local DM density [GeV/cm**3] multiplied by the effective   ***
***              distance [kpc] out to which the DM flux is assumed to      ***
***              originate from:                                            ***
***                Deff = (\int dV \rho_DM \rho_CR /(4pi d^2)               ***
***                       / (\rho_DM \rho_CR)_local                         ***
***    CRtype  - Specifies which incident CR particle, see                  ***
***              src/cr_aux/dscrISRflux.f for details                       ***
***                                                                         ***
***  Output: d\Phi/dT, units: [cm**2 s GeV]^-1                              ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-06-22                                                         ***
***  mod 2019-03-10 (diff. scattering x-section now in dsddDMCRsigCR)       ***
***  mod tb, hk 2022-01-14 (added further elements)                         ***
***  mod 2025-05-26 edited for Co-SIMP CRDM                                 ***
*******************************************************************************
      real*8 function dsddDMCRflux(Tdm, rhoDeff)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h !
      include 'dsddcom.h'

      real*8 Tdm, Tcrmin, Tcrmax, Tdmcom, lnTcrmin, lnTcrmax, mcr, mdm
      real*8 norm, rhoDeff, intres, res
      integer CRtype, CRtypecom

      common/dsddDMCR_flux/ Tdmcom, mdm, CRtypecom  ! Pass to integrand_CRflux

c... functions
      external integrand_CRflux
      real*8 dsmwimp

      mdm = dsmwimp()
      Tdmcom = Tdm        ! required for passing Tdm to aux routine
      norm = kpc*1.d21*rhoDeff/mdm   ! cm/kpc ~ 3.086*1.d21
      dsddDMCRflux = 0.d0
      res = 0.d0

c... Loop over all CR types and sum up the contributions
      do CRtype = 1, NCRelements
        if (CR_inc(CRtype)) then
          mcr = mNCRau(CRtype)*atomicmassunit
          CRtypecom = CRtype  ! required for passing CRtype to aux routine

          Tcrmin = 1.d0/4.d0*((mdm - 2.d0*mcr + 2.d0*Tdm)
     &            + sqrt((mdm**2 + 4.d0*mdm*mcr + 4.d0*mcr**2 + 2.d0*mdm*Tdm)
     &            * (2.d0*mdm**2 + 5.d0*mdm*Tdm + 2.d0*Tdm**2)/(mdm*Tdm)))

          Tcrmax = 6.d1*Tcrmin  
          lnTcrmin = log(Tcrmin)
          lnTcrmax = log(Tcrmax)
          if (lnTcrmax.lt.3) lnTcrmax = 3.d0 
          call dgadap(lnTcrmin, lnTcrmax, integrand_CRflux, CRDM_acc/5.d0, intres)

          if (intres.ne.intres) then
            write(*,*) 'WARNING: dsddDMCRflux integral returned NaN for CRtype =', CRtype
            cycle
          endif
          if (intres.lt.0.d0) then
            write(*,*) 'WARNING: dsddDMCRflux integral was negative for CRtype =', CRtype
            cycle
          endif

          res = res + intres
        endif
      enddo

      dsddDMCRflux = norm*res
     &             * 2.d0     ! Factor of 2 from symmetry.
     &             * 1.d-20   ! Undo rescaling in integrand.

      return
      end

*******************************************************************************
*** auxiliary routine for integration
*******************************************************************************
      real*8 function integrand_CRflux(lnTcr)
      implicit none
      include 'dsddcom.h'
      real*8 lnTcr, Tcr, Tdmcom, sigCR, mdm
      real*8 dscrISRflux, dsddDMCRsigCRff, dsddDMCRsigCR
      integer CRtypecom
      common/dsddDMCR_flux/ Tdmcom, mdm, CRtypecom

      Tcr = exp(lnTcr)
      if (CRDM_form_factor) then
c... If we include form factors we must perform all the integrals numerically.
c... The above masses are the only ones which have been tabulated.
        sigCR = dsddDMCRsigCRff(Tdmcom, Tcr, CRtypecom)
      else
c... Otherwise we use the (semi-)analytical result.
        sigCR = dsddDMCRsigCR(Tdmcom, Tcr, CRtypecom)
      endif

      integrand_CRflux = dscrISRflux(Tcr, CRtypecom) ! dPhi/dTcr [cm^2 s GeV]^-1
     &                 * sigCR                       ! dsigma/dTdm [cm^2/GeV]
     &                 * Tcr                         ! From conversion to ln(Tcr).
     &                 * 1.d20

      return
      end