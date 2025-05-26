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
***                                                                         ***
***  Output: d\Phi/dT, units: [cm**2 s GeV]^-1                              ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-06-22                                                         ***
***  mod 2019-03-10 (diff. scattering x-section now in dsddDMCRsigCR)       ***
***  mod tb, hk 2022-01-14 (added further elements)                         ***
***  mod 2025-05-26 edited for 2-loop Co-SIMP CRDM                          ***
*******************************************************************************
      real*8 function dsddDMCRflux_2loop(Tkin, rhoDeff)
      implicit none
      include 'dsnuclides.h' ! also includes dsmpconst.h !
      include 'dsddcom.h'

      real*8 Tkin, rhoDeff
      real*8 norm, Tcrmin, Tcrmax, intres, res
      real*8 mcr, mdm, Tdm
      integer CRtype
      common/DMCRflux/ Tdm, CRtype

c... TB IMPROVE: Tabulate this (once per model), for code speedup !

c... functions
      real*8 dsmwimp
      external dsddDMCRflux_2loop_aux
      
      dsddDMCRflux_2loop=0.0d0

c... debug
c      if (Tkin.gt.0.1d0) return
      
      mdm  = dsmwimp()
      Tdm  = Tkin
      norm = kpc*1.0d21*rhoDeff/mdm
      res  = 0.0d0

      do CRtype = 1, NCRelements
        if (CR_inc(CRtype)) then
          mcr = mNCRau(CRtype)*atomicmassunit
          Tcrmin = sqrt((mcr-Tkin/2.)**2 + (mcr+mdm)**2*Tkin/2./mdm) - mcr + Tkin/2.
          Tcrmin = log(Tcrmin) ! change to log to better sample fast falling integrand
          Tcrmax = Tcrmin+4.0d0
          if (Tcrmax.lt.3) Tcrmax = 3.0 ! make sure to sample peak of spectrum
          call dgadap(Tcrmin,Tcrmax,dsddDMCRflux_2loop_aux,CRDM_acc/5.,intres)
          res = res + intres
        endif
      enddo
      
      dsddDMCRflux_2loop = res*norm*1.0d-30 ! undo spurious normalization 
                                            ! from integration
      return
      end


*******************************************************************************
*** auxiliary routine just for integration
*******************************************************************************
      real*8 function dsddDMCRflux_2loop_aux(Tcr)
      implicit none
      real*8 Tcr, Tdm
      real*8 dscrISRflux, dsddDMCRsigCR_2loop
      integer CRtype
      common/DMCRflux/ Tdm, CRtype

c... add large normalization to improve integral conversion   
      dsddDMCRflux_2loop_aux = 1.0d30*dscrISRflux(exp(Tcr), CRtype)
     &                   *exp(Tcr) ! from conversion to log(Tcr) 
     &                   *dsddDMCRsigCR_2loop(Tdm,exp(Tcr),CRtype)     

      return
      end
