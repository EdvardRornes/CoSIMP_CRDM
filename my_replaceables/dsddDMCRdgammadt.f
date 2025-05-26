*******************************************************************************
*** Function dsddDMCRdgammadt returns the differential rate of recoils      ***
*** for DM hitting a nucleus at rest. Here, the DM flux is assumed to be    ***
*** the one that results from CR interactions with DM particles in the      ***
*** halo, see Bringmann & Pospelov (2018).                                  ***
***                                                                         ***
***  Input:                                                                 ***
***    TN       - recoil energy of nucleus [GeV]                            ***
***    mN       - mass of target nucleus [GeV]                              ***
***    depth    - penetration depth (detector location) [cm]                ***
***    rhoDeff  - local DM density [GeV/cm**3] multiplied by the effective  ***
***               distance [kpc] out to which the DM flux is assumed to     ***
***               originate from:                                           ***
***               Deff = (\int dV \rho_DM \rho_CR /(4pi d^2)                ***
***                       / (\rho_DM \rho_CR)_local                         ***
***                                                                         ***
***  Output: d\Gamma/dTN, units: [1/(s GeV)/nucleus]                        ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-06-24                                                         ***
*** mod  2019-03-13 removed scattering cross section as argument            ***
***                 (DMCR routines now use dsddDMCRsigtarget and            ***
***                  dsddsigmanucleon)                                      ***
*** mod  2019-10-04 changed argument zlfree -> depth                        ***
*** mod  2025-05-26 edited for Co-SIMP CRDM                                 ***
*******************************************************************************
      real*8 function dsddDMCRdgammadt(TN, mN, depth)
      implicit none
      include 'dsddcom.h'   ! CRDM_acc
      include 'dsnuclides.h' ! also includes dsmpconst.h

      real*8 TNcom, TN, Tdmmin, Tdmmax, lnTdmmin, lnTdmmax, mdm, mN, mNcom
      real*8 intres, mFe56, Tzmin, depth, depthcom
      common/dsddDMCR_dgammadt/ TNcom, mNcom, mdm, depthcom  ! Pass to aux integrand_dgammadt

c... functions
      external integrand_dgammadt
      real*8 dsmwimp

      dsddDMCRdgammadt = 1.d-60   ! Improve convergence in case of early return
      mdm = dsmwimp()
      TNcom = TN   ! Required for passing TN to aux routine
      mNcom = mN   ! Required for passing mN to aux routine
      depthcom = depth ! Required for passing depth to aux routine

c... What is done below is equivalent to switching between + and - sol
c... since x*sqrt(1/x^2)=sgn(x) so it accounts for both solutions
c... Note that when x switches sign is when we should switch solutions.
      if (CRDM_2loop) then
c... 2-loop elastic scattering
        Tzmin = TN/2.d0-mdm+sqrt((TN/2.d0-mdm)**2+TN/(2.d0*mN)*(mN+mdm)**2)
        mFe56 = 55.93494d0*atomicmassunit ! Use Fe-56 as heaviest abundant element
c... We assume that the 2->3 process occurs past this energy scale and exponentially
c... cut-off the Tdm corresponding to this range.
        CRDM_absorption_Ecut = 3.d0*mdm**2/(2.d0*mFe56)+mdm

        if (CRDM_absorption_Ecut.gt.1.d-1 ! Approximate inelastic effects
     &      .and.CRDM_attenuation         ! when including attenuation
     &      .and.attenuation_how.eq.2     ! but not the full treatment.
     &      .and.(.not.CRDM_inelastic)) CRDM_absorption_Ecut=1.d-1
        if (CRDM_attenuation) then
          call dsddTDMattenuation(Tdmmin, Tzmin, depth, 2) ! get Tmin from Tzmin
c... If the entire range of the integral is above the energies where we assume
c... that the 2->3 process will diminish the flux and/or the approximated
c... inelastic effect will remove all of the flux we simply return 0.
          if (Tdmmin.gt.1.2d0*CRDM_absorption_Ecut) return
        else
          Tdmmin = Tzmin
        endif
        ! Tdmmax = 6.d1*Tdmmin            ! This is a mild shortcut, but seems to
                                          ! affect the results slightly.
        Tdmmax = 1.2*CRDM_absorption_Ecut ! This is the physically accurate upper bound
        if (Tdmmin.ge.Tdmmax) return
      else
c... Use the 2->3 result
        Tdmmin = 1.d0/(4.d0*mN)*(2.d0*mN*(TN - 2.d0*mdm) + 3.d0*mdm**2
     &       + sqrt(1.d0/TN*(2.d0*mN + TN)*(2.d0*mN*TN + mdm**2)
     &       * (2.d0*mN*TN + 9.d0*mdm**2)))

        Tdmmax = 6.d1*Tdmmin  ! This should be inf, but the result does not 
                              ! change much with this choice. Increasing it 
                              ! simply makes this function slower.
      endif
      lnTdmmin = log(Tdmmin)
      lnTdmmax = log(Tdmmax)

      call dgadap(lnTdmmin, lnTdmmax, integrand_dgammadt, CRDM_acc/5.d0, intres)

      if (intres.ne.intres) then
        write(*,*) 'WARNING: dsddDMCRdgammadt returned NaN!'
        write(*,*) 'Returning dsddDMCRdgammadt = 0'
        return
      endif
      if (intres.lt.0.d0) then
        write(*,*) 'WARNING: dsddDMCRdgammadt integral was found to be negative!'
        write(*,*) 'Returning dsddDMCRdgammadt = 0'
        return
      endif

      dsddDMCRdgammadt = intres
     &                 * 1.d-30    ! undo rescaling in integrand
     
      return
      end


*******************************************************************************
*** auxiliary routine for integration                                       ***
*******************************************************************************
      real*8 function integrand_dgammadt(lnTdm)
      implicit none
      include 'dsddcom.h'   ! for access to Deff, rho0
      real*8 mdm, mNcom, Tdm, lnTdm, TNcom, Tz, depthcom, cut_smooth
      real*8 Co_SIMP_Flux, dsigtarget
      common/dsddDMCR_dgammadt/ TNcom, mNcom, mdm, depthcom
      integer mdmerr
c... functions
      real*8 dsddDMCRflux, dsddDMCRsigtarget, dsddDMCRflux_Tab
      real*8 dsddDMCRsigtarget_2loop, dsddDMCRflux_2loop

      Tdm = exp(lnTdm)

      integrand_dgammadt = 1.d-20 ! >0 to improve convergence in case of early return

      dsigtarget = 0.d0
c... Compute the differential cross section dsig/dTN
      if (CRDM_2loop) then
c... Use 2-loop elastic scattering cross section
        cut_smooth = 1.d0
        if (Tdm.gt.1.2*CRDM_absorption_Ecut) return
        if (Tdm.gt.CRDM_absorption_Ecut) then ! smoothen cut starting at CRDM_absorption_Ecut
c... This seems more realistic for Co-SIMPs due to the 
c... cross section being suppressed around the absorption cut
          cut_smooth = exp(-2d2*(Tdm/CRDM_absorption_Ecut-1d0)**2)
        endif
        if (CRDM_attenuation) then
          call dsddTDMattenuation(Tdm, Tz, depthcom, 1) ! get Tz from Tdm
          ! This is the Tdm which reaches the detector
          dsigtarget = cut_smooth*dsddDMCRsigtarget_2loop(TNcom, mNcom, Tz)
        else
          dsigtarget = cut_smooth*dsddDMCRsigtarget_2loop(TNcom, mNcom, Tdm)
        endif
      else
c... Use tree-level 2->3 cross section
        dsigtarget = dsddDMCRsigtarget(TNcom, mNcom, Tdm)
      endif

      Co_SIMP_Flux = 0.d0
      mdmerr = 1  ! dsddDMCRflux_Tab returns 0 if no error, 1 otherwise
c... Retrieve Co_SIMP_flux after dealing with attenuation 
      if (CRDM_form_factor.and.CRDM_tab) then
c... Use tabulated flux for form factor
        Co_SIMP_Flux = dsddDMCRflux_Tab(Tdm, mdmerr)
      endif
c... If tabulated flux returns an error, compute the flux numerically
      if (mdmerr.ne.0) then
        Co_SIMP_Flux = dsddDMCRflux(Tdm, rholocal*Deff)
      endif
c... Additionally, include the subdominant flux from the 2-loop elastic process
      if (CRDM_both) Co_SIMP_Flux = Co_SIMP_Flux + dsddDMCRflux_2loop(Tdm, rholocal*Deff)

      integrand_dgammadt = Tdm
     &                   * Co_SIMP_Flux
     &                   * dsigtarget
     &                   * 1.d30 ! rescale for better numerical result

      return
      end
