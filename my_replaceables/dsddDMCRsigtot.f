*******************************************************************************
*** This program calculates the total cross section for a 2->3 process with ***
*** a proton as the target. The user must make a few changes to switch      ***
*** between using kinetic energy and CM energy squared.                     ***
***                                                                         ***
***  Input (only 1 should be used):                                         ***
***    Tdm     - kinetic energy of proton in the lab frame [GeV]            ***
***              (recoil energy after scattering)                           ***
***    s       - CM energy [GeV^2]                                          ***
***                                                                         ***
***  Output:  \sigma_{\chi p}, units: cm**2                                 ***
***           Note: These units depend on how sigsip and sigsin in          ***
***                 dsgivemodel_generic_wimp are setup.                     ***
***                 When used with DDCR_limits cm^2 is default.             ***
***                                                                         ***
*** author: edvard-rornes@hotmail.com                                       ***
*** date 2025-05-23                                                         ***
*******************************************************************************

      ! real*8 function dsddDMCRsigtot(Tkin)
      real*8 function dsddDMCRsigtot(s)
      implicit none
      include 'dsmpconst.h' ! For m_p
      include 'dsddcom.h'
      real*8 mN, mdm, Tdm, TNmin, TNmax, intres, s, Tkin, Term1, Term2
      common/dsddDMCRsigtot_T/ Tdm, mN, mdm ! Pass to integrand_sigtot

c... functions
      external integrand_sigtot
      real*8 dsmwimp

      dsddDMCRsigtot = 0.d0


      mdm = dsmwimp()
      mN = m_p    ! use proton mass instead to compare to reported limits

      if (mdm.lt.1.d-50) return

      ! Tdm = Tkin                         ! Needed when Tkin is the argument
      ! s = (mN+mdm)**2 + 2.d0*mN*Tdm      ! Needed when Tkin is the argument
      Tdm = (s-(mN+mdm)**2)/(2.d0*mN)    ! Needed when s is the argument

      if (Tdm.le.mdm+3.d0*mdm**2/(2.d0*mN)) then
    !     write(*,*) 'WARNING: Input energy in dsddDMCRsigtot is ',
    !  &             'too low to produce another DM particle!'
    !     write(*,*) 'Returning sigtot = 0.'
        return
      endif
      
      Term1 = (s + mN**2 - mdm**2)*(s + mN**2 - 4.d0*mdm**2)/(4.d0*s*mN) - mN
      Term2 = sqrt((s**2 - 2.d0*s*(mN**2 + mdm**2)
     &      + (mN**2 - mdm**2)**2)*(s**2 - 2.d0*s*(mN**2 + 4.d0*mdm**2)
     &      + (mN**2 - 4.d0*mdm**2)**2))/(4.d0*s*mN)

      TNmin = Term1 - Term2
      TNmax = Term1 + Term2
      
      if (TNmin.lt.0.d0) then
      !   write(*,*) 'WARNING: TNmin < 0 in dsddDMCRsigtot! TNmin =', TNmin
      !   write(*,*) 'Setting TNmin = 0'
        TNmin = 1.d-5
      endif
      
c... Doing this to TNmin is problematic as its usually zero (or slightly negative).
c TNmin = log(TNmin)
c TNmax = log(TNmax)

      call dgadap(TNmin, TNmax, integrand_sigtot, CRDM_acc/5.d0, intres)

      if (intres.ne.intres) then
        write(*,*) 'WARNING: dsddDMCRsigtot integral returned NaN!'//
     &             ' Returning dsddDMCRsigtot = 0'
        write(*,*) 's=',s,'mdm=',mdm,'mN=',mN,'TNmin=',TNmin,'TNmax=',TNmax
        return
      endif
      if (intres.lt.0.d0) then
        write(*,*) 'WARNING: dsddDMCRsigtot integral < 0!'//
     &             ' Returning dsddDMCRsigtot = 0'
        return
      endif
      
      dsddDMCRsigtot = intres
     &               * 1.d-30   ! Undo normalization in integrand

      return
      end


c Define a helper function for the integrand
      real*8 function integrand_sigtot(TN)
      implicit none
      real*8 dsddDMCRsigtarget
      real*8 TN, Tdm, mdm, mN!, lnTN
      common/dsddDMCRsigtot_T/ Tdm, mN, mdm

      integrand_sigtot = dsddDMCRsigtarget(TN, mN, Tdm)
     &                 * 1.d30  ! For better numerical result      

      return
      end

