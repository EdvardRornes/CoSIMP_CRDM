*******************************************************************************
*** Provides the tabulated results for the DMCR flux when form factors are  ***
*** included. The results are read from a file and interpolated for the     ***
*** requested values.                                                       ***
***                                                                         ***
***  Input:                                                                 ***
***    Tkin    - kinetic energy of outgoing DM particle [GeV]               ***
***    mdmerr  - error flag, 0 if no error, 1 if mdm out of range           ***
***              (mdm < 10^-4 GeV or mdm > 10^1 GeV)                        ***
***                                                                         ***
***  Output: d\Phi/dT, units: [cm**2 s GeV]^-1                              ***
***                                                                         ***
*** author: edvard-rornes@hotmail.com                                       ***
*** date 2025-04-13                                                         ***
*******************************************************************************
      real*8 function dsddDMCRflux_Tab(Tdm, mdmerr)
      implicit none
      include 'dsddcom.h'
      include 'dsnuclides.h'
      real*8 Tdm
      real*8 dm_masses(100), flux(100,100), Tdm_values(100)  ! Flux indices: flux(Tdm, mdm)
      integer i_mlow, i_mhigh, i_Tlow, i_Thigh, mdmerr
      real*8 dsmwimp, mdm, flux_interp, alpha
      real*8 sip, sin, sdp, sdn, C
      integer ierr, dsidnumber
      integer idold
      data idold/-123456789/
      save idold, sip, sin, sdp, sdn
    
c... functions
      real*8 dsddDMCRinvert_cpsi

      mdm = dsmwimp()
      mdmerr = 0

      dsddDMCRflux_Tab = 0.d0

      if (mdm.lt.1.d-4 .or. mdm.gt.1.d1) then
          write(*,*) "dsddDMCRflux_Tab: mdm out of range:", mdm
          write(*,*) "Performing the numerical calculation of the flux instead"
          mdmerr = 1
        return
      endif

      if (idold.ne.dsidnumber()) then
         call dsddsigmanucleon(0.0d0, 0.0d0, sip, sin, sdp, sdn, ierr)
      endif

      call read_flux_data(dm_masses, flux, Tdm_values)
      
      call find_closest_indices(mdm, dm_masses, 100, i_mlow, i_mhigh)

      ! Extrapolation if Tdm is outside the range
      if (Tdm .lt. Tdm_values(1)) then
         alpha = (log(flux(2, i_mlow))-log(flux(1, i_mlow)))
     &         / (log(Tdm_values(2))-log(Tdm_values(1)))
         flux_interp = flux(1, i_mlow)*(Tdm/Tdm_values(1))**alpha
      else if (Tdm .gt. Tdm_values(100)) then
         alpha = (log(flux(100, i_mhigh))-log(flux(99, i_mhigh)))
     &         / (log(Tdm_values(100))-log(Tdm_values(99)))
         flux_interp = flux(100, i_mhigh)*(Tdm/Tdm_values(100))**alpha
      else
        call find_closest_indices(Tdm, Tdm_values, 100, i_Tlow, i_Thigh)

        call interpolate_flux_2D(Tdm, mdm, dm_masses, Tdm_values, flux,
     &                         i_mlow, i_mhigh, i_Tlow, i_Thigh, flux_interp)
      endif
      if (flux_interp.ne.flux_interp) return


      if (CRDM_2loop) then
        C = dsddDMCRinvert_cpsi()  ! This function provides 2->3 constant C based 
                                   ! on the sigmaNR the program is currently using
      else
        C = (sip+sin)/2.d0
        if (C.eq.0.d0) C = sdp
      endif

c... We multiply by 1.d20 to remove the normalization used when computing the flux
      dsddDMCRflux_Tab = flux_interp
     &                 * C
     &                 * 1.d20


      mdmerr = 0

      return
      end


      subroutine find_closest_indices(value, array, n, i_low, i_high)
      implicit none
      integer n, i_low, i_high, k
      real*8 value, array(n)

      if (value < array(1) .or. value > array(n)) then
         write(*,*) "Value outside array bounds:", value, array(1), array(n)
         stop
      endif

      do k = 1, n-1
         if (value >= array(k) .and. value <= array(k+1)) then
            i_low = k
            i_high = k+1
            return
         endif
      end do
      end


      subroutine interpolate_flux_2D(Tdm, mdm, dm_masses, Tdm_values, flux,
     &                               i_mlow, i_mhigh, i_Tlow, i_Thigh, flux_interp)
      implicit none
      real*8 Tdm, mdm, dm_masses(100), Tdm_values(100), flux(100,100)
      integer i_mlow, i_mhigh, i_Tlow, i_Thigh
      real*8 flux_interp
      real*8 m1, m2, T1, T2, x, y
      real*8 f11, f21, f12, f22, logf

      ! If both mdm and Tdm match exactly in the table, return the flux
      if (mdm.eq.dm_masses(i_mlow) .and. Tdm.eq.Tdm_values(i_Tlow)) then
        flux_interp = flux(i_Tlow, i_mlow)
        return
      endif

      ! If mdm matches exactly, do 1D interpolation for Tdm
      if (mdm.eq.dm_masses(i_mlow)) then
        ! Interpolate 1D along Tdm (log-linear)
        T1 = Tdm_values(i_Tlow)
        T2 = Tdm_values(i_Thigh)
        y = (log(Tdm) - log(T1)) / (log(T2) - log(T1))
        
        f11 = log(flux(i_Tlow, i_mlow))
        f12 = log(flux(i_Thigh, i_mlow))

        logf = (1.d0 - y)*f11 + y*f12
        flux_interp = exp(logf)
        return
      endif

      ! If Tdm matches exactly, do 1D interpolation for mdm
      if (Tdm.eq.Tdm_values(i_Tlow)) then
        m1 = dm_masses(i_mlow)
        m2 = dm_masses(i_mhigh)
        x = (log(mdm) - log(m1)) / (log(m2) - log(m1))
        
        f11 = log(flux(i_Tlow, i_mlow))
        f21 = log(flux(i_Tlow, i_mhigh))

        logf = (1.d0-x)*f11 + x*f21
        flux_interp = exp(logf)
        return
      endif

      ! If neither mdm nor Tdm match exactly, perform 2D interpolation
      m1 = dm_masses(i_mlow)
      m2 = dm_masses(i_mhigh)
      T1 = Tdm_values(i_Tlow)
      T2 = Tdm_values(i_Thigh)

      ! Logarithmic interpolation for mdm and Tdm
      x = (log(mdm) - log(m1)) / (log(m2) - log(m1))
      y = (log(Tdm) - log(T1)) / (log(T2) - log(T1))

      ! Interpolate flux values based on log(mdm) and log(Tdm)
      f11 = log(flux(i_Tlow, i_mlow))
      f12 = log(flux(i_Tlow, i_mhigh))
      f21 = log(flux(i_Thigh, i_mlow))
      f22 = log(flux(i_Thigh, i_mhigh))

      logf = (1.d0-x)*(1.d0-y)*f11 + x*(1.d0-y)*f12 + (1.d0-x)*y*f21 + x*y*f22

      flux_interp = exp(logf)

      end


      subroutine read_flux_data(dm_masses, flux, Tdm_values)
      implicit none
      include 'dsddcom.h'
      real*8 dm_masses(100), flux(100,100), Tdm_values(100)
      integer i, j, ios

      if (CRDM_cs) then
         open(unit=11, file='data/Co-SIMP_flux_ff_cs.dat', status='old', action='read')
      else
         open(unit=11, file='data/Co-SIMP_flux_ff.dat', status='old', action='read')
      endif

      read(11, *, iostat=ios) dm_masses
      if (ios.ne.0) then
         write(*,*) "Error reading DM masses", dm_masses
         stop
      end if

      do i=1, 100
         read(11, *, iostat=ios) Tdm_values(i), (flux(i, j), j=1, 100)
         if (ios.ne.0) then
            write(*,*) "Error reading flux and Tdm values for line", i
            stop
         end if
      end do

      close(11)
      return
      end
