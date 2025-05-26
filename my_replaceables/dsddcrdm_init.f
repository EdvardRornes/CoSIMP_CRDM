      subroutine dsddcrdm_init
c...  initialize the crdm routines
      implicit none
      include 'dsddcom.h'
      include 'dsnuclides.h' ! also includes dsmpconst.h
      integer CRtype

c... astro parameters, so far only used in DDCR routines
      rholocal = 0.3d0 ! GeV/cm^3; eventually, this should be extracted from halo model instead 
      vlocal = 2./sqrt(pi)*220*1.d5 ! mean local DM velocity in cm/s
      Deff = 5.d0 ! kpc; eff distance out to which source term is integrated
                  ! Deff = 0.997 (8.02) kpc corresponds to 1 (10) kpc in real
                  ! distance
                  ! (new) default, based on 2111.05559

c... further settings for CRDM routines
      attenuation_how = 2 ! 1: use analytic expressions from 1810.10543,
                          !    assuming a constant scattering cross section 
                          ! 2: [default] solve differential equation for soil attenuation
                          !    numerically (important for q-dependent scattering)
      
      do CRtype=1,NCRelements    ! H. Kolesova: choose which elements to include in CR flux
        CR_inc(CRtype) = .true.  ! Default: Include all CR types
      enddo
c... Set the above to false and uncomment the lines below to only include the given CR species
      ! CR_inc(1) = .true.  ! include protons in CR flux
      ! CR_inc(2) = .true.  ! include He in CR flux
      ! CR_inc(3) = .true.  ! include C in CR flux
      ! CR_inc(4) = .true.  ! include O in CR flux
      CRDM_form_factor = .true.  ! include nuclear form factors
      CRDM_cs = .true.            ! Set Msq=C/s, if false then Msq=const
                                  ! Note that Msq=const violates unitarity!
      CRDM_high_acc = .false.     ! use high accuracy computing the 2->3 cross section
      CRDM_tab = .true.           ! use tabulated cross section
                                  ! (only available for 1.d-4<mdm<1.d1 GeV)
      CRDM_2loop = .true.         ! use 2-loop cross section in the detector
      CRDM_EnergyLoss = .false.   ! DO NOT CHANGE THIS! This is only used to calculate 
                                  ! the average energy loss of the outgoing DM particle
      CRDM_attenuation = .false.   ! include attenuation of the DM flux in the soil
      CRDM_inelastic = .true.     ! include inelastic scattering in CRDM soil absorption
      CRDM_both = .true.          ! Additionally includes 2-loop contribution to the CRDM flux.
                                  ! Recommended to set to false when CRDM_2loop is false.
      CRDM_absorption_Ecut = 1.d1 ! Initial DM energy [GeV] above which the TOA CRDM flux
                                  ! is assumed to be fully absorbed
      CRDM_inel_Eref = 9.5d0 ! reference energy [GeV] beyond which we extrapolate tabulated
                             ! inelastic scattering results by simple rescaling
                             ! results should be insensitive to this value in the range ~7-9.99
      CRDM_acc = 5d-3 ! accuracy of final integration in DMCR rate (@ detector)
      n_att1   = 40  ! # tabulation points for soil attenuation
      n_att2   = 40  ! # additional tabulation points for soil attenuation in
                     !   region of large gradients
                     !   default: (40,40); increase for higher accuracy (n_att_max=3000)
                     !   for strong Q2-dependence, start to test increasing n_att1 !

c..   make sure that settings above are not inconsistent...
      n_att = n_att1+n_att2
      if (n_att.gt.n_att_max) then
        n_att  = n_att_max/5
        n_att2 = n_att_max - n_att -1
        n_att  = n_att+n_att2
      endif
      if (CRDM_acc.lt.1.d-5) CRDM_acc=1.d-5
      if (CRDM_EnergyLoss) CRDM_EnergyLoss = .false. ! I REPEAT, DO NOT CHANGE THIS !


      return
      end
