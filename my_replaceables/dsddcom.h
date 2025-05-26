*     -*- mode: fortran -*-
      integer ddng,ddnsf,ddnme
      parameter (ddng=27)
      parameter (ddnsf=6)
      parameter (ddnme=9)
! sf : 1_m, 2_sigma, 3_delta, 4_phi, 5_deltasigma, 6_mphi
      character*16 ddsf(6)
! me : 1_s, 2_v, 3_vm, 4_t, 5_a, 6_am, 7_p, 8_t2e, 9_t2o
      character*16 ddme(9)
      real*8 fme(9,2,4) ! me ;  p,n ; u,d,s,g
      common /ddcom/
     &     fme,
     &     ddsf,
     &     ddme
      save /ddcom/

! added by TB -- quenching for Borexino etc
      real*8 lntqdat(300),lnTrecdat(300),lndTdTdat(300)
      integer nqdat,khi,klo,quenchhow
      logical quenching_set
      common /ddquench/ lntqdat,lnTrecdat,lndTdTdat,
     &                  nqdat,khi,klo,quenchhow,quenching_set
      save /ddquench/

! added by TB -- astro parameters, so far only used by dsddCR routines
      real*8 rholocal, vlocal, Deff
      common /DDCRastro/ rholocal, vlocal, Deff
      save /DDCRastro/

c... needed by CRDM routines
      real*8 CRDM_acc, CRDM_absorption_Ecut, CRDM_inel_Eref  ! set in dsddinit
      real*8 lnT0min,lnT0max,lnTzmin,lnTzmax                 ! min/max (+ number of) DM energies
                                                             ! to tabulate for attenuation
      integer n_att, n_att1, n_att2, n_att_max
      parameter (n_att_max = 3000)
      parameter (lnT0min=log(1.d-8), lnT0max=log(1.d11))
      parameter (lnTzmin=log(1.d-8), lnTzmax=log(1.d11))
      integer targetoption
      integer attenuation_how, exp_loc
! ER -- included additional options (CRDM_form_factor, CRDM_cs, CRDM_EnergyLoss)
!       for Co-SIMP CRDM
      logical CRDM_inelastic, CRDM_form_factor, CRDM_cs, CRDM_EnergyLoss,
     &        CRDM_high_acc, CRDM_tab, CRDM_2loop, CRDM_attenuation, CRDM_both
      common /crdmtargets/ CRDM_acc, CRDM_absorption_Ecut, CRDM_inel_Eref, n_att, n_att1, n_att2,
     &                     targetoption, attenuation_how, exp_loc, CRDM_inelastic, CRDM_form_factor,
     &                     CRDM_cs, CRDM_EnergyLoss, CRDM_high_acc, CRDM_tab, CRDM_2loop, CRDM_attenuation, CRDM_both
      save /crdmtargets/

c... tabulated 'inelastic neutrino cross sections', c.f. Alvey+ (2022) and data/sigma_inel_*.txt
c...     grid spacing in mDM and tdm everywhere
      integer nel_inel, nch_inel, ntdmmax_inel, nomegamax_inel
      parameter (nel_inel       = 7,   ! number of simulated elements
     &           nch_inel       = 5,   ! number of simulated channels (QE,delta,res,DIS,tot)
     &           ntdmmax_inel   = 200, ! maximal number of tabulated tdm
     &           nomegamax_inel = 200  ! maximal number of tabulated omega-values (per tdm)
     &          )
      integer ntdm_inel(nel_inel,nch_inel), nomega_inel(ntdmmax_inel,nel_inel,nch_inel),
     &        nomegamax_el_inel(nel_inel,nch_inel) ! actual # of tab. bins
      real*8 sig_ineltab(ntdmmax_inel,nomegamax_inel,nel_inel,nch_inel)
      real*8 tdmtab_inel(ntdmmax_inel,nel_inel,nch_inel),
     &       omegatab_inel(ntdmmax_inel,nomegamax_inel,nel_inel,nch_inel)
      common /nu_inel/ sig_ineltab,omegatab_inel,tdmtab_inel,
     &                 nomega_inel,nomegamax_el_inel,ntdm_inel
      save /nu_inel/


! obsolescent
      real*8 ftp(7:12),ftn(7:12),delu,deld,dels
      common /ddcomlegacy/
     &     ftp,ftn,delu,deld,dels
      save /ddcomlegacy/
