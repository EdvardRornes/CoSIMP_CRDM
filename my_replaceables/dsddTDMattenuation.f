*******************************************************************************
*** Subroutine dsddTDMattenuation relates average DM kinetic energies       *** 
*** before and after propagating through a dense medium                     ***
***                                                                         ***
***  Input:                                                                 ***
***    Tin    - initial kinetic energy  [GeV]                               ***
***    Tz     - average kinetic energy after scattering in medium [GeV]     ***
***    depth  - penetration depth (detector location) [cm]                  ***
***    how    - convert from Tin to Tz (how=1)                              ***
***             or from Tz to Tin (how=2)                                   ***
***                                                                         ***
***  Output:                                                                ***
***    Tin  - initial kinetic energy  [GeV]                                 ***
***    Tz    - average kinetic energy after scattering in medium [GeV]      ***
***                                                                         ***
***  NB: i) For how=1, the input value of Tz is overwritten on output       ***
***         For how=2, the input value of Tin is overwritten on output      ***
***      ii) If this routine is called directly, you first need to set      ***
***          'targetoption' (typically set in dsddDMCRcountrate) for the    ***
***          correct experimental location                                  ***
***                                                                         ***
*** WARNING: While this now should work for 'arbitrary' Q-dependence, it    ***
*** has only been tested for single mediators that are neither too light    ***
*** (<<MeV) nor too heavy (>>10 GeV). Since the numerical integration is    ***
*** potentially unstable, the behaviour of this routine must be carefully   ***
*** tested outside this regime and/or for more complicated energy           ***
*** dependendences.                                                         ***
***                                                                         ***
*** The default behaviour (full Q2-dependence, including inelastic          ***
*** scattering) can be changed by setting the common block flag             ***
*** 'attenuation_how' to 1 in dsdd_init. In that case, dsddTDMattenuation   ***
*** will NOT take into account a possible Q-dependence of the scattering    ***
*** cross section. Effectively assuming the cross section at zero momentum  ***
*** transfer, the code runs much faster in this case -- but typically       ***
*** over-estimates the stopping power for e.g. light mediators.             ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-07-07                                                         ***
*** mod tb, 31/10/2019: added full Q2 dependence; argument zlfree -> depth  ***
*** mod tb, 12/01/2021: added target location dependence                    ***
*******************************************************************************
      subroutine dsddTDMattenuation(Tin, Tz, depth, how)
      implicit none
      include 'dsddcom.h'
      include 'dsnuclides.h'

      real*8 Tin, Tz, depth
      integer how, ierr
      real*8 dsmwimp, mdm, zlfree, sip,sin,sdp,sdn, sigsi
      save zlfree

      real*8 dsddlfreesimp, dsgetTzratio

ccc to make sure that we only tabulate once per model / constraint
      integer dsidnumber
      integer idold
      real*8 depthold
      data idold/-123456789/
      data depthold/1.d100/
      save idold
      save depthold

      mdm = dsmwimp()
      if (how.ne.1.and.how.ne.2) goto 20
      if (idold.ne.dsidnumber().or.depth.ne.depthold) then ! new model / constraint
        if (attenuation_how.eq.2) then ! fully take into account Q2-dependence
          call dstabulateTinTz(depth)
        else ! take the Q2->0 limit and assume only elastic scattering independent of Q2
          call dsddsigmanucleon(0.0d0,0.0d0,sip,sin,sdp,sdn,ierr)
          if ((targetoption/10).ne.2) then ! assume that DM scatters on n and p
            sigsi = (sip+sin)/2.           ! take *average* for simplified mean free path       
            zlfree = depth/dsddlfreesimp(sigsi, depth, 1)
            ! write(*,*) 'Free length [m]: ', dsddlfreesimp(sigsi, depth, 1)/1.d2
          else ! Borexino SD limits -> only scattering on protons 
            zlfree = depth/dsddlfreesimp(sdp, depth, 2) ! NB: using sigsi->sdp assumes 
                                                        ! an unrealistically high
                                                        ! stopping power!
            ! write(*,*) sdp, dsddlfreesimp(sdp, depth, 2), zlfree
            if (targetoption.eq.21) zlfree = 1.0d-5*zlfree ! In reality, spin-dependent 
                                                           ! scattering should have 
                                                           ! 'MUCH' larger mean free path
          endif
        endif
        idold=dsidnumber() 
        depthold=depth
      endif

c... this typically indicates that exp() has taken too large arguments
c... (but should be captured now)
      if (Tin.ne.Tin.or.Tz.ne.Tz) then
        write(*,*) 'WARNING in dsddTDMattenuation: Tin, Tz = ',Tin, Tz
      endif

      if (how.eq.1) then ! convert Tin to Tz
        if (Tin.lt.0.0d0) goto 10
        if (attenuation_how.eq.2) then 
          Tz = Tin*dsgetTzratio(Tin,1)
          ! write(*,*) 'Tz = ',Tz, dsgetTzratio(Tin,1)
        else
          if (zlfree.gt.30.0d0) then ! nothing arrives that far...
             Tz = 1.0d-50
          else
             Tz = 2.*mdm*Tin / ((2.*mdm+Tin)*exp(zlfree) - Tin)
          endif 
        endif          
      elseif (how.eq.2) then ! convert Tz to Tin
        if (Tz.lt.0.0d0) goto 10
        if (attenuation_how.eq.2) then 
          Tin = Tz*dsgetTzratio(Tz,2)        
        else
          if (zlfree.gt.30.0d0) then ! nothing arrives that far...
            Tin = 1.0d25 ! 'infinitely' high initial energy needed
          elseif (Tz.gt.(2.*mdm/(exp(zlfree)-1.0d0))) then
            Tin = 1.0d25 ! 'infinitely' high initial energy needed
          else
            Tin = 2.*mdm*Tz / ((2.*mdm+Tz)/exp(zlfree) - Tz)
          endif
        endif          
      endif      

      return

 10   write(*,*) 'FATAL ERROR in dsddTDMattenuation:'
      write(*,*) 'You have supplied a non-positive input energy!'
      stop

 20   write(*,*) 'FATAL ERROR in dsddTDMattenuation:'
      write(*,*) 'unknown option how = ',how
      stop 
      end


******************************************************************************
*** auxiliary routines to tabulate and read out conversion between Tin and Tz
******************************************************************************

ccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dstabulateTinTz(depth)
c... tabulates ratios (T0,Tz/T0) and (Tz,Tz/T0)      
ccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'dsddcom.h'
      real*8 depth
      
      real*8 lnT0dm, lnTzdm, zi, zf, stepguess,stepmin, logslope
      real*8 dTdzrescale, lnTzmax2, lnTzmin2
      parameter (dTdzrescale=1d2) ! NB: must be the same as further down!
      integer ier, i, imax,imin
      real*8 LogTinTzdat(n_att_max,2,2)
      common /soilcom/ LogTinTzdat
      save /soilcom/
      external dsdlogTdxrhs
      logical dsisnan

      imax = 1
      imin = n_att1

      do i = n_att1,1,-1
c... sampling in Tz is much better than sampling in T0..
        lnTzdm = lnTzmin + (i-1.)/(1.*n_att1-1.)*(lnTzmax-lnTzmin)
        zf = 1.0d0 ! integrate up to 1cm depth to avoid zero density at surface
        zi = depth*dTdzrescale
        stepguess = zi/2d1
        stepmin = CRDM_acc*depth*1d-3 ! TB debug 1d-5
        lnT0dm = lnTzdm
        ! NB: we need a rather high accuracy here because of large derivatives
        call dskdinty(lnT0dm,zi,zf,CRDM_acc/1d1,stepguess,stepmin,dsdlogTdxrhs,ier)

        if (lnT0dm.lt.lnTzdm) lnT0dm=lnTzdm
        if (lnT0dm.gt.7d1.or.ier.ne.0.or.dsisnan(lnT0dm)) then
c          if (ier.ne.0) write(*,*) 'error', ier, lnTzdm,lnT0dm, zi,zf,depth
          lnT0dm = 7d1
        endif
        LogTinTzdat(i,1,1) = lnT0dm   ! T0 -> Tz/T0
        LogTinTzdat(i,1,2) = lnTzdm - lnT0dm
        LogTinTzdat(i,2,1) = lnTzdm   ! Tz -> T0/Tz
        LogTinTzdat(i,2,2) = lnT0dm - lnTzdm
c... now replace T0 values that look wrong with previous, well-behaved T0
c... (-1d-5 to avoid infinite derivatives):
        if (i.lt.n_att1) then
          logslope=(LogTinTzdat(i+1,2,1)-LogTinTzdat(i,2,1))/
     &             (LogTinTzdat(i+1,1,1)-LogTinTzdat(i,1,1))
c          write(*,*) logslope, LogTinTzdat(i+1,1,1), LogTinTzdat(i,1,1)
          if (.not.(logslope.gt.0d0)              ! dTz/dT0 < 0, incl. NaN/inf
     &        .or.(1d0/logslope.gt.1d4)) then     ! dT0/dTz >>> 0
            LogTinTzdat(i,1,1)=LogTinTzdat(i+1,1,1) - 1d-5
            LogTinTzdat(i,1,2)=LogTinTzdat(i,2,1)-LogTinTzdat(i+1,1,1) + 1d-5
            LogTinTzdat(i,2,2)=LogTinTzdat(i+1,1,1)-LogTinTzdat(i,2,1) - 1d-5
          else
            if (logslope.lt.0.95) then ! we want a better resolution between imin and imax
                                       ! IMPROVE (speed): Actually, we probably only need a good resolution
                                       ! for something like 0.95>loglope>0.5
                                       ! -> loop over all occurences like this !?
              if (i.gt.imax) imax = i
              if (i.lt.imin) imin = i
            endif
          endif
        endif
      enddo
c      write(*,*)

c      goto 100

      lnTzmax2 = LogTinTzdat(imax+1,2,1)
      lnTzmin2 = LogTinTzdat(imin,2,1)
      n_att    = n_att1+n_att2-imax+imin
      do i = n_att,n_att2+imin+1,-1
        LogTinTzdat(i,1,1) = LogTinTzdat(i-n_att2+imax-imin,1,1)
        LogTinTzdat(i,1,2) = LogTinTzdat(i-n_att2+imax-imin,1,2)
        LogTinTzdat(i,2,1) = LogTinTzdat(i-n_att2+imax-imin,2,1)
        LogTinTzdat(i,2,2) = LogTinTzdat(i-n_att2+imax-imin,2,2)
      enddo
      do i = n_att2,1,-1
        lnTzdm = lnTzmin2 + (i-1.)/(1.*n_att2-1.)*(lnTzmax2-lnTzmin2)
        zf = 1.0d0
        zi = depth*dTdzrescale
        stepguess = zi/2d1
        stepmin = CRDM_acc*depth*1d-3
        lnT0dm = lnTzdm
        call dskdinty(lnT0dm,zi,zf,CRDM_acc/1d1,stepguess,stepmin,dsdlogTdxrhs,ier)
        if (lnT0dm.lt.lnTzdm) lnT0dm=lnTzdm
        if (lnT0dm.gt.7d1.or.ier.ne.0.or.dsisnan(lnT0dm)) then
          lnT0dm = 7d1
        endif
        LogTinTzdat(i+imin,1,1) = lnT0dm   ! T0 -> Tz/T0
        LogTinTzdat(i+imin,1,2) = lnTzdm - lnT0dm
        LogTinTzdat(i+imin,2,1) = lnTzdm   ! Tz -> T0/Tz
        LogTinTzdat(i+imin,2,2) = lnT0dm - lnTzdm
        logslope=(LogTinTzdat(i+imin+1,2,1)-LogTinTzdat(i+imin,2,1))/
     &           (LogTinTzdat(i+imin+1,1,1)-LogTinTzdat(i+imin,1,1))
        if (.not.(logslope.gt.0d0)              ! dTz/dT0 < 0, incl. NaN/inf
     &      .or.(1d0/logslope.gt.1d4)) then     ! dT0/dTz >>> 0
          LogTinTzdat(i+imin,1,1)=LogTinTzdat(i+imin+1,1,1) - 1d-5
          LogTinTzdat(i+imin,1,2)=LogTinTzdat(i+imin,2,1)-LogTinTzdat(i+imin+1,1,1) + 1d-5
          LogTinTzdat(i+imin,2,2)=LogTinTzdat(i+imin+1,1,1)-LogTinTzdat(i+imin,2,1) - 1d-5
        endif
      enddo
      
c      do i = n_att,1,-1
c        write(*,*) i, exp(LogTinTzdat(i,1,1)), exp(LogTinTzdat(i,2,1))
c      enddo
c      write(*,*) imax,imin,n_att1,n_att2
c      write(*,*)
c      read(*,*)

      goto 100
 100  continue !    read(*,*)

      end ! dstabulateTinTz
  
      
ccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function dsgetTzratio(T,how) 
ccccccccccccccccccccccccccccccccccccccccccccc
c interpolate. assuming that soilcom is ordered
      implicit none
      include 'dsddcom.h'
    
      real*8 T, lnt
      integer how  ! 1: T=T0->Tz; 2: T=Tz->T0
      
      integer k, dk, klo_att(2), khi_att(2)
      real*8 res, LogTinTzdat(n_att_max,2,2)
      common /soilcom/ LogTinTzdat
      save /soilcom/
      logical initialized
      data initialized /.false./
      save initialized, klo_att, khi_att

      lnt = log(T)  
      res = 0.0d0
      if (.not.initialized) then
        klo_att(1) = 1
        klo_att(2) = 1
        khi_att(1) = n_att
        khi_att(2) = n_att
        initialized = .true.
      endif

      if (lnt.le.LogTinTzdat(1,how,1)) then
         klo_att(how) = 1
         res = LogTinTzdat(1,how,2)
      elseif (lnt.ge.LogTinTzdat(n_att,how,1)) then
         khi_att(how) = n_att
         res = LogTinTzdat(n_att,how,2)
      else
        dk=n_att/40+1
 100    if (LogTinTzdat(khi_att(how),how,1).lt.lnt) then
           khi_att(how) = khi_att(how) + dk
           dk = 2*dk
           if (khi_att(how).lt.n_att) goto 100
           khi_att(how)=n_att
        endif
        dk=n_att/40+1
 110    if (LogTinTzdat(klo_att(how),how,1).gt.lnt) then
           klo_att(how) = klo_att(how) - dk
           dk = 2*dk
           if (klo_att(how).gt.1) goto 110
           klo_att(how)=1
        endif
 120    if (khi_att(how)-klo_att(how).gt.1) then
           k=(khi_att(how)+klo_att(how))/2
           if (LogTinTzdat(k,how,1).gt.lnt) then
              khi_att(how)=k
           else
              klo_att(how)=k
           endif
           goto 120
        endif
        
c... quad interpolation -> no better...
c        if (khi_att(how).lt.n_att) then
c           x(1) = LogTinTzdat(klo_att(how),how,1)
c           x(2) = LogTinTzdat(khi_att(how),how,1)
c           x(3) = LogTinTzdat(khi_att(how)+1,how,1)
c           y(1) = LogTinTzdat(klo_att(how),how,2)
c           y(2) = LogTinTzdat(khi_att(how),how,2)
c           y(3) = LogTinTzdat(khi_att(how)+1,how,2)
c        else
c           x(1) = LogTinTzdat(klo_att(how)-1,how,1)
c           x(2) = LogTinTzdat(klo_att(how),how,1)
c           x(3) = LogTinTzdat(khi_att(how),how,1)
c           y(1) = LogTinTzdat(klo_att(how)-1,how,2)
c           y(2) = LogTinTzdat(klo_att(how),how,2)
c           y(3) = LogTinTzdat(khi_att(how),how,2)
c        endif
 
c        res =   y(1)*(lnt-x(2))*(lnt-x(3))/(x(1)-x(2))/(x(1)-x(3))
c     &        + y(2)*(lnt-x(1))*(lnt-x(3))/(x(2)-x(1))/(x(2)-x(3))
c     &        + y(3)*(lnt-x(1))*(lnt-x(2))/(x(3)-x(1))/(x(3)-x(1))
c... this is for linear interpolation
         res = LogTinTzdat(klo_att(how),how,2)+
     &               (LogTinTzdat(khi_att(how),how,2)-LogTinTzdat(klo_att(how),how,2))
     &               *(lnt-LogTinTzdat(klo_att(how),how,1))/
     &               (LogTinTzdat(khi_att(how),how,1)-LogTinTzdat(klo_att(how),how,1))
      endif
      dsgetTzratio = 1.0d0
      dsgetTzratio = exp(res)

      return
      end ! dsgetTzratio
      

ccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dsdlogTdxrhs(z,lnTdm,dlogTdx)   
ccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'dsddcom.h'
      include 'dsnuclides.h' ! also includes dsmpconst.h !

      real*8 z,lnTdm,dlogTdx

      integer i
      real*8 res, tmp, omegamid, omegamax, mdm
      real*8 lnomegaa, lnomegab, lnTdmsav
      real*8 dsmwimp, dsddTrmax, dsddsoilcomp, dsTdsigdomega
      real*8 dTdzrescale
      parameter (dTdzrescale=1d2) ! NB: must be the same as further up!
      external dsTdsigdomega
      character*30 exploc
      integer tosav
      data tosav /-100/
      save tosav
      real*8 n(Nelements), mN(Nelements), Tdmcom
      common /Tdxrhscom/ n, mN, Tdmcom
      save /Tdxrhscom/
      logical initialized
      data initialized /.false./
      save initialized

c... these parameters are needed by dqagse
      real*8 abserr
      integer neval,ier, limit
      parameter (limit=40)
      real*8 alist(limit),blist(limit),rlist(limit),elist(limit)
      integer iord(limit),last
      
      if (.not.initialized.or.targetoption.ne.tosav) then
        do i=1, Nelements
          exploc = 'default'
          if (targetoption.eq.10.or.targetoption.eq.20) exploc='Gran_Sasso'
          exploc=exploc(1:index(exploc,' ')-1)
          n(i) = dsddsoilcomp(exploc,zn(i)) ! density in cm*w*-3
          mN(i) = mNaU(i)*atomicmassunit    ! mass in GeV
        enddo
        initialized = .true.
        tosav = targetoption
      endif

      res=00d0
      dlogTdx=0d0
      lnTdmsav=lnTdm
      if (lnTdm.lt.lnTzmin) lnTdm=lnTzmin
      if (lnTdm.gt.7d1) lnTdm=7d1
      Tdmcom = exp(lnTdm)
      mdm = dsmwimp()
      omegamid=0d0
      omegamax=0d0
      do i=1, Nelements
        tmp=0d0
        if (n(i).gt.1.0d10) tmp=dsddTrmax(Tdmcom, mdm, mN(i))
        if (i.eq.1.or.tmp.gt.0d0.and.tmp.lt.omegamid) omegamid=tmp
        if (tmp.gt.omegamax) omegamax=tmp
      enddo
      ! omegamid = omegamax !DEBUG
      lnomegaa = log(1d-10*omegamax)
      if (lnomegaa.lt.-50d0) lnomegaa=-50d0
      lnomegab = Log(omegamid)
      if (lnomegab.gt.70d0) lnomegab=70.
      call ! choose large abserr -> only rel. error relevant
     &    dqagse(dsTdsigdomega,lnomegaa,lnomegab,1.0d15,CRDM_acc/15.,limit,res,abserr,
     &           neval,ier,alist,blist,rlist,elist,iord,last)
      dlogTdx = res
c      goto 50 ! DEBUG
      if (lnomegab.ge.70d0) goto 100
      lnomegaa = lnomegab
      lnomegab = Log(omegamax)
      if (lnomegab.gt.70d0) lnomegab=70.
      call ! separate integral to better capture kinematic steps due to different nuclei
     &    dqagse(dsTdsigdomega,lnomegaa,lnomegab,1.0d15,CRDM_acc/5.,limit,res,abserr,
     &           neval,ier,alist,blist,rlist,elist,iord,last)
      dlogTdx = dlogTdx + res
      goto 50
 50   if (lnomegab.ge.70d0) goto 100
     
      if (CRDM_inelastic) then ! add integral up to Tdmcom in case of inelastic scattering
        lnomegaa = lnomegab
        lnomegab = Log(Tdmcom)
        if (lnomegab.gt.lnomegaa) then
          call
     &    dqagse(dsTdsigdomega,lnomegaa,lnomegab,1.0d15,CRDM_acc/5.,limit,res,
     &           abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          dlogTdx = dlogTdx + res
        endif
      endif
      

100   lnTdm = lnTdmsav
      
      dlogTdx = -dlogTdx/Tdmcom/dTdzrescale ! return dlogTdx rather than dTdm/dx
      
      end ! dsdlogTdxrhs

ccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function dsTdsigdomega(lnomega)
c... auxiliary routine for integration
ccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'dsddcom.h'
      include 'dsnuclides.h' ! also includes dsmpconst.h !

      real*8 omega, lnomega, dsddDMCRsigsoil, res
      integer i
      real*8 n(Nelements), mN(Nelements), Tdmcom
      common /Tdxrhscom/ n, mN, Tdmcom

      dsTdsigdomega = 0.0d0
      if (lnomega.lt.-60.0.or.lnomega.gt.70.0) then
        return
      endif
      
      omega = exp(lnomega)

      res=0d0
      do i=1, Nelements
        if (n(i).gt.1.0d10) then ! Helena Kolesova: neglect underabundant elements
          res = res + n(i)*dsddDMCRsigsoil(Tdmcom,omega,i)
          ! write(*,*) res, n(i), dsddDMCRsigsoil(Tdmcom,omega,i), Tdmcom, omega
          if (res.ne.res) write(*,*) 'ERROR in dsTdsigdomega:', omega,Tdmcom,dsTdsigdomega
        endif
      enddo
      
c... additional factor of omega from log-integration
      dsTdsigdomega = omega**2*res
      
      return
      end ! dsTdsigdomega
