*******************************************************************************
*** Function dsddlfreesimp returns a simple estimate for the mean free      ***
*** path of strongly interacting DM through the Earth crust. We include     ***
*** contributions from the 11 most abundant elements, with densities as     ***
*** returned from dsddsoilcomp.                                             ***
***                                                                         ***
***  Input:                                                                 ***
***    sigsi - spin-independent cross section per nucleon [cm**2]           ***
***    depth - detector location below surface [cm]                         ***
***    how   - assume equal couplings of DM to protons and neutrons (how=1) ***
***            or assume DM only couples to protons (how=2)                 ***
***                                                                         ***
***  Output: *average* mean free path in cm, between surface and depth      ***
***          specified as input.                                            ***
***                                                                         ***
***  NB: If this routine is called directly, you first need to set          ***
***          'targetoption' (typically set in dsddDMCRcountrate) for the    ***
***          correct experimental location                                  ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2019-02-06                                                         ***
*** mod tb, 12/01/2021: added target location dependence                    ***
*******************************************************************************
      real*8 function dsddlfreesimp(sigsi, depth, how)
      implicit none
      include 'dsio.h'
      include 'dsddcom.h'
      include 'dsnuclides.h' ! also includes dsmpconst.h !

      real*8 sigsi, depth
      integer how
      real*8 linv, dsmwimp, mdm
      integer i
      real*8 nav(Nelements), mN(Nelements)
      character*30 exploc
      integer acom
      common /earthdensaux/ acom 
      
      real*8 dsddsoilcomp
 
      linv=1.0d-50
      if (sigsi.lt.0.0d0.or.depth.lt.0.0) goto 100
      if (sigsi.lt.1.d-35) goto 50 ! no attenuation

      mdm = dsmwimp()
      do i=1,Nelements
         mn(i)  = mNaU(i)*atomicmassunit ! convert masses to GeV
         acom   = an(i)
         exploc = 'default'
         if (targetoption.eq.10.or.targetoption.eq.20) exploc='Gran_Sasso'
         exploc=exploc(1:index(exploc,' ')-1)
         nav(i) = dsddsoilcomp(exploc,zn(i)) ! density in cm**-3
      enddo

      if (how.eq.1) then
        linv=0.0
        do i=1,Nelements
           linv = linv + nav(i)*AN(i)**2*mN(i)**3/(mn(i)+mdm)**4
        enddo
        linv = linv*2*sigsi*(mdm+m_p)*(mdm+m_n)*mdm/m_p/m_n

      elseif (how.eq.2) then
        linv=0.0
        do i=1,Nelements
           linv = linv + nav(i)*ZN(i)**2*mN(i)**3/(mn(i)+mdm)**4
        enddo
        linv = linv*2*sigsi*(mdm+m_p)**2*mdm/m_p**2

      else
        if (prtlevel.gt.1) then
          write(*,*) 'ERROR in dsddlfreesimp:'
          write(*,*) 'unknown option how = ',how
          write(*,*) 'Will return 1.0d40 cm for mean free path...'
        endif  
      endif      

 50   dsddlfreesimp=1.0/linv            
      return

 100  write(*,*) 'FATAL ERROR in dsddlfreesimp:'
      write(*,*) 'You have supplied a negative cross section or depth!'
      stop

      end
