*******************************************************************************
*** Function dsddsigmarel returns the full differential cross section for   ***
*** elastic scattering of DM on a nucleus, or vice versa, expressed in      ***
*** units of the same quantity in the highly non-relativistic limit.        ***
*** Everything is expressed in terms of Lorentz-invariant quantities, so    ***
*** no assumption is made about the frame in which the cross section is     ***
*** is evaluated.                                                           ***
***                                                                         ***
***  type : replaceable                                                     ***
***  desc : relativistic DM - nucleus scattering cross section              ***
***                                                                         ***
***  Input:                                                                 ***
***    Q2    - momentum transfer squared [GeV**2]                           ***
***    s     - CMS energy squared [GeV**2]                                  ***
***    mN    - mass of nucleus [GeV]                                        ***
***    nuc2S - twice the nuclear spin                                       ***
***            (i.e. 0 for scalar nuclei, 1 for spin 1/2 nuclei etc.)       ***
***                                                                         ***
***  Output:  \sigma_N(Q2) / \sigma_N(Q2)_NR, dimensionless                 ***
***           = (mN+mDM)**2/s * |M|^2 / |M|^2_NR                            ***
***           [The output thus does not contain the nuclear form factor]    ***
***                                                                         ***
***  NB: This is the default function that corresponds to the asumption of  ***
***      a *constant* cross section, where the return value by definition   ***
***      is 1d0. A particle module or main program can replace this         ***
***      function with any expression for the relativistic amplitude. For   ***
***      convenience, a few examples for contact interactions and light     ***
***      mediators are provided below (just comment in the corresponding    ***
***      block).                                                            ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date  : 2022-07-29                                                      ***
***    mod: added nuclear spin as argument (H. Kolesova)                    ***
***    mod: 2025-05-26 (E.Rørnes): adapted to effective 2->2 el scatt       ***
*******************************************************************************
      real*8 function dsddsigmarel(Q2,s,mN,nuc2S)
      implicit none
      include 'dsio.h'
      
      real*8 Q2, s, mN
      integer nuc2S

      real*8 mdm, mred, res

c... comment in for use with light mediators
c      real*8 mmed
c      common /mymed/ mmed

c... functions
      real*8 dsmwimp
      

      mdm  = dsmwimp() ! needed because it may appear in the matrix element,
                       ! not needed to calculate kinematical factors
      mred = mN*mdm/(mN+mdm)

c... this is the result for a 'constant' scattering cross section (independent of nuclear spin)
      res = 1.0d0
      if (nuc2S.eq.1) then
        res = (mN+mdm)**2*(1.d0+Q2/(4.d0*mN**2))/s
      else
        res = 1.d0 ! We have not calculated the 2-loop diagram for four scalars.
                   ! Thus, the program is not accurate in this case...
        res = (mdm+mN)**2/s/mN**2 ! However, we can approximate it as follows.
                                  ! The mass dimension is different due to how
                                  ! we treat cpsi. Since [cpsi]=-1 when N is a
                                  ! fermion then [cpsi]=0 when N is a scalar.
                                  ! Thus, we must include an additional factor of
                                  ! mN in cpsi which is the correction done here
                                  ! compared to the usual 4-point interaction below.
      endif


c... this is the result for a constant amplitude, e.g. a scalar 4-point interaction
c      res = (mdm+mN)**2/s


c... this is the result for a SCALAR mediator (Dirac DM)
c      res = mmed**4/(mmed**2 + Q2)**2
c     &               * mN**2*(Q2+4*mdm**2) / (4.*s*mred**2) ! spin 0 nuclei
c      if (nuc2S.eq.1) then ! extra factor spin 1/2 nuclei
c        res = res * (1d0 + Q2/4./mN**2)
c      else
c        if (prtlevel.gt.0) then
c           write(*,*) 'No available cross section for nuclear spin 2S=',nuc2S
c           write(*,*) 'Will use cross section for scalar nuclei instead'
c        endif
c      endif


c... this is the result for a VECTOR mediator (Dirac DM)
c      res = mmed**4/(mmed**2 + Q2)**2
c     &               / (4.*s*mred**2)
c      if (nuc2S.eq.0) then ! spin 0 nuclei
c        res = res * (Q2*mdm**2 - Q2*s + (s-mN**2-mdm**2)**2)
c      elseif (nuc2S.eq.1) then !spin 1/2 nuclei
c        res = res * (Q2**2/2d0 - Q2*s + (s-mN**2-mdm**2)**2)
c      else
c        res = res * (Q2*mdm**2 - Q2*s + (s-mN**2-mdm**2)**2)
c        if (prtlevel.gt.0) then
c           write(*,*) 'No available cross section for nuclear spin 2S=',nuc2S
c           write(*,*) 'Will use cross section for scalar nuclei instead'
c        endif
c      endif

      dsddsigmarel = res

      return
      end

