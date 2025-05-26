*******************************************************************************
*** Function dsddDMCRinvert_cpsi uses the bisection method to find the C in ***
*** msq=C or msq=C/s such that the computed cpsi matches a target value.    ***
***                                                                         ***
***  Hidden input:                                                          ***
***    cpsi_target - target cpsi value                                      ***
***                                                                         ***
***  Output:  C corresponding to target cpsi                                ***
***                                                                         ***
*** author: E.Rørnes                                                        ***
*** date 2025-04-26                                                         ***
*******************************************************************************

      real*8 function dsddDMCRinvert_cpsi()
      implicit none
      include 'dsmpconst.h' ! for pi and m_p
      include 'dsddcom.h'

      real*8 cpsi_target, msq_min, msq_max, msq_mid
      real*8 cpsi_mid, tol
      integer max_iter, iter
      
      real*8 sip, sin, sdp, sdn, mdm

      integer ierr, dsidnumber ! to make sure cross section is only called once per model
      integer idold
      data idold/-123456789/
      save idold, sip, sin, sdp, sdn

c... Helper function
      real*8 calculate_cpsi, dsmwimp
      external calculate_cpsi

c... Parameters
      tol = 1.d-8
      max_iter = 100

      dsddDMCRinvert_cpsi = 1.d0

      if (idold.ne.dsidnumber()) then
        call dsddsigmanucleon(0.0d0, 0.0d0, sip, sin, sdp, sdn, ierr)
      endif

      mdm = dsmwimp()
c... In the non-rel limit sigma0=cpsi**2*mN**2/(4*pi*(mdm+mN)**2)
      cpsi_target = 2*sqrt(pi*((sip+sin)/2.d0+sdp))*(mdm+m_p)/m_p  ! The program uses sigma_NR

      if (cpsi_target.ge.1.d-4 .or. cpsi_target.lt.1.d-50) return 
                                              ! The function is not invertible
                                              ! at the higher end.
                                              ! Corresponds to unphysical msq.
                                              ! Lower end is to avoid not finding
                                              ! a solution.

      msq_min = 1.d-50
      msq_max = 1.d10

c... Start bisection
      do iter = 1, max_iter
        msq_mid = sqrt(msq_min*msq_max)
        cpsi_mid = calculate_cpsi(msq_mid)
c... Check relative convergence
        if (abs(cpsi_target/cpsi_mid-1.d0).lt.tol) then
          dsddDMCRinvert_cpsi = msq_mid
          ! Remove the non-rel factor of 1/s in Msq to get C
          if (CRDM_cs) dsddDMCRinvert_cpsi = dsddDMCRinvert_cpsi*(2.d0*mdm+m_p)**2
          return
        endif

c... Update interval
        if (cpsi_mid .gt. cpsi_target) then
          msq_max = msq_mid
        else
          msq_min = msq_mid
        endif
      enddo
      if (iter.ge.max_iter) then
        write(*,*) "Warning: dsddDMCRinvert_cpsi did not converge!"
        write(*,*) 'Estimated relative error on msq:', 1.d0+sqrt(abs(cpsi_target/cpsi_mid - 1.d0)**3)
      endif

c... If maximum iterations reached, return best estimate
      dsddDMCRinvert_cpsi = msq_mid
      return
      end


*******************************************************************************
*** Function calculate_cpsi computes cpsi given msq                         ***
*******************************************************************************
      real*8 function calculate_cpsi(msq)
      implicit none
      include 'dsmpconst.h' ! for pi and m_p
      include 'dsddcom.h'

      real*8 mdm, msq, s_NR, sigv32, cut_off_sq, term1, term2, term3
      real*8 dsmwimp

      mdm = dsmwimp()

      if (mdm.le.0.d0) then
        calculate_cpsi = 0.d0
        return
      endif
      
      s_NR = (2.d0*mdm+m_p)**2

      sigv32 = (1.d0/(2.d0*8.d0*mdm**2*m_p))
     &       * sqrt(s_NR**2-2.d0*s_NR*(m_p**2+mdm**2)+(m_p**2-mdm**2)**2)
     &       / (8.d0*pi*s_NR)
     &       * msq

      ! if (CRDM_cs) sigv32=sigv32/s_NR

      cut_off_sq = sqrt(sqrt(3.d0)/(4.d0*pi*mdm*sigv32))

      term1 = m_p/((4.d0*pi)**4)/cut_off_sq
      term2 = 1.d0 - (mdm**2/cut_off_sq)
      term3 = log((cut_off_sq + m_p**2)/(4.d0*mdm**2))

      calculate_cpsi = term1*term2*term3

      return
      end
