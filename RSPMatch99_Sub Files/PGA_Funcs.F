      subroutine FindPgaIndex ( nQ, w, pgaIndex, scalePeriod )
      
c     Passed Variables      
      integer pgaIndex, nQ
      real w(1), scaleperiod
      
c     Local Variables      
      real twoPi, w1, delta, period1
      integer i
      
      twoPi = 2. * 3.141592653589790
      
c     Find the index for pga or other period speficfied by scalePeriod
c     Pga is defined as largest freq and smallest damp
      if ( scalePeriod .eq. 0. ) then
        w1 = 0.
        do i=1,nQ
          if ( w(i) .gt. w1  ) then
            pgaIndex = i
            w1 = w(i)
          endif
        enddo
      else
        delta = 1.0e30
        do i=1,nQ
          period1 = w(i) / twoPi
          if ( abs(period1 - scalePeriod) .lt. delta ) then
            delta = abs(period1 - scalePeriod)
            pgaIndex = i
          endif
        enddo
      endif
      
      return
      end subroutine FindPgaIndex

c --------------------------------------------------------------

      subroutine FindPga ( acc, nPts, pga )
      
c     Passed Variables      
      integer nPts
      real*8 acc(1)
      real pga
      
c     Local Variables
      integer i 
      real a1    
      
c     Find the pga 
      pga = 0.
      do i=1,nPts
        a1 = abs(acc(i))
        if ( a1 .gt. pga  ) then
          pga = a1
        endif
      enddo
      
      return
      end subroutine FindPga

c --------------------------------------------------------------

      subroutine ScaleToPga ( pgaIndex, acc, npts, R, nQ, Q )
      
      implicit none     
      integer pgaIndex, npts, nQ
      real*8 acc(1), R(1)
	real Q(1), scale

      scale = (Q(pgaIndex) / R(pgaIndex)) 
      write (*,'( 2x,''Scale to PGA, scale factor ='',f10.4)') scale
      
c     Apply scale factor to accelerogram and spectrum
      call ApplyScale ( acc, nPts, scale)
      call ApplyScale ( R, nQ, scale)

      return
      end subroutine ScaleToPga
      
c ----------------------------------------------------------------------