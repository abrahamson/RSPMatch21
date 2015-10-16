	subroutine PadAcc ( acc, nPts, dt, w_all, nQ_all, tpeak,  
	1      t_before, nPts_before )

	! Added by Linda on Feb 23, 2009
	implicit none
	include 'rspMatch.h'
c	Passed Variables
	real*8 acc(1)
	real dt, t_before, w_all(1), tPeak(1)
 	integer nPts, nQ_all
c	Local Variables
	real*8 acc_temp(MAXPTS)
	real t_temp_b, tmin 
	integer i, nPts_before

c	  Padding the beginning of the acc time series if needed
	  t_temp_b = 0
	  nPts_before = 0
	  do i=1, nQ_all
	      tmin = 18.54 * ( w_all(i)**(-0.84) )
	      if ( tpeak(i) < tmin ) then
	          if ( t_temp_b < (tmin - tpeak(i)) ) then
	              t_temp_b = tmin - tpeak(i)
	          endif
	      endif
	  enddo

	  if ( t_temp_b > t_before) then
		    nPts_before = NINT( (t_temp_b - t_before)/dt ) + 1
		    do i=1, nPts_before
			    acc_temp(i) = 0
		    enddo
		    do i=nPts_before+1, nPts_before + nPts  
			    acc_temp(i) = acc(i-nPts_before)
		    enddo
		    nPts = nPts_before + nPts
		    do i=1, nPts	
			    acc(i) = acc_temp(i)
		    enddo
		    t_before = t_temp_b 
	  endif

       end subroutine PadAcc
      
c ----------------------------------------------------------------------

      subroutine ApplyScale ( y, n, scale )
      
      implicit none
c     Passed Variables
      integer n
      real scale
      real*8 y(1)
      
c     Local Variables
      integer i
      
      do i=1,n
        y(i) = y(i) * scale
      enddo
      
      return
      end subroutine ApplyScale
      
c ----------------------------------------------------------------------

      subroutine AddPga ( acc, npts, nQ, Q, maxFreq, dt, pgaAcc, 
     1      iModPGA )
     
c     Passed Variables      
      real Q(1), maxFreq, pgaAcc(1), dt
      real*8 acc(1)
      integer npts, nQ, iModPGA, n1
      
c     Passed Variables
      real twoPi, a1, polarity, gamma, scale, omega, t1, t0, t, pga,
     1  Adjust
      integer i, iTime      
      
      twoPi = 2. * 3.1415926
      
      if ( iModPga .eq. 0 ) then
        do i=1,npts
          pgaacc(i) = 0.
        enddo
      endif
      
c     Find time and polarity of pga
      pga = 0.
      do i=1,nPts
        pgaacc(i) = 0.
        a1 = abs(acc(i))
        if ( a1 .gt. pga  ) then
          pga = a1
          iTime = i
          if (acc(i) .lt. 0. ) then
            polarity = -1.
          else
            polarity = 1.
          endif
        endif
      enddo
 
c     Add a small bump to acc at time of pga
c     Use cosine bell at max freq  (max freq in radians)
c     gamma is a convergence factor
c      gamma = 0.5
      gamma = 1.5      
      scale = ( pga - Q(nQ) ) * polarity * gamma 
      omega = maxFreq
      t1 = twoPi / omega
      n1 = int( t1 / dt )
      t0 = iTime * dt
      write (*,'( 6f10.4,i5)') pga, Q(nQ), polarity, gamma, scale, t0,n1
      do i=iTime-n1,iTime+n1
        t = i * dt
        adjust = scale * 0.5 * ( cos( omega*(t-t0) ) + 1. )  
        pgaacc(i) = - adjust
      enddo
      return
      end subroutine AddPga
          
c ----------------------------------------------------------------------