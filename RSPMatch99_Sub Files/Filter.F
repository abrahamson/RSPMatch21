      subroutine BandPass ( acc, nPts, dt, maxFreq, minFreq, nPole  )
      
      implicit none
      include 'rspMatch.h'

      real*8 acc(1)
      real dt, df, minFreq, maxFreq, fc1, fc2
      integer tb, te, npts, nmin, mpad, i, j, npts1,
     1        nPole
      complex cu(MAXPTS)

      if ( maxFreq .eq. 0. .and. minFreq .eq. 0 ) return

c     Set tapers and min padding
      tb = 1
      te = 2

c     REMOVE THE DC        
      call rdc ( acc, npts, 0 )

c     TAPER
      call CosTaper ( acc, npts, tb, te )

c     PAD TO POWER OF 2
      npts1 = npts
	  nmin = 0
      call pad ( acc, npts1, nmin, mpad )

c     FILL COMPLEX ARRAY
      do j=1,npts1
        cu(j) = cmplx(acc(j),0.0)
      enddo

c     CALCULATE FORWARD FFT
      call cool ( -1., mpad, cu )
	  df = 1. / ( npts1*dt )
	  
c     Filter
	  fc1 = minFreq 
	  fc2 = maxFreq
      if ( fc1 .gt. 0. ) then
c	    write (*,'( 2x,''highpass'')') 
        call hipass  ( fc1, nPole, df, npts1, cu )
      endif
	  if ( fc2 .gt. 0. ) then
	    call lowpass ( fc2, nPole, df, npts1, cu )
	  endif
	  
c     Compute Inverse FFT
      call cool ( 1., mpad, cu )
	  
c     Apply scale factor (only original length)
      do i=1,npts
	    acc(i) = real(cu(i)) / npts1 
	  enddo

      return
      end subroutine BandPass

	  
c ----------------------------------------------------------------------

      subroutine hipass ( fc, nPole, df, npts, cx )

c     This subroutine applies a butterworth high-pass filter
c     in the frequency domain

      real fc, freq, fBar, scale, df
      integer nPole, npts
	complex cx(1)
	
c     Local Variabled
      integer i	

	    cx(1) = cmplx( 0.,0. )
		freq = df
        do i=2,npts/2+1
          fbar = ( freq / fc )**( 2 * nPole )
          scale = sqrt( fBar / (1. + fBar) )
          cx(i) = cx(i) * scale
          cx(npts-i+2) = cx(npts-i+2) * scale
          freq = freq + df
        enddo
      return
      end subroutine hipass

c ----------------------------------------------------------------------

      subroutine lowpass ( fc, nPole, df, npts, cx )

c     This subroutine applies a butterworth low-pass filter
c     in the frequency domain

      complex cx(1)
      real fc, freq, fBar, amp, df
      integer nPole, npts
      
c     Local Variables
      integer i      

      freq = df
      do i=2,npts/2 + 1
          fBar = ( freq / fc )**( 2 * npole )
          amp = sqrt( 1. / (1.+fBar) )
          cx(i) = cx(i) * amp
          cx(npts-i+2) = cx(npts-i+2) * amp
          freq = freq + df
      enddo
      
      return
      end subroutine lowpass

c ----------------------------------------------------------------------

      subroutine cool ( signi, n, cx )
c     FFT subroutine.
c     signi = -1.  forward transform
c           =  1.  inverse transform
c     n = log base 2 (npts)

      complex cx(1), carg, temp, cw
      integer istep, i, j, n
      real signi
      
      real pi
      integer lx, m, l
      
      pi = 4. * atan(1.) * signi
      lx = 2**n
      j = 1
      do 30 i=1,lx
        if (i .gt. j) goto 10
        temp = cx(j)
        cx(j) = cx(i)
        cx(i) = temp
  10    m = lx/2
  20    if (j .le. m) goto 25
        j = j-m
        m = m/2
        if (m .ge. 1) goto 20
  25    j = j+m
  30  continue
      l = 1
  40  istep = l+l
      do 50 m=1,l
        carg = cmplx( 0., pi * float(m-1) / float(l) )
        cw = cexp(carg)
        do 45 i=m,lx,istep
          temp = cw * cx(i+l) 
          cx(i+l) = cx(i) - temp
          cx(i) = cx(i) + temp
  45    continue
  50  continue
      l = istep
      if (l .lt. lx) goto 40
      
      return
      end subroutine cool
c ----------------------------------------------------------------------

      subroutine pad (x,npts,nmin,m)

c     This subroutine pads the x array to a power of 2 for FFT

      real*8 x(1)
      integer test,npts,m,nmin
      
      integer i

      do 10 i=1,20
        test = 2**i
        if (test .ge. npts .and. test .ge. nmin) goto 20
  10  continue
  20  m = i
c       write (*,1000) test
 1000 format( 2x,'Pad to ',i6,' points')

      do 30 i=npts,test
        x(i) = 0.0
  30  continue
      npts = test

      return
      end subroutine pad

c ----------------------------------------------------------------------

      subroutine CosTaper (x,npts,tb,te)

c     This subroutine tapers the x array
    
      real*8 x(1)
      real arg, pi
      integer npts,tb,te, i, n
      
      pi = 3.1415926

      if (tb .ne. 0.) then
        n = (npts*tb)/100
        do 10 i=1,n
          arg = pi*float(i-1)/float(n) + pi
          x(i) = x(i)*(1.+cos(arg))/2.
  10    continue
c         write (*,1000) tb
 1000   format( 2x,'Taper beginning ',i2,' percent')
      endif

      if (te .ne. 0.) then
        n = (npts*te)/100
        do 20 i=1,n
          arg = pi*float(i-1)/float(n) + pi
          x(npts-i+1) = x(npts-i+1) * (1.+cos(arg))/2.
  20    continue
c         write (*,1001) te
 1001   format( 2x,'Taper end ',i2,' percent')
      endif
      
      return
      end subroutine CosTaper

c ----------------------------------------------------------------------
      subroutine rdc (x,npts,iflag)

c     This subroutine removes a dc shift from the data

c     IFLAG = 0    remove the mean
c     IFLAG = 1    remove the mean value of the first 10 points
c     IFLAG = 2    manual set of DC value to be removed

      real*8 x(1)
      real sum,mean
      integer npts,iflag, i

      if (iflag .eq. 0) then
        sum = 0.0
        do 10 i=1,npts
          sum = x(i) + sum
  10    continue
        mean = sum/float(npts)

      elseif (iflag .eq. 1) then
        sum = 0.0
        do 20 i=1,10
          sum = x(i) + sum
  20    continue
        mean = sum / float(npts)

      else
        write (*,1000)
        read (*,1001) mean

      endif

      do 100 i=1,npts
        x(i) = x(i) - mean
  100 continue

c       write (*,1010) mean
c 1010 format( 2x,'Remove DC of ',f12.8)

      return
 1000 format( 2x,'Enter mean to be removed')
 1001 format( f12.8)
      end subroutine rdc

c ----------------------------------------------------------------------

      subroutine taper ( x, n )
      
      real x(1), arg, pi
      integer n, i
      pi = 3.1415926
      
      if ( n .eq. 0 ) return
      do i=1,n
        arg = pi*float(i-1)/float(n)
        x(i) = x(i) * 0.5 * ( 1. - cos(arg) )
      enddo
      
      return
      end subroutine taper

c ----------------------------------------------------------------------