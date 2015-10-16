      subroutine RdRunParam ( tol, maxIter, gamma, scaleFlag, dtFlag,
     1      evMin, groupSize, maxFreq, freqBand, nPole1, scalePeriod, 
     2      iModPGA, iSeed, ranFactor, freqMatch, flag2, scale1 )
     
C      implicit none
c     Passed Variables
      integer maxIter, scaleFlag, dtFlag,
     1        groupSize, nPole1, iModPGA, iSeed, flag2
      real tol, gamma, a1, a2, f1, f2, evMin, maxFreq, freqBand(1), 
     1          scalePeriod
      real ranFactor, freqMatch(2), scale1

c     Commmon
      integer iModel
      common / Model / iModel
      common /alphaModel/ a1, a2, f1, f2
      
c     Enter run parameters
      write (*,'( 2x,''Enter maximum number of iterations (eg 30)'')')
      read (5,*) maxIter
      write (*,'( 2x,''Enter tolerance for spectral match (eg 0.02)'')')
      read (5,*) tol
      write (*,'( 2x,''Enter convergence damping (eg 0.5)'')')
      read (5,*) gamma
      write (*,'( 2x,''Select model:  1=reverse impulse response'',/,
     2           17x,''6=tapered cosine'',/, 
     3           17x,''7=improved tapered cosine'')')
      read (5,*) iModel
      write (*,'( 2x,''Enter alpha model (a1, a2, f1, f2)'')')
      read (5,*) a1, a2, f1, f2
      write (*,'( 2x,''Enter scaleflag, scalePeriod'')')
      read (5,*) scaleFlag, scalePeriod
      write (*,'( 2x,''Enter interpolation flag?'')')
      read (5,*) dtFlag
      write (*,'( 2x,''Enter minimum eigenvalue (e.g. 1.0e-5)'')')
      read (5,*) evMin
      write (*,'( 2x,''Enter group size (eg 20)'')')
      read (5,*) groupSize
      write (*,'( 2x,''Enter maxFreq (HZ)'')')
      read (5,*) maxFreq
      write (*,'( 2x,''Enter Initial filters (min, max in Hz) and nPole'')')
      read (5,*) freqBand(1), freqBand(2), nPole1
      write (*,'( 2x,''Mod PGA?'')')
      read (5,*) iModPGA
      write (*,'( 2x,''Randomize target? (0=no), factor'')')
      read (5,*) iSeed, ranFactor
      write (*,'( 2x,''Enter freq bands for matching (fmin, fmax in Hz)'')')
      read (5,*) freqMatch(1), freqMatch(2)
      write (*,'( 2x,''Enter baseline cor flag (1=yes)'')')
      read (5,*) flag2
      write (*,'( 2x,''Enter scale fact '')')
      read (5,*) scale1
      
c     Convert maxFreq to radians
      maxFreq = maxFreq * 2. * 3.1415926
      
      return
      end subroutine RdRunParam
      
c ----------------------------------------------------------------------

      subroutine RdTarget ( w, damp, fileName, minFreq, minTime, 
     1      minTime_all, beta_all, Q_all, nQ_all, nFreq, nDamp, target2,
     2      w_all, maxTime, maxTime_all, freqmatch )
     
      implicit none
      include 'rspMatch.h'

c     Passed Variables
      integer nQ_all, nFreq, nDamp
      real w(1), damp(1), minFreq, minTime(1), target2(MAX_FREQ,1)
      real Q_all(1), beta_all(1), minTime_all(1), w_all(1), maxTime(1),
     1     maxTime_all(1), freqmatch(2)
      character*80 fileName
      
c     Local Variables
      character*80 dummy
      integer i, j, k, m
	real freq, pi, minTime_tmp, maxTime_tmp, 
     1      Target_tmp(MAX_FREQ,MAX_DAMP)
      pi = 3.1415926
     
c     Open Target file
      write (*,'( 2x,''Enter target spectrum file'')')
      read (5,'( a60)') fileName
      write (*,'( 2x,a80)') fileName
      open ( 10, file=fileName, status='old' )
      
c     Read Header
      do 20 i=1,1
        read (10,'( a1)') dummy
  20  continue
      read (10,*) nFreq, nDamp
      read (10,*) (damp(k),k=1,nDamp)
  
c     Check array dimensions
      if ( nFreq*nDamp .gt. MAX_MATCH ) then
        write (*,'( 2x,''Increase MAX_MATCH to'',i5)') nFreq*nDamp
        stop 99
      endif
      if ( nFreq .gt. MAX_FREQ ) then
        write (*,'( 2x,''Increase MAX_FREQ to'',i5)') nFreq
        stop 99
      endif
      
      m = 1
      minFreq = 100.     
c     Read Target and find minimum frequency removing frequencies that 
c     are outside the specified range in the input file - LA 05/08/09
      do 30 i=1,nFreq
        read (10,*) freq, minTime_tmp, maxTime_tmp, 
     1      (target_tmp(i,j),j=1,nDamp) 
        if ( freq .ge. freqMatch(1).and.freq .le. freqMatch(2) ) then
            w(m) = 2.*pi * freq
            minTime(m) = minTime_tmp
            maxTime(m) = maxTime_tmp
            do j=1, nDamp
                target2(m,j) = target_tmp(i,j)
            enddo
            m = m+1
            if ( freq .lt. minFreq ) minFreq = freq
        endif
  30  continue
      nFreq = m-1

c     Load Full Set of damping and Freq pairs
      k = 1
      do i=1,nFreq
        do j=1,nDamp
            Q_all(k) = target2(i,j)
            beta_all(k) = damp(j)
            minTime_all(k) = minTime(i)
            maxTime_all(k) = maxTime(i)
            w_all(k) = w(i)
            k = k + 1
        enddo
      enddo
      nQ_all = k-1     
      
      close (10)
      
      return
      end subroutine RdTarget 
      
c ----------------------------------------------------------------------

      subroutine RdAcc ( acc, nPts, dt, dtFlag, fileName, nAdded )
      
      implicit none
      include 'rspMatch.h'
c     Passed Variables
      integer nPts, dtFlag, nAdded
      character*80 fileName
      real*8 acc(1)
      real dt
      
c     Local Variables
      character*80 dummy
      integer i, i1, j, ifactor
     
c     Open Accelerogram file
      write (*,'( 2x,''Enter accelerogram file'')')
      read (5,'( a60)') fileName
      open ( 10, file=fileName, status='old' )
      
c     Read Header 
C      do 10 i=1,4 
        read (10,'( a1)') dummy
C  10  continue
      read (10,*) nPts, dt, nAdded
      write (*,'( 2x,''nPts ='',i5)') nPts
      if ( npts .gt. MAXPTS ) then
        write (*,'( 2x,''*** truncating number of points to'',i5)')
     1        MAXPTS
        npts = MAXPTS
      endif
      
c     Read Data
      read (10,*) (acc(i),i=1,nPts)
      close (10)
      
c     Interpolate to smaller time step
      ifactor = dtFlag
      if ( dtFlag .gt. 1 ) then
        i1 = npts
        npts = iFactor*(nPts-1) + 1
        if ( npts .gt. MAXPTS) then
          write (*,'( 2x,''ERROR: Increase MAXPTS to '',i5)') npts
          stop 99
        endif
        dt = dt / ifactor
        do 20 i=npts,1,-ifactor
          acc(i) = acc(i1)
          i1 = i1 -1
  20    continue
        do 30 i=1,nPts,ifactor
          do 25 j=1,ifactor-1
            acc(i+j) = acc(i) + ( acc(i+ifactor) - acc(i) ) *
     1                 float(j) / float(ifactor)
  25      continue
  30    continue
  
        write (*,'( 2x,''interpolated to smaller time step'')')
        write (*,'( 2x,''dt ='',f10.4)') dt
        write (*,'( 2x,''New nPts ='',i5)') nPts
      endif
      
      return
      end  subroutine RdAcc
      
c ----------------------------------------------------------------------
