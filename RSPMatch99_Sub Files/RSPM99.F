      program RSPM99
c
c     Version: 99
c     Date: July 14, 1999
c
c     NOTE: This spectral matching program is based on the computer
c           program spectralMatch developed by N. Abrahamson.
c
c     Spectral matching program by N. Abrahamson.Copyright July 7,1993.
c     Based on time domain method of Tseng and Lilanand (1988).
c     Modified to preserve non-stationarity at long periods by using 
c     different functional forms for the adjustment time history.
c
c     Note: This version of the code was taken from the RSPM22B version 
c           of the code. Only non-critical input/output changes were 
c           made to this new version of the code. 
c
c     This source code is concatenated from the following subprograms:
c    
c          match_main
c          BASELINE
c          calcnewth
c          calc_c
c          fftsub
c          ft_th
c          match_io
c          matrix
c          median
c          misfit
c          response
c
c     If this were run on the Macintosh computer, then the 
c     following subprograms could also be included for plotting:
c          updateplot
c          tek
c
c
c              Also, there may be tabs in this code that may confuse
c              your compiler.  Edit as appropriate.
c
      
      implicit none
      include 'rspMatch.h'
      
      real acc(MAXPTS), Q(MAX_MATCH), R(MAX_MATCH), tPeak(MAX_MATCH), 
     1     w(MAX_MATCH), beta(MAX_MATCH), fc, a1, a2, f1, f2, 
     2     R0(MAX_MATCH)
      real w_unmatched(MAX_MATCH), Q_unmatched(MAX_MATCH)
      real Q1(MAXGROUP,MAX_MATCH), w1(MAXGROUP,MAX_MATCH), 
     1     beta1(MAXGROUP,MAX_MATCH), accSave(MAXPTS), bestMisfit
      real Q_all(MAX_MATCH), beta_all(MAX_MATCH),Q_all_noRand(MAX_MATCH)
      real Q_all_noRand_unmatched(MAX_MATCH)
      real target2(MAX_FREQ,MAX_DAMP), damp(MAX_DAMP), 
	1     minTime_all(MAX_MATCH),
     1     w_all(MAX_MATCH), pgaacc(MAXPTS)
      real*8 deltaR(MAX_MATCH) 
      real tol, gamma, dt, misfit, scale, tol1, aveMisfit, minFreq,
     1     evMin, misFit1, lastMisfit, minTime(MAX_MATCH), 
     2     minTime1(MAXGROUP,MAX_MATCH), minTime2(MAX_MATCH),
     3     maxfreq, pga, freqBand(2), scalePeriod, ranFactor
      real freqMatch(2), scale1
      real maxTime(MAX_MATCH), maxTime_all(MAX_MATCH),
	1     Time_array(MAX_MATCH),
     2     maxTime1(MAXGROUP,MAX_MATCH), maxTime2(MAX_MATCH), 
     3     zero(MAX_MATCH)
      integer iter, maxIter, nQ, P(MAX_MATCH), nPts, i, j, nPad, iFreq,
     1        scaleFlag, dtFlag, pgaIndex, iDamp, iModPga
      integer iGroup, nGroup, nQ1(MAXGROUP), groupSize, flag1, nPole1, 
	1        flag2
      character*80 targetFile, thFilein, thFileout, version, file1
      integer iModel, nQ_all, nFreq, nDamp, iSeed
      common / Model / iModel
      common /alphaModel/ a1, a2, f1, f2
      data version / 'Version:99, Date:7_14_99'/

      flag1 = 0

      write (*,'(2x,''********************************'')')
      write (*,'(2x,''RSPM - Spectral Matching Program'')')
      write (*,'( 2x,''Program Version: '',a25)') version
      write (*,'(2x, ''Based on Code by, N. Abrahamson (1993)'')')
      write (*,'(2x,''********************************'')')

      write (*,'( 2x,''Enter run file'')')
      read (*,'( a80)') file1
      open (5,file=file1,status='old')
      do i=1,MAX_MATCH
        zero(i) = 0.
      enddo
      
c     Initialize BestMisfit
      bestMisfit = 1.0e30
           
c     Read Run parameters, target spectrum, and reference accelerogram
      call RdRunParam ( tol, maxIter, gamma, scaleFlag, dtFlag, evMin,
     1                groupSize, maxFreq, freqBand, nPole1, scalePeriod, 
     1             iModPGA, iSeed, ranFactor, freqMatch, flag2, scale1 )
      call RdTarget (w, damp, targetFile, minFreq, minTime, minTime_all,
     1             beta_all, Q_all, nQ_all,nFreq, nDamp, target2, w_all,
     2                maxTime, maxTime_all)
      call RdAcc ( acc, nPts, dt, dtFlag, thFilein )

      do i=1,MAX_MATCH
        Time_array(i) = npts*dt
      enddo

c     Randomize Target?
      do i=1,nQ_all
        Q_all_norand(i) = Q_all(i)
      enddo
      if ( iSeed .ne. 0 ) then
        call RandomizeTarget ( iSeed, Q_all, nQ_all, ranFactor, target2,
	1                              nDamp, nFreq )
      endif

c     Initial Band-pass filter
      freqBand(2) = freqBand(2) * 6.28
      call BandPass ( acc, npts, dt, freqBand(2), freqBand(1), nPole1 )
      
c     Calculate response spectrum of reference time history
      call CalcRsp ( acc, nPts, dt, w_all, beta_all, nQ_all, tPeak, R0,
	2               P, zero,
     1               time_array)
      
c     Calulate initial misfit from target 
      call CalcMisfit ( R0, Q_all, nQ_all, P, deltaR, misfit, iFreq, 
	1          aveMisfit, freqMatch, w_all )
      write (*,'( /,2x,''Initial Solution    AveMisfit  MaxMisfit  Freq(
	1Hz)  Damping  Pga(g)'')')
      write (*,'( 8x,12x,f10.4,1x,f10.4,2x,f8.3,1x,2f8.3)')
     1           aveMisfit, (R0(iFreq)-Q_all(iFreq))/Q_all(iFreq), 
     2           w_all(iFreq)/6.28, beta_all(iFreq)

c     Make initial adjustment to time history by scaling to PGA 
c          (or other period)
      call FindPgaIndex ( nQ_all, w_all, pgaIndex, scalePeriod )
      call FindPga ( acc, npts, pga )
      write (*,'( 2x,''scaleIndex ='',i5,2x,''pga ='',f10.4)') pgaIndex,
	1            pga
      if ( scaleFlag .ge. 1 ) then
        call ScaleToPga ( pgaIndex, acc, npts, R0, nQ_all, Q_all )
        call FindPga ( acc, npts, pga )
        call CalcMisfit ( R0, Q_all, nQ_all, P, deltaR, misfit, iFreq, 
	1           aveMisfit, freqMatch, w_all)
        write (*,'( /,2x,''Scaled Solution     AveMisfit  MaxMisfit  Fre
	1q(Hz)  Damping  Pga(g)'')')
        write (*,'( 8x,12x,f10.4,1x,f10.4,2x,f8.3,1x,3f8.3)')
     1           aveMisfit, (R0(iFreq)-Q_all(iFreq))/Q_all(iFreq), 
     2           w_all(iFreq)/6.28, beta_all(iFreq), pga, R0(pgaIndex)
      endif
      do i=1,npts
        acc(i) = acc(i) * scale1
      enddo

c     Initialize Graphics
c      call plots
c      call factor (1.5)
c      call UpdatePlot ( w_all, Q_all, beta_all, R0, nFreq, nDamp, 0,0,
c    1                aveMisfit, misfit, iFreq )  
      
c     Iterate until misfit is less than tolerance
      iter = 1
      lastMisfit = 1.0e10
      doWhile ( misfit .gt. tol .and. iter .le. maxIter )
        
        do iDamp=1,nDamp
          call LoadQ_array ( nFreq, nDamp, nQ, Q, beta, target2, damp, 
	1           iDamp, w )

c         Subdivide Target into subgroups
          call SubdivideTarget ( nQ, Q, w, beta, nQ1, Q1, w1, beta1, 
	1             nGroup, 
     1                       groupSize, minTime, minTime1, maxfreq, 
     2             maxTime, maxTime1 )

          call DoMatch ( nGroup, Q1, w1, beta1, acc, nPts, dt, misfit, 
     1                 R, Q_all, nQ1, evMin, gamma,lastMisfit,scaleFlag,
     2                 w_all, beta_all, nQ_all, iFreq, flag1, aveMisfit,
     2         minTime1,
     3                 minTime_all, maxFreq, minFreq, iter, iDamp, 
     3           scalePeriod, freqMatch, flag2,
     4                 maxTime1, maxTime_all, pgaacc)
          call FindPga ( acc, npts, pga )
          write (*,'( 2x,i6,'' full set'',3x,f10.4,1x,f10.4,2x,f8.3,1x,
	1        	4f8.3)')
     1           iter, aveMisfit, (R(iFreq)-Q_all(iFreq))/Q_all(iFreq), 
     2           w_all(iFreq)/6.28, beta_all(iFreq), pga, R(pgaIndex)
c          call UpdatePlot ( w_all, Q_all, beta_all, R, nFreq, nDamp, 
c    1                  iter, iDamp, aveMisfit, misfit, iFreq )  

c         If this solution is best so far, then save time history
c         This allows you to keep the best fit in case the solution 
c         starts to diverge.
          if ( misfit .lt. bestMisfit ) then
            call SaveSoln ( acc, nPts, accSave, bestMisfit, misfit )
          endif
          
          lastMisfit = misfit
        enddo

c       Add pga modification (except if finished)
         if ( iter .ne. MaxIter .and. misfit .gt. tol ) then
             call AddPga ( acc, npts, nQ_all, Q, maxFreq, dt, pgaAcc,
	1          	    iModPga )
         endif
        
        iter = iter + 1
      enddo
c      call exitgp (.false.)      
      do i=1,MAX_MATCH
        maxTime_all(i) = npts*dt
      enddo

c     Recompute spectrum for best solution and write output
      call CalcRsp ( accSave, nPts, dt, w_all, beta_all, nQ_all, tPeak, 
	2               R, P, 
     1               zero, maxTime_all )
      call CalcMisfit ( R, Q_all, nQ_all, P, deltaR, misfit, iFreq, 
	1            aveMisfit, freqMatch, w_all )
      call FindPga ( accSave, npts, pga )
      call WrOutTH ( accSave, nPts, dt, misfit, targetFile, thFileout)
      call WrOutSpect ( R, Q_all, nQ_all, w_all, beta_all, iter, misfit,
	1                  maxiter,
     1                  nDamp, damp, version,tol, gamma, iModel,
     2                  a1, a2, f1, f2, scaleFlag, dtFlag, evmin,
     3                  groupSize, maxFreq, thFilein, thFileout, 
     1                  targetFile, pga, R0,
     4                  Q_all_noRand, iSeed, ranFactor, tPeak, file1, 
     5                  scaleperiod,freqBand,npole1,
     6                  iModPGA,FreqMatch,scale1 )

c     Check unmatched frequencies
      do i=1,nQ_all-1
        w_unmatched(i) = (w_all(i)+w_all(i+1))/2.
        Q_unmatched(i) = sqrt( Q_all(i)*Q_all(i+1) )
       Q_all_noRand_unmatched(i)=sqrt(Q_all_noRand(i)*Q_all_noRand(i+1))
      enddo
      nQ_all = nQ_all-1
      call CalcRsp ( accSave, nPts, dt, w_unmatched, beta_all, nQ_all, 
	1               tPeak, R, P, 
     1               zero, maxTime_all )
      call CalcMisfit ( R, Q_all, nQ_all, P, deltaR, misfit, iFreq, 
	1                aveMisfit, freqMatch, w_all )
      call WrOutSpect ( R, Q_unmatched, nQ_all, w_unmatched, beta_all,
	1                  iter, misfit, maxiter,
     1                  nDamp, damp, version,tol, gamma, iModel,
     2                  a1, a2, f1, f2, scaleFlag, dtFlag, evmin,
     3                  groupSize, maxFreq, thFilein, thFileout, 
     1                  targetFile, pga, R0,
     4                  Q_all_noRand, iSeed, ranFactor, tPeak, file1,
     5                  scaleperiod,freqBand,npole1,
     6                  iModPGA,FreqMatch,scale1 )
      
      stop 
      end

c --------------------------------------------------------------

      subroutine DoMatch ( nGroup, Q1, w1, beta1, acc, nPts, dt, misfit,
     1                  R, Q, nQ1, evMin, gamma, lastMisfit, scaleFlag,
     2                  w, beta, nQ, iFreq, flag1, aveMisfit, minTime1,
     3              minTime, maxFreq, minFreq, iter, iDamp, scalePeriod, 
     4                    freqMatch, flag2, maxTime1, maxTime, pgaacc )

      implicit none
      include 'rspMatch.h'
      integer nGroup, nPts, nQ1(1), scaleFlag, iFreq, nQ, flag1, 
	1           pgaIndex, iter, iDamp, flag2
      real Q1(MAXGROUP,1), w1(MAXGROUP,1), beta1(MAXGROUP,1), w(1), 
	1       beta(1)
      real acc(1), dt, R(1), Q(1), evMin, gamma, lastMisfit, misfit,
     1     minTime1(MAXGROUP,1), minTime(1), minFreq, maxFreq, 
     1       scalePeriod, freqMatch(1)
      real freq1, freq2, maxTime1(MAXGROUP,1), maxTime(1), pgaacc(1)

c     Local Variables      
      real Q2(MAX_MATCH), w2(MAX_MATCH), beta2(MAX_MATCH), 
	1         minTime2(MAX_MATCH), maxTime2(MAX_MATCH)
      real tPeak(MAX_MATCH), aveMisfit, scale
      real*8 deltaR(MAX_MATCH)
      integer iGroup, P(MAX_MATCH), i
      
c     Write Header to screen
      write (*,'( /,2x,''Iteration SubGroup  AveMisfit  MaxMisfit  Freq(
	1Hz)  Damping  Pga(g)'')')

c     Loop over each subgroup
      do iGroup=1,nGroup
          call LoadSubgroup ( Q1, w1, beta1, Q2, w2, beta2, iGroup, nQ1,
	1                   minTime1, minTime2, maxTime1, maxTime2 )
	 
c         Calculate misfit for subgroup
          call CalcRsp ( acc, nPts, dt, w2, beta2, nQ1(iGroup), tPeak,
     1                   R, P, minTime2, maxTime2 )
          call CalcMisfit ( R, Q2, nQ1(iGroup), P, deltaR, misfit, 
     1                      iFreq, aveMisfit, freqMatch, w2 )
          write (*,'( 2x,i6,''.'',i1,i7,3x,f10.4,1x,f10.4,2x,f8.3,1x,
	1           3f8.3)') 
     1           iter, iDamp, iGroup, aveMisfit, 
     2      (R(iFreq)-Q2(iFreq))/Q2(iFreq), w2(iFreq)/6.28, beta2(iFreq)

c         Modify time history 
          call CalcNewTH ( acc, nPts, dt, Q2, w2, beta2, nQ1(iGroup), 
	1                 gamma,
     1          R, P, tPeak, deltaR, evMin, freqMatch, pgaacc, iGroup )  
      enddo
      
c     Apply Baseline Correction
      if ( flag2 .eq. 1 ) then
        call Baseline ( acc, npts, dt, flag1 )
      endif
      
c     Band-pass Filter to restrict energy to specifed freq band
      freq1 = maxFreq/6.28
      freq2 = minFreq * 0.6
      call BandPass ( acc, npts, dt, freq1, freq2, 5 )
c      call BandPass ( acc, npts, dt, freq1, 0., 5 )

c     Calculate response spectrum of new time history
      call CalcRsp ( acc, nPts, dt, w, beta, nQ, tPeak, R, P, minTime, 
	1          maxTime )

c     Scale to PGA
      if ( scaleFlag .eq. 1 ) then
        call FindPgaIndex ( nQ, w, pgaIndex, scalePeriod )
        call ScaleToPga ( pgaIndex, acc, npts, R, nQ, Q )
      endif

c     Calulate misfit from target (deltaR and maximum misfit)
      call CalcMisfit ( R, Q, nQ, P, deltaR, misfit, iFreq, aveMisfit, 
	1            freqMatch, w )

      return
      end

c --------------------------------------------------------------

      subroutine LoadQ_array ( nFreq, nDamp, nQ, Q, beta, target2, damp,
	1          iDamp, w )
     
      include 'rspMatch.h'
      integer nFreq, nDamp, nQ, iDamp, pgaIndex
      real Q(1), beta(1), target2(MAX_FREQ,1), damp(1) , w(MAX_MATCH)
      
      j = iDamp
      do i=1,nFreq
        Q(i) = target2(i,j)
        beta(i) = damp(j)
        j = j + 1
        if ( j .gt. nDamp )  j = 1
      enddo
      nQ = nFreq
      
c     Find the index for pga 
c      call FindPgaIndex ( nQ, w, pgaIndex, scalePeriod )
      return
      end

c ----------------------------------------

      subroutine FindPgaIndex ( nQ, w, pgaIndex, scalePeriod )
      integer pgaIndex, nQ
      real w(1), twoPi
      twoPi = 2. * 3.1415926
      
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
      end


c --------------------------------------------------------------

      subroutine FindPga ( acc, nPts, pga )
      integer nPts
      real acc(1)
      
c     Find the pga 
      pga = 0.
      do i=1,nPts
        a1 = abs(acc(i))
        if ( a1 .gt. pga  ) then
          pga = a1
        endif
      enddo
      return
      end

c --------------------------------------------------------------

      subroutine AddPga ( acc, npts, nQ, Q, maxFreq, dt, pgaAcc, 
	1          iModPGA )
      real acc(1), Q(1), maxFreq, pgaAcc(1)
      integer npts, nQ, iModPGA
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
      end
          
c ----------------------------------------------------------------------

      subroutine RandomizeTarget ( iSeed, Q, nQ, ranFactor, target2, 
	1          nDamp, nFreq )
      include 'rspMatch.h'
      real ranFactor, Q(1), x, y
      integer iSeed, nQ
      real target2(MAX_FREQ,MAX_DAMP)
      
      integer i
      
      i1 = 1
      do i=1,nFreq-5
        x = RAN0(iSeed)
        y = (x-0.5) * 2. * ranFactor
        do j=1,nDamp
          target2(i,j) = target2(i,j) * ( 1 + y )
          Q(i1) = Q(i1) * ( 1 + y )
          i1 = i1 + 1
        enddo
      enddo
      return
      end

c ----------------------------------------------------------------------

      FUNCTION RAN0(IDUM)
C     From numerical recipes
C     This function returns a uniform random deviate between 0.0 and 1.0
C     using a function called RANN(ISEED).
C     Note:  For example, to generate an integer number between 1 and 10
C            use j = 1 +int(10. * RAN0(IDUM))

      DIMENSION V(97)
      DATA IFF /0/
      common /temp/iff, iseed
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        ISEED=ABS(IDUM)
        IDUM=1
        DO 11 J=1,97
          DUM=RANN(ISEED)
11      CONTINUE
        DO 12 J=1,97
          V(J)=RANN(ISEED)
12      CONTINUE
        Y=RANN(ISEED)
      ENDIF
      J=1+INT(97.*Y)
C     IF(J.GT.97.OR.J.LT.1)PAUSE  original statement
      IF(J.GT.97.OR.J.LT.1) j=97
      Y=V(J)
      RAN0=Y
      V(J)=RANN(ISEED)
      RETURN
      END

C **********************************************************************
      FUNCTION RANN(ISEED)
C     From numerical recipes
      PARAMETER(IA=7141,IC=54773,IM=259200)
      ISEED=MOD(ISEED*IA+IC,IM)
      RANN=FLOAT(ISEED)/FLOAT(IM)
      RETURN
      END

C-----------------------------------------------------------------------
      subroutine baseline ( acc, nPts, dt, flag1 )
      implicit none
      include 'rspMatch.h'

c     Passed Variables
      real acc(1), dt
      integer nPts, flag1
      
c     Local Variables
      integer i, j, nParam
      real t
      real*8 b(MAXPTS,1), A(MAXPTS,MAXPARAM)
      real*8 xhat(MAXPARAM,1)
      real*8 ATAIAT(MAXPARAM,MAXPTS)

c      common /Temp1/ nParam, ATAIAT

c     Compute ATAIAT (first time only)
      if ( flag1 .eq. 0 ) then 
      
c       Load ATA and AT matrices (Ax=b) for baseline correction
        call LoadA ( A, npts, dt, nParam, MAXPTS )

c       Compute ATAIAT ( part of OLS solution)   
        call ols ( MAXPTS, MAXPARAM, A, npts, nParam, ATAIAT  )
        flag1 = 1
      endif
      
c     Load b vector
      do i=1,nPts
        b(i,1) = acc(i)
      enddo

c     Compute xhat vector
      call mult ( ATAIAT, MAXPARAM, nParam, nPts, b, MAXPTS,
     1              nPts, 1, xhat, MAXPARAM )
     
c     Remove baseline
      t = 0
      do i=1,nPts
        acc(i) = acc(i) - ( xhat(1,1)*t + xhat(2,1)*t**2 )
        t = t + dt
      enddo

         
      return
      end

c ---------------------------------------------------------------------- 

      subroutine LoadA ( A, npts, dt, nParam, MAXPTS )
      implicit none
      integer npts, MAXPTS
      real dt
      real*8 A (MAXPTS,1)

c     Local Variables
      integer i, nParam
      real sumI, t

c     Load A Matrix
      t = 0
      do i=1,npts
        A(i,1) = t
        A(i,2) = t**2
        t = t + dt
      enddo
      
      nParam = 2
      return
      end

c ----------------------------------------------------------------------
      
      subroutine ols ( arow, acol, A, npts, nterm, ATAIAT  )
      include 'rspMatch.h'
C	  implicit none
      integer arow, acol, i, j, nterm, npts
      real*8 A(arow,acol)
      real*8 AT(MAXPARAM,MAXPTS), ATA(MAXPARAM,MAXPARAM), 
	1             ATAI(MAXPARAM,MAXPARAM)
      real*8 ATAIAT(MAXPARAM,MAXPTS), eps, deter, work(MAXPTS)
c     real*8 xhat(acol,1), b(1,1)
      eps = 1.0e-11

c     Form (A transpose A) inverse
      call transpose ( A, arow, npts, nterm, AT, acol )
      call mult ( AT, acol, nterm, npts, A, arow, npts, nterm, ATA,acol)
      do i=1,nterm
        do j=1,nterm
          ATAI(i,j) = ATA(i,j)
        enddo
      enddo
      call simul ( nterm, ATA, work, eps, -1, acol, deter )
      
c     Form ATAIAT
      call mult ( ATA, acol, nterm, nterm, AT, acol, nterm, npts, 
     1            ATAIAT, acol )
c      call mult ( ATAIAT, acol, nterm, npts, b, arow, npts, 1, 
C    1             xhat, acol )
      return
      end 
      
c ----------------------------------------------------------------------
c ---------------------------------------------------------------------- 

      subroutine Calc_delta_a ( delta_a, nPts, dt, b, w, damp, 
     1                          tPeak, alpha, nQ )
      implicit none
      include 'rspMatch.h'
      
c     Passed Variables
      real dt, w(1), damp(1), tPeak(1), alpha(1)
      real*8 b(1), delta_a(1)
      integer nPts, nQ 

c     External Functions
      real f_func

c     Local Variables
      real tau, temp_a(MAXPTS), twoPi
      integer i, j, n1
      twoPi = 2.*3.1415926

c     Initialize delta_A
      do i=1,nPts
        delta_a(i) = 0.0
      enddo

c     Compute delta_a      
      do 20 j=1,nQ
        tau = 0.0 
        do 10 i=1,nPts
          temp_a(i) = b(j) * f_func( tau, tPeak(j), w(j), damp(j),
     1                 alpha(j) )
          tau = tau + dt 
  10    continue
  
c       Apply taper to begining of adjustment
        n1 = int( 0.1 / (dt*(w(j)/twoPi)) )
        call taper ( temp_a, n1 )
        
c       Add to delta_a
        do 15 i=1,nPts
          delta_a(i) = delta_a(i) + temp_a(i)
  15    continue
  20  continue
  
      return
      end
    
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
      end

c ----------------------------------------------------------------------

      subroutine AdjustTH ( acc, nPts, delta_a, gamma )
      implicit none
c     Passed Variables
      real acc(1), gamma 
      real*8 delta_a(1)
      integer nPts 
      
c     Local Variables
      integer i 

c     Compute new acc  
      do 10 i=1,nPts
        acc(i) = acc(i) + delta_a(i) * gamma 
  10  continue
      
      return
      end

c ----------------------------------------------------------------------

      subroutine Calc_C_Matrix ( c, damp, w, tPeak, nQ, alpha )
      implicit none
      include 'rspMatch.h'
c     Passed Variables
      integer nQ
      real damp(1), w(1), tPeak(1), alpha(1)
      real*8 c(MAX_MATCH,1)      

c     Local Variables
      integer i, j
      
      do 20 i=1,nQ
        do 10 j=1,nQ
          call Calc_C ( tPeak(i), tPeak(j), damp(i), damp(j),
     1                  w(i), w(j), alpha(j), c(i,j) )
  10    continue
  20  continue
 
      return
      end

c ----------------------------------------------------------------------

      subroutine LoadDeltaR ( nQ, c, deltaR, MAX_MATCH )
      implicit none
      integer nQ, MAX_MATCH
      real*8 c(MAX_MATCH,1)
      real deltaR(1)
      
      integer i

      do 10 i=1,nQ
        c(nQ+1,i) = deltaR(i)
  10  continue
      
      return
      end
      
      
      
c ----------------------------------------------------------------------

      subroutine Calc_alpha ( w, alpha, nQ )
      implicit none
c     Passed Variables
      real w(1), alpha(1)
      integer nQ

c     External func
      real alpha_func

c     Local Variables      
      integer i
      
      do 10 i=1,nQ
        alpha(i) = alpha_func( w(i) )
  10  continue
      return
      end
      
      
      
      
c ----------------------------------------------------------------------

      real function h_func( w, beta, t )
      implicit none
c     Passed Variables
      real w, beta, t
      
c     Local Variables
      real tmp1, theta
      
      if ( t .lt. 0. ) then
          h_func = 0.0
          return
      endif

      tmp1 = sqrt( 1. - beta**2 )
      theta = w * tmp1 * t
       
      h_func = - w / tmp1 * 
     1        exp( -w*beta*t ) * ( (2.*beta**2 - 1. ) * sin ( theta )
     2          - 2.*beta * tmp1 * cos ( theta ) )
      return
      end

c ----------------------------------------------------------------------

      real function f_func ( tau, tj, w, beta, alpha )
      implicit none
c     Passed Variables
      real tau, tj, w, beta, alpha

c     External Functions and common
      real h_func, func_deltaT, DiracDelta
      integer iModel
      common /Model/ iModel

c     Local Variables
      real t, deltaT, t1, dt, wPrime
      
      deltaT = func_deltaT ( beta, w )
      wPrime = w * sqrt( 1. - beta**2 )
      
c     Check for model (iModel=1 is for Lilhanand and Tseng model)      
      if ( iModel .eq. 1 .or. iModel .eq. 3 ) then
        t = tj - tau 
      elseif ( iModel .eq. 2 .or. iModel .eq. 4 ) then
        t = abs( tj - tau )
      elseif ( iModel .eq. 6 ) then
        t = tau - tj + deltaT
      endif

c     Use Lilhanand and Tseng f function except for pga (f>=1000Hz)
      if ( w .lt. 6280. ) then
        if ( iModel .le. 4 ) then
c         f_func = h_func( w, beta, t )*exp( -abs(tau-tj+deltaT)*alpha)
          f_func = h_func( w, beta, t ) * exp( -abs(tau-tj)*alpha )
        elseif ( iModel .eq. 6 ) then
          f_func = cos( wPrime*t ) * exp( -abs(t)*alpha )
        endif
      else
        f_func = DiracDelta ( t )
      endif
      return
      end
c ----------------------------------------------------------------------

      real function DiracDelta ( x )
      real x
      
      if ( abs(x) .lt. 1.0e-4 ) then
        DiracDelta = 1.
      else
        DiracDelta = 0.
      endif
      return
      end
      
c ----------------------------------------------------------------------

      real function alpha_func ( w )
      implicit none
      real w, alpha, freq, a1, a2, f1, f2
      integer iModel
      common /Model/ iModel
      common /alphaModel/ a1, a2, f1, f2
      freq = w / ( 6.28)

c     Check for Lilhanand and Tseng model
      if ( iModel .eq. 1 .or. iModel .eq. 4 ) then
        alpha_func = 0.
        return
      endif

c     Use Abrahamson model      
      if ( freq .lt. f1 ) then
        alpha = a1 * freq
      elseif ( freq .lt. f2 ) then
        alpha = ( a1 + (freq-f1) * (a2-a1)/(f2-f1) ) * freq
      else
        alpha = a2 * freq
      endif
      alpha_func = alpha  
      return
      end

c ----------------------------------------------------------------------
       
      subroutine Calc_C ( ti, tj, beta_i, beta_j,
     1                    wi, wj, alpha, C )
      
      implicit none
c     Passed Varaibles
      real ti, tj, beta_i, beta_j, wi, wj, alpha
      real*8 C
      
c     External Functions and common
      real Integ1, Integ2, Integ3, func_deltaT, h_func, f_func
      integer iModel
      common /Model/ iModel

c     Local Variables
      real deltaT, t1, tPeak

      
c     Check if this is the pga value (e.g. f>=1000 Hz)
      if ( wj .ge. 6280. ) then
        if ( wi .eq. wj ) then
          C = 1.0
        else
          C = h_func ( wi, beta_i, ti-tj )
        endif
        return
      endif
      if ( wi .ge. 6280. ) then
        C = f_func ( ti, tj, wj, beta_j, alpha ) 
        return
      endif

      deltaT = func_deltaT ( beta_j, wj )
      t1 = tj - deltaT

c     Check which model to use
      if ( iModel .eq. 1 .or. iModel .eq. 3 ) then

c       Evaluate integral for model 1
        if ( ti .ge. tj ) then
          C = Integ1(t1,1., 1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha)
     1       -Integ1(0.,1., 1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha)
     2       +Integ1(tj,1.,-1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha)
     3       -Integ1(t1,1.,-1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha)
        elseif ( ti .gt. t1 ) then
          C = Integ1(t1,1., 1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha)
     1       -Integ1(0.,1., 1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha)
     2       +Integ1(ti,1.,-1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha)
     3       -Integ1(t1,1.,-1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha)
        else
          C = Integ1(ti, 1., 1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha)
     1       -Integ1(0., 1., 1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha)
        endif

      elseif ( iModel .eq. 2 .or. iModel .eq. 4 ) then
      
c       Evaluate integral for model 2
        if ( ti .gt. tj ) then
          C = Integ2(t1, 1., 1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha )
     1       -Integ2(0., 1., 1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha )
     2       +Integ2(tj, 1.,-1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha )
     3       -Integ2(t1, 1.,-1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha )
     4       +Integ2(ti,-1.,-1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha )
     5       -Integ2(tj,-1.,-1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha )
        elseif ( ti .gt. t1 ) then
          C = Integ2(t1, 1., 1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha )
     1       -Integ2(0., 1., 1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha )
     2       +Integ2(ti, 1.,-1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha )
     3       -Integ2(t1, 1.,-1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha )
        else
          C = Integ2(ti, 1., 1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha )
     1       -Integ2(0., 1., 1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha )
        endif
        
        
        elseif ( iModel .eq. 6 ) then
c         Evaluate integral for model 3 (tapered cosine)
          tPeak = tj - deltaT
          if ( ti .gt. tPeak ) then
        C =Integ3( tPeak,  1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha ) -
     1   Integ3( 0.,     1.,deltaT, beta_i, beta_j,wi,wj,ti,tj,alpha ) +
     2   Integ3( ti,    -1.,deltaT, beta_i, beta_j,wi,wj,ti,tj,alpha ) -
     3     Integ3( tPeak, -1.,deltaT,beta_i,beta_j,wi,wj,ti,tj,alpha )
          else
        C = Integ3( ti, 1., deltaT, beta_i, beta_j,wi,wj,ti,tj,alpha ) -
     1      Integ3( 0., 1., deltaT, beta_i, beta_j,wi,wj,ti,tj,alpha )
          endif
      endif   
     
      return
      end
      
c ----------------------------------------------------------------------

      real function Integ2 ( tau, s, s1, deltaT, beta_i, beta_j, 
     1                      wi, wj, ti, tj, alpha )
c     This function computes the integral for c(i,j) for the
c     modified form of Lihanand and Tseng.
c     Includes alpha term and symmetry  (iModel=2)

c     Passed Variables
      real tau, s, s1, deltaT, beta_i, beta_j, wi, wj, ti, tj,
     1     alpha, xi, xj, theta_i, theta_j, B_Prime, C_Prime,
     2     D_Prime

c     Local Variables
      real Ai, Aj, B, C, D, term1 , E, F
      complex G, H, ctmp1, ctmp2, ctmp3, ctmp4, 
     1        carg1, carg2, cInteg
      
c     Set constants
      Ai = sqrt( 1.-beta_i**2 )
      Aj = sqrt( 1.-beta_j**2 )
      theta_i = wi * Ai
      theta_j = wj * Aj
      C = theta_j * abs(tau-tj)
      B = theta_i * ( ti - tau )
c      D = abs( tau - tj )*(wj*beta_j) + wi*beta_i*(ti-tau) 
c     1    + alpha * abs( tau - tj + deltaT ) 
      D = abs( tau - tj )*(wj*beta_j) + wi*beta_i*(ti-tau) 
     1    + alpha * abs( tau - tj ) 
      E = 2.*beta_i**2 - 1.
      F = 2. * beta_i * Ai
      G = cmplx ( -F, -E ) / 2.
      H = cmplx ( F, -E ) / 2.
      
      B_Prime = -theta_i
      C_Prime = -s*theta_j
      D_Prime = -s*wj*beta_j - wi*beta_i - s1*alpha

      carg1 = cmplx ( 0., B+C )
      carg2 = cmplx ( 0., B-C )
      
c     Set up integral terms
      term1 = wi*wj * exp(-D) / (Ai*Aj)
      ctmp1 = cexp( carg1 ) / cmplx ( -D_Prime, B_Prime + C_Prime )
      ctmp2 = cexp( carg2 ) / cmplx ( -D_Prime, B_Prime - C_Prime )
      ctmp3 = cexp( -carg2 ) / cmplx ( -D_Prime, -B_Prime + C_Prime )
      ctmp4 = cexp( -carg1 ) / cmplx ( -D_Prime, -B_Prime - C_Prime )

c     Evaluate integral
      cInteg = term1 * ( G*G*ctmp1 - G*H*ctmp2 - G*H*ctmp3 + H*H*ctmp4 )
      Integ2 = real( cInteg )
      return
      end
c ----------------------------------------------------------------------

      real function func_deltaT ( beta, w )

      real beta, w 

c     Local Variables
      real tmp1      
      
      tmp1 = sqrt( 1. - beta**2 )
      func_deltaT = atan( tmp1/beta ) / ( w * tmp1 )
      return
      end
      
c ----------------------------------------------------------------------
      real function Integ1 ( tau, s, s1, deltaT, beta_i, beta_j, 
     1                      wi, wj, ti, tj, alpha )
c     This function computes the integral for c(i,j) for the
c     modified form of Lihanand and Tseng.
c     Includes alpha term only (no symmetry) (iModel=1)

c     Passed Variables
      real tau, s, s1, deltaT, beta_i, beta_j, wi, wj, ti, tj,
     1     alpha, xi, xj, theta_i, theta_j, B_Prime, C_Prime,
     2     D_Prime

c     Local Variables
      real Ai, Aj, B, C, D, term1 , E, F
      complex G, H, ctmp1, ctmp2, ctmp3, ctmp4, 
     1        carg1, carg2, cInteg
      
c     Set constants
      Ai = sqrt( 1.-beta_i**2 )
      Aj = sqrt( 1.-beta_j**2 )
      theta_i = wi * Ai
      theta_j = wj * Aj
      C = theta_j * (tj-tau)
      B = theta_i * ( ti - tau )
c      D = (tj-tau)*wj*beta_j + wi*beta_i*(ti-tau) 
c     1    + alpha * abs( tau - tj + deltaT ) 
      D = (tj-tau)*wj*beta_j + wi*beta_i*(ti-tau) 
     1    + alpha * abs( tau - tj ) 
      E = 2.*beta_i**2 - 1.
      F = 2. * beta_i * Ai
      G = cmplx ( -F, -E ) / 2.
      H = cmplx ( F, -E ) / 2.
      
      B_Prime = -theta_i
      C_Prime = -theta_j
      D_Prime = -wj*beta_j - wi*beta_i - s1*alpha

      carg1 = cmplx ( 0., B+C )
      carg2 = cmplx ( 0., B-C )
      
c     Set up integral terms
      term1 = wi*wj * exp(-D) / (Ai*Aj)
      ctmp1 = cexp( carg1 ) / cmplx ( -D_Prime, B_Prime + C_Prime )
      ctmp2 = cexp( carg2 ) / cmplx ( -D_Prime, B_Prime - C_Prime )
      ctmp3 = cexp( -carg2 ) / cmplx ( -D_Prime, -B_Prime + C_Prime )
      ctmp4 = cexp( -carg1 ) / cmplx ( -D_Prime, -B_Prime - C_Prime )

c     Evaluate integral
      cInteg = term1 * ( G*G*ctmp1 - G*H*ctmp2 - G*H*ctmp3 + H*H*ctmp4 )
      Integ1 = real( cInteg )
      return
      end
c ----------------------------------------------------------------------

      real function Integ3 ( tau, s, deltaT, beta_i, beta_j, 
     1                      wi, wj, ti, tj, alpha_j )

c     Passed Variables
      real tau, s, deltaT, beta_i, beta_j, wi, wj, ti, tj,
     1     alpha_j, xi, xj, theta_i, theta_j, B_Prime, C_Prime,
     2     D_Prime

c     Local Variables
      real Ai, Aj, B, C, D, term1 , E, F
      complex G, H, ctmp1, ctmp2, ctmp3, ctmp4, 
     1        carg1, carg2, cInteg
            
c     Set constants
      Ai = sqrt( 1.-beta_i**2 )
      Aj = sqrt( 1.-beta_j**2 )
      theta_i = wi * Ai
      theta_j = wj * Aj
      B = theta_j * (tau-tj+deltaT)
      C = theta_i * ( ti - tau )
      D = abs(tau-tj+deltaT)*alpha_j + wi*beta_i*(ti-tau)
      E = 2.*beta_i**2 - 1.
      F = 2. * beta_i * Ai
      G = cmplx ( F, E )
      H = cmplx ( F, -E )
      
      B_Prime = theta_j
      C_Prime = -theta_i
      D_Prime = -s*alpha_j - wi*beta_i 

      carg1 = cmplx ( 0., B+C )
      carg2 = cmplx ( 0., B-C )
      
c     Set up integral terms
      term1 = wi * exp(-D) / (4*Ai)
      ctmp1 = cexp( carg1 ) / cmplx ( -D_Prime, B_Prime + C_Prime )
      ctmp2 = cexp( carg2 ) / cmplx ( -D_Prime, B_Prime - C_Prime )
      ctmp3 = cexp( -carg2 ) / cmplx ( -D_Prime, -B_Prime + C_Prime )
      ctmp4 = cexp( -carg1 ) / cmplx ( -D_Prime, -B_Prime - C_Prime )

c     Evaluate integral
      cInteg = term1 * ( G*ctmp1 + H*ctmp2 + G*ctmp3 + H*ctmp4 )
      Integ3 = real( cInteg )
      return
      end
      
c ----------------------------------------------------------------------
      subroutine rdc (x,npts,iflag)

c     This subroutine removes a dc shift from the data

c     IFLAG = 0    remove the mean
c     IFLAG = 1    remove the mean value of the first 10 points
c     IFLAG = 2    manual set of DC value to be removed

      real x(1),sum,mean
      integer npts,iflag

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
      end

c ----------------------------------------------------------------------

      subroutine CosTaper (x,npts,tb,te)

c     This subroutine tapers the x array
    
      real x(1), arg
      integer npts,tb,te
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
      end

c ----------------------------------------------------------------------

      subroutine pad (x,npts,nmin,m)

c     This subroutine pads the x array to a power of 2 for FFT

      real x(1)
      integer test,npts,m,nmin

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
      end

c ----------------------------------------------------------------------

      subroutine cool ( signi, n, cx )
c     FFT subroutine.
c     signi = -1.  forward transform
c           =  1.  inverse transform
c     n = log base 2 (npts)

      complex cx(1), carg, temp, cw
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
      end
c ----------------------------------------------------------------------

      subroutine BandPass ( acc, nPts, dt, maxFreq, minFreq, nPole  )
      implicit none
      include 'rspMatch.h'

      real acc(1), dt, df, minFreq, maxFreq, fc1, fc2
      integer tb, te, npts, nmin, mpad, i, j, icomp, npts1,
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
	    write (*,'( 2x,''highpass'')') 
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
      end

	  
c ----------------------------------------------------------------------

      subroutine hipass ( fc, nPole, df, npts, cx )

c     This subroutine applies a butterworth high-pass filter
c     in the frequency domain

      real fc, freq, fBar, scale, df
      integer nPole, npts
	  complex cx(1)

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
      end

c ----------------------------------------------------------------------

      subroutine lowpass ( fc, nPole, df, npts, cx )

c     This subroutine applies a butterworth low-pass filter
c     in the frequency domain

      complex cx(1)
      real fc, freq, fBar, amp, df
      integer nPole, npts

      freq = df
      do i=2,npts/2 + 1
          fBar = ( freq / fc )**( 2 * npole )
          amp = sqrt( 1. / (1.+fBar) )
          cx(i) = cx(i) * amp
          cx(npts-i+2) = cx(npts-i+2) * amp
          freq = freq + df
      enddo
      return
      end

c ----------------------------------------------------------------------
      subroutine RdRunParam ( tol, maxIter, gamma, scaleFlag, dtFlag, 
	1       evMin, 
     1      groupSize, maxFreq, freqBand, nPole1, scalePeriod, iModPGA, 
     2                    iSeed, ranFactor, freqMatch, flag2, scale1 )
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
      write (*,'( 2x,''Select model:  1=no taper, no sym'',/,
     1           17x,''2=taper, sym'',/,17x,''3=taper, no sym'',/,
     2           17x,''4=no taper, sym'',/,17x,''6=tapered cosine'')')
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
      write (*,'( 2x,''Enter Initial filters (min, max in Hz) and nPole'
	1')')
      read (5,*) freqBand(1), freqBand(2), nPole1
      write (*,'( 2x,''Mod PGA?'')')
      read (5,*) iModPGA
      write (*,'( 2x,''Randomize target? (0=no), factor'')')
      read (5,*) iSeed, ranFactor
      write (*,'( 2x,''Enter freq bands for matching (fmin, fmax in Hz)'
	1')')
      read (5,*) freqMatch(1), freqMatch(2)
      write (*,'( 2x,''Enter baseline cor flag (1=yes)'')')
      read (5,*) flag2
      write (*,'( 2x,''Enter scale fact '')')
      read (5,*) scale1
      
c     Convert maxFreq to radians
      maxFreq = maxFreq * 2. * 3.1415926
      
      return
      end
      
c ----------------------------------------------------------------------

      subroutine RdTarget ( w, damp, fileName, minFreq, minTime, 
	1     minTime_all,
     1                  beta_all, Q_all, nQ_all, nFreq, nDamp, target2,
     2                      w_all, maxTime, maxTime_all )
      implicit none
      include 'rspMatch.h'

c     Passed Variables
      integer nQ_all, nFreq, nDamp
      real Q(1), w(1), damp(1), minFreq, minTime(1), target2(MAX_FREQ,1)
      real Q_all(1), beta_all(1), minTime_all(1), w_all(1), maxTime(1),
     1     maxTime_all(1)
      character*80 fileName
      
c     Local Variables
      character*80 dummy
      integer i, j, k
	real freq, pi, w1, damp1
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
            
c     Read Target and find minimum frequency
      minFreq = 100.
      do 30 i=1,nFreq
        read (10,*) freq, minTime(i),maxTime(i),(target2(i,j),j=1,nDamp) 
        w(i) = 2.*pi * freq
        if ( freq .lt. minFreq ) minFreq = freq
  30  continue
  
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
      nQ_all = nDamp * nFreq
      
      close (10)
      return
      end
      
c ----------------------------------------------------------------------

      subroutine RdAcc ( acc, nPts, dt, dtFlag, fileName )
      implicit none
      include 'rspMatch.h'
c     Passed Variables
      integer nPts, dtFlag
      character*80 fileName
      real acc(1), dt
      
c     Local Variables
      character*80 dummy
      integer i, i1, j, step, ifactor
     
c     Open Accelerogram file
      write (*,'( 2x,''Enter accelerogram file'')')
      read (5,'( a60)') fileName
      open ( 10, file=fileName, status='old' )
      
c     Read Header 
C      do 10 i=1,4 
        read (10,'( a1)') dummy
C  10  continue
      read (10,*) nPts, dt
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
      endif
      
      return
      end
      
c ----------------------------------------------------------------------

      subroutine WrOutTH ( acc, nPts, dt, misfit, spectFile, thFileout)
      implicit none
c     Passed Variables
      real acc(1), dt, misfit
      integer nPts
      character*80 spectFile, thFileout
      
c     Local Variables
      integer i
      
c     Open output file
      write (*,'( 2x,''Enter accelerogram output file'')')
      read (5,'( a60)') thFileout
      open ( 20,file=thFileout, status='new')
      
c     Write header
      write (20,'( 2x,''Time history matched to spectrum:'',a60)')
     1       spectFile
      write (20,'( i5,f10.4)') nPts, dt 
     
c     Write data
      write (20,'( 5e14.5)') (acc(i),i=1,nPts)
      close (20)
      return
      end
      
c ----------------------------------------------------------------------

      subroutine WrOutSpect ( R, Q, nQ, w, beta, iter, misfit, maxiter, 
     1                      nDamp, damp, version, tol, gamma, iModel,
     2                      a1, a2, f1, f2, scaleFlag, dtFlag, evmin,
     3                      groupSize, maxFreq, thFilein, thFileout, 
     3                      targetFile, pga, R0,
     4                      Q_noRand, iSeed, ranFactor, tPeak, file1,
     5                      scaleperiod,freqBand,npole1,
     6                      iModPGA,FreqMatch,scale1 )
      implicit none
c     Passed Variables
      real R(1), Q(1), w(1), damp(1), misfit, beta(1), R0(1), tPeak(1),
     1     tol, gamma, a1, a2, f1, f2, evmin, maxFreq, pga, Q_noRand(1),
     1      ranFactor, scaleperiod, freqBand(1), freqmatch(1), scale1
      integer nQ, iter, maxiter, nDamp, dtFlag, iModel, scaleFlag, 
	1     groupSize, iSeed, npole1, iModPGA
      character*80 version, thFilein, thFileout, targetFile, file1
      
c     Local Variables
      character*80 fileName
      integer i, j
      real twoPi
      twoPi = 2. * 3.1415926
      
c     Open output file
      write (*,'( 2x,''Enter spectrum output file'')')
      read (5,'( a60)') fileName
      open ( 20,file=fileName, status='new')
      
c     Write header
      write (20,'(2x,''********************************'')')
      write (20,'(2x,''RSPM - Spectral Matching Program'')')
      write (20,'( 2x,''Program Version: '',a25)') version
      write (20,'(2x, ''Based on Code by, N. Abrahamson (1993)'')')
      write (20,'(2x,''********************************'')')

      if ( iter .le. maxiter) then
        write (20,'( 2x,''Solution converged in'',i5,'' iterations'')')
     1        iter
      else
        write (20,'( 2x,''Solution did not converge in maximum'',
     1         '' number of iterations'')')
      endif  
      write (20,'( 2x,''Maximum misfit ='',f8.2,'' percent'')') 
     1       misfit*100.
      write (20,'( 2x,''Input Parameter File:'',a80)') file1
      write (20,'( /,2x,''Run Parameters:'')')
      write (20,'( 2x,''Input time history:  '',a80)') thfilein
      write (20,'( 2x,''Output time history: '',a80)') thfileout
      write (20,'( 2x,''Target spectrum:     '',a80)') targetfile
      write (20,'( 2x,''Maximum number of iterations ='',i5)') maxIter
      write (20,'( 2x,''Convergence tolerance ='',f10.4)') tol
      write (20,'( 2x,''Convergence Damping factor ='',f10.4)') gamma
      write (20,'( 2x,''Adjustment time history model ='',i5)') iModel
      write (20,'( 2x,''Adjustment TH tapers (a1,a2,f1,f2) ='',4f10.4)')
	1      a1, a2, f1, f2
      write (20,'( 2x,''Interpolate TH Flag ='',i5)') dtFlag
      write (20,'( 2x,''Scale TH Flag ='',i5)') scaleFlag
      if (scaleFlag .eq. 0) then
         write (20,'(2x,''No Scaling Applied'')')
      elseif (scaleFlag .eq. 1) then
         write (20,'(2x,''Scaling Applied'')')
         write (20,'(2x,''Period of Scaling = '',f8.3)') scaleperiod
      elseif (scaleFlag .eq. 2) then
         write (20,'(2x,''Scaling Applied Only Once'')')
         write (20,'(2x,''Period of Scaling = '',f8.3)') scaleperiod
      endif
      write (20,'( 2x,''minimum eigenvalue ='',e12.4)') evmin
      write (20,'( 2x,''Group Size ='',i5)') groupSize
      write (20,'( 2x,''Max Freq (Hz) ='',f10.4)') maxFreq
      write (20,'( 2x,''Filter Parameters: '')')
      write (20,'( 2x,''F1, F2, Npole: '',f8.2,f8.2,i5)') freqBand(1),
     1             freqBand(2),nPole1
      if (iModPGA .eq. 0) then
         write (20,'( 2x,''PGA Not Modified'')')
      elseif (iModPGA .eq. 1) then   
         write (20,'( 2x,''PGA Modified'')')
      endif
      write (20,'(2x,''Frequency Range of Spectral Matching: '',2f8.2)')
     1       freqMatch(1), freqMatch(2)
      write (20,'(2x,''Scaling Factor for Initial Time History = '',
     1           f8.2)') scale1
      write (20,'( 2x,''PGA ='',f10.4)') pga
      write (20,'( /,2x,''Matched spectrum:'')')
      write (20,'( i5,'' = Number of spectral frequencies'')') nQ
      write (20,'( i5,2x,''= Number of Damping Values'')') nDamp
      write (20,'( 10f10.3,2x,''= Damping Values'')')(damp(i),i=1,nDamp)
      write (20,'( 2x,''Random flags'',i5,f10.4)') iSeed, ranFactor
      write (20,'( /,4x,'' Freq     Damping   Target   Computed   Initia
	1l  Randomized  tPeak'')')

c     Write Target and computed spectra
      do j=1,nDamp
        do i=j,nQ,nDamp
          write (20,'( 7f10.4)') w(i)/twoPi, beta(i), Q_noRand(i), R(i),
	1          R0(i), Q(i), tPeak(i)
        enddo
      enddo
      close (20)
      return
      end

c ----------------------------------------------------------------------

      subroutine transpose ( A, arow, m, n, At, atrow )
      integer arow, atrow
      real*8 A(arow,1), At(atrow,1)

      do i=1,m
        do j=1,n
          At(j,i) = A(i,j)
        enddo
      enddo
      return
      end

c ----------------------------------------------------------------------
      subroutine ztilde (zrow,m,n,z)
c     This subroutine converts a matrix from partitioned format
c     to tilde format.

      integer zrow
      real*8 z(zrow,1)
   
      do 100 i=1,m/2
        do 90 j=1,n
          z(i,j+n) = -z(i+m/2,j)
          z(i+m/2,j+n) = z(i,j)
  90    continue
 100  continue
      n = 2*n
      return
      end
c ----------------------------------------------------------------------

      subroutine mult (a,arow,m1,n1,b,brow,m2,n2,c,crow)
c     This subroutine multiplies two matrices.
c       C = A B
c     where A is m1 x n1  and  B is m2 x n2

      INTEGER AROW, BROW, CROW
      double precision a(arow,1), b(brow,1), c(crow,1)

      if (n1 .ne. m2) stop 99
      do 100 i=1,m1
        do 90 j=1,n2
          c(i,j) = 0.0
          do 80 ii=1,n1
            c(i,j) = c(i,j) + a(i,ii)*b(ii,j)
  80      continue
  90    continue
  100 continue
      return
      end

c ----------------------------------------------------------------------

      subroutine svdcmp ( A, m, n, mp, np, W, V )

      include 'rspMatch.h'      
c      parameter (NMAX=100)
C      implicit none
      integer m,n,mp,np, nmax
      real*8 A(mp,np), W(np), V(np,np), rv1(MAX_MATCH)
      integer i, j, k, l, jj, ii, its, nm
      real*8 anorm, c, f, g, h, s, x, y, z, scale
      if (m .lt. n ) pause 'You must augment A with extra zero rows.'
      NMAX = MAX_MATCH
	  
      g = 0.0
      scale = 0.0
      anorm = 0.0
      do 25 i=1,n
        l = i + 1
        rv1(i) = scale*g
        g = 0.0
        s = 0.0
        scale = 0.0
        if ( i .le. m) then
          do 11 k=i,m
            scale = scale + abs(A(k,i))
  11      continue
          if ( scale .ne. 0.0 ) then
            do 12 k=i,m
              A(k,i) = A(k,i) / scale
              s = s + A(k,i)*A(k,i)
  12        continue
            f = A(i,i)
            g = -sign(sqrt(s),f)
            h = f*g - s
            A(i,i) = f-g
            if ( i .ne. n ) then
              do 15 j=l,n
                s = 0.0
                do 13 k=i,m
                  s = s + A(k,i)*A(k,j)
  13            continue
                f = s/h
                do 14 k=i,m
                  A(k,j) = A(k,j) + f*A(k,i)
  14            continue
  15          continue
            endif
            do 16 k=i,m
              A(k,i) = scale * A(k,i)
  16        continue
          endif
        endif
        
        W(i) = scale * g
        g = 0.0
        s = 0.0
        scale = 0.0
        if ( (i .le. m) .and. (i .ne. n) ) then
          do 17 k=l,n
            scale = scale + abs(A(i,k))
  17      continue
          if ( scale .ne. 0.0) then
            do 18 k=l,n
              A(i,k) = A(i,k) / scale
              s = s + A(i,k) * A(i,k)
  18        continue
            f = A(i,l)
            g = -sign( sqrt(s), f)
            h = f*g - s
            A(i,l) = f-g
            do 19 k=l,n
              rv1(k) = A(i,k) / h
  19        continue
            if ( i .ne. m ) then
              do 23 j=l,m
                s = 0.0
                do 21 k=l,n
                  s = s + A(j,k)*A(i,k)
  21            continue
                do 22 k=l,n
                  A(j,k) = A(j,k) + s*rv1(k)
  22            continue
  23          continue
            endif
            do 24 k=l,n
              A(i,k) = scale * A(i,k)
  24        continue
          endif
        endif
        anorm = max( anorm, abs(W(i)) + abs(rv1(i)) )
  25  continue
  
c     Accumulation of right-hand tranforms
      do 32 i=n,1,-1
        if ( i .lt. n ) then
          if ( g .ne. 0.0 ) then
            do 26 j=l,n
              V(j,i) = (A(i,j) / A(i,l)) / g
  26        continue
            do 29 j=l,n
              s = 0.0
              do 27 k=l,n
                s = s + A(i,k)*V(k,j)
  27          continue
              do 28 k=l,n
                V(k,j) = V(k,j) + s*V(k,i)
  28          continue
  29        continue
          endif
          do 31 j=l,n
            V(i,j) = 0.0
            V(j,i) = 0.0
  31      continue
        endif
        V(i,i) = 1.0
        g = rv1(i)
        l = i
  32  continue
      
c     Accumulation of left-hand transformations.
      do 39 i=n,1,-1
        l = i+1
        g = W(i)
        if ( i .lt. n ) then
          do 33 j=l,n
            A(i,j) = 0.0
  33      continue
        endif
        
        if ( g .ne. 0.0 ) then
          g = 1.0 / g
          if ( i .ne. n ) then
            do 36 j=l,n
              s = 0.0
              do 34 k=l,m
                s = s + A(k,i)*A(k,j)
  34          continue
              f = (s/A(i,I)) * g
              do 35 k=i,m
                A(k,j) = A(k,j) + f*A(k,i)
  35          continue
  36        continue
          endif
          do 37 j=i,m
            A(j,i) = A(j,i) * g
  37      continue
        else
          do 38 j=i,m
            A(j,i) = 0.0
  38      continue
        endif
        A(i,i) = A(i,i) + 1.0
  39  continue
            
c     Diagonalization of the bidiagonal form
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm = l-1
            if ( (abs(rv1(l)) + anorm) .eq. anorm ) goto 2
            if ( (abs(W(nm)) + anorm) .eq. anorm ) goto 1
 41       continue
  1       c = 0.0
          s = 1.0
          do 43 i=l,k
            f = s*rv1(i)
            rv1(i) = c * rv1(i)
            if ( abs(f)+anorm .eq. anorm ) goto 2
            g = W(i)
            h = sqrt( f*f + g*g)
            W(i) = h
            h = 1.0 / h
            c = g*h
            s = -(f*g)
            do 42 j=1,m
              y = A(j,nm)
              z = A(j,i)
              A(j,nm) = (y*c) + (z*s)
              A(j,i) = -(y*s) + (z*c)
  42        continue
  43      continue
   2      z = W(k)
          if ( l .eq. k ) then
            if ( z .lt. 0.0 ) then
              W(k) = -z
              do 44 j=1,n
                V(j,k) = -V(j,k)
  44          continue
            endif
            goto 3
          endif
          if ( its .eq. 30 ) pause 'no convergence in 30 iterations'
          x = W(l)
          nm = k-1
          y = W(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ( (Y-z)*(y+z) + (g-h)*(g+h) ) / ( 2.0*h*y)
          g = sqrt( f*f + 1.0 )
          f = ( (x-z)*(x+z) + h*((y/(f+sign(g,f)))-h))/x
          
c         Next QR transformation
          c = 1.0
          s = 1.0
          do 47 j=l,nm
            i = j+1
            g = rv1(i)
            y = W(i)
            h = s*g
            g = c*g
            z = sqrt( f*f + h*h )
            rv1(j) = z
            c = f / z
            s = h / z
            f = (x*c) + (g*s)
            g = -(x*s) + (g*c)
            h = y*s
            y = y*c
            do 45 jj=1,n
              x = V(jj,j)
              z = V(jj,i)
              V(jj,j) = (x*c) + (z*s)
              V(jj,i) = -(x*s) + (z*c)
  45        continue
            
            z = sqrt(f*f + h*h)
            W(j) = z
            if ( z .ne. 0.0 ) then
              z = 1.0 / z
              c = f*z
              s = h*z
	      endif
            f = (c*g) + (s*y)
            x = -(s*g) + (c*y)
            do 46 jj=1,m
              y = A(jj,j)
              z = A(jj,i)
              A(jj,j) = (y*c) + (z*s)
              A(jj,i) = -(y*s) + (z*c)
  46        continue
  47      continue
          rv1(l) = 0.0
          rv1(k) = f
          W(k) = x
  48    continue
   3    continue
  49  continue
      return
      end
	  
c ----------------------------------------------------------------------

      subroutine SvdInv ( W, U, V, n, AInv )
      include 'rspMatch.h'	  
C      implicit none
      integer n
      real*8 W(1), U(MAX_MATCH,1), V(MAX_MATCH,1), AInv(MAX_MATCH,1)
      integer i, j
	  real*8 precission, wMax, W_Inv(MAX_MATCH,MAX_MATCH),
     1       Ut(MAX_MATCH,MAX_MATCH), VW(MAX_MATCH,MAX_MATCH)
	  precission = 1.0e-10
	  
c     Find largest eigenvalue
      wMax = 0.0
      do i=1,n
	    if ( abs(W(i)) .gt. wMax) wMax = abs(W(i))
      enddo
	  
c     Set W inverse
      do i=1,n
	    do j=1,n
		  W_Inv(i,j) = 0.0
		enddo
	    if ( abs(W(i)/wMax) .gt. precission ) then
	      W_Inv(i,i) = 1./W(i)
		else
		  W_Inv(i,i) = 0.0
		endif
      enddo

c     Find U transpose
      call Transpose ( U, MAX_MATCH, n, n, Ut, MAX_MATCH )

c     Compute inverse
      call mult ( V, MAX_MATCH, n, n, W_Inv, MAX_MATCH, n, n, VW, 
	1      MAX_MATCH )
      call mult ( VW, MAX_MATCH, n, n, Ut, MAX_MATCH, n, n, AInv, 
	1      MAX_MATCH )
	  
	  return
	  end
	  
c ----------------------------------------------------------------------

      subroutine Svbksb ( U, W, V, m, n, mp, np, b, x, evMin )
      include 'rspMatch.h'
c	  parameter ( NMAX=100 )
C      implicit none
      integer m, n, mp, np, i, j, jj, NMAX
	  real*8 U(mp,np), W(np), V(np,np), b(mp), x(np), tmp(MAX_MATCH)
	  real*8 s, precission, wMax, evMIn
	  NMAX = MAX_MATCH

c     This routine solves Ax=b using U, W, and V from SVDCMP.
	  precission = evMin
	  
c     Find largest eigenvalue
      wMax = 0.0
      do i=1,n
	    if ( abs(W(i)) .gt. wMax) wMax = abs(W(i))
      enddo
c	  write (*,'( 2x,''wMax ='',e12.4)') wMax
	  
c     Zero out small values
      do i=1,n
	    if ( abs(W(i)/wMax) .lt. precission ) then
		  W(i) = 0.0
c		  write (*,'( 2x,''Zero eigenvalue'',i5)') i
		endif
      enddo

c     Calculate UT B	  
	  do 12 j=1,n
	    s = 0.0
		if ( W(j) .ne. 0. ) then
		  do 11 i=1,m
		    s = s + U(i,j)*b(i)
  11      continue
          s = s / W(j)
		endif
		tmp(j) = s
  12  continue

c     Multiply by V
      do 14 j=1,n
        s = 0.
		do 13 jj=1,n
		  s = s + V(j,jj)*tmp(jj)
  13    continue
        x(j) = s
  14  continue
      return
	  end
      
c ----------------------------------------------------------------------	  
	  
       SUBROUTINE SIMUL( N, A, X, EPS, INDIC, NRC , DETER )
C
C        WHEN INDIC IS NEGATIVE, SIMUL COMPUTES THE INVERSE OF THE N BY
C        N MATRIX A IN PLACE.  WHEN INDIC IS ZERO, SIMUL COMPUTES THE
C        N SOLUTIONS X(1)...X(N) CORRESPONDING TO THE SET OF LINEAR
C        EQUATIONS WITH AUGMENTED MATRIX OF COEFFICIENTS IN THE N BY
C        N+1 ARRAY A AND IN ADDITION COMPUTES THE INVERSE OF THE
C        COEFFICIENT MATRIX IN PLACE AS ABOVE.  IF INDIC IS POSITIVE,
C        THE SET OF LINEAR EQUATIONS IS SOLVED BUT THE INVERSE IS NOT
C        COMPUTED IN PLACE. THE GAUSS-JORDAN COMPLETE ELIMINATION METHOD
C        IS EMPLOYED WITH THE MAXIMUM PIVOT STRATEGY.  ROW AND COLUMN
C        SUBSCRIPTS OF SUCCESSIVE PIVOT ELEMENTS ARE SAVED IN ORDER IN
C        THE IROW AND JCOL ARRAYS RESPECTIVELY.  K IS THE PIVOT COUNTER,
C        PIVOT THE ALGEBRAIC VALUE OF THE PIVOT ELEMENT, MAX
C        THE NUMBER OF COLUMNS IN A AND DETER THE DETERMINANT OF THE
C        COEFFICIENTS MATRIX.  THE SOLUTIONS ARE COMPUTED IN THE (N+1)TH
C        COLUMN OF A AND THEN UNSCRAMBLED AND PUT IN PROPER ORDER IN
C        X(1)...X(N) USING THE PIVOT SUBSCRIPT INFORMATION AVAILABLE
C        IN THE IROW AND JCOL ARRAYS.  THE SIGN OF THE DETERMINANT IS
C        ADJUSTED, IF NECESSARY, BY DETERMINING IF AN ODD OR EVEN NUMBER
C        OF PAIRWISE INTERCHANGES IS REQUIRED TO PUT THE ELEMENTS OF THE
C        JORD ARRAY IN ASCENDING SEQUENCE WHERE JORD(IROW(I)) = JCOL(I).
C        IF THE INVERSE IS REQUIRED, IT IS UNSCRAMBLED IN PLACE USING
C        Y(1)...Y(N) AS TEMPORARY STORAGE.  THE VALUE OF THE DETERMINANT
C        IS RETURNED AS THE VALUE OF THE FUNCTION.  SHOULD THE POTENTIAL
C        PIVOT OF LARGEST MAGNITUDE BE SMALLER IN MAGNITUDE THAN EPS,
C        THE MATRIX IS CONSIDERED TO BE SINGULAR AND A TRUE ZERO IS
C        RETURNED AS THE VALUE OF THE FUNCTION.
C
      IMPLICIT REAL*8(A-H, O-Z)
      real*8 A(nrc,nrc), x(N), Y(200), eps, deter
      integer irow(200), jcol(200), jord(200)
C
      MAX = N
      IF ( INDIC.GE.0 )  MAX = N + 1
C
C     ..... IS N LARGER THAN 78 .....
      IF ( N.LE.78 )  GO TO 5
      WRITE (6,200)
      DETER = 0.
      RETURN
C
C     ..... BEGIN ELIMINATION PROCEDURE .....
 5    DETER = 1.
      DO 18 K = 1, N
      KM1 = K - 1
C     ..... SEARCH FOR THE PIVOT ELEMENT .....
      PIVOT = 0.
      DO 11 I = 1, N
      DO 11 J = 1, N
C     ..... SCAN IROW AND JCOL ARRAYS FOR INVALID PIVOT SUBSCRIPTS .....
      IF ( K.EQ.1 ) GO TO 9
      DO 8 ISCAN = 1, KM1
      DO 8 JSCAN = 1, KM1
      IF ( I.EQ.IROW(ISCAN) ) GO TO 11
      IF ( J.EQ.JCOL(JSCAN) ) GO TO 11
 8    CONTINUE
 9    CONTINUE
      PIVOT = A(I,J)
      IROW(K) = I
      JCOL(K) = J
 11   CONTINUE
C
C     ..... INSURE THAT SELECTED PIVOT IS LARGER THAN EPS .....
      IF ( DABS(PIVOT).GT.EPS ) GO TO 13
      DETER = 0.
      RETURN
C
C     ..... UPDATE THE DETERMINANT VALUE .....
 13   IROWK = IROW(K)
      JCOLK = JCOL(K)
      DETER = DETER*PIVOT
C
C     ..... NORMALIZE PIVOT ROW ELEMENTS .....
      DO 14 J = 1, MAX
 14   A(IROWK,J) = A(IROWK,J)/PIVOT
C
C     ..... CARRY OUT ELIMINATION AND  DEVELOP INVERSE .....
      A(IROWK,JCOLK) = 1./PIVOT
      DO 18 I = 1, N
      AIJCK = A(I,JCOLK)
      IF ( I.EQ.IROWK ) GO TO 18
      A(I,JCOLK) = - AIJCK/PIVOT
      DO 17 J = 1, MAX
 17   IF ( J.NE.JCOLK ) A(I,J) = A(I,J) - AIJCK*A(IROWK,J)
 18   CONTINUE
C
C     ..... ORDER SOLUTION VALUES (IF ANY) AND CREATE JORD ARRAY .....
      DO 20 I = 1, N
      IROWI = IROW(I)
      JCOLI = JCOL(I)
      JORD(IROWI) = JCOLI
 20   IF ( INDIC.GE.0 ) X(JCOLI) = A(IROWI,MAX)
C
C     ..... ADJUST SIGN OF DETERMINANT .....
      INTCH = 0
      NM1 = N - 1
      DO 22 I = 1, NM1
      IP1 = I + 1
      DO 22 J = IP1,N
      IF ( JORD(J).GE.JORD(I) ) GO TO 22
      JTEMP = JORD(J)
      JORD(J) = JORD(I)
      JORD(I) = JTEMP
      INTCH = INTCH + 1
 22   CONTINUE
      IF( INTCH/2*2.NE.INTCH ) DETER = - DETER
C
C     ..... IF INDIC IS POSITIVE RETURN WITH RESULTS .....
      IF ( INDIC.LE.0 ) GO TO 26
c     DETER = DETER
      RETURN
C
C     ..... IF INDIC IS NEGATIVE OR ZERO, UNSCRAMBLE THE INVERSE
C           FIRST BY ROWS .....
 26   DO 28 J = 1, N
      DO 27 I = 1, N
      IROWI = IROW(I)
      JCOLI = JCOL(I)
 27   Y(JCOLI) = A(IROWI,J)
      DO 28 I = 1, N
 28   A(I,J) = Y(I)
C     ..... THEN BY COLUMNS .....
      DO 30 I = 1, N
      DO 29 J = 1, N
      IROWJ = IROW(J)
      JCOLJ = JCOL(J)
 29   Y(IROWJ) = A(I,JCOLJ)
      DO 30 J = 1, N
 30   A(I,J) = Y(J)
C
C     ..... RETURN FOR INDIC NEGATIVE OR ZERO .....
c     DETER = DETER
      RETURN
C
C     ..... FORMAT FOR OUTPUT STATEMENT .....
 200  FORMAT( 10H0N TOO BIG )
C
      END
	  
c ----------------------------------------------------------------------

      subroutine mdian1 ( x, n, xmed )
c     Routine to find the median of a distribution
c     Passed Variables
      integer n
      real x(1), xmed

c     Local Variables      
      integer i, icount
      real xmin, xmax
      
c     Find min and Max
      xmin = 1.0e30
      xmax = -1.0e30
      do  20 i=1,n
        if ( x(i) .gt. xmax ) xmax = x(i)
        if ( x(i) .lt. xmin ) xmin = x(i)
  20  continue
      
c     Find median by iteration
  10  xmed = (xmax + xmin) / 2.
      icount = 0
      do 30 i=1,n
        if ( x(i) .lt. xmed ) icount = icount + 1
  30  continue
      if ( icount .eq. n/2 ) return
      if ( icount .gt. n/2 ) then
          xmax = xmed
      else
          xmin = xmed
      endif
      goto 10          
      
      end
      
c ----------------------------------------------------------------------

      subroutine CalcScale ( R, Q, nQ, median )
      implicit none
      include 'rspMatch.h'

c     Passed Variables
      real R(1), Q(1), median
      integer nQ
      
c     Local Variables
      real ratio(MAX_MATCH)
      integer i
      
c     Convert to ratio
      do i=1,nQ
        ratio(i) = Q(i) / R(i)
      enddo
      
c     Find Median
      call mdian1 ( ratio, nQ, median )
      
      return
      end
      
c ----------------------------------------------------------------------

      subroutine CalcMisfit ( R, Q, nQ, polarity, deltaR, misfit, 
     1                        iFreq, aveMisfit, freqMatch, w )
      implicit none
c     Passed Variables
      integer nQ, polarity(1)
      real R(1), Q(1), misfit, aveMisfit, freqMatch(1), w(1)
      real*8 deltaR(1)
      
c     Local Variables
      integer i, iFreq
      real temp, sum, maxDeltaR, w1, w2
      maxDeltaR = 0.1
      w1 = freqMatch(1) * 2. * 3.1415926
      w2 = freqMatch(2) * 2. * 3.1415926
      
      misfit = 0
      sum = 0.0
      do i=1,nQ
        deltaR(i) = ( Q(i) - R(i) ) * polarity(i)
        temp = abs( deltaR(i) / Q(i) )
        sum = sum + temp 
        if ( temp .gt. misfit .and. w(i) .ge. w1 .and. w(i) .le. w2 ) 
     1          then
            misfit = temp
            iFreq = i
        endif
        
        
c       Restrict DeltaR to be less than "maxdeltaR" of R
        if ( abs(deltaR(i)) .gt. R(i)*maxDeltaR ) then
          deltaR(i) = R(i)*maxDeltaR * deltaR(i) / abs(deltaR(i))
        endif
      enddo
      aveMisfit = sum / nQ
      
      return
      end
      
c ----------------------------------------------------------------------

      subroutine ApplyScale ( y, n, scale )
      implicit none
c     Passed Variables
      integer n
      real y(1), scale
      
c     Local Variables
      integer i
      
      do i=1,n
        y(i) = y(i) * scale
      enddo
      
      return
      end
      
c ----------------------------------------------------------------------
      
      subroutine InterpFreq ( Q, freq, damp, nQ, n1 )
      include 'rspMatch.h'
c     Passed Variables
      real Q(1), freq(1), damp(1)
      integer nQ
      
c     Local Variables
      integer i, j, n1, k
      real freq1(MAX_MATCH), Q1(MAX_MATCH)

c     Interpolate (linear on log-log)
      k = 1
      do i=1,nQ-1
        do j=0,n1-1
          freq1(k) = log(freq(i)) + (log(freq(i+1))-log(freq(i))) * 
	1                float(j)/n1
          Q1(k) = log(Q(i)) + (log(Q(i+1))-log(Q(i))) * float(j)/n1
          k = k+1
        enddo
      enddo
      freq1(k) = log(freq(nQ))
      Q1(k) = log(Q(nQ))

c     Copy back to original array      
      nQ = k
      do i=1,nQ
        freq(i) = exp( freq1(i) )
        Q(i) = exp( Q1(i) )
        damp(i) = damp(1)
      enddo
      
      return
      end

c ----------------------------------------------------------------------


      subroutine SubdivideTarget ( nQ, Q, w, beta, nQ1, Q1, w1, beta1, 
	1     nGroup, 
     1     groupSize, minTime, minTime1, maxFreq, maxTime, maxTime1 )
      implicit none
      include 'rspMatch.h'
c     Passed Variables
      integer nQ, nQ1(1), nGroup, groupSize
      real Q(1), w(1), beta(1), Q1(MAXGROUP,1), w1(MAXGROUP,1),
     1     beta1(MAXGROUP,1), minTime(1), minTime1(MAXGROUP,1),
     2     maxFreq, maxTime(1), maxTime1(MAXGROUP,1)

c     Local Variables
      integer i, j, k, i1

      nGroup = int( float(nQ) / groupSize + 0.5 )
c      write (*,'( 2x,''nGroup ='',i5)') nGroup
      
c     Initialize nQ1
      do i=1,nGroup
        nQ1(i) = 0
      enddo
      
      i1 = 0
      do i=0,nQ-1
        if ( w(i+1) .lt. maxFreq .or. w(i+1) .gt. 6280. ) then
          j = i1 - (i1/nGroup)*nGroup + 1
          k = i1/nGroup + 1
          w1(j,k)    = w(i+1)
          beta1(j,k) = beta(i+1)
          Q1(j,k)    = Q(i+1)
          minTime1(j,k) = minTime(i+1)
          maxTime1(j,k) = maxTime(i+1)
          nQ1(j) = nQ1(j) + 1
          i1 = i1 + 1
        endif
      enddo
      
      return
      end
      
c ---------------------------------------------------------------------- 

      subroutine LoadSubgroup (Q1, w1, beta1, Q2, w2, beta2, iGroup,nQ1,
     1                          minTime1, minTime2, maxTime1, maxTime2 )
      implicit none
      include 'rspMatch.h'
c     Passed Variables
      integer nQ1(1), iGroup
      real Q2(1), w2(1), beta2(1), Q1(MAXGROUP,1), w1(MAXGROUP,1),
     1     beta1(MAXGROUP,1), minTime1(MAXGROUP,1), minTime2(1),
     2     maxTime1(MAXGROUP,1), maxTime2(1)

c     Local Variables
      integer i, j, k
      
      do i=1,nQ1(iGroup)
        w2(i)    = w1(iGroup,i)
        beta2(i) = beta1(iGroup,i)
        Q2(i)    = Q1(iGroup,i)
        minTime2(i) = minTime1(iGroup,i)
        maxTime2(i) = maxTime1(iGroup,i)
      enddo
      
      return
      end
c ----------------------------------------------------------------------
 
      subroutine SaveSoln ( acc, nPts, accSave, bestMisfit, misfit )
      implicit none
      real acc(1), accSave(1), bestMisfit, misfit
      integer nPts, i
      
      write (*,'( 2x,''Saving Solution'')')

c     Copy acc array and reset bestMisfit
      do i=1,nPts
        accSave(i) = acc(i)
      enddo
      bestMisfit = misfit
      
      return
      end

c ----------------------------------------------------------------------

      subroutine ScaleToPga ( pgaIndex, acc, npts, R, nQ, Q )
      implicit none
      
      integer pgaIndex, npts, nQ
	real acc(1), R(1), Q(1), scale

      scale = (Q(pgaIndex) / R(pgaIndex)) 
      write (*,'( 2x,''Scale to PGA, scale factor ='',f10.4)') scale
      
c     Apply scale factor to accelerogram and spectrum
      call ApplyScale ( acc, nPts, scale)
      call ApplyScale ( R, nQ, scale)

      return
      end
      
c ----------------------------------------------------------------------
      
      subroutine CalcRsp ( acc, npts, dt, w, damping, 
     +                     nQ, time, SA, polarity, minTime, maxTime)

c     This subroutine calculates the response spectrum 

      implicit none

      real w(1), damping(1), acc(1), dt, SA(1), time(1), minTime(1), 
	1           maxTime(1)
      integer npts, nQ, i, j , polarity(1), timeIndex


c     LOOP FOR EACH FREQUENCY
      do i=1,nQ
	  
	    if ( w(i) .gt. 628. ) then
		  call CalcPGA ( acc, npts, sa(i), timeIndex, polarity(i) )
		else
 	      call CalcOneRsp ( acc, npts, dt, w(i), damping(i), sa(i),
     1                    timeIndex, polarity(i), minTime(i),maxTime(i))
        endif
		time(i) = (timeIndex-1) * dt
	  enddo
	  
	  return
	  end
	  
c ----------------------------------------------------------------------

      subroutine CalcPGA ( acc, npts, sa, timeIndex, polarity  )

      real acc(1), sa
	  integer npts, polarity, timeIndex
	  
	  accMax = 0.0
	  do i=1,npts
	    if ( abs(acc(i)) .gt. accMax ) then
          accMax = abs(acc(i) )
		  sa = accMax
		  if ( acc(i) .ge. 0. ) then
  		    polarity = 1.
		  else 
		    polarity = -1
		  endif
		  timeIndex = i
		endif
      enddo
	  write (*,'( 2x,''pga ='',f10.4)') sa
	  return
	  end
	    
c ----------------------------------------------------------------------

      subroutine CalcOneRsp ( acc, npts, dt, w, damping, 
     +                        sa, timeIndex, polarity, minTime, maxTime)
      implicit none
      include 'rspMatch.h'

      real rsp(MAXPTS), w, damping, acc(1), dt, SA, minTime, maxTime
      integer npts, i, j, polarity, timeIndex, iTime, iTime2
	  
c     Compute coeff 
      call coeff ( w, damping, dt )

c     CALCULATE THE RESPONSES
      call brs ( acc, w, damping, npts, rsp )

c     FIND THE MAXIMUM OF THE RESPONSE
      SA = 0.0
	  iTime = int( minTime/dt) + 1
	  iTime2 = int( maxTime/dt) -1
      do j=iTime,iTime2
          if (abs(rsp(j)) .gt. SA) then
		     SA = abs(rsp(j))
			 timeIndex = j
			 if ( rsp(j) .ge. 0. ) then
			    polarity = 1
			 else
			    polarity = -1
			 endif
		  endif
      enddo

      return
      end


c ----------------------------------------------------------------------

      subroutine coeff ( w, beta1, dt1 )
      implicit none
      real beta1, dt1, w
      real*8 a11, a12, a21, a22, b11, b12, b21, b22
      real*8 beta, dt, t1, t2, t3, t4, s1, s2
      common /coef/a11,a12,a21,a22,b11,b12,b21,b22

      beta = dble( beta1 )
      dt = dble( dt1 )

c
c     Set up repeated terms
      t1 = sqrt(1.-beta**2)
      t2 = sin (w*t1*dt)
      t3 = cos (w*t1*dt)
      t4 = exp (-beta*w*dt)
      s1 = (2.*beta**2-1.) / (w**2*dt)
      s2 = 2.*beta / (w**3*dt)
c	  write (*,'( 10f10.5)') t1, t2, t3, t4, s1, s2
c
c     calculate the a's
      a11 = t4*(beta*t2/t1+t3)
      a12 = t4*t2 / (w*t1)
      a21 = -t4*w*t2 / t1
      a22 = t4*(t3-beta*t2/t1)
c
c     calculate the b's
      b11 = t4*((s1+beta/w)*t2 / (w*t1) + (s2+1./w**2)*t3) - s2
      b12 = -t4*(s1*t2/(w*t1)+s2*t3) - 1./w**2 + s2
      b21 = (s1+beta/w) * (t3-beta*t2/t1)
      b21 = t4*(b21 - (s2+1./w**2)*(w*t1*t2+beta*w*t3)) + 1./(w**2*dt)
      b22 = s1*(t3-beta*t2/t1)
      b22 = -t4*(b22 - s2*(w*t1*t2+beta*w*t3)) - 1./(w**2*dt)
      return
      end

c ----------------------------------------------------------------------

      subroutine brs (x,w,beta,npts,rsp)
      real x(1), rsp(1), beta
      real w
      real*8 d, v, a, z, ap1, dp1, vp1, t1, t2
      real*8 a11, a12, a21, a22, b11, b12, b21, b22
      common /coef/ a11,a12,a21,a22,b11,b12,b21,b22
	  
c
c     initialize
      t1 = 2.*beta*w
      t2 = w**2
      d = 0.
      v = 0.
      a = 0.
c
c     calculate the response
      do 10 i=1,npts
        ap1 = dble( x(i) )
        dp1 = a11*d + a12*v + b11*a + b12*ap1
        vp1 = a21*d + a22*v + b21*a + b22*ap1
        z = -(t1*vp1 + t2*dp1)
        rsp(i) = sngl( z )
        a = ap1
        v = vp1
        d = dp1
  10  continue
      return
      end 

c ----------------------------------------------------------------------

      subroutine CheckAdjust ( delta_a, nPts, dt, tPeak0, nQ, 
     +                         deltaR0, w, damp, P )
      include 'rspMatch.h'
      real*8 delta_a(1), deltaR0(1)
      real dt, tPeak0(1), w(1), damp(1), acc(MAXPTS), rsp(MAXPTS)
      integer nPts, nQ, P(1)

c     Copy array
      do i=1,nPts
         acc(i) = sngl( delta_a(i) )
      enddo
	 
c     LOOP FOR EACH FREQUENCY
      do i=1,nQ
	  
c       Compute coeff 
        call coeff ( w(i), damp(i), dt )

c       CALCULATE THE RESPONSE
        call brs ( acc, w(i), damp(i), npts, rsp )
		
c       Find response at time tPeak0
        i1 = int( tPeak0(i) / dt )
        rsp1 = rsp(i1)

c        write (15,'( 4f10.5,i5)') w(i)/6.28, rsp1, deltaR0(i), 
C    1                 tPeak0(i), P(i)

      enddo
	  
      return
      end
	 
      
      
