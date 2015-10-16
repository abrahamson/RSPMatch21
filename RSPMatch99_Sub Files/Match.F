      subroutine DoMatch ( nGroup, Q1, w1, beta1, acc, nPts, dt, misfit,
     1           R, Q, nQ1, evMin, gamma, lastMisfit, scaleFlag, w, 
     2           beta, nQ, iFreq, flag1, aveMisfit, minTime1, minTime, 
     3           maxFreq, minFreq, iter, iDamp, scalePeriod, freqMatch, 
     4           flag2, maxTime1, maxTime, pgaacc, n_b, t_before )

      implicit none
      include 'rspMatch.h'
      
c     Passed Variables
      integer nGroup, nPts, nQ1(1), scaleFlag, iFreq, nQ, flag1, 
     1      iter, iDamp, flag2, n_b
      real Q1(MAXGROUP,1), w1(MAXGROUP,1), beta1(MAXGROUP,1), w(1), 
     1      beta(1), dt, Q(1), evMin, gamma, lastMisfit, misfit,
     1      minTime1(MAXGROUP,1), minTime(1), minFreq, maxFreq, 
     1      scalePeriod, freqMatch(2), maxTime1(MAXGROUP,1), maxTime(1),
     1      pgaacc(1), t_before     
      real*8 acc(1), R(1)
      integer iModel
      common / Model / iModel

c     Local Variables      
      real Q2(MAX_MATCH), w2(MAX_MATCH), beta2(MAX_MATCH), 
     1      minTime2(MAX_MATCH), maxTime2(MAX_MATCH),tPeak(MAX_MATCH)
     1      , aveMisfit, freq1, freq2
      real*8 deltaR(MAX_MATCH)
      integer iGroup, P(MAX_MATCH), nPts_before, pgaIndex
      
c     Write Header to screen
      write (*,'( /,2x,''Iteration SubGroup  AveMisfit  MaxMisfit  Freq(
     1Hz)  Damping  Pga(g)'')')

c     Loop over each subgroup
      do iGroup=1,nGroup
          call LoadSubgroup ( Q1, w1, beta1, Q2, w2, beta2, iGroup, nQ1,
     1                   minTime1, minTime2, maxTime1, maxTime2 )
     
          call CalcRsp ( acc, nPts, dt, w2, beta2, nQ1(iGroup), tPeak,
     1          R, P, minTime2, maxTime2 ) 
     
         if (iModel == 7) then
c		    Zero-pad the acceleration time series if needed
		    call PadAcc ( acc, nPts, dt, w2, nQ1(iGroup), tPeak, 
	1		    t_before, nPts_before )
		    n_b = n_b + nPts_before
		    call CalcRsp ( acc, nPts, dt, w2, beta2, nQ1(iGroup), tPeak,
     1          R, P, minTime2, maxTime2 )
		 endif	
		 	 
c         Calculate misfit for subgroup       
          call CalcMisfit ( R, Q2, nQ1(iGroup), P, deltaR, misfit, 
     1                      iFreq, aveMisfit, freqMatch, w2 )
          write (*,'( 2x,i6,''.'',i1,i7,3x,f10.4,1x,f10.4,2x,f8.3,1x,
     1           3f8.3)') iter, iDamp, iGroup, aveMisfit, 
     2      (R(iFreq)-Q2(iFreq))/Q2(iFreq), w2(iFreq)/6.28, beta2(iFreq)

c         Modify time history 
          call CalcNewTH ( acc, nPts, dt, Q2, w2, beta2, nQ1(iGroup), 
     1         gamma, R, P, tPeak, deltaR, evMin, freqMatch, pgaacc, 
     2         iGroup )  
      enddo
    
      if (iModel .ne. 7) then
c       Apply Baseline Correction
        if ( flag2 .eq. 1 ) then
           call Baseline ( acc, npts, dt, flag1 )
        endif      
c       Band-pass Filter to restrict energy to specifed freq band
        freq1 = maxFreq/6.28
        freq2 = minFreq * 0.6
        call BandPass ( acc, npts, dt, freq1, freq2, 5 )
      endif

c     Calculate response spectrum of new time history
      call CalcRsp ( acc, nPts, dt, w, beta, nQ, tPeak, R, P, minTime, 
     1     maxTime )

c     Scale to PGA
      if ( scaleFlag .eq. 1 ) then
        call FindPgaIndex ( nQ, w, pgaIndex, scalePeriod )
        call ScaleToPga ( pgaIndex, acc, npts, R, nQ, Q )
      endif

c     Calulate misfit from target (deltaR and maximum misfit)
      call CalcMisfit ( R, Q, nQ, P, deltaR, misfit, iFreq, aveMisfit, 
     1     freqMatch, w )

      return
      end subroutine DoMatch

c --------------------------------------------------------------

      subroutine CalcNewTH ( acc, nPts, dt, Q, w, damp, nQ, gamma,
     1           R, P, tPeak, deltaR, evMin1, freqMatch,pgaacc, iGroup )
     
      implicit none
      include 'rspMatch.h'
      
c     Passed Variables
      real Q(1), dt, damp(1), w(1), gamma, 
     1     tPeak(1), evMin1, freqMatch(2), pgaacc(1)
      real*8 deltaR(1), R(1), acc(1)
      integer P(1), nPts, nQ, iGroup, iModel
      common / Model / iModel
      
c     Local Variables
      integer i
      real alpha(MAX_MATCH), freqMin, freqMax
      real*8 b(MAX_MATCH), c(MAX_MATCH,MAX_MATCH), evMin, D(MAX_MATCH), 
     1      V(MAX_MATCH, MAX_MATCH), delta_a(MAXPTS), 
     2      delta_a1(MAXPTS)
     
      evMin = dble( evMin1 )

c     Load the Alpha matrix
      call Calc_Alpha ( w, alpha, nQ )

c     Calulate C matrix
      call Calc_C_Matrix ( c, damp, w, tPeak, nQ, alpha )

c     Singular value decomposition of C matrix        
      call SvdCmp ( c, nQ, nQ, MAX_MATCH, MAX_MATCH, D, V )

c     Compute b vector
      call Svbksb ( c, D, V, nQ, nQ, MAX_MATCH, MAX_MATCH, deltaR, 
     1       b, evMin )
   
c     Compute delta a(t)
      call Calc_delta_a ( delta_a, nPts, dt, b, w, damp, 
     1      tPeak, alpha, nQ )
     
c      write(100,*) (delta_a(k),k=1,npts)
     
      do i=1,npts
        delta_a1(i) = delta_a(i)
      enddo
      
      if (iModel .ne. 7) then
        if ( iGroup .eq. 1 ) then
            do i=1,npts
                delta_a1(i) = delta_a1(i) + pgaacc(i)
            enddo
        endif
        freqMax = freqMatch(2)*1.1
        freqMin = freqMatch(1)/ 1.1
        call BandPass ( delta_a1, npts, dt, freqMax, freqMin, 4 )
      endif

      
      do i=1,npts
        delta_a(i) = delta_a1(i)
      enddo
    
c     Adjust time history 
      call AdjustTH ( acc, nPts, delta_a, gamma )
            
      return
      end subroutine CalcNewTH

c ---------------------------------------------------------------------- 

      subroutine Calc_delta_a ( delta_a, nPts, dt, b, w, damp, 
     1      tPeak, alpha, nQ )
     
      implicit none
      include 'rspMatch.h'
      
c     Passed Variables
      real dt, w(1), damp(1), tPeak(1), alpha(1)
      real*8 b(1), delta_a(1)
      integer nPts, nQ, iModel
      common / Model / iModel

c     External Functions
      real f_func, f_funcXX

c     Local Variables
      real tau, temp_a(MAXPTS), twoPi
      integer i, j, n1
      twoPi = 2.*3.141592653589790

c     Initialize delta_A
      do i=1,nPts
        delta_a(i) = 0.0
      enddo

c     Compute delta_a      
      do 20 j=1,nQ
        tau = 0.0 
        do 10 i=1,nPts
          if ( iModel .eq. -6 ) then
                     temp_a(i) = b(j) * f_funcXX( tau, tPeak(j), w(j), damp(j),alpha(j) )

          ELSE
            temp_a(i) = b(j) * f_func( tau, tPeak(j), w(j), damp(j),
     1          alpha(j) )
          ENDIF
          tau = tau + dt 
  10    continue
      
      if (iModel .ne. 7) then
c       Apply taper to begining of adjustment
        n1 = int( 0.1 / (dt*(w(j)/twoPi)) )
        call taper ( temp_a, n1 )
      endif
        
c       Add to delta_a
        do 15 i=1,nPts
          delta_a(i) = delta_a(i) + temp_a(i)
  15    continue
  20  continue
  
      return
      end subroutine Calc_delta_a 
    
c ----------------------------------------------------------------------

      subroutine AdjustTH ( acc, nPts, delta_a, gamma )
      
      implicit none
c     Passed Variables
      real gamma 
      real*8 delta_a(1), acc(1)
      integer nPts 
      
c     Local Variables
      integer i 

c     Compute new acc  
      do 10 i=1,nPts
        acc(i) = acc(i) + delta_a(i) * gamma 
  10  continue
      
      return
      end subroutine AdjustTH

c ----------------------------------------------------------------------

      subroutine CalcMisfit ( R, Q, nQ, polarity, deltaR, misfit, 
     1      iFreq, aveMisfit, freqMatch, w )
     
      implicit none
c     Passed Variables
      integer nQ, polarity(1)
      real Q(1), misfit, aveMisfit, freqMatch(2), w(1)
      real*8 deltaR(1), R(1)
      
c     Local Variables
      integer i, iFreq
      real temp, sum, maxDeltaR, w1, w2
      
      maxDeltaR = 0.1
      w1 = freqMatch(1) * 2. * 3.141592653589790
      w2 = freqMatch(2) * 2. * 3.141592653589790
      
      misfit = 0
      sum = 0.0
      do i=1,nQ
        deltaR(i) = ( Q(i) - R(i) ) * polarity(i)
        temp = abs( deltaR(i) / Q(i) )
        sum = sum + temp 
        if ( temp .gt. misfit .and. w(i) .ge. w1 .and. w(i) .le. w2 ) 
     1      then
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
      end subroutine CalcMisfit
      
c ----------------------------------------------------------------------
