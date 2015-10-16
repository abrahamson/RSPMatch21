      program RSPM09
c
c     Version: 2009
c     Date: May 13, 2009
c
c     Spectral matching program developed by N. Abrahamson.Copyright July 7,1993.
c     Based on time domain method of Tseng and Lilanand (1988).
c     Modified in 1999 to preserve non-stationarity at long periods by using 
c     different functional forms for the adjustment time history.
c
c     Version 2009: by Linda Al Atik - UC Berkeley
c     email: l_atik@berkeley.edu
c     Modified by adding a new adjustment function that ensures stability
c     and efficinecy of the numerical solution while preventing drift in 
c     the resulting velocity and displacement time series.
c     Refer to user manual for more informtaion.
c
c     If this were run on the Macintosh computer, then the 
c     following subprograms could also be included for plotting:
c          updateplot
c          tek
c
c     Also, there may be tabs in this code that may confuse your compiler. 
c     Edit as appropriate.
      
      implicit none
      include 'rspMatch.h'
      
      real*8 acc(MAXPTS), deltaR(MAX_MATCH), accSave(MAXPTS), 
     1      R0(MAX_MATCH), R(MAX_MATCH)   
      real Q(MAX_MATCH), tPeak(MAX_MATCH), w(MAX_MATCH), beta(MAX_MATCH)
     1      , w_unmatched(MAX_MATCH), Q_unmatched(MAX_MATCH), 
     2      Q1(MAXGROUP,MAX_MATCH), w1(MAXGROUP,MAX_MATCH), 
     4      beta1(MAXGROUP,MAX_MATCH), Q_all(MAX_MATCH), 
     5      beta_all(MAX_MATCH),Q_all_noRand(MAX_MATCH), 
     6      Q_all_noRand_unmatched(MAX_MATCH), damp(MAX_DAMP), 
     7      target2(MAX_FREQ,MAX_DAMP), minTime_all(MAX_MATCH),
     8      w_all(MAX_MATCH), pgaacc(MAXPTS), minTime(MAX_MATCH), 
     9      minTime1(MAXGROUP,MAX_MATCH), 
     1      freqBand(2), freqMatch(2), maxTime(MAX_MATCH), 
     2      maxTime_all(MAX_MATCH), Time_array(MAX_MATCH),
     3      maxTime1(MAXGROUP,MAX_MATCH), 
     4      zero(MAX_MATCH) 
      real a1, a2, f1, f2, bestMisfit, tol, gamma, dt, misfit, 
     1      aveMisfit, minFreq, evMin, lastMisfit,   
     3      maxfreq, pga, scalePeriod, ranFactor, scale1, t_before  
      integer P(MAX_MATCH), nQ1(MAXGROUP) 
      integer iter, maxIter, nQ, nPts, i, iFreq, scaleFlag, 
     1      dtFlag, pgaIndex, iDamp, iModPga, nGroup, groupSize,
     2      flag1, nPole1, flag2, iModel, nQ_all, nFreq, nDamp, iSeed, 
     3      nAdded, n_b, irun, nb_runs
      character*80 targetFile, thFilein, thFileout, version, file1,
     1      input_file
      common / Model / iModel
      common /alphaModel/ a1, a2, f1, f2
      data version / 'Version:09, Date:5_13_09'/ 

      pause 'test'
      write (*,'(2x,''********************************'')')
      write (*,'(2x,''RSPM - Spectral Matching Program'')')
      write (*,'( 2x,''Program Version: '',a25)') version
      write (*,'(2x, ''Based on Code by, N. Abrahamson (1993)'')')
      write (*,'(2x,''********************************'')')
      
      write(*,'( 2x,''Enter input file'')')
      read (*,'( a80)') input_file
      open (2, file=input_file, status='old')
      read (2,*) nb_runs

c     Loop over the run files      
      do irun = 1, nb_runs
        write(*,*)
        read (2,'( a80)') file1
        write (*,1) 'Run File: ', file1
1       format('  ',A10, A80)
        open (5,file=file1,status='old')
      
c       Initialize      
        flag1 = 0
        do i=1,MAX_MATCH
            zero(i) = 0.
        enddo
        bestMisfit = 1.0e30
        n_b = 0   !Initialize number of points added to the acceleration 
                  !record
        t_before = 0. 
           
c       Read Run parameters, target spectrum, and reference accelerogram
        call RdRunParam ( tol, maxIter, gamma, scaleFlag, dtFlag, evMin,
     1      groupSize, maxFreq, freqBand, nPole1, scalePeriod, iModPGA, 
     2      iSeed, ranFactor, freqMatch, flag2, scale1 )
        call RdTarget (w, damp, targetFile, minFreq, minTime, minTime_all,
     1      beta_all, Q_all, nQ_all,nFreq, nDamp, target2, w_all, 
     2      maxTime, maxTime_all, freqmatch)
        call RdAcc ( acc, nPts, dt, dtFlag, thFilein, nAdded )

        do i=1,MAX_MATCH
            Time_array(i) = (npts-1)*dt     
        enddo

c       Randomize Target?
        do i=1,nQ_all
            Q_all_norand(i) = Q_all(i)
        enddo
        if ( iSeed .ne. 0 ) then
            call RandomizeTarget ( iSeed, Q_all, nQ_all, ranFactor, target2,
     1          nDamp, nFreq )
        endif

c       Initial Band-pass filter
        freqBand(2) = freqBand(2) * 6.28
        call BandPass ( acc, npts, dt, freqBand(2), freqBand(1), nPole1 )
      
c       Calculate response spectrum of reference time history
        call CalcRsp ( acc, nPts, dt, w_all, beta_all, nQ_all, tPeak, R0,
     2      P, zero, time_array )
      
c       Calulate initial misfit from target 
        call CalcMisfit ( R0, Q_all, nQ_all, P, deltaR, misfit, iFreq, 
     1      aveMisfit, freqMatch, w_all )
        write (*,'( /,2x,''Initial Solution    AveMisfit  MaxMisfit  
     1      Freq(Hz)  Damping'')')
        write (*,'( 8x,12x,f10.4,1x,f10.4,2x,f8.3,1x,2f8.3)') aveMisfit, 
     1      (R0(iFreq)-Q_all(iFreq))/Q_all(iFreq), w_all(iFreq)/6.28, 
     2      beta_all(iFreq)

c       Make initial adjustment to time history by scaling to PGA (or other period)
        call FindPgaIndex ( nQ_all, w_all, pgaIndex, scalePeriod )
        call FindPga ( acc, npts, pga )
        write (*,'( 2x,''scaleIndex ='',i5,2x,''pga ='',f10.4)') pgaIndex,
     1      pga
     
        if ( scaleFlag .ge. 1 ) then
            call ScaleToPga ( pgaIndex, acc, npts, R0, nQ_all, Q_all )
            call FindPga ( acc, npts, pga )
            call CalcMisfit ( R0, Q_all, nQ_all, P, deltaR, misfit, iFreq, 
     1          aveMisfit, freqMatch, w_all )
            write (*,'( /,2x,''Scaled Solution     AveMisfit  MaxMisfit  Fre
     1      q(Hz)  Damping  Pga(g)'')')
            write (*,'( 8x,12x,f10.4,1x,f10.4,2x,f8.3,1x,3f8.3)')aveMisfit, 
     1          (R0(iFreq)-Q_all(iFreq))/Q_all(iFreq), w_all(iFreq)/6.28, 
     2          beta_all(iFreq), pga
        endif
      
c       Apply scale factor to the acceleration time series 
        do i=1,npts
            acc(i) = acc(i) * scale1
        enddo
        
c       Initialize Graphics
c       call plots
c       call factor (1.5)
c       call UpdatePlot ( w_all, Q_all, beta_all, R0, nFreq, nDamp, 0,0,
c    1      aveMisfit, misfit, iFreq ) 
      
c       Iterate until misfit is less than tolerance
        iter = 1
        lastMisfit = 1.0e10
c	  Save initial solution just in case it's good enough and no iterations required
        call SaveSoln ( acc, nPts, accSave, bestMisfit, misfit )     
        doWhile ( misfit .gt. tol .and. iter .le. maxIter )        
            do iDamp=1,nDamp
                call LoadQ_array ( nFreq, nDamp, nQ, Q, beta, target2, damp, 
     1              iDamp, w )
c               Subdivide Target into subgroups
                call SubdivideTarget ( nQ, Q, w, beta, nQ1, Q1, w1, beta1, 
     1              nGroup, groupSize, minTime, minTime1, maxfreq, 
     2              maxTime, maxTime1 )
                call DoMatch ( nGroup, Q1, w1, beta1, acc, nPts, dt, misfit, 
     1              R, Q_all, nQ1, evMin, gamma,lastMisfit,scaleFlag,
     2              w_all, beta_all, nQ_all, iFreq, flag1, aveMisfit,
     2              minTime1, minTime_all, maxFreq, minFreq, iter, iDamp, 
     3              scalePeriod, freqMatch, flag2, maxTime1, maxTime_all, 
     4              pgaacc, n_b, t_before)
                call FindPga ( acc, npts, pga )
                write (*,'( 2x,i6,'' full set'',3x,f10.4,1x,f10.4,2x,f8.3,1x,
     1        	    4f8.3)') iter, aveMisfit, 
     2              (R(iFreq)-Q_all(iFreq))/Q_all(iFreq), w_all(iFreq)/6.28,
     3              beta_all(iFreq), pga
c               If this solution is best so far, then save time history
c               This allows you to keep the best fit in case the solution 
c               starts to diverge.
                if ( misfit .lt. bestMisfit ) then
                    call SaveSoln ( acc, nPts, accSave, bestMisfit, misfit )
                endif          
                lastMisfit = misfit
            enddo

c           Add pga modification (except if finished)
            if (iModel .ne. 7) then
               if ( iter .ne. MaxIter .and. misfit .gt. tol ) then
                   call AddPga ( acc, npts, nQ, Q, maxFreq, dt, pgaAcc,
     1                 iModPga )
               endif 
            endif      
            iter = iter + 1
        enddo
   
        do i=1,MAX_MATCH
            maxTime_all(i) = (npts-1)*dt        ! was npts*dt
        enddo

c       Recompute spectrum for best solution and write output
        call CalcRsp ( accSave, nPts, dt, w_all, beta_all, nQ_all, tPeak, 
     2      R, P, zero, maxTime_all )
        call CalcMisfit ( R, Q_all, nQ_all, P, deltaR, misfit, iFreq, 
     1      aveMisfit, freqMatch, w_all )
        call FindPga ( accSave, npts, pga )
        call WrOutTH ( accSave, nPts, dt, targetFile, thFileout,
     1      n_b, nAdded)
        call WrOutSpect ( R, Q_all, nQ_all, w_all, beta_all, iter, misfit,
     1      maxiter, nDamp, damp, version,tol, gamma, iModel, a1, a2, 
     2      f1, f2, scaleFlag, dtFlag, evmin, groupSize, maxFreq, 
     3      thFilein, thFileout, targetFile, pga, R0, Q_all_noRand, 
     4      iSeed, ranFactor, tPeak, file1, scaleperiod,freqBand,npole1,
     6      iModPGA,FreqMatch,scale1 )

c       Check unmatched frequencies
        do i=1,nQ_all-1
            w_unmatched(i) = (w_all(i)+w_all(i+1))/2.
            Q_unmatched(i) = sqrt( Q_all(i)*Q_all(i+1) )
            Q_all_noRand_unmatched(i)=sqrt(Q_all_noRand(i)*
     1          Q_all_noRand(i+1))
        enddo

        call CalcRsp ( accSave, nPts, dt, w_unmatched, beta_all, nQ_all, 
     1      tPeak, R, P, zero, maxTime_all )
        call CalcMisfit ( R, Q_all, nQ_all, P, deltaR, misfit, iFreq, 
     1      aveMisfit, freqMatch, w_all )
        call WrOutSpect ( R, Q_unmatched, nQ_all, w_unmatched, beta_all,
     1      iter, misfit, maxiter, nDamp, damp, version,tol, gamma, 
     2      iModel, a1, a2, f1, f2, scaleFlag, dtFlag, evmin, groupSize,
     3      maxFreq, thFilein, thFileout, targetFile, pga, R0, 
     4      Q_all_noRand, iSeed, ranFactor, tPeak, file1, scaleperiod,
     5      freqBand,npole1, iModPGA,FreqMatch,scale1 )
      enddo
      stop 
      end program RSPM09

c --------------------------------------------------------------
      