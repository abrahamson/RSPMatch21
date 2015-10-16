      subroutine WrOutTH ( acc, nPts, dt, spectFile, thFileout,
     1      n_b, nAdded )
     
      implicit none
c     Passed Variables
      real dt
      real*8 acc(1)
      integer nPts, n_b, nAdded
      character*80 spectFile, thFileout
      
c     Local Variables
      integer i
      
      write(*,*) 'N added in previous run(s): ', nAdded
	write(*,*) 'N_before ',n_b
	write(*,*) 'nPts', nPts
      
c     Open output file
      write (*,'( 2x,''Enter accelerogram output file'')')
      read (5,'( a60)') thFileout
      open ( 20,file=thFileout, status='new')
      
c     Write header
      write (20,'( 2x,''Time history matched to spectrum:'',a60)')
     1       spectFile
      write (20,'( i5,f12.4,i8)') nPts, dt, n_b+nAdded
           
c     Write data
      write (20,'( 5e14.5)') (acc(i),i=1,nPts)
      close (20)
      
      return
      end subroutine WrOutTH
      
c ----------------------------------------------------------------------

      subroutine WrOutSpect ( R, Q, nQ, w, beta, iter, misfit, maxiter, 
     1      nDamp, damp, version, tol, gamma, iModel, a1, a2, f1, f2, 
     2      scaleFlag, dtFlag, evmin, groupSize, maxFreq, thFilein, 
     3      thFileout, targetFile, pga, R0, Q_noRand, iSeed, ranFactor,
     4      tPeak, file1, scaleperiod,freqBand,npole1, iModPGA, 
     5      FreqMatch,scale1 )
     
      implicit none
c     Passed Variables
      real*8 R(1), R0(1) 
      real Q(1), w(1), damp(1), misfit, beta(1), tPeak(1),
     1     tol, gamma, a1, a2, f1, f2, evmin, maxFreq, pga, Q_noRand(1),
     1      ranFactor, scaleperiod, freqBand(1), freqmatch(2), scale1
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
      end subroutine WrOutSpect

c ----------------------------------------------------------------------

      subroutine SaveSoln ( acc, nPts, accSave, bestMisfit, misfit )
      
      implicit none
      real bestMisfit, misfit
      real*8 acc(1), accSave(1)
      integer nPts, i
      
      write (*,'( 2x,''Saving Solution'')')

c     Copy acc array and reset bestMisfit
      do i=1,nPts
        accSave(i) = acc(i)
      enddo
      bestMisfit = misfit
      
      return
      end subroutine SaveSoln

c ----------------------------------------------------------------------
