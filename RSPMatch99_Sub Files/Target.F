      subroutine RandomizeTarget ( iSeed, Q, nQ, ranFactor, target2, 
     1          nDamp, nFreq )
     
      include 'rspMatch.h'
      
c     Passed Variables      
      real ranFactor, Q(1), target2(MAX_FREQ,MAX_DAMP) 
      integer iSeed, nQ, nFreq, nDamp
  
c     Local Variables      
      integer i, i1, j
      real x, y
      
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
      end subroutine RandomizeTarget

c ----------------------------------------------------------------------

      subroutine SubdivideTarget ( nQ, Q, w, beta, nQ1, Q1, w1, beta1, 
     1      nGroup, groupSize, minTime, minTime1, maxFreq, maxTime, 
     2      maxTime1 )
     
      implicit none
      include 'rspMatch.h'
c     Passed Variables
      integer nQ, nQ1(1), nGroup, groupSize
      real Q(1), w(1), beta(1), Q1(MAXGROUP,1), w1(MAXGROUP,1),
     1     beta1(MAXGROUP,1), minTime(1), minTime1(MAXGROUP,1),
     2     maxFreq, maxTime(1), maxTime1(MAXGROUP,1)

c     Local Variables
      integer i, j, k, i1

      nGroup = nint( float(nQ) / groupSize + 0.5 )
      if ( nGroup .eq. 0 ) then     ! small number of spectral spoints
		nGroup = 1
	end if
      
c     Initialize nQ1
      do i=1,nGroup
        nQ1(i) = 0
      enddo
      
      i1 = 0
      do i = 0, nQ-1
        if ( w(i+1) .lt. maxFreq .or. w(i+1) .gt. 6280. ) then
c         j cycles between subgroups and k keeps track of the number of values in each subgroup          
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
      end subroutine SubdivideTarget
      
c ---------------------------------------------------------------------- 

      subroutine LoadQ_array ( nFreq, nDamp, nQ, Q, beta, target2, damp,
     1          iDamp, w )
     
      include 'rspMatch.h'
 
c     Passed Variables 
      integer nFreq, nDamp, nQ, iDamp
      real Q(1), beta(1), target2(MAX_FREQ,1), damp(1) , w(MAX_MATCH)
      
c     Local Variables        
      integer j, i
      
      j = iDamp
      do i=1,nFreq
        Q(i) = target2(i,j)
        beta(i) = damp(j)
        j = j + 1
        if ( j .gt. nDamp )  j = 1
      enddo
      nQ = nFreq

      return
      end subroutine LoadQ_array

c ----------------------------------------

      subroutine LoadSubgroup (Q1, w1, beta1, Q2, w2, beta2, iGroup,nQ1,
     1      minTime1, minTime2, maxTime1, maxTime2 )
      implicit none
      include 'rspMatch.h'
c     Passed Variables
      integer nQ1(1), iGroup
      real Q2(1), w2(1), beta2(1), Q1(MAXGROUP,1), w1(MAXGROUP,1),
     1     beta1(MAXGROUP,1), minTime1(MAXGROUP,1), minTime2(1),
     2     maxTime1(MAXGROUP,1), maxTime2(1)

c     Local Variables
      integer i
      
      do i=1,nQ1(iGroup)
        w2(i)    = w1(iGroup,i)
        beta2(i) = beta1(iGroup,i)
        Q2(i)    = Q1(iGroup,i)
        minTime2(i) = minTime1(iGroup,i)
        maxTime2(i) = maxTime1(iGroup,i)
      enddo
      
      return
      end subroutine LoadSubgroup
c ----------------------------------------------------------------------