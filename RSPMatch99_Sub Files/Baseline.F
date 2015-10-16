      subroutine baseline ( acc, nPts, dt, flag1 )
      
      implicit none
      include 'rspMatch.h'

c     Passed Variables
      real*8 acc(1)
      real dt
      integer nPts, flag1
      
c     Local Variables
      integer i, nParam
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
     1      nPts, 1, xhat, MAXPARAM )
     
c     Remove baseline
      t = 0
      do i=1,nPts
        acc(i) = acc(i) - ( xhat(1,1)*t + xhat(2,1)*t**2 )
        t = t + dt
      enddo
     
      return
      end subroutine baseline

C-----------------------------------------------------------------------

      subroutine LoadA ( A, npts, dt, nParam, MAXPTS )
      
      implicit none
 
c     Passed Variables
      integer npts, MAXPTS
      real dt
      real*8 A (MAXPTS,1)

c     Local Variables
      integer i, nParam
      real t

c     Load A Matrix
      t = 0
      do i=1,npts
        A(i,1) = t
        A(i,2) = t**2
        t = t + dt
      enddo     
      nParam = 2
      
      return
      end subroutine LoadA

c ----------------------------------------------------------------------