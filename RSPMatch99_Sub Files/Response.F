      subroutine CalcRsp ( acc, npts, dt, w, damping, nQ, time, SA, 
     1      polarity, minTime, maxTime)

c     This subroutine calculates the response spectrum 
      implicit none
      include 'rspMatch.h'      
      real*8 acc(1), SA(1)
      real w(1), damping(1), dt, time(1), minTime(1), 
     1           maxTime(1)
      integer npts, nQ, i, polarity(1), timeIndex, iModel
      common / Model / iModel

c	Local variables
	real*8 acc_mod(MAXPts_RSP)
	real dt_mod
	integer npts_mod

c	Initialize
	dt_mod = dt
	nPts_mod = nPts
	do i=1, nPts_mod
		acc_mod(i) = acc(i)
	enddo
	
c     LOOP FOR EACH FREQUENCY
      do i=1,nQ 
        ! Compute the pseudo-acceleration response spectrum using Newmark's method
        if ( iModel == 7 ) then
            call set_dt ( w(i), acc_mod, npts_mod, dt_mod ) 
            call CalcOneRsp_Newmark( acc_mod, npts_mod, dt_mod, w(i), 
     1          damping(i), sa(i),timeIndex, polarity(i), minTime(i),
     2          maxTime(i) )
            time(i) = (timeIndex-1) * dt_mod 
        else                
c         if ( w(i) .gt. 628. ) then
c           call CalcPGA ( acc, npts, sa(i), timeIndex, polarity(i) )
c          else           
            call CalcOneRsp( acc, npts, dt, w(i), damping(i), sa(i),
     1          timeIndex, polarity(i), minTime(i),maxTime(i) )    
            time(i) = (timeIndex-1) * dt
c          endif
        endif
      enddo
	  
      return
      end subroutine CalcRsp 

c ----------------------------------------------------------------------

	subroutine set_dt (w, acc_mod, npts_mod, dt_mod) 

c	This subroutine finds the time-step required to ensure accuracy of the  
c	numerical solution for the response spectrum. If input dt is greater than the  
c	required time step, this subroutine applies linear interpolation to get the 
c	required time step.

	implicit none
	include 'rspMatch.h'

c	Passed variables
	real*8 acc_mod(1)
	real dt_mod, w
	integer nPts_mod
c	Local varaibles
	real F, freq, period, pi
	integer iFactor, i1, i, j

	pi = 3.1415926
	freq = w / (2*pi)
	period = 1 / freq
	F  = dt_mod/period

	if ( F > 0.01 ) then
		i1 = npts_mod
		iFactor = INT( dt_mod / (0.01 * period) )
		nPts_mod = 	iFactor*(nPts_mod-1) + 1
	    if ( npts_mod .gt. MAXPTS_RSP) then
			write (*,'( 2x,''ERROR: Increase MAXPTS_RSP to '',i10)') npts_mod
			stop 99
		endif
		dt_mod = dt_mod / iFactor
          do 20 i = npts_mod, 1, -ifactor
			acc_mod(i) = acc_mod(i1)
			i1 = i1 -1
  20		continue
          do 30 i = 1, nPts_mod, ifactor
			do 25 j = 1, ifactor-1
				acc_mod(i+j) = acc_mod(i) + ( acc_mod(i+ifactor) - acc_mod(i) ) *
     1                 float(j) / float(ifactor)
  25			continue
  30		continue

      endif

	return
	end subroutine set_dt
c-----------------------------------------------------------------------
	subroutine CalcOneRsp_Newmark( acc, npts, dt, w, damping, sa, timeIndex, 
     +      polarity, minTime, maxTime)

	implicit none
	include 'rspMatch.h'

c	Passed variables	
      real*8 acc(1), Sa
	real dt, w, damping, minTime, maxTime
	integer nPts, timeindex, polarity

c	Local variables
	real gamma, beta
	real*8 Rsp(MAXPts_RSP), a, v, d, k, a_term, b_term, delta_p, delta_acc,
	1	delta_a, delta_v, delta_d
	integer i, iTime, iTime2, j

c	Initial Calculations
	gamma = 0.5
	beta = 0.25
	k = w**2 + gamma/(beta*dt)*2*damping*w + 1/(beta * dt**2)
	a_term = 1/(beta*dt) + gamma/beta*2*damping*w
	b_term = 1/(2*beta) + dt * (gamma/(2*beta)-1) * 2*damping*w

c	Calculation of the response for each time step, i
	do i=1, nPts
		if (i == 1) then
			d = 0.
			v = 0.
			a = - acc(1) - 2 * damping * w * v - w**2 * d
		else
			delta_acc = acc(i) - acc(i-1)
			delta_p = - delta_acc + a_term * v + b_term * a
			delta_d	= delta_p / k
			delta_v	= gamma/(beta*dt) * delta_d - gamma/beta * v + 
	1			dt * ( 1 - gamma/(2*beta) ) * a
			delta_a = 1/(beta * dt**2) * delta_d - 1/(beta*dt) * v -
	1			1/(2*beta) * a

			d = d + delta_d
			v = v + delta_v
			a = a + delta_a
		endif
		Rsp(i) = w**2 * d   !Pseudo-spectral acceleration
	enddo

c	Finding the maximum of the response
      Sa = 0.0
	iTime = int( minTime/dt) + 1
	iTime2 = int( maxTime/dt) -1
	if ( iTime2 .gt. npts ) iTime2=npts
	if ( maxTime .eq. 0. ) iTime2=npts
      do j=iTime,iTime2
          if (abs(Rsp(j)) .gt. Sa) then
		     Sa = abs(Rsp(j))
			 timeIndex = j
			 if ( rsp(j) .ge. 0. ) then
			    polarity = 1
			 else
			    polarity = -1
			 endif
		  endif
      enddo

	return
      end subroutine CalcOneRsp_Newmark
c--------------------------------------------------------------------------------

       subroutine CalcOneRsp ( acc, npts, dt, w, damping, sa, timeIndex,
     1      polarity, minTime, maxTime)
     
      implicit none
      include 'rspMatch.h'
      real*8 acc(1), SA
      real rsp(MAXPTS), w, damping, dt, minTime, maxTime
      integer npts, i, j, polarity, timeIndex, iTime, iTime2
	  
c     Compute coeff 
      call coeff ( w, damping, dt )

c     CALCULATE THE RESPONSES
      call brs ( acc, w, damping, npts, rsp )

c     FIND THE MAXIMUM OF THE RESPONSE
      SA = 0.0
	  iTime = int( minTime/dt) + 1
	  iTime2 = int( maxTime/dt) -1
	  if ( iTime2 .gt. npts ) iTime2=npts
	  if ( maxTime .eq. 0. ) iTime2=npts
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
      end subroutine CalcOneRsp


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
      end subroutine coeff

c ----------------------------------------------------------------------

      subroutine brs (x,w,beta,npts,rsp)
      
      real rsp(1), beta
      real w
      real*8 d, v, a, z, ap1, dp1, vp1, t1, t2, x(1)
      real*8 a11, a12, a21, a22, b11, b12, b21, b22
      common /coef/ a11,a12,a21,a22,b11,b12,b21,b22

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
        z = -(t1*vp1 + t2*dp1)  !spectral acceleration
	!  z = w**2 * dp1         !pseudo-spectral acceleration
        rsp(i) = sngl( z )
        a = ap1
        v = vp1
        d = dp1
c        write (23,'( 5e12.6)') d
  10  continue
  
      return
      end subroutine brs 

c ----------------------------------------------------------------------  