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
     1          w(i), w(j), alpha(j), c(i,j) )
  10    continue
  20  continue
 
      return
      end subroutine Calc_C_Matrix
 
c ----------------------------------------------------------------------

      subroutine Calc_C ( ti, tj, beta_i, beta_j, wi, wj, alpha, C )
      
      implicit none
c     Passed Varaibles
      real ti, tj, beta_i, beta_j, wi, wj, alpha
      real*8 C
      
c     External Functions and common
      real Integ1, Integ2, Integ3, func_deltaT, h_func, f_func
      integer iModel
      double precision integ4
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
      if ( iModel .eq. 1 ) then

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
        
        elseif ( iModel .eq. 6 .or. iModel .eq. -6 ) then
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
          
       elseif ( iModel .eq. 7 ) then
c         Evaluate integral for model 4 (New Model) - Linda, 02/03/09
        C = Integ4( ti, deltaT, beta_i, beta_j,wi,wj,ti,tj,alpha ) -
     1      Integ4( 0., deltaT, beta_i, beta_j,wi,wj,ti,tj,alpha ) 
      endif   
           
      return
      end subroutine Calc_C
      
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
      end subroutine Calc_alpha      
      
c ----------------------------------------------------------------------

      real function alpha_func ( w )
      
      implicit none
      real w, alpha, freq, a1, a2, f1, f2
      integer iModel
      common /Model/ iModel
      common /alphaModel/ a1, a2, f1, f2
      
      freq = w / ( 6.28)

c     Check for Lilhanand and Tseng model
      if ( iModel .eq. 1 ) then
        alpha_func = 0.
        return
        
      elseif ( iModel .eq. 6 .or. imodel .eq. -6 ) then
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
        
	elseif ( iModel .eq. 7 ) then 
        alpha = 1.178 * ( freq**(-0.93) )
	  alpha_func = alpha
        return
	endif
        
      return
      end function alpha_func

c ----------------------------------------------------------------------

      real function func_deltaT ( beta, w )

      real beta, w 

c     Local Variables
      real tmp1      
      
      tmp1 = sqrt( 1. - beta**2 )
      func_deltaT = atan( tmp1/beta ) / ( w * tmp1 )
      return
      end function func_deltaT 
      
c ----------------------------------------------------------------------

      real function h_func( w, beta, t )
      
      implicit none
c     Passed Variables
      real w, beta, t
      integer iModel
      common /Model/ iModel
      
c     Local Variables
      real tmp1, theta
      
      if ( t .lt. 0. ) then
          h_func = 0.0
          return
      endif

      tmp1 = sqrt( 1. - beta**2 )
      theta = w * tmp1 * t
      
      if (iModel == 7) then
c       Pseudo-acceleration impulse response      
        h_func = - w / tmp1 * exp( -w*beta*t ) * sin ( theta )
      else
c       Acceleration impulse response
        h_func = - w / tmp1 * exp( -w*beta*t ) * ( (2.*beta**2 - 1. ) 
     1      * sin ( theta ) - 2.*beta * tmp1 * cos ( theta ) )
      endif
     
      return
      end function h_func

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
      real t, deltaT, wPrime
      
      deltaT = func_deltaT ( beta, w )
      wPrime = w * sqrt( 1. - beta**2 )
      
c     Check for model (iModel=1 is for Lilhanand and Tseng model)      
      if ( iModel .eq. 1 ) then
        t = tj - tau 
      elseif ( iModel .eq. 6 .or. iModel .eq. 7 .or. imodel .eq. -6 ) then
        t = tau - tj + deltaT
      endif

c     Use Lilhanand and Tseng f function except for pga (f>=1000Hz)
      if ( w .lt. 6280. ) then
        if ( iModel .le. 4 ) then
c         f_func = h_func( w, beta, t )*exp( -abs(tau-tj+deltaT)*alpha)
          f_func = h_func( w, beta, t ) * exp( -abs(tau-tj)*alpha )
        elseif ( iModel .eq. 6 .or. imodel .eq. -6 ) then
          f_func = cos( wPrime*t ) * exp( -abs(t)*alpha )
        elseif ( iModel .eq. 7 ) then
          f_func =  cos( wPrime*t ) * exp( - (t/alpha)**2 ) 
        endif
      else
        f_func = DiracDelta ( t )
      endif
      
      return
      end function f_func
c ----------------------------------------------------------------------

      real function DiracDelta ( x )
      
      real x
      
      if ( abs(x) .lt. 1.0e-4 ) then
        DiracDelta = 1.
      else
        DiracDelta = 0.
      endif
      return
      
      end function DiracDelta
      
c ----------------------------------------------------------------------

      real function f_funcXX ( tau, tj, w, beta, alpha )
      
      implicit none
c     Passed Variables
      real tau, tj, w, beta, alpha

c     External Functions and common
      real h_func, func_deltaT, DiracDelta
      integer iModel
      common /Model/ iModel

c     Local Variables
      real t, deltaT, wPrime
      
      deltaT = func_deltaT ( beta, w )
      wPrime = w * sqrt( 1. - beta**2 )
      
c     Check for model (iModel=1 is for Lilhanand and Tseng model)      
      if ( iModel .eq. 1 ) then
        t = tj - tau 
      elseif ( iModel .eq. 6 .or. iModel .eq. 7 .or. imodel .eq. -6 ) then
        t = tau - tj + deltaT
      endif

c     Use Lilhanand and Tseng f function except for pga (f>=1000Hz)
      if ( w .lt. 6280. ) then
        if ( iModel .le. 4 ) then
c         f_func = h_func( w, beta, t )*exp( -abs(tau-tj+deltaT)*alpha)
          f_funcXX = h_func( w, beta, t ) * exp( -abs(tau-tj)*alpha )
        elseif ( iModel .eq. 6  ) then
          f_funcXX = cos( wPrime*t ) * exp( -abs(t)*alpha )
        elseif ( iModel .eq. 7 ) then
          f_funcXX =  cos( wPrime*t ) * exp( - (t/alpha)**2 ) 
        elseif ( imodel .eq. -6 ) then
          f_funcXX = cos( wPrime*t ) * exp( -abs(t)*alpha )
          if (f_funcXX .lt. 0.3 .and. abs(t) .lt. 20. ) f_funcXX=0.3
         endif
      else
        f_funcXX = DiracDelta ( t )
      endif
      
      return
      end function f_funcXX

