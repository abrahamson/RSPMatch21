      real function Integ1 ( tau, s, s1, deltaT, beta_i, beta_j, 
     1      wi, wj, ti, tj, alpha )
c     This function computes the integral for c(i,j) for the
c     modified form of Lihanand and Tseng.
c     Includes alpha term only (no symmetry) (iModel=1)

c     Passed Variables
      real tau, s, s1, deltaT, beta_i, beta_j, wi, wj, ti, tj,
     1     alpha, theta_i, theta_j, B_Prime, C_Prime,
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
      end function Integ1
c ----------------------------------------------------------------------

      real function Integ3 ( tau, s, deltaT, beta_i, beta_j, 
     1      wi, wj, ti, tj, alpha_j )

c     Passed Variables
      real tau, s, deltaT, beta_i, beta_j, wi, wj, ti, tj,
     1     alpha_j, theta_i, theta_j, B_Prime, C_Prime,
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
      end function Integ3
      
c ----------------------------------------------------------------------
      double precision function Integ4 ( tau, deltaT, beta_i, beta_j, 
     1      wi, wj, ti, tj, alpha )
	! New Function - Linda, 02/03/09 

	implicit none

c     Passed Variables
      real tau, deltaT, beta_i, beta_j, wi, wj, ti, tj,
     1     alpha 

c     Local Variables
      real wj1, wi1, tc, t1
	double complex Exp1_a, Exp1_b, Exp2_a, Exp2_b, Arg1, Arg2, Arg3, Arg4, 
	1	CArg1, term1, term2, CArg2, Exp3_a, Exp3_b, term3, term4, cInteg, 
     2	F, A_1, expo_1, B_1, A_2, expo_2, B_2, A_3, expo_3, B_3, A_4,
     3	expo_4, B_4, Exp4_a, Exp4_b, C1, C2, C3        
	double precision A, B, term1_r, PI

c     Set constants
	  PI = 3.14159265358979323846264338328
	  wj1 = wj*SQRT(1-beta_j**2)
	  wi1 = wi*SQRT(1-beta_i**2)
	  deltaT = ATAN( SQRT(1-beta_j**2)/beta_j ) / wj1
	  tc = deltaT + ti - tj

	  C1 = dcmplx( -sin(tc*wi1), cos(tc*wi1) )
	  C2 = dcmplx( cos(tc*wi1), -sin(tc*wi1) )
	  C3 = dcmplx( sin(tc*wi1), cos(tc*wi1) )

	  A = (alpha**2 * beta_i**2 * wi**2 - 4*beta_i*wi*tc -2* alpha**2
	1	*(wi1**2+wj1**2))/4
	  B = -2*beta_i*wi*alpha**2 * (wi1 + wj1)/4
	  F = dcmplx( A, B)

	  t1 = tau - tj + deltaT

c	  Set up integral terms 
	  term1_r = SQRT(PI) * alpha * wi  / (8*SQRT(1-beta_i**2))
	  term1 = dcmplx( 0, term1_r)

	  Arg1 = dcmplx ( alpha*(wi1-wj1)/2, (-2*t1+alpha**2*beta_i*wi)/
     1      (2*alpha) )
	  Arg2 = dcmplx ( alpha*(wi1+wj1)/2, (-2*t1+alpha**2*beta_i*wi)/
     1      (2*alpha) )
	  Arg3 = dcmplx ( -t1/alpha+alpha*beta_i*wi/2, alpha*(wi1+wj1)/2 )
	  Arg4 = dcmplx ( alpha*(wi1-wj1)/2, t1/alpha-alpha*beta_i*wi/2 )

	  call erfi(Arg1, A_1, expo_1, B_1)
	  call erfi(Arg2, A_2, expo_2, B_2)
	  call erf(Arg3, A_3, expo_3, B_3)
	  call erfi(Arg4, A_4, expo_4, B_4)

	  Exp1_a = (alpha**2)/4 * dcmplx( (wi1+wj1)**2, 4*beta_i*wi*wj1 ) + F
	  Exp1_b = (alpha**2)/4 * dcmplx( (wi1+wj1)**2, 4*beta_i*wi*wj1 ) + F 
	1	    + expo_1

 	  Exp2_a = (alpha**2)/4 * (wi1-wj1)**2 + F 
	  Exp2_b = (alpha**2)/4 * (wi1-wj1)**2 + F + expo_2

	  Exp3_a = (alpha**2)/4 * dcmplx( (wi1-wj1)**2, 4*beta_i*wi*(wi1+wj1)) 
     1      + F 
	  Exp3_b = (alpha**2)/4 * dcmplx( (wi1-wj1)**2, 4*beta_i*wi*(wi1+wj1)) 
	2	    + F + expo_3

	  Exp4_a = (alpha**2)/4 * dcmplx( (wi1+wj1)**2, 4*beta_i*wi*wi1 ) + F 
	  Exp4_b = (alpha**2)/4 * dcmplx( (wi1+wj1)**2, 4*beta_i*wi*wi1 ) + F 
	3	    + expo_4
  
	  CArg1 = CDEXP(Exp1_a) * A_1 + CDEXP(Exp1_b) * B_1  		
	  CArg2 = CDEXP(Exp2_a) * A_2 + CDEXP(Exp2_b) * B_2
	  term2 = ( CArg1 + CArg2 ) * C1  
		
	  term3 = CDEXP(Exp3_a) * A_3 * C2 + CDEXP(Exp3_b) * B_3 * C2

	  term4 = CDEXP(Exp4_a) * A_4 * C3 + CDEXP(Exp4_b) * B_4 * C3
		 
c	  Evaluate integral
	  cInteg =  term1 * (term2 + term3 + term4)
        Integ4 = dreal( cInteg )
        
c        write(102,*) wi, ti, wj, tj, deltaT, alpha, tau, Integ4

      return
      end function Integ4

c ----------------------------------------------------------------------
