      subroutine ols ( arow, acol, A, npts, nterm, ATAIAT  )
      
      include 'rspMatch.h'
C	  implicit none
      integer arow, acol, i, j, nterm, npts
      real*8 A(arow,acol), AT(MAXPARAM,MAXPTS), ATA(MAXPARAM,MAXPARAM),
     1      ATAI(MAXPARAM,MAXPARAM), ATAIAT(MAXPARAM,MAXPTS), eps, 
     2      deter, work(MAXPTS)

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

      return
      end subroutine ols 
      
c ----------------------------------------------------------------------

      subroutine transpose ( A, arow, m, n, At, atrow )
      
c     Passed Variables            
      integer arow, atrow, m, n
      real*8 A(arow,1), At(atrow,1)
      
c     Local Variables
      integer i, j      

      do i=1,m
        do j=1,n
          At(j,i) = A(i,j)
        enddo
      enddo
      
      return
      end subroutine transpose

c ----------------------------------------------------------------------

      subroutine mult (a,arow,m1,n1,b,brow,m2,n2,c,crow)
c     This subroutine multiplies two matrices.
c       C = A B
c     where A is m1 x n1  and  B is m2 x n2

c     Passed Variabled      
      INTEGER AROW, BROW, CROW, m1, n1, m2, n2
      double precision a(arow,1), b(brow,1), c(crow,1)
      
c     Local Variabled      
      integer i, j, ii

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
      end subroutine mult

c ----------------------------------------------------------------------

      subroutine svdcmp ( A, m, n, mp, np, W, V )

      include 'rspMatch.h'      
c      parameter (NMAX=100)
C      implicit none
      integer m,n,mp,np, nmax
      real*8 A(mp,np), W(np), V(np,np), rv1(MAX_MATCH)
      integer i, j, k, l, jj, its, nm
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
      end subroutine svdcmp
	  
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
	  end subroutine Svbksb
      
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
      END SUBROUTINE SIMUL
       
c ----------------------------------------------------------------------

	subroutine erf(z, A, expo, B)

	implicit double complex (c,z)
	double complex z, A, B, expo
	real a0,pi
	integer k
	
	a0=cdabs(z)
	expo=-z*z
	pi=3.141592653589793d0
	z1=z
	if (real(z).lt.0.0) then
	   z1=-z
	endif
	if (a0.le.5.8d0) then    
	   cs=z1
	   cr=z1
	   do 10 k=1,120
	      cr=cr*z1*z1/(k+0.5d0)
	      cs=cs+cr
	      if (cdabs(cr/cs).lt.1.0d-15) go to 15
10         continue
15		 A = 0
           B =2.0d0*cs/sqrt(pi)
	else                              
	   cl=1.0d0/z1              
	   cr=cl
	   do 20 k=1,13
	      cr=-cr*(k-0.5d0)/(z1*z1)
	      cl=cl+cr
	      if (cdabs(cr/cl).lt.1.0d-15) go to 25
20         continue
25		 A = 1
		 if ( real(expo) > 307d0 ) then
			A = 0
		 else
			A = 1 
		 endif
           B = -cl/sqrt(pi)
	endif
	if (real(z).lt.0.0) then
	   A = -A
	   B = -B
	endif	
     
      return
	end subroutine erf
c-------------------------------------------------------------------------------------
	subroutine erfi(z, A, expo, B )

	implicit none
	double complex z, z1, A, expo, B, m, c, n
	double precision z_r, z_i

	z_r = dreal(z)
	z_i = dimag(z)
	z1 = dcmplx( -z_i, z_r )   ! z1 = i*z
	call erf(z1, m, c, n)	 

	expo = c
	A = -m * cmplx(0, 1)
	B = -n * cmplx(0, 1)

	return
	end subroutine erfi
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
      END FUNCTION RAN0

C ---------------------------------------------------------------------
      FUNCTION RANN(ISEED)
C     From numerical recipes
      PARAMETER(IA=7141,IC=54773,IM=259200)
      ISEED=MOD(ISEED*IA+IC,IM)
      RANN=FLOAT(ISEED)/FLOAT(IM)
      RETURN
      END FUNCTION RANN

c ---------------------------------------------------------------------- 