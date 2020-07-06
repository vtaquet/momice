MODULE INTLIB

USE VARIABLES

CONTAINS

subroutine avint ( ftab, xtab, ntab, a, b, result )
!
!***********************************************************************
!
!! AVINT estimates the integral of unevenly spaced data.
!
!
!  Discussion:
!
!    The method uses overlapping parabolas and smoothing.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    P E Hennion,
!    Algorithm 77,
!    Interpolation, Differentiation and Integration,
!    Communications of the Association for Computing Machinery,
!    Volume 5, page 96, 1962.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real FTAB(NTAB), the function values,
!    FTAB(I) = F(XTAB(I)).
!
!    Input, real XTAB(NTAB), the abscissas at which the
!    function values are given.  The XTAB's must be distinct
!    and in ascending order.
!
!    Input, integer NTAB, the number of entries in FTAB and
!    XTAB.  NTAB must be at least 3.
!
!    Input, real A, the lower limit of integration.  A should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Input, real B, the upper limit of integration.  B should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Output, real RESULT, the approximate value of the integral.
!
  implicit none
!
  integer(kind=long) ntab
!
  real(kind=dp) a
  real atemp
  real(kind=dp) b
  real btemp
  real ca
  real cb
  real cc
  real ctemp
  real(kind=dp) ftab(ntab)
  integer i
  integer ihi
  integer ilo
  integer ind
  real(kind=dp) result
  real sum1
  real syl
  real term1
  real term2
  real term3
  real x1
  real x2
  real x3
  real(kind=dp) xtab(ntab)
!
  if ( ntab < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVINT - Fatal error!'
    write ( *, '(a,i6)' ) '  NTAB is less than 3.  NTAB = ', ntab
    stop
  end if
 
  do i = 2, ntab
 
    if ( xtab(i) <= xtab(i-1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'AVINT - Fatal error!'
      write ( *, '(a)' ) '  XTAB(I) is not greater than XTAB(I-1).'
      write ( *, '(a,i6)' ) '  Here, I = ', I
      write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ', xtab(i-1)
      write ( *, '(a,g14.6)' ) '  XTAB(I) =   ', xtab(i)
      stop
    end if
 
  end do
 
  result = 0.0E+00
 
  if ( a == b ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVINT - Warning!'
    write ( *, '(a)' ) '  A = B, integral=0.'
    return
  end if
!
!  If A > B, temporarily switch A and B, and store sign.
!
  if ( a > b ) then
    syl = b
    b = a
    a = syl
    ind = -1
  else
    syl = a
    ind = 1
  end if
!
!  Bracket A and B between XTAB(ILO) and XTAB(IHI).
!
  ilo = 1
  ihi = ntab

  do i = 1, ntab
    if ( xtab(i) >= a ) then
      exit
    end if
    ilo = ilo + 1
  end do

  ilo = max ( 2, ilo )
  ilo = min ( ilo, ntab-1 )

  do i = 1, ntab
    if ( b >= xtab(i) ) then
      exit
    end if
    ihi = ihi - 1
  end do
  
  ihi = min ( ihi, ntab-1 )
  ihi = max ( ilo, ihi-1 )
!
!  Carry out approximate integration from XTAB(ILO) to XTAB(IHI).
!
  sum1 = 0.0E+00
 
  do i = ilo, ihi
 
    x1 = xtab(i-1)
    x2 = xtab(i)
    x3 = xtab(i+1)
 
    term1 = ftab(i-1) / ((x1-x2)*(x1-x3))
    term2 = ftab(i) / ((x2-x1)*(x2-x3))
    term3 = ftab(i+1) / ((x3-x1)*(x3-x2))
 
    atemp = term1 + term2 + term3
    btemp = -(x2+x3)*term1-(x1+x3)*term2-(x1+x2)*term3
    ctemp = x2*x3*term1+x1*x3*term2+x1*x2*term3
 
    if ( i <= ilo ) then
      ca = atemp
      cb = btemp
      cc = ctemp
    else
      ca = 0.5E+00 * ( atemp + ca )
      cb = 0.5E+00 * ( btemp + cb )
      cc = 0.5E+00 * ( ctemp + cc )
    end if
 
    sum1 = sum1 &
          + ca * ( x2**3 - syl**3 ) / 3.0E+00 &
          + cb * 0.5E+00 * ( x2**2 - syl**2 ) &
          + cc * ( x2 - syl )
 
    ca = atemp
    cb = btemp
    cc = ctemp
 
    syl = x2
 
  end do
 
  result = sum1 &
        + ca * ( b**3 - syl**3 ) / 3.0E+00 &
        + cb * 0.5E+00 * ( b**2 - syl**2 ) &
        + cc * ( b - syl )
!
!  Restore original values of A and B, reverse sign of integral
!  because of earlier switch.
!
  if ( ind /= 1 ) then
    ind = 1
    syl = b
    b = a
    a = syl
    result = -result
  end if
 
  return
end subroutine avint
!!$subroutine cadre ( func, a, b, abserr, relerr, error, result, ind )
!!$!
!!$!***********************************************************************
!!$!
!!$!! CADRE estimates the integral of F(X) from A to B.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    CADRE is the Cautious Adaptive Romberg Extrapolator.
!!$!
!!$!  Reference:
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!    Carl DeBoor and J R Rice,
!!$!    CADRE: An algorithm for numerical quadrature,
!!$!    Mathematic Software, pages 417-449,
!!$!    Academic Press, New York, 1971.
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real, external FUNC, the name of the function to be integrated.
!!$!    The user must declare the name an external parameter in the calling
!!$!    program, write a function routine of the form FUNCTION FUNC(X) which
!!$!    evaluates the function at X, and pass the name of the function
!!$!    in FUNC.
!!$!
!!$!    Input, real A, the lower limit of integration.
!!$!
!!$!    Input, real B, the upper limit of integration.
!!$!
!!$!    Input, real ABSERR, the absolute error tolerance.
!!$!
!!$!    Input, real RELERR, the relative error tolerance.
!!$!
!!$!    Output, real ERROR, an estimate of the absolute error.
!!$!
!!$!    Output, real RESULT, the approximate value of the integral.
!!$!
!!$!    Output, integer IND, reliability indicator.
!!$!    If IND <= 2, RESULT is very reliable.  Higher values of
!!$!    IND indicate less reliable values of RESULT.
!!$!
!!$  implicit none
!!$!
!!$  integer, parameter :: mxstge = 30
!!$  integer, parameter :: maxtbl = 10
!!$  integer, parameter :: maxts = 2049
!!$!
!!$  real a
!!$  real abserr
!!$  real ait(maxtbl)
!!$  logical aitken
!!$  real, parameter :: aitlow = 1.1E+00
!!$  real, parameter :: aittol = 0.1E+00
!!$  real astep
!!$  real b
!!$  real beg
!!$  real begin(mxstge)
!!$  real bma
!!$  real curest
!!$  real dif(maxtbl)
!!$  real diff
!!$  real end
!!$  real ergoal
!!$  real erra
!!$  real errer
!!$  real error
!!$  real errr
!!$  real est(mxstge)
!!$  real fbeg
!!$  real fbeg2
!!$  real fend
!!$  real fextm1
!!$  real fextrp
!!$  real finis(mxstge)
!!$  real fn
!!$  real fnsize
!!$  real, external :: func
!!$  logical h2conv
!!$  real h2next
!!$  real h2tfex
!!$  real, parameter :: h2tol = 0.15E+00
!!$  real hovn
!!$  integer i
!!$  integer ibeg
!!$  integer ibegs(mxstge)
!!$  integer iend
!!$  integer ii
!!$  integer iii
!!$  integer ind
!!$  integer istage
!!$  integer istep
!!$  integer istep2
!!$  integer it
!!$  integer l
!!$  integer lm1
!!$  integer n
!!$  integer n2
!!$  integer nnleft
!!$  real prever
!!$  real r(maxtbl)
!!$  logical reglar
!!$  logical reglsv(mxstge)
!!$  real relerr
!!$  real result
!!$  logical right
!!$  real rn(4)
!!$  real rnderr
!!$  real sing
!!$  real singnx
!!$  real slope
!!$  real stage
!!$  real step
!!$  real stepmn
!!$  real sum1
!!$  real sumabs
!!$  real t(maxtbl,maxtbl)
!!$  real tabs
!!$  real tabtlm
!!$  real, parameter :: tljump = 0.01E+00
!!$  real ts(2049)
!!$  real vint
!!$!
!!$  if ( a == b ) then
!!$    result = 0.0E+00
!!$    return
!!$  end if
!!$ 
!!$  begin(1:mxstge) = 0.0E+00
!!$  est(1:mxstge) = 0.0E+00
!!$  finis(1:mxstge) = 0.0E+00
!!$  ibegs(1:mxstge) = 0
!!$  reglsv(1:mxstge) = .false.
!!$ 
!!$  vint = 0.0E+00
!!$ 
!!$  rn(1:4) = (/ 0.7142005E+00, 0.3466282E+00, 0.8437510E+00, 0.1263305E+00 /)
!!$ 
!!$  rnderr = epsilon ( rnderr )
!!$  result = 0.0E+00
!!$  error = 0.0E+00
!!$  ind = 1
!!$  bma = abs ( b - a )
!!$  errr = min ( 0.1E+00, max ( abs ( relerr ), 10.0E+00*rnderr) )
!!$  erra = abs ( abserr )
!!$  stepmn = max ( bma / 2**mxstge, max ( bma, abs ( a ), abs ( b ) ) * rnderr )
!!$  stage = 0.5
!!$  istage = 1
!!$  curest = 0.0E+00
!!$  fnsize = 0.0E+00
!!$  prever = 0.0E+00
!!$  reglar = .false.
!!$  beg = a
!!$  fbeg = func(beg) / 2.0E+00
!!$  ts(1) = fbeg
!!$  ibeg = 1
!!$  end = b
!!$  fend = func(end) / 2.0E+00
!!$  ts(2) = fend
!!$  iend = 2
!!$ 
!!$10 continue
!!$ 
!!$  right = .false.
!!$ 
!!$20 continue
!!$
!!$  step = end - beg
!!$  astep = abs ( step )
!!$ 
!!$  if ( astep < stepmn ) then
!!$    ind = 5
!!$    result = curest + vint
!!$    return
!!$  end if
!!$ 
!!$  t(1,1) = fbeg+fend
!!$  tabs = abs ( fbeg ) + abs ( fend )
!!$  l = 1
!!$  n = 1
!!$  h2conv = .false.
!!$  aitken = .false.
!!$  go to 40
!!$ 
!!$30 continue
!!$ 
!!$40 continue
!!$ 
!!$  lm1 = l
!!$  l = l+1
!!$  n2 = n*2
!!$  fn = n2
!!$  istep = (iend-ibeg)/n
!!$
!!$  if ( istep > 1 ) then
!!$    go to 60
!!$  end if
!!$
!!$  ii = iend
!!$  iend = iend + n
!!$
!!$  if ( iend > maxts ) then
!!$    go to 440
!!$  end if
!!$
!!$  hovn = step / fn
!!$ 
!!$  iii = iend
!!$  do i = 1, n2, 2
!!$    ts(iii) = ts(ii)
!!$    ts(iii-1) = func(end-i*hovn)
!!$    iii = iii-2
!!$    ii = ii-1
!!$  end do
!!$ 
!!$  istep = 2
!!$ 
!!$60 continue
!!$ 
!!$  istep2 = ibeg+istep/2
!!$ 
!!$  sum1 = 0.0E+00
!!$  sumabs = 0.0E+00
!!$  do i = istep2, iend, istep
!!$    sum1 = sum1 + ts(i)
!!$    sumabs = sumabs + abs ( ts(i) )
!!$  end do
!!$ 
!!$  t(l,1) = t(l-1,1) / 2.0E+00 + sum1 / fn
!!$  tabs = tabs / 2.0E+00 + sumabs / fn
!!$ 
!!$  n = n2
!!$  it = 1
!!$  vint = step * t(l,1)
!!$  tabtlm = tabs * rnderr
!!$  fnsize = max ( fnsize, abs ( t(l,1) ) )
!!$  ergoal = max ( astep * rnderr * fnsize, &
!!$    stage * max ( erra , errr * abs ( curest+vint ) ) )
!!$  fextrp = 1.0E+00
!!$  do i = 1, lm1
!!$    fextrp = fextrp * 4.0E+00
!!$    t(i,l) = t(l,i) - t(l-1,i)
!!$    t(l,i+1) = t(l,i) + t(i,l) / ( fextrp - 1.0 )
!!$  end do
!!$ 
!!$  errer = astep * abs ( t(1,l) )
!!$  if ( l > 2 ) go to 90
!!$  if ( abs ( t(1,2) ) <= tabtlm ) go to 290
!!$  go to 40
!!$ 
!!$90 continue
!!$ 
!!$  do i = 2, lm1
!!$
!!$    if ( abs ( t(i-1,l) ) > tabtlm ) then
!!$      diff = t(i-1,lm1) / t(i-1,l)
!!$    else
!!$      diff = 0.0E+00
!!$    end if
!!$
!!$    t(i-1,lm1) = diff
!!$
!!$  end do
!!$ 
!!$  if ( abs ( 4.0 - t(1,lm1) ) <= h2tol ) go to 130
!!$  if ( t(1,lm1) == 0.0 ) go to 120
!!$  if ( abs ( 2.0 - abs ( t(1,lm1) ) ) < tljump ) go to 280
!!$  if (l==3) go to 30
!!$  h2conv = .false.
!!$  if ( abs ( ( t(1,lm1) - t(1,l-2) ) / t(1,lm1) ) <= aittol ) go to 160
!!$ 
!!$  if ( .not. reglar .and. l == 4 ) go to 30
!!$ 
!!$120 continue
!!$ 
!!$  if ( errer <= ergoal ) go to 310
!!$  go to 380
!!$
!!$130 continue
!!$
!!$  if ( .not. h2conv ) then
!!$    aitken = .false.
!!$    h2conv = .true.
!!$  end if
!!$
!!$140 continue
!!$
!!$  fextrp = 4.0E+00
!!$
!!$150 continue
!!$
!!$  it = it+1
!!$  vint = step * t(l,it)
!!$  errer = abs ( step / (fextrp-1.0) * t(it-1,l))
!!$  if ( errer <= ergoal ) go to 340
!!$  if ( it == lm1 ) go to 270
!!$  if ( t(it,lm1) == 0.0 ) go to 150
!!$  if ( t(it,lm1) <= fextrp ) go to 270
!!$
!!$  if ( abs ( t(it,lm1) / 4.0 - fextrp ) / fextrp < aittol ) then
!!$    fextrp = fextrp*4.0E+00
!!$  end if
!!$
!!$  go to 150
!!$ 
!!$160 continue
!!$
!!$  if ( t(1,lm1) < aitlow ) then
!!$    go to 380
!!$  end if
!!$ 
!!$  if ( .not. aitken ) then
!!$    h2conv = .false.
!!$    aitken = .true.
!!$  end if
!!$ 
!!$170 continue
!!$
!!$  fextrp = t(l-2,lm1)
!!$  if ( fextrp > 4.5 ) go to 140
!!$  if ( fextrp < aitlow ) go to 380
!!$
!!$  if ( abs ( fextrp - t(l-3,lm1) ) / t(1,lm1) > h2tol ) then
!!$    go to 380
!!$  end if
!!$
!!$  sing = fextrp
!!$  fextm1 = fextrp - 1.0E+00
!!$
!!$  ait(1) = 0.0E+00
!!$  do i = 2, l
!!$    ait(i) = t(i,1) + (t(i,1)-t(i-1,1)) / fextm1
!!$    r(i) = t(1,i-1)
!!$    dif(i) = ait(i) - ait(i-1)
!!$  end do
!!$
!!$  it = 2
!!$
!!$190 continue
!!$
!!$  vint = step*ait(l)
!!$
!!$200 continue
!!$
!!$  errer = errer / fextm1
!!$ 
!!$  if ( errer <= ergoal ) then
!!$    ind = max ( ind, 2 )
!!$    go to 340
!!$  end if
!!$ 
!!$210 continue
!!$
!!$  it = it+1
!!$  if ( it == lm1 ) go to 270
!!$
!!$  if ( it <= 3 ) then
!!$    h2next = 4.0E+00
!!$    singnx = 2.0E+00 * sing
!!$  end if
!!$
!!$  if ( h2next < singnx ) go to 230
!!$  fextrp = singnx
!!$  singnx = 2.0E+00 * singnx
!!$  go to 240
!!$
!!$230 continue
!!$
!!$  fextrp = h2next
!!$  h2next = 4.0E+00 * h2next
!!$
!!$240 continue
!!$ 
!!$  do i = it, lm1
!!$    if ( abs ( dif(i+1) ) > tabtlm ) then
!!$      r(i+1) = dif(i) / dif(i+1)
!!$    else
!!$      r(i+1) = 0.0E+00
!!$    end if
!!$  end do
!!$ 
!!$  h2tfex = -h2tol*fextrp
!!$  if ( r(l) - fextrp < h2tfex ) go to 270
!!$  if ( r(l-1) - fextrp < h2tfex ) go to 270
!!$  errer = astep * abs ( dif(l) )
!!$  fextm1 = fextrp - 1.0E+00
!!$  do i = it, l
!!$    ait(i) = ait(i)+dif(i) / fextm1
!!$    dif(i) = ait(i)-ait(i-1)
!!$  end do
!!$ 
!!$  go to 190
!!$ 
!!$270 continue
!!$
!!$  fextrp = max(prever/errer,aitlow)
!!$  prever = errer
!!$  if (l<5) go to 40
!!$  if (l-it>2.and.istage<mxstge) go to 370
!!$  if (errer/fextrp**(maxtbl-l)<ergoal) go to 40
!!$  go to 370
!!$ 
!!$280 continue
!!$
!!$  if ( errer > ergoal ) go to 370
!!$  diff = abs ( t(1,l) ) * 2.0E+00 * fn
!!$  go to 340
!!$ 
!!$290 continue
!!$
!!$  slope = (fend-fbeg) * 2.0E+00
!!$  fbeg2 = fbeg * 2.0E+00
!!$ 
!!$  do i = 1, 4
!!$    diff = abs ( func ( beg + rn(i) * step ) - fbeg2 - rn(i) * slope )
!!$    if ( diff > tabtlm ) go to 330
!!$  end do
!!$ 
!!$  go to 340
!!$ 
!!$310 continue
!!$
!!$  slope = (fend-fbeg)*2.0E+00
!!$  fbeg2 = fbeg*2.0E+00
!!$  i = 1
!!$ 
!!$320 continue
!!$
!!$  diff = abs ( func(beg+rn(i)*step) - fbeg2 - rn(i) * slope )
!!$ 
!!$330 continue
!!$
!!$  errer = max ( errer, astep * diff )
!!$  if (errer > ergoal) go to 380
!!$  i = i+1
!!$  if ( i <= 4 ) go to 320
!!$  ind = 3
!!$ 
!!$340 continue
!!$
!!$  result = result+vint
!!$  error = error+errer
!!$ 
!!$350 continue
!!$
!!$  if (right) go to 360
!!$  istage = istage-1
!!$  if (istage==0) return
!!$  reglar = reglsv(istage)
!!$  beg = begin(istage)
!!$  end = finis(istage)
!!$  curest = curest-est(istage+1)+vint
!!$  iend = ibeg-1
!!$  fend = ts(iend)
!!$  ibeg = ibegs(istage)
!!$  go to 400
!!$ 
!!$360 continue
!!$
!!$  curest = curest+vint
!!$  stage = stage*2.0E+00
!!$  iend = ibeg
!!$  ibeg = ibegs(istage)
!!$  end = beg
!!$  beg = begin(istage)
!!$  fend = fbeg
!!$  fbeg = ts(ibeg)
!!$  go to 10
!!$ 
!!$370 continue
!!$
!!$  reglar = .true.
!!$ 
!!$380 continue
!!$ 
!!$  if ( istage == mxstge ) then
!!$    ind = 5
!!$    result = curest+vint
!!$    return
!!$  end if
!!$ 
!!$390 continue
!!$
!!$  if (right) go to 410
!!$  reglsv(istage+1) = reglar
!!$  begin(istage) = beg
!!$  ibegs(istage) = ibeg
!!$  stage = stage/2.0E+00
!!$
!!$400 continue
!!$
!!$  right = .true.
!!$  beg = (beg+end)/2.0E+00
!!$  ibeg = (ibeg+iend)/2
!!$  ts(ibeg) = ts(ibeg) / 2.0E+00
!!$  fbeg = ts(ibeg)
!!$  go to 20
!!$
!!$410 continue
!!$
!!$  nnleft = ibeg-ibegs(istage)
!!$  if (end+nnleft>=maxts) go to 440
!!$  iii = ibegs(istage)
!!$  ii = iend
!!$  do i = iii, ibeg
!!$    ii = ii+1
!!$    ts(ii) = ts(i)
!!$  end do
!!$ 
!!$  do i = ibeg, ii
!!$    ts(iii) = ts(i)
!!$    iii = iii+1
!!$  end do
!!$ 
!!$  iend = iend+1
!!$  ibeg = iend-nnleft
!!$  fend = fbeg
!!$  fbeg = ts(ibeg)
!!$  finis(istage) = end
!!$  end = beg
!!$  beg = begin(istage)
!!$  begin(istage) = end
!!$  reglsv(istage) = reglar
!!$  istage = istage+1
!!$  reglar = reglsv(istage)
!!$  est(istage) = vint
!!$  curest = curest+est(istage)
!!$  go to 10
!!$
!!$440 continue
!!$
!!$  ind = 4
!!$
!!$460 continue
!!$
!!$  result = curest+vint
!!$
!!$  return
!!$end
!!$subroutine chinsp ( func, a, b, epsin, epsout, result )
!!$!
!!$!***********************************************************************
!!$!
!!$!! CHINSP estimates an integral using a modified Clenshaw-Curtis scheme.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    The integral is approximated by Chebyshev polyonomials over each
!!$!    subinterval.  These are integrated to give the approximate integral.
!!$!    If the error estimate is unsatisfactory, the integration is repeated
!!$!    with smaller intervals.
!!$!
!!$!    The internal parameter NUPPER is currently set to 9,
!!$!    corresponding to 1024 subintervals for the unfolded integral,
!!$!    and 1025 function evaluations.  This parameter may be changed
!!$!    if necessary.
!!$!
!!$!  Reference:
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!    T Havie
!!$!    BIT 9 (1969), pages 338-350.
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real, external FUNC, the name of the function to be
!!$!    integrated.  The user must declare the name an external
!!$!    parameter in the calling program, pass the name of the
!!$!    function in FUNC, and write a function of the form
!!$!
!!$!      FUNCTION FUNC(X)
!!$!
!!$!    which evaluates the function at the point X.
!!$!
!!$!    Input, real A, the lower limit of integration.
!!$!
!!$!    Input, real B, the upper limit of integration.
!!$!
!!$!    Input, real EPSIN, the relative error tolerance.
!!$!
!!$!    Output, real EPSOUT, estimated integration error.
!!$!
!!$!    Output, real RESULT, the approximate value of the integral.
!!$!
!!$  implicit none
!!$!
!!$  integer, parameter :: nupper = 9
!!$!
!!$  real a
!!$  real a0
!!$  real a1
!!$  real a2
!!$  real acof(257)
!!$  real alf
!!$  real alfnj
!!$  real alfno
!!$  real b
!!$  real bcof(257)
!!$  real bet
!!$  real betnj
!!$  real betno
!!$  real bounds
!!$  real ccof(513)
!!$  real cof
!!$  real cofmax
!!$  real const1
!!$  real const2
!!$  real deln
!!$  real deltan
!!$  real epsin
!!$  real epsout
!!$  real error
!!$  real etank
!!$  real, external :: func
!!$  real gamman
!!$  real hnstep
!!$  integer i
!!$  integer index
!!$  integer j
!!$  integer k
!!$  integer ksign
!!$  integer n
!!$  integer ncof
!!$  integer nhalf
!!$  integer nn
!!$  real, parameter :: one = 1.0E+00
!!$  real r1
!!$  real r2
!!$  real result
!!$  real rk
!!$  real rn
!!$  real rnderr
!!$  real rounde
!!$  real tend
!!$  real tnew
!!$  real triarg
!!$  real umid
!!$  real wmean
!!$  real xmin
!!$  real xplus
!!$  real xsink
!!$!
!!$  if ( a == b ) then
!!$    result = 0.0E+00
!!$    return
!!$  end if
!!$!
!!$!  ROUNDE = RNDERR*(R1+R2*N), where R1, R2 are two empirical constants.
!!$!
!!$!  Set coefficients in formula for accumulated roundoff error.
!!$!  N is the current number of function values used.
!!$!
!!$  rnderr = epsilon ( 1.0E+00 )
!!$ 
!!$  r1 = 1.0E+00
!!$  r2 = 2.0E+00
!!$  error = epsin
!!$!
!!$!  Integration interval parameters.
!!$!
!!$  alf = 0.5E+00 * ( b - a )
!!$  bet = 0.5E+00 * ( b + a )
!!$!
!!$!  Parameters for trigonometric recurrence relations.
!!$!
!!$  triarg = atan ( 1.0E+00 )
!!$  alfno = -1.0E+00
!!$!
!!$!  Parameters for integration stepsize and loops.
!!$!
!!$  rn = 2.0E+00
!!$  n = 2
!!$  nhalf = 1
!!$  hnstep = 1.0E+00
!!$!
!!$!  Initial calculation for the end-point approximation.
!!$!
!!$  const1 = 0.5E+00 * ( func(a) + func(b) )
!!$  const2 = func(bet)
!!$  acof(1) = 0.5E+00 * (const1+const2)
!!$  acof(2) = 0.5E+00 * (const1-const2)
!!$  bcof(2) = acof(2)
!!$  tend = 2.0E+00 * ( acof(1) - acof(2) / 3.0E+00 )
!!$!
!!$!  Start actual calculations.
!!$!
!!$  do i = 1, nupper
!!$!
!!$!  Compute function values.
!!$!
!!$    const1 = -sin(triarg)
!!$    const2 = 0.5E+00 * alfno / const1
!!$    alfno = const1
!!$    betno = const2
!!$    gamman = 1.0E+00 - 2.0E+00 * alfno**2
!!$    deltan = -2.0E+00 * alfno * betno
!!$    bcof(1) = 0.0E+00
!!$ 
!!$    do j = 1, nhalf
!!$      alfnj = gamman * const1 + deltan*const2
!!$      betnj = gamman * const2 - deltan*const1
!!$      xplus = alf * alfnj+bet
!!$      xmin = -alf * alfnj+bet
!!$      ccof(j) = func(xplus) + func(xmin)
!!$      bcof(1) = bcof(1) + ccof(j)
!!$      const1 = alfnj
!!$      const2 = betnj
!!$    end do
!!$ 
!!$    bcof(1) = 0.5E+00 * hnstep * bcof(1)
!!$!
!!$!  Calculation of first B-coefficient finished compute the higher
!!$!  coefficients if NHALF greater than one.
!!$!
!!$    if ( nhalf <= 1 ) go to 60
!!$    const1 = one
!!$    const2 = 0.0E+00
!!$    ncof = nhalf-1
!!$    ksign = -1
!!$ 
!!$    do k = 1, ncof
!!$!
!!$!  Compute trigonometric sum for B-coefficient.
!!$!
!!$      etank = gamman * const1 - deltan*const2
!!$      xsink = gamman * const2 + deltan*const1
!!$      cof = 2.0E+00 * ( 2.0E+00 * etank**2 - 1.0E+00 )
!!$      a2 = 0.0E+00
!!$      a1 = 0.0E+00
!!$      a0 = ccof(nhalf)
!!$ 
!!$      do j = 1, ncof
!!$        a2 = a1
!!$        a1 = a0
!!$        index = nhalf-j
!!$        a0 = ccof(index) + cof * a1 - a2
!!$      end do
!!$ 
!!$      bcof(k+1) = hnstep * (a0-a1) * etank
!!$      bcof(k+1) = ksign * bcof(k+1)
!!$      ksign = -ksign
!!$      const1 = etank
!!$      const2 = xsink
!!$ 
!!$    end do
!!$!
!!$!  Calculation of B-coefficients finished.
!!$!
!!$!  Compute new modified mid-point approximation when the interval
!!$!  of integration is divided in N equal sub intervals.
!!$!
!!$60  continue
!!$ 
!!$    umid = 0.0E+00
!!$    rk = rn
!!$    nn = nhalf+1
!!$    do k = 1, nn
!!$      index = nn+1-k
!!$      umid = umid+bcof(index)/(rk**2-one)
!!$      rk = rk-2.0E+00
!!$    end do
!!$ 
!!$    umid = -2.0E+00 * umid
!!$!
!!$!  Compute new C-coefficients for end-point approximation and largest
!!$!  absolute value of coefficients.
!!$!
!!$    nn = n+2
!!$    cofmax = 0.0E+00
!!$ 
!!$    do j = 1, nhalf
!!$      index = nn-j
!!$      ccof(j) = 0.5E+00 * (acof(j)+bcof(j))
!!$      ccof(index) = 0.5E+00 * (acof(j)-bcof(j))
!!$      const1 = abs ( ccof(j) )
!!$      cofmax = max ( cofmax, const1 )
!!$      const1 = abs ( ccof(index) )
!!$      cofmax = max ( cofmax, const1 )
!!$    end do
!!$ 
!!$    ccof(nhalf+1) = acof(nhalf+1)
!!$!
!!$!  Calculation of new coefficients finished.
!!$!
!!$!  Compute new end-point approximation when the interval of
!!$!  integration is divided in 2N equal sub intervals.
!!$!
!!$    wmean = 0.5E+00 * (tend+umid)
!!$    bounds = 0.5E+00 * (tend-umid)
!!$    deln = 0.0E+00
!!$    rk = 2.0E+00 * rn
!!$    do j = 1, nhalf
!!$      index = n+2-j
!!$      deln = deln+ccof(index) / (rk**2-one)
!!$      rk = rk-2.0E+00
!!$    end do
!!$ 
!!$    deln = -2.0E+00 * deln
!!$    tnew = wmean+deln
!!$    epsout = abs ( bounds / tnew )
!!$
!!$    if ( cofmax < rnderr ) then
!!$      go to 160
!!$    end if
!!$
!!$    rounde = rnderr*(r1+r2*rn)
!!$    if ( epsout < rounde ) epsout = rounde
!!$    if ( error < rounde ) error = rounde
!!$    if ( epsout > error ) go to 160
!!$!
!!$!  Required accuracy obtained or the maximum number of function
!!$!  values used without obtaining the required accuracy.
!!$!
!!$120 continue
!!$ 
!!$    n = 2*n+1
!!$    tend = alf*(tend+deln)
!!$    umid = alf*(umid+deln)
!!$    deln = alf*deln
!!$    result = alf*tnew
!!$    return
!!$!
!!$!  If I = NUPPER then the required accuracy is not obtained.
!!$!
!!$160 continue
!!$ 
!!$    if ( i == nupper ) go to 120
!!$ 
!!$    acof(1:n) = ccof(1:n)
!!$    acof(n+1) = ccof(n+1)
!!$    bcof(n+1) = ccof(n+1)
!!$    tend = tnew
!!$    nhalf = n
!!$    n = 2 * n
!!$    rn = 2.0E+00 * rn
!!$    hnstep = 0.5E+00 * hnstep
!!$    triarg = 0.5E+00 * triarg
!!$ 
!!$  end do
!!$ 
!!$  return
!!$end
!!$subroutine class ( kind, n, alpha, beta, b, a, muzero )
!!$!
!!$!***********************************************************************
!!$!
!!$!! CLASS sets recurrence coeeficients for various orthogonal polynomials.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    CLASS supplies the coefficients A(J), B(J) of the recurrence relation
!!$!
!!$!      B(J)*P(J) (X) = (X-A(J))*P(J-1)(X) - B(J-1)*P(J-2)(X)
!!$!
!!$!    for the various classical (normalized) orthogonal polynomials,
!!$!    and the zero-th moment
!!$!
!!$!      MUZERO = Integral W(X) DX
!!$!
!!$!    of the given polynomial's weight function W(X).  Since the
!!$!    polynomials are orthonormalized, the tridiagonal matrix is
!!$!    guaranteed to be symmetric.
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real ALPHA, BETA, parameters needed for Laguerre and Jacobi
!!$!    polynomials.
!!$!
!!$!    Input, integer KIND, specifies which polynomial is to be handled:
!!$!
!!$!    1: Legendre polynomials P(X) on (-1, +1),
!!$!    W(X) = 1.
!!$!
!!$!    2: Chebyshev polynomials of the first kind T(X) on (-1, +1),
!!$!    W(X) = 1 / SQRT(1 - X*X)
!!$!
!!$!    3: Chebyshev polynomials of the second kind U(X) on (-1, +1),
!!$!    W(X) = SQRT(1 - X*X)
!!$!
!!$!    4: Hermite polynomials H(X) on (-infinity,+infinity),
!!$!    W(X) = EXP(-X**2)
!!$!
!!$!    5: Jacobi polynomials P(ALPHA,BETA)(X) on (-1, +1),
!!$!    W(X) = (1-X)**ALPHA + (1+X)**BETA,
!!$!    ALPHA and BETA greater than -1.
!!$!
!!$!    6: Laguerre polynomials, L(ALPHA)(X) on (0, +infinity),
!!$!    W(X) = EXP(-X) * X**ALPHA,
!!$!    ALPHA greater than -1.
!!$!
!!$!    Input, integer N, specifies the number of coefficients to
!!$!    calculate.
!!$!
!!$!    Input, real ALPHA, the value of the ALPHA parameter,
!!$!    required only for Jacobi or Laguerre polynomials.
!!$!
!!$!    Input, real BETA, the value of the BETA parameter,
!!$!    required only for Jacobi polynomials.
!!$!
!!$!    Output, real B(N-1), the offdiagonal coefficients.
!!$!
!!$!    Output, real A(N), the diagonal coefficients.
!!$!
!!$!    Output, real MUZERO, the zero-th moment, Integral W(X) DX,
!!$!    of the polynomial's weight function over its interval of
!!$!    definition.
!!$!
!!$  implicit none
!!$!
!!$  integer n
!!$!
!!$  real a(n)
!!$  real abi
!!$  real alpha
!!$  real b(n-1)
!!$  real beta
!!$  real gamma
!!$  integer i
!!$  integer kind
!!$  real muzero
!!$  real pi
!!$!
!!$!  KIND = 1:
!!$!
!!$!  Legendre polynomials P(X) on (-1, +1),
!!$!  W(X) = 1.
!!$!
!!$  if ( kind == 1 ) then
!!$ 
!!$    muzero = 2.0E+00
!!$ 
!!$    a(1:n) = 0.0E+00
!!$ 
!!$    do i = 1, n-1
!!$      b(i) = real(i) / sqrt(4.0*real(i*i) - 1.0)
!!$    end do
!!$!
!!$!  KIND = 2:
!!$!
!!$!  Chebyshev polynomials of the first kind T(X) on (-1, +1),
!!$!  W(X) = 1 / SQRT(1 - X*X)
!!$!
!!$  else if ( kind == 2 ) then
!!$ 
!!$    muzero = pi()
!!$    a(1:n) = 0.0E+00
!!$    b(1) = sqrt ( 0.5E+00 )
!!$    b(1:n-1) = 0.5E+00
!!$!
!!$!  KIND = 3:
!!$!
!!$!  Chebyshev polynomials of the second kind U(X) on (-1, +1),
!!$!  W(X) = SQRT(1 - X*X)
!!$!
!!$  else if ( kind == 3 ) then
!!$ 
!!$    muzero = pi() / 2.0E+00
!!$    a(1:n) = 0.0E+00
!!$    b(1:n-1) = 0.5E+00
!!$!
!!$!  KIND = 4:
!!$!
!!$!  Hermite polynomials H(X) on (-infinity,+infinity),
!!$!  W(X) = EXP(-X**2)
!!$!
!!$  else if ( kind == 4 ) then
!!$ 
!!$    muzero = sqrt ( pi() )
!!$    a(1:n) = 0.0E+00
!!$    do i = 1, n-1
!!$      b(i) = sqrt ( real(i) / 2.0E+00 )
!!$    end do
!!$!
!!$!  KIND = 5:
!!$!
!!$!  Jacobi polynomials P(ALPHA,BETA)(X) on (-1, +1),
!!$!  W(X) = (1-X)**ALPHA + (1+X)**BETA,
!!$!  ALPHA and BETA greater than -1
!!$!
!!$  else if ( kind == 5 ) then
!!$ 
!!$    muzero = 2.0**( alpha + beta + 1.0E+00 ) * gamma ( alpha + 1.0 ) &
!!$      * gamma ( beta + 1.0 ) / gamma ( 2.0 + alpha + beta )
!!$ 
!!$    do i = 1, n
!!$      a(i) = ( beta**2 - alpha**2 )/ &
!!$        ((2.0*(i-1)+alpha+beta) * ( 2.0 * real ( i ) + alpha + beta ) )
!!$    end do
!!$ 
!!$    abi = 2.0E+00 + alpha+beta
!!$    b(1) = sqrt ( 4.0*(1.0E+00 + alpha)*(1.0 + beta)/((abi + 1.0)*abi*abi))
!!$ 
!!$    do i = 2, n-1
!!$      abi = real ( 2 * i ) + alpha + beta
!!$      b(i) = sqrt ( 4.0E+00 * real ( i ) * ( real ( i ) + alpha ) &
!!$        * ( i + beta ) * ( i + alpha + beta) / &
!!$        ( ( abi*abi - 1.0E+00 ) * abi * abi ) )
!!$    end do
!!$!
!!$!  KIND = 6:
!!$!
!!$!  Laguerre polynomials
!!$!
!!$!  L(ALPHA)(X) on (0, +infinity),
!!$!  W(X) = EXP(-X) * X**ALPHA,
!!$!  ALPHA greater than -1.
!!$!
!!$  else if ( kind == 6 ) then
!!$ 
!!$    muzero = gamma ( alpha + 1.0E+00 )
!!$ 
!!$    do i = 1, n
!!$      a(i) = 2.0E+00 * real ( i ) - 1.0E+00 + alpha
!!$    end do
!!$ 
!!$    do i = 1, n-1
!!$      b(i) = sqrt ( real ( i ) * ( real ( i ) + alpha ) )
!!$    end do
!!$ 
!!$  end if
!!$ 
!!$  return
!!$end
!!$subroutine cspint ( ftab, xtab, ntab, a, b, y, e, work, result )
!!$!
!!$!***********************************************************************
!!$!
!!$!! CSPINT estimates the integral of a tabulated function.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    The routine is given the value of a function F(X) at a set of 
!!$!    nodes XTAB, and estimates
!!$!
!!$!      INTEGRAL (A to B) F(X) DX
!!$!
!!$!    by computing the cubic natural spline S(X) that interpolates
!!$!    F(X) at the nodes, and then computing
!!$!
!!$!      INTEGRAL (A to B) S(X) DX
!!$!
!!$!    exactly.
!!$!
!!$!    Other output from the program includes the definite integral
!!$!    from X(1) to X(I) of S(X), and the coefficients necessary for
!!$!    the user to evaluate the spline S(X) at any point.
!!$!
!!$!  Reference:
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real FTAB(NTAB), contains the tabulated values of
!!$!    the function, FTAB(I) = F(XTAB(I)).
!!$!
!!$!    Input, real XTAB(NTAB), contains the points at which the
!!$!    function was evaluated.  The XTAB's must be distinct and
!!$!    in ascending order.
!!$!
!!$!    Input, integer NTAB, the number of entries in FTAB and
!!$!    XTAB.  NTAB must be at least 3.
!!$!
!!$!    Input, real A, lower limit of integration.
!!$!
!!$!    Input, real B, upper limit of integration.
!!$!
!!$!    Output, real Y(3,NTAB), will contain the coefficients
!!$!    of the interpolating natural spline over each subinterval.
!!$!
!!$!    For XTAB(I) <= X <= XTAB(I+1),
!!$!
!!$!      S(X) = FTAB(I) + Y(1,I)*(X-XTAB(I))
!!$!                   + Y(2,I)*(X-XTAB(I))**2
!!$!                   + Y(3,I)*(X-XTAB(I))**3
!!$!
!!$!    Output, real E(NTAB), E(I) = the definite integral from
!!$!    XTAB(1) to XTAB(I) of S(X).
!!$!
!!$!    Workspace, real WORK(NTAB).
!!$!
!!$!    Output, real RESULT, the estimated value of the integral.
!!$!
!!$  implicit none
!!$!
!!$  integer ntab
!!$!
!!$  real a
!!$  real b
!!$  real e(ntab)
!!$  real ftab(ntab)
!!$  integer i
!!$  integer j
!!$  real r
!!$  real result
!!$  real s
!!$  real term
!!$  real u
!!$  real work(ntab)
!!$  real xtab(ntab)
!!$  real y(3,ntab)
!!$!
!!$  if ( ntab < 3 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'CSPINT - Fatal error!'
!!$    write ( *, '(a,i6)' ) '  NTAB must be at least 3, but input NTAB = ',ntab
!!$    stop
!!$  end if
!!$ 
!!$  do i = 1, ntab-1
!!$ 
!!$    if ( xtab(i+1) <= xtab(i) ) then
!!$      write ( *, '(a)' ) ' '
!!$      write ( *, '(a)' ) 'CSPINT - Fatal error!'
!!$      write ( *, '(a)' ) '  Nodes not in strict increasing order.'
!!$      write ( *, '(a,i6)' ) '  XTAB(I) <= XTAB(I-1) for I=',i
!!$      write ( *, '(a,g14.6)' ) '  XTAB(I) = ',xtab(i)
!!$      write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ',xtab(i-1)
!!$      stop
!!$    end if
!!$ 
!!$  end do
!!$ 
!!$  s = 0.0E+00
!!$  do i = 1, ntab-1
!!$    r = ( ftab(i+1) - ftab(i) ) / ( xtab(i+1) - xtab(i) )
!!$    y(2,i) = r - s
!!$    s = r
!!$  end do
!!$ 
!!$  result = 0.0E+00
!!$  s = 0.0E+00
!!$  r = 0.0E+00
!!$  y(2,1) = 0.0E+00
!!$  y(2,ntab) = 0.0E+00
!!$ 
!!$  do i = 2, ntab-1
!!$    y(2,i) = y(2,i)+r*y(2,i-1)
!!$    work(i) = 2.0E+00 * ( xtab(i-1) - xtab(i+1) ) - r * s
!!$    s = xtab(i+1) - xtab(i)
!!$    r = s / work(i)
!!$  end do
!!$ 
!!$  do j = 2, ntab-1
!!$    i = ntab+1-j
!!$    y(2,i) = ((xtab(i+1)-xtab(i))*y(2,i+1)-y(2,i)) / work(i)
!!$  end do
!!$ 
!!$  do i = 1, ntab-1
!!$    s = xtab(i+1)-xtab(i)
!!$    r = y(2,i+1)-y(2,i)
!!$    y(3,i) = r / s
!!$    y(2,i) = 3.0E+00 * y(2,i)
!!$    y(1,i) = (ftab(i+1)-ftab(i)) / s-(y(2,i)+r)*s
!!$  end do
!!$ 
!!$  e(1) = 0.0E+00
!!$  do i = 1, ntab-1
!!$    s = xtab(i+1)-xtab(i)
!!$    term = (((y(3,i)* 0.25E+00 *s+y(2,i) / 3.0 ) *s+y(1,i)* 0.5 )*s+ftab(i))*s
!!$    e(i+1) = e(i) + term
!!$  end do
!!$!
!!$!  Determine where the endpoints A and B lie in the mesh of XTAB's.
!!$!
!!$  r = a
!!$  u = 1.0E+00
!!$ 
!!$  do j = 1, 2
!!$!
!!$!  The endpoint is less than or equal to XTAB(1).
!!$!
!!$    if ( r <= xtab(1) ) then
!!$      result = result-u*((r-xtab(1))*y(1,1)*0.5E+00 +ftab(1))*(r-xtab(1))
!!$!
!!$!  The endpoint is greater than or equal to XTAB(NTAB).
!!$!
!!$    else if ( r >= xtab(ntab) ) then
!!$
!!$      result = result-u*(e(ntab)+(r-xtab(ntab))*(ftab(ntab)+ &
!!$        0.5E+00 *(ftab(ntab-1)+(xtab(ntab)-xtab(ntab-1))*y(1,ntab-1)) &
!!$        *(r-xtab(ntab))))
!!$!
!!$!  The endpoint is strictly between XTAB(1) and XTAB(NTAB).
!!$!
!!$    else
!!$      do i = 1,ntab-1
!!$ 
!!$        if ( r <= xtab(i+1) ) then
!!$          r = r-xtab(i)
!!$          result = result-u*(e(i)+(((y(3,i)*0.25*r+y(2,i)/3.0)*r &
!!$            +y(1,i)*0.5E+00 )*r+ftab(i))*r)
!!$          go to 120
!!$        end if
!!$ 
!!$      end do
!!$ 
!!$    end if
!!$ 
!!$  120   continue
!!$ 
!!$    u = -1.0E+00
!!$    r = b
!!$ 
!!$  end do
!!$ 
!!$  return
!!$end
!!$subroutine cubint ( ftab, xtab, ntab, ia, ib, result, error )
!!$!
!!$!***********************************************************************
!!$!
!!$!! CUBINT approximates an integral using cubic interpolation of data.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    The integral to be approximated is
!!$! 
!!$!      INTEGRAL (XTAB(IB) to XTAB(IA)) F(X) DX
!!$!
!!$!    The routine estimates the error in integration.
!!$!
!!$!  Reference:
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!    P E Gill and G F Miller
!!$!    An algorithm for the integration of unequally spaced data,
!!$!    Comput J, Number 15, 1972, pages 80-83.
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real FTAB(NTAB), contains the tabulated function
!!$!    values, FTAB(I) = F(XTAB(I)).
!!$!
!!$!    Input, real XTAB(NTAB), contains the points at which the
!!$!    function was tabulated.  XTAB should contain distinct
!!$!    values, given in ascending order.
!!$!
!!$!    Input, integer NTAB, the number of tabulated points.
!!$!    NTAB must be at least 4.
!!$!
!!$!    Input, integer IA, the entry of XTAB at which integration
!!$!    is to begin.  IA must be no less than 1 and no greater
!!$!    than NTAB.
!!$!
!!$!    Input, integer IB, the entry of XTAB at which integration
!!$!    is to end.  IB must be no less than 1 and no greater than
!!$!    NTAB.
!!$!
!!$!    Output, real RESULT, the approximate value of the
!!$!    integral from XTAB(IA) to XTAB(IB) of the function.
!!$!
!!$!    Output, real ERROR, an estimate of the error in
!!$!    integration.
!!$!
!!$  implicit none
!!$!
!!$  integer ntab
!!$!
!!$  real c
!!$  real d1
!!$  real d2
!!$  real d3
!!$  real error
!!$  real ftab(ntab)
!!$  real h1
!!$  real h2
!!$  real h3
!!$  real h4
!!$  integer i
!!$  integer ia
!!$  integer ib
!!$  integer ind
!!$  integer it
!!$  integer j
!!$  integer k
!!$  real r1
!!$  real r2
!!$  real r3
!!$  real r4
!!$  real result
!!$  real s
!!$  real term
!!$  real xtab(ntab)
!!$!
!!$  result = 0.0E+00
!!$  error = 0.0E+00
!!$ 
!!$  if ( ia == ib ) then
!!$    return
!!$  end if
!!$ 
!!$  if ( ntab < 4 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'CUBINT - Fatal error!'
!!$    write ( *, '(a,i6)' ) '  NTAB must be at least 4, but input NTAB = ',ntab
!!$    stop
!!$  end if
!!$ 
!!$  if ( ia < 1 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'CUBINT - Fatal error!'
!!$    write ( *, '(a,i6)' ) '  IA must be at least 1, but input IA = ',ia
!!$    stop
!!$  end if
!!$ 
!!$  if ( ia > ntab ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'CUBINT - Fatal error!'
!!$    write ( *, '(a,i6)' ) '  IA must be <= NTAB, but input IA=',ia
!!$    stop
!!$  end if
!!$ 
!!$  if ( ib < 1 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'CUBINT - Fatal error!'
!!$    write ( *, '(a,i6)' ) '  IB must be at least 1, but input IB = ',ib
!!$    stop
!!$  end if
!!$ 
!!$  if ( ib > ntab ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'CUBINT - Fatal error!'
!!$    write ( *, '(a,i6)' ) '  IB must be <= NTAB, but input IB=',ib
!!$    stop
!!$  end if
!!$!
!!$!  Temporarily switch IA and IB, and store minus sign in IND
!!$!  so that, while integration is carried out from low X's
!!$!  to high ones, the sense of the integral is preserved.
!!$!
!!$  if ( ia > ib ) then
!!$    ind = -1
!!$    it = ib
!!$    ib = ia
!!$    ia = it
!!$  else
!!$    ind = 1
!!$  end if
!!$ 
!!$  s = 0.0E+00
!!$  c = 0.0E+00
!!$  r4 = 0.0E+00
!!$  j = ntab-2
!!$  if ( ia < ntab-1 .or. ntab == 4 ) then
!!$    j=max(3,ia)
!!$  end if
!!$
!!$  k = 4
!!$  if ( ib > 2 .or. ntab == 4 ) then
!!$    k=min(ntab,ib+2)-1
!!$  end if
!!$ 
!!$  do i = j, k
!!$ 
!!$    if ( i <= j ) then
!!$ 
!!$      h2 = xtab(j-1)-xtab(j-2)
!!$      d3 = (ftab(j-1)-ftab(j-2)) / h2
!!$      h3 = xtab(j)-xtab(j-1)
!!$      d1 = (ftab(j)-ftab(j-1)) / h3
!!$      h1 = h2+h3
!!$      d2 = (d1-d3)/h1
!!$      h4 = xtab(j+1)-xtab(j)
!!$      r1 = (ftab(j+1)-ftab(j)) / h4
!!$      r2 = (r1-d1) / (h4+h3)
!!$      h1 = h1+h4
!!$      r3 = (r2-d2) / h1
!!$ 
!!$      if ( ia <= 1 ) then
!!$        result = h2 * (ftab(1)+h2*(0.5*d3-h2*(d2/6.0-(h2+h3+h3)*r3/12.)))
!!$        s = -h2**3 * (h2*(3.0*h2+5.0*h4)+10.0*h3*h1)/60.0
!!$      end if
!!$ 
!!$    else
!!$ 
!!$      h4 = xtab(i+1)-xtab(i)
!!$      r1 = (ftab(i+1)-ftab(i))/h4
!!$      r4 = h4+h3
!!$      r2 = (r1-d1)/r4
!!$      r4 = r4+h2
!!$      r3 = (r2-d2)/r4
!!$      r4 = (r3-d3)/(r4+h1)
!!$ 
!!$    end if
!!$ 
!!$    if ( i > ia .and. i <= ib ) then
!!$ 
!!$      term = h3*((ftab(i)+ftab(i-1))*0.5-h3*h3*(d2+r2+(h2-h4)*r3) / 12.0 )
!!$      result = result+term
!!$      c = h3**3*(2.0E+00 *h3*h3+5.*(h3*(h4+h2) + 2.0 * h2 * h4 ) ) / 120.0E+00
!!$      error = error+(c+s)*r4
!!$ 
!!$      if ( i /= j ) then
!!$        s = c
!!$      else
!!$        s = s+c+c
!!$      end if
!!$ 
!!$    else
!!$ 
!!$      error = error+r4*s
!!$ 
!!$    end if
!!$ 
!!$    if ( i >= k ) then
!!$ 
!!$      if ( ib >= ntab ) then
!!$        term = h4*(ftab(ntab) - h4*(0.5*r1+h4*(r2/6.0 +(h3+h3+h4)*r3/12.)))
!!$        result = result + term
!!$        error = error - h4**3 * r4 * &
!!$          ( h4 * ( 3.0 * h4 + 5.0 * h2 ) &
!!$          + 10.0 * h3 * ( h2 + h3 + h4 ) ) / 60.0E+00
!!$      end if
!!$ 
!!$      if ( ib >= ntab-1 ) error=error+s*r4
!!$    else
!!$      h1 = h2
!!$      h2 = h3
!!$      h3 = h4
!!$      d1 = r1
!!$      d2 = r2
!!$      d3 = r3
!!$    end if
!!$ 
!!$  end do
!!$!
!!$!  Restore original values of IA and IB, reverse signs
!!$!  of RESULT and ERROR, to account for integration
!!$!  that proceeded from high X to low X.
!!$!
!!$  if ( ind /= 1 ) then
!!$    it = ib
!!$    ib = ia
!!$    ia = it
!!$    result = -result
!!$    error = -error
!!$  end if
!!$ 
!!$  return
!!$end
!!$subroutine filon_cos ( ftab, a, b, ntab, t, result )
!!$!
!!$!***********************************************************************
!!$!
!!$!! FILON_COS uses Filon's method on integrals with a cosine factor.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    The integral to be approximated has the form:
!!$!
!!$!      Integral ( A <= X <= B ) F(X) * COS(T*X) dX
!!$!
!!$!    where T is user specified.
!!$!
!!$!    The function is interpolated over each subinterval by
!!$!    a parabolic arc.
!!$!
!!$!  Reference:
!!$!
!!$!    Abramowitz and Stegun,
!!$!    Handbook of Mathematical Functions,
!!$!    pages 890-891.
!!$!
!!$!    S M Chase and L D Fosdick,
!!$!    Algorithm 353, Filon Quadrature,
!!$!    Communications of the Association for Computing Machinery,
!!$!    Volume 12, 1969, pages 457-458.
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!  Modified:
!!$!
!!$!    19 February 2002
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real FTAB(NTAB), contains the value of the function
!!$!    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(NTAB-1).
!!$!
!!$!    Input, real A, B, the limits of integration.
!!$!
!!$!    Input, integer NTAB, the number of data points.
!!$!    NTAB must be odd, and greater than 1.
!!$!
!!$!    Input, real T, the multiplier of the X argument of the cosine.
!!$!
!!$!    Output, real RESULT, the approximate value of the integral.
!!$!
!!$  implicit none
!!$!
!!$  integer ntab
!!$!
!!$  real a
!!$  real alpha
!!$  real b
!!$  real beta
!!$  real c2n
!!$  real c2nm1
!!$  real cost
!!$  real ftab(ntab)
!!$  real gamma
!!$  real h
!!$  real result
!!$  real sint
!!$  real t
!!$  real theta
!!$  real xtab(ntab)
!!$!
!!$  if ( a == b ) then
!!$    result = 0.0E+00
!!$    return
!!$  end if
!!$ 
!!$  if ( ntab <= 1 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'FILON_COS - Fatal error!'
!!$    write ( *, '(a)' ) '  NTAB < 2'
!!$    write ( *, '(a,i6)' ) '  NTAB = ', ntab
!!$    stop
!!$  end if
!!$ 
!!$  if ( mod ( ntab, 2 ) /= 1 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'FILON_COS - Fatal error!'
!!$    write ( *, '(a)' ) '  NTAB must be odd.'
!!$    write ( *, '(a,i6)' ) '  NTAB = ', ntab
!!$    stop
!!$  end if
!!$!
!!$!  Set up a vector of the NTAB X values.
!!$! 
!!$  call rvec_even ( a, b, ntab, xtab )
!!$
!!$  h = ( b - a ) / real ( ntab - 1 )
!!$  theta = t * h
!!$  sint = sin ( theta )
!!$  cost = cos ( theta )
!!$
!!$  alpha = ( theta**2 + theta * sint * cost &
!!$    - 2.0E+00 * sint**2 ) / theta**3
!!$
!!$  beta = ( 2.0E+00 * theta + 2.0E+00 * theta * cost**2 &
!!$    - 4.0E+00 * sint * cost ) / theta**3
!!$
!!$  gamma = 4.0E+00 * ( sint - theta * cost ) / theta**3
!!$  
!!$  c2n = sum ( ftab(1:ntab:2) * cos ( t * xtab(1:ntab:2) ) ) &
!!$    - 0.5E+00 * ( ftab(ntab) * cos ( t * xtab(ntab) ) &
!!$                + ftab(1) * cos ( t * xtab(1) ) )
!!$
!!$  c2nm1 = sum ( ftab(2:ntab-1:2) * cos ( t * xtab(2:ntab-1:2) ) )
!!$ 
!!$  result = h * ( &
!!$      alpha * ( ftab(ntab) * sin ( t * xtab(ntab) ) & 
!!$              - ftab(1) * sin ( t * xtab(1) ) ) &
!!$    + beta * c2n &
!!$    + gamma * c2nm1 )
!!$
!!$  return
!!$end
!!$subroutine filon_sin ( ftab, a, b, ntab, t, result )
!!$!
!!$!***********************************************************************
!!$!
!!$!! FILON_SIN uses Filon's method on integrals with a sine factor.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    The integral to be approximated has the form
!!$!
!!$!      Integral ( A <= X <= B ) F(X) * SIN(T*X) dX
!!$!
!!$!    where T is user specified.
!!$!
!!$!    The function is interpolated over each subinterval by
!!$!    a parabolic arc.
!!$!
!!$!  Reference:
!!$!
!!$!    Abramowitz and Stegun,
!!$!    Handbook of Mathematical Functions,
!!$!    pages 890-891.
!!$!
!!$!    S M Chase and L D Fosdick,
!!$!    Algorithm 353, Filon Quadrature,
!!$!    Communications of the Association for Computing Machinery,
!!$!    Volume 12, 1969, pages 457-458.
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!  Modified:
!!$!
!!$!    19 February 2002
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real FTAB(NTAB), contains the value of the function
!!$!    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(NTAB-1).
!!$!
!!$!    Input, real A, B, the limits of integration.
!!$!
!!$!    Input, integer NTAB, the number of data points, including the
!!$!    endpoints.  NTAB must be odd, and greater than 1.
!!$!
!!$!    Input, real T, multiplier of the X argument of the sine.
!!$!
!!$!    Output, real RESULT, the approximate value of the integral.
!!$!
!!$  implicit none
!!$!
!!$  integer ntab
!!$!
!!$  real a
!!$  real alpha
!!$  real b
!!$  real beta
!!$  real cost
!!$  real ftab(ntab)
!!$  real gamma
!!$  real h
!!$  real result
!!$  real s2n
!!$  real s2nm1
!!$  real sint
!!$  real t
!!$  real theta
!!$  real xtab(ntab)
!!$!
!!$  if ( a == b ) then
!!$    result = 0.0E+00
!!$    return
!!$  end if
!!$ 
!!$  if ( ntab <= 1 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'FILON_SIN - Fatal error!'
!!$    write ( *, '(a)' ) '  NTAB < 2'
!!$    write ( *, '(a,i6)' ) '  NTAB = ',ntab
!!$    stop
!!$  end if
!!$ 
!!$  if ( mod ( ntab, 2 ) /= 1 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'FILON_SIN - Fatal error!'
!!$    write ( *, '(a)' ) '  NTAB must be odd.'
!!$    write ( *, '(a,i6)' ) '  NTAB = ',ntab
!!$    stop
!!$  end if
!!$!
!!$!  Set up a vector of the NTAB X values.
!!$! 
!!$  call rvec_even ( a, b, ntab, xtab )
!!$
!!$  h = ( b - a ) / real ( ntab - 1 )
!!$  theta = t * h
!!$  sint = sin ( theta )
!!$  cost = cos ( theta )
!!$ 
!!$  alpha = ( theta**2 + theta * sint * cost &
!!$    - 2.0E+00 * sint**2 ) / theta**3
!!$
!!$  beta = ( 2.0E+00 * theta + 2.0E+00 * theta * cost**2 &
!!$    - 4.0E+00 * sint * cost ) / theta**3
!!$
!!$  gamma = 4.0E+00 * ( sint - theta * cost ) / theta**3
!!$   
!!$  s2n = sum ( ftab(1:ntab:2) * sin ( t * xtab(1:ntab:2) ) ) &
!!$    - 0.5E+00 * ( ftab(ntab) * sin ( t * xtab(ntab) ) &
!!$                + ftab(1) * sin ( t * xtab(1) ) )
!!$
!!$  s2nm1 = sum ( ftab(2:ntab-1:2) * sin ( t * xtab(2:ntab-1:2) ) )
!!$
!!$  result = h * ( &
!!$      alpha * ( ftab(1) * cos ( t * xtab(1) ) &
!!$              - ftab(ntab) * cos ( t * xtab(ntab) ) ) &
!!$    + beta * s2n &
!!$    + gamma * s2nm1 )
!!$ 
!!$  return
!!$end
!!$function gamma ( x )
!!$!
!!$!*******************************************************************************
!!$!
!!$!! GAMMA calculates the Gamma function for a real argument X.
!!$!
!!$!
!!$!  Definition:
!!$!
!!$!    GAMMA(X) = Integral ( 0 <= T <= Infinity ) T**(X-1) EXP(-T) DT
!!$!
!!$!  Recursion:
!!$!
!!$!    GAMMA(X+1) = X * GAMMA(X)
!!$!
!!$!  Special values:
!!$!
!!$!    GAMMA(0.5) = SQRT(PI)
!!$!    If N is a positive integer, GAMMA(N+1) = N!, the standard factorial.
!!$!
!!$!  Discussion:
!!$!
!!$!    Computation is based on an algorithm outlined in reference 1.
!!$!    The program uses rational functions that approximate the GAMMA
!!$!    function to at least 20 significant decimal digits.  Coefficients
!!$!    for the approximation over the interval (1,2) are unpublished.
!!$!    Those for the approximation for X .GE. 12 are from reference 2.
!!$!    The accuracy achieved depends on the arithmetic system, the
!!$!    compiler, the intrinsic functions, and proper selection of the
!!$!    machine-dependent constants.
!!$!
!!$!  Machine-dependent constants:
!!$!
!!$!    BETA: radix for the floating-point representation.
!!$!    MAXEXP: the smallest positive power of BETA that overflows.
!!$!    XBIG: the largest argument for which GAMMA(X) is representable
!!$!      in the machine, i.e., the solution to the equation
!!$!      GAMMA(XBIG) = BETA**MAXEXP.
!!$!    XINF: the largest machine representable floating-point number;
!!$!      approximately BETA**MAXEXP.
!!$!    EPS: the smallest positive floating-point number such that
!!$!      1.0+EPS .GT. 1.0.
!!$!    XMININ: the smallest positive floating-point number such that
!!$!      1/XMININ is machine representable.
!!$!
!!$!    Approximate values for some important machines are:
!!$!
!!$!                               BETA       MAXEXP        XBIG
!!$!
!!$!    CRAY-1         (S.P.)        2         8191        966.961
!!$!    Cyber 180/855
!!$!      under NOS    (S.P.)        2         1070        177.803
!!$!    IEEE (IBM/XT,
!!$!      SUN, etc.)   (S.P.)        2          128        35.040
!!$!    IEEE (IBM/XT,
!!$!      SUN, etc.)   (D.P.)        2         1024        171.624
!!$!    IBM 3033       (D.P.)       16           63        57.574
!!$!    VAX D-Format   (D.P.)        2          127        34.844
!!$!    VAX G-Format   (D.P.)        2         1023        171.489
!!$!
!!$!                               XINF         EPS        XMININ
!!$!
!!$!    CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
!!$!    Cyber 180/855
!!$!      under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
!!$!    IEEE (IBM/XT,
!!$!      SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
!!$!    IEEE (IBM/XT,
!!$!      SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
!!$!    IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
!!$!    VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
!!$!    VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!!$!
!!$!  Reference:
!!$!
!!$!    W J Cody,
!!$!    "An Overview of Software Development for Special Functions",
!!$!    Lecture Notes in Mathematics, 506,
!!$!    Numerical Analysis Dundee, 1975,
!!$!    G. A. Watson (ed.),
!!$!    Springer Verlag, Berlin, 1976.
!!$!
!!$!    Hart et al,
!!$!    Computer Approximations,
!!$!    Wiley and sons, New York, 1968.
!!$!
!!$!  Author:
!!$!
!!$!    W. J. Cody and L. Stoltz,
!!$!    Applied Mathematics Division,
!!$!    Argonne National Laboratory,
!!$!    Argonne, Illinois, 60439.
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real X, the argument of the function.
!!$!
!!$!    Output, real GAMMA, the value of the function.  The program
!!$!    returns the value XINF for singularities or when overflow would occur.
!!$!    The computation is believed to be free of underflow and overflow.
!!$!
!!$  implicit none
!!$!
!!$  real, parameter, dimension ( 7 ) :: c = (/ &
!!$    -1.910444077728E-03, &
!!$     8.4171387781295E-04, &
!!$    -5.952379913043012E-04, &
!!$     7.93650793500350248E-04, &
!!$    -2.777777777777681622553E-03, &
!!$     8.333333333333333331554247E-02, &
!!$     5.7083835261E-03 /)
!!$  real, parameter :: EPS = 1.19E-07
!!$  real fact
!!$  real gamma
!!$  integer i
!!$  integer n
!!$  real, parameter, dimension ( 8 ) :: p = (/ &
!!$    -1.71618513886549492533811E+00, &
!!$     2.47656508055759199108314E+01, &
!!$    -3.79804256470945635097577E+02, &
!!$     6.29331155312818442661052E+02, &
!!$     8.66966202790413211295064E+02, &
!!$    -3.14512729688483675254357E+04, &
!!$    -3.61444134186911729807069E+04, &
!!$     6.64561438202405440627855E+04 /)
!!$  logical parity
!!$  real, parameter :: PI = &
!!$    3.14159265358979323846264338327950288419716939937510E+00
!!$  real, parameter, dimension ( 8 ) :: q = (/ &
!!$    -3.08402300119738975254353e+01, &
!!$     3.15350626979604161529144e+02, &
!!$    -1.01515636749021914166146e+03, &
!!$    -3.10777167157231109440444e+03, &
!!$     2.25381184209801510330112e+04, &
!!$     4.75584627752788110767815e+03, &
!!$    -1.34659959864969306392456e+05, &
!!$    -1.15132259675553483497211e+05 /)
!!$  real, parameter :: SQRTPI = 0.9189385332046727417803297E+00
!!$  real sum1
!!$  real x
!!$  real, parameter :: XBIG = 35.040E+00
!!$  real xden
!!$  real, parameter :: XINF = 3.4E+38
!!$  real, parameter :: XMININ = 1.18E-38
!!$  real xnum
!!$  real y
!!$  real y1
!!$  real ysq
!!$  real z
!!$!
!!$  parity = .false.
!!$  fact = 1.0E+00
!!$  n = 0
!!$  y = x
!!$!
!!$!  Argument is negative.
!!$!
!!$  if ( y <= 0.0E+00 ) then
!!$
!!$    y = - x
!!$    y1 = aint ( y )
!!$    gamma = y - y1
!!$
!!$    if ( gamma /= 0.0E+00 ) then
!!$
!!$      if ( y1 /= aint ( y1 * 0.5E+00 ) * 2.0E+00 ) then
!!$        parity = .true.
!!$      end if
!!$
!!$      fact = - PI / sin ( PI * gamma )
!!$      y = y + 1.0E+00
!!$
!!$    else
!!$
!!$      gamma = XINF
!!$      return
!!$
!!$    end if
!!$
!!$  end if
!!$!
!!$!  Argument < EPS
!!$!
!!$  if ( y < EPS ) then
!!$
!!$    if (y >= XMININ ) then
!!$      gamma = 1.0E+00 / y
!!$    else
!!$      gamma = XINF
!!$      return
!!$    end if
!!$
!!$  else if ( y < 12.0E+00 ) then
!!$
!!$    y1 = y
!!$!
!!$!  0.0E+00 < argument < 1.0E+00
!!$!
!!$    if ( y < 1.0E+00 ) then
!!$      z = y
!!$      y = y + 1.0E+00
!!$!
!!$!  1.0E+00 < argument < 12.0, reduce argument if necessary.
!!$!
!!$    else
!!$      n = int ( y ) - 1
!!$      y = y - real ( n )
!!$      z = y - 1.0E+00
!!$    end if
!!$!
!!$!  Evaluate approximation for 1.0E+00 < argument < 2.0.
!!$!
!!$    xnum = 0.0E+00
!!$    xden = 1.0E+00
!!$    do i = 1, 8
!!$      xnum = ( xnum + p(i) ) * z
!!$      xden = xden * z + q(i)
!!$    end do
!!$
!!$    gamma = xnum / xden + 1.0E+00
!!$!
!!$!  Adjust result for case  0.0E+00 < argument < 1.0.
!!$!
!!$    if ( y1 < y ) then
!!$      gamma = gamma / y1
!!$!
!!$!  Adjust result for case  2.0E+00 < argument < 12.0.
!!$!
!!$    else if ( y1 > y ) then
!!$
!!$      do i = 1, n
!!$        gamma = gamma * y
!!$        y = y + 1.0E+00
!!$      end do
!!$
!!$    end if
!!$!
!!$!  Evaluate for 12 <= argument.
!!$!
!!$  else
!!$
!!$    if ( y <= XBIG ) then
!!$
!!$      ysq = y**2
!!$      sum1 = c(7)
!!$      do i = 1, 6
!!$        sum1 = sum1 / ysq + c(i)
!!$      end do
!!$      sum1 = sum1 / y - y + SQRTPI
!!$      sum1 = sum1 + ( y - 0.5E+00 ) * log ( y )
!!$      gamma = exp ( sum1 )
!!$
!!$    else
!!$
!!$      gamma = XINF
!!$      return
!!$
!!$    end if
!!$
!!$  end if
!!$!
!!$!  Final adjustments and return.
!!$!
!!$  if ( parity ) then
!!$    gamma = - gamma
!!$  end if
!!$
!!$  if ( fact /= 1.0E+00 ) then
!!$    gamma = fact / gamma
!!$  end if
!!$
!!$  return
!!$end
!!$subroutine gaus8 ( func, a, b, err, result, ierr )
!!$!
!!$!***********************************************************************
!!$!
!!$!! GAUS8 estimates the integral of a function.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    GAUS8 integrates real functions of one variable over finite
!!$!    intervals using an adaptive 8-point Legendre-Gauss
!!$!    algorithm.
!!$!
!!$!    GAUS8 is intended primarily for high accuracy integration or
!!$!    integration of smooth functions.
!!$!
!!$!  Reference:
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!  Author:
!!$!
!!$!    R E Jones,
!!$!    Sandia National Laboratory,
!!$!    Los Alamos, New Mexico
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real, external FUNC, name of external function to be integrated.
!!$!    This name must be in an external statement in the calling program.
!!$!    FUNC must be a function of one real argument.  The value
!!$!    of the argument to FUNC is the variable of integration
!!$!    which ranges from A to B.
!!$!
!!$!    Input, real A, the lower limit of integration.
!!$!
!!$!    Input, real B, the upper limit of integration.
!!$!
!!$!    Input/output, real ERR.
!!$!    On input, ERR is a requested pseudorelative error
!!$!    tolerance.  Normally pick a value of ABS ( ERR ) so that
!!$!    STOL < ABS ( ERR ) <= 1.0E-3 where STOL is the single precision
!!$!    unit roundoff  = R1MACH(4).
!!$!
!!$!    RESULT will normally have no more error than
!!$!    ABS ( ERR ) times the integral of the absolute value of
!!$!    FUN(X).  Usually, smaller values for ERR yield more
!!$!    accuracy and require more function evaluations.
!!$!
!!$!    A negative value for ERR causes an estimate of the
!!$!    absolute error in RESULT to be returned in ERR.  Note that
!!$!    ERR must be a variable (not a constant) in this case.
!!$!    Note also that the user must reset the value of ERR
!!$!    before making any more calls that use the variable ERR.
!!$!
!!$!    On output, ERR will be an estimate of the absolute error
!!$!    in RESULT if the input value of ERR was negative.  ERR is
!!$!    unchanged if the input value of ERR was non-negative.
!!$!
!!$!    The estimated error is solely for information to the user
!!$!    and should not be used as a correction to the computed integral.
!!$!
!!$!    Output, real RESULT, the computed value of the integral.
!!$!
!!$!    Output, integer IERR, a status code.
!!$!
!!$!    Normal Codes:
!!$!
!!$!     1 RESULT most likely meets requested error tolerance, or A = B.
!!$!    -1 A and B are too nearly equal to allow normal
!!$!        integration.  RESULT is set to zero.
!!$!
!!$!     Abnormal Code:
!!$!
!!$!     2 RESULT probably does not meet requested error tolerance.
!!$!
!!$  implicit none
!!$!
!!$  real a
!!$  real aa(30)
!!$  real ae
!!$  real anib
!!$  real area
!!$  real b
!!$  real c
!!$  real ce
!!$  real ee
!!$  real ef
!!$  real eps
!!$  real err
!!$  real est
!!$  real, external :: func
!!$  real g8
!!$  real gl
!!$  real glr
!!$  real gr(30)
!!$  real h
!!$  real hh(30)
!!$  integer i1mach
!!$  integer, save :: icall = 0
!!$  integer ierr
!!$  integer k
!!$  integer, save :: kml = 6
!!$  integer, save :: kmx = 5000
!!$  integer l
!!$  integer lmn
!!$  integer lmx
!!$  integer lr(30)
!!$  integer mxl
!!$  integer nbits
!!$  integer nib
!!$  integer, save :: nlmn = 1
!!$  integer nlmx
!!$  real r1mach
!!$  real result
!!$  real tol
!!$  real vl(30)
!!$  real vr
!!$  real, save :: w1
!!$  real, save :: w2
!!$  real, save :: w3
!!$  real, save :: w4
!!$  real x
!!$  real, save :: x1
!!$  real, save :: x2
!!$  real, save :: x3
!!$  real, save :: x4
!!$!
!!$  data x1, x2, x3, x4/ &
!!$          1.83434642495649805e-01,     5.25532409916328986e-01, &
!!$          7.96666477413626740e-01,     9.60289856497536232e-01/
!!$!
!!$  data w1, w2, w3, w4/ &
!!$          3.62683783378361983e-01,     3.13706645877887287e-01, &
!!$          2.22381034453374471e-01,     1.01228536290376259e-01/
!!$!
!!$!  Warning!  Statement function!
!!$!
!!$  g8(x,h) = h*((w1*(func(x-x1*h) + func(x+x1*h)) &
!!$             +w2*(func(x-x2*h) + func(x+x2*h))) &
!!$            +(w3*(func(x-x3*h) + func(x+x3*h)) &
!!$             +w4*(func(x-x4*h) + func(x+x4*h))))
!!$!
!!$  if ( a == b ) then
!!$    err = 0.0E+00
!!$    result = 0.0E+00
!!$    return
!!$  end if
!!$ 
!!$  if ( icall /= 0 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'GAUS8 - Fatal error!'
!!$    write ( *, '(a)' ) '  GAUS8 was called recursively.'
!!$    stop
!!$  end if
!!$
!!$  icall = 1
!!$  k = i1mach(11)
!!$  anib = r1mach(5) * real(i1mach(11)) / 0.30102000
!!$  nbits = int(anib)
!!$  nlmx = min ( 30, (nbits*5)/8 )
!!$  result = 0.0E+00
!!$  ierr = 1
!!$  ce = 0.0E+00
!!$  result = 0.0E+00
!!$  lmx = nlmx
!!$  lmn = nlmn
!!$ 
!!$  if ( b /= 0.0 ) then
!!$
!!$    if ( sign ( 1.0, b ) * a <= 0.0E+00 ) then
!!$      go to 10
!!$    end if
!!$
!!$    c = abs ( 1.0 - a / b )
!!$    if ( c > 0.1 ) go to 10
!!$    if ( c <= 0.0E+00 ) go to 140
!!$    anib = 0.5 - log(c) / 0.69314718E+00
!!$    nib = int(anib)
!!$    lmx = min(nlmx,nbits-nib-7)
!!$    if ( lmx < 1 ) go to 130
!!$    lmn = min ( lmn, lmx )
!!$
!!$  end if
!!$ 
!!$10    continue
!!$ 
!!$  tol = max ( abs ( err ), 2.0**(5-nbits)) / 2.0E+00
!!$  if ( err == 0.0 ) then
!!$    tol=sqrt( epsilon ( 1.0E+00 ) )
!!$  end if
!!$
!!$  eps = tol
!!$  hh(1) = (b-a) / 4.0E+00
!!$  aa(1) = a
!!$  lr(1) = 1
!!$  l = 1
!!$  est = g8 ( aa(l) + 2.0*hh(l),2.0*hh(l) )
!!$  k = 8
!!$  area = abs ( est )
!!$  ef = 0.5
!!$  mxl = 0
!!$!
!!$!  Compute refined estimates, estimate the error, etc.
!!$!
!!$20 continue
!!$ 
!!$  gl = g8 ( aa(l)+hh(l),hh(l) )
!!$  gr(l) = g8(aa(l)+3.0*hh(l),hh(l))
!!$  k = k + 16
!!$  area = area + ( abs ( gl ) + abs ( gr(l) ) - abs ( est ) )
!!$ 
!!$  glr = gl + gr(l)
!!$  ee = abs ( est - glr ) * ef
!!$  ae = max ( eps * area, tol * abs ( glr ) )
!!$
!!$  if ( ee - ae <= 0.0E+00 ) then
!!$    go to 40
!!$  else
!!$    go to 50
!!$  end if
!!$ 
!!$30 continue
!!$ 
!!$  mxl = 1
!!$ 
!!$40 continue
!!$ 
!!$  ce = ce + (est-glr)
!!$ 
!!$  if ( lr(l) <= 0 ) then
!!$    go to 60
!!$  else
!!$    go to 80
!!$  end if
!!$!
!!$!  Consider the left half of this level
!!$!
!!$50 continue
!!$
!!$  if ( k > kmx ) lmx = kml
!!$  if ( l >= lmx ) go to 30
!!$  l = l + 1
!!$  eps = eps * 0.5E+00
!!$  ef = ef / sqrt ( 2.0E+00 )
!!$  hh(l) = hh(l-1) * 0.5E+00
!!$  lr(l) = -1
!!$  aa(l) = aa(l-1)
!!$  est = gl
!!$  go to 20
!!$!
!!$!  Proceed to right half at this level
!!$!
!!$60 continue
!!$
!!$  vl(l) = glr
!!$
!!$70 continue
!!$
!!$  est = gr(l-1)
!!$  lr(l) = 1
!!$  aa(l) = aa(l) + 4.0E+00 * hh(l)
!!$  go to 20
!!$!
!!$!  Return one level
!!$!
!!$80 continue
!!$
!!$  vr = glr
!!$
!!$90 continue
!!$
!!$  if ( l <= 1 ) go to 120
!!$  l = l - 1
!!$  eps = eps * 2.0E+00
!!$  ef = ef * sqrt ( 2.0E+00 )
!!$ 
!!$  if ( lr(l) <= 0 ) then
!!$    vl(l) = vl(l+1) + vr
!!$    go to 70
!!$  else
!!$    vr = vl(l+1) + vr
!!$    go to 90
!!$  end if
!!$!
!!$!  Exit
!!$!
!!$120   continue
!!$ 
!!$  result = vr
!!$  if ( mxl == 0 .or. abs ( ce ) <= 2.0E+00 * tol * area ) go to 140
!!$  ierr = 2
!!$  write ( *, '(a)' ) ' '
!!$  write ( *, '(a)' ) 'GAUS8 - Warning!'
!!$  write ( *, '(a)' ) '  RESULT is probably insufficiently accurate.'
!!$  icall = 0
!!$
!!$  if ( err < 0.0E+00 ) then
!!$    err = ce
!!$  end if
!!$
!!$  return
!!$ 
!!$130   continue
!!$ 
!!$  ierr = -1
!!$  write ( *, '(a)' ) ' '
!!$  write ( *, '(a)' ) 'GAUS8 - Warning!'
!!$  write ( *, '(a)' ) '  A and B are too close to carry out integration.'
!!$ 
!!$140   continue
!!$ 
!!$  icall = 0
!!$ 
!!$  if ( err < 0.0E+00 ) then
!!$    err = ce
!!$  end if
!!$ 
!!$  return
!!$end
!!$subroutine gausq2 ( n, d, e, z, ierr )
!!$!
!!$!***********************************************************************
!!$!
!!$!! GAUSQ2 finds the eigenvalues of a symmetric tridiagonal matrix.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    GAUSQ2 finds the eigenvalues and first components of the
!!$!    eigenvectors of a symmetric tridiagonal matrix by the implicit QL
!!$!    method.
!!$!
!!$!    GAUSQ2 is a translation of an ALGOL procedure,
!!$!
!!$!      Martin and Wilkinson,
!!$!      Numerische Mathematik,
!!$!      Volume 12, pages 377-383, 1968,
!!$!
!!$!    as modified by
!!$!
!!$!      Dubrulle,
!!$!      Numerische Mathematik,
!!$!      Volume 15, page 450, 1970.
!!$!
!!$!    Handbook for Automatic Computation,
!!$!    vol.ii-linear algebra,
!!$!    pages 241-248, 1971.
!!$!
!!$!    GAUSQ2 is a modified version of the EISPACK routine IMTQL2.
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, integer N, is the order of the matrix.
!!$!
!!$!    Input/output, real D(N).
!!$!    On input, D contains the diagonal elements of the matrix.
!!$!    On output, D contains the eigenvalues in ascending order.
!!$!    If an error exit is made, the eigenvalues are correct but
!!$!    unordered for indices 1, 2, ..., IERR-1;
!!$!
!!$!    Input/output, real E(N).
!!$!    On input, E contains the subdiagonal elements of the input matrix
!!$!    in its first N-1 positions.  E(N) is arbitrary.
!!$!    On output, E has been destroyed.
!!$!
!!$!    Input/output, real Z(N).
!!$!    On input, Z contains the first row of the identity matrix.
!!$!    On output, Z contains the first components of the orthonormal
!!$!    eigenvectors of the symmetric tridiagonal matrix.  If an error exit is
!!$!    made, Z contains the eigenvectors associated with the stored
!!$!    eigenvalues.
!!$!
!!$!    Output, integer IERR.
!!$!    0, for normal return,
!!$!    J, if the j-th eigenvalue has not been determined after 30 iterations.
!!$!
!!$  implicit none
!!$!
!!$  integer n
!!$!
!!$  real b
!!$  real c
!!$  real d(n)
!!$  real e(n)
!!$  real epmach
!!$  real f
!!$  real g
!!$  integer i
!!$  integer ierr
!!$  integer ii
!!$  integer j
!!$  integer k
!!$  integer l
!!$  integer m
!!$  integer mml
!!$  real p
!!$  real r
!!$  real s
!!$  real z(n)
!!$!
!!$  epmach = epsilon ( epmach )
!!$ 
!!$  ierr = 0
!!$
!!$  if ( n == 1 ) then
!!$    return
!!$  end if
!!$ 
!!$  e(n) = 0.0E+00
!!$ 
!!$  do l = 1, n
!!$ 
!!$    j = 0
!!$!
!!$!  Look for a small sub-diagonal element
!!$!
!!$    do m = l, n
!!$
!!$      if ( m == n ) then
!!$        exit
!!$      end if
!!$
!!$      if ( abs ( e(m) ) <= epmach * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
!!$        exit
!!$      end if
!!$
!!$    end do
!!$ 
!!$10  continue
!!$ 
!!$    p = d(l)
!!$    if ( m == l ) go to 20
!!$ 
!!$    if ( j == 30 ) then
!!$      ierr = l
!!$      return
!!$    end if
!!$ 
!!$    j = j + 1
!!$!
!!$!  Form shift
!!$!
!!$    g = (d(l+1) - p) / (2.0E+00 * e(l))
!!$    r = sqrt ( g*g + 1.0E+00 )
!!$    g = d(m) - p + e(l) / (g + sign(r, g))
!!$    s = 1.0E+00
!!$    c = 1.0E+00
!!$    p = 0.0E+00
!!$    mml = m - l
!!$ 
!!$    do ii = 1, mml
!!$ 
!!$      i = m - ii
!!$      f = s * e(i)
!!$      b = c * e(i)
!!$ 
!!$      if ( abs ( f ) >= abs ( g ) ) then
!!$ 
!!$        c = g / f
!!$        r = sqrt ( c * c + 1.0E+00 )
!!$        e(i+1) = f * r
!!$        s = 1.0E+00 / r
!!$        c = c * s
!!$ 
!!$      else
!!$ 
!!$        s = f / g
!!$        r = sqrt ( s*s + 1.0E+00 )
!!$        e(i+1) = g * r
!!$        c = 1.0E+00 / r
!!$        s = s * c
!!$ 
!!$      end if
!!$ 
!!$      g = d(i+1) - p
!!$      r = (d(i) - g) * s + 2.0E+00 * c * b
!!$      p = s * r
!!$      d(i+1) = g + p
!!$      g = c * r - b
!!$!
!!$!  Form the first component of the vector.
!!$!
!!$      f = z(i+1)
!!$      z(i+1) = s * z(i) + c * f
!!$      z(i) = f * z(i) - s * f
!!$    end do
!!$ 
!!$    d(l) = d(l) - p
!!$    e(l) = g
!!$    e(m) = 0.0E+00
!!$    go to 10
!!$ 
!!$20  continue
!!$ 
!!$  end do
!!$!
!!$!  Order the eigenvalues and eigenvectors.
!!$!
!!$  do ii = 2, n
!!$ 
!!$    i = ii - 1
!!$    k = i
!!$    p = d(i)
!!$ 
!!$    do j = ii, n
!!$      if ( d(j) < p ) then
!!$        k = j
!!$        p = d(j)
!!$      end if
!!$    end do
!!$ 
!!$    if ( k /= i ) then
!!$      d(k) = d(i)
!!$      d(i) = p
!!$      call r_swap ( z(i), z(k) )
!!$    end if
!!$ 
!!$  end do
!!$ 
!!$  return
!!$end
!!$subroutine gaussq ( kind, norder, alpha, beta, kpts, endpts, b, xtab, &
!!$  weight )
!!$!
!!$!***********************************************************************
!!$!
!!$!! GAUSSQ computes a Gauss quadrature rule.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    GAUSSQ computes the nodes and weights for Gaussian-type quadrature
!!$!    rules with pre-assigned nodes.
!!$!
!!$!    These are used when one wishes to approximate
!!$!
!!$!      Integral (from A to B)  F(X) W(X) DX
!!$!
!!$!    by
!!$!
!!$!      Sum (J = 1 to NORDER) WEIGHT(I)*F(XTAB(I))
!!$!
!!$!
!!$!    GAUSSQ includes six integration rules that are applicable
!!$!    to this problem, for particular weight functions and particular
!!$!    intervals, including infinite and semi-infinite intervals.
!!$!
!!$!    Associated with each weight function W(X) is a set of
!!$!    orthogonal polynomials.  The nodes XTAB are just the zeroes
!!$!    of the proper NORDER-th degree polynomial.
!!$!
!!$!    GAUSSQ allows the user to modify the rule to require that
!!$!    one or both of the endpoints of the interval are to be
!!$!    included as quadrature nodes.
!!$!
!!$!  References:
!!$!
!!$!    Golub, G H, and Welsch, J H,
!!$!    Calculation of Gaussian Quadrature Rules,
!!$!    Mathematics of Computation
!!$!    Volume 23, April 1969, pages 221-230.
!!$!
!!$!    Golub, G H,
!!$!    Some Modified Matrix Eigenvalue Problems,
!!$!    SIAM Review
!!$!    Volume 15, April 1973, pages 318-334, section 7.
!!$!
!!$!    Stroud and Secrest,
!!$!    Gaussian Quadrature Formulas,
!!$!    Prentice-Hall,
!!$!    Englewood Cliffs, New Jersey, 1966.
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, integer KIND, chooses the quadrature rule to be calculated.
!!$!    1:  Legendre quadrature,
!!$!         W(X) = 1
!!$!         on (-1, 1)
!!$!
!!$!    2:  Chebyshev quadrature of the first kind
!!$!         W(X) = 1/SQRT(1 - X*X)
!!$!         on (-1, +1)
!!$!
!!$!    3:  Chebyshev quadrature of the second kind
!!$!         W(X) = SQRT(1 - X*X)
!!$!         on (-1, 1)
!!$!
!!$!    4:  Hermite quadrature,
!!$!         W(X) = EXP(-X*X)
!!$!         on (-infinity, +infinity)
!!$!
!!$!    5:  Jacobi quadrature,
!!$!         W(X) = (1-X)**ALPHA * (1+X)**BETA
!!$!         on (-1, 1),
!!$!         ALPHA > -1, BETA > -1.
!!$!
!!$!         Note that KIND = 2 and 3 are a special case of this.
!!$!
!!$!    6:  Generalized Laguerre quadrature,
!!$!         W(X) = EXP(-X)*X**ALPHA
!!$!         on (0, +infinity),
!!$!         ALPHA > -1
!!$!
!!$!    Input, integer NORDER, the number of points used for the quadrature rule.
!!$!
!!$!    Input, real ALPHA.
!!$!    ALPHA is only required for Gauss-Jacobi and Gauss-Laguerre
!!$!    quadrature.  Its value is ignored in other cases.
!!$!
!!$!    Input, real BETA.
!!$!    BETA is only required for Gauss-Jacobi quadrature.
!!$!    Its value is ignored in other cases.
!!$!
!!$!    Input, integer KPTS.
!!$!    KPTS is normally zero.
!!$!
!!$!    If KPTS is nonzero, it signals that one or both of the
!!$!    endpoints of the interval is required to be a node.
!!$!    This is called Gauss-Radau or Gauss-Lobatto quadrature.
!!$!    Then KPTS is the number of endpoints that must be
!!$!    included, either 1 or 2.
!!$!
!!$!    Input, real ENDPTS(2).
!!$!    If KPTS is 1 or 2, ENDPTS contains the locations of the
!!$!    endpoints to be fixed.
!!$!
!!$!    Workspace, real B(NORDER).
!!$!
!!$!    Output, real XTAB(NORDER).
!!$!    XTAB contains the nodes for the quadrature rule.
!!$!
!!$!    Output, real WEIGHT(NORDER).
!!$!    WEIGHT contains the weights for the quadrature rule.
!!$!
!!$  implicit none
!!$!
!!$  integer norder
!!$!
!!$  real alpha
!!$  real b(norder)
!!$  real beta
!!$  real endpts(2)
!!$  real gam
!!$  integer i
!!$  integer ierr
!!$  integer kind
!!$  integer kpts
!!$  real muzero
!!$  real solve
!!$  real t1
!!$  real weight(norder)
!!$  real xtab(norder)
!!$!
!!$!  Get the diagonal coefficients XTAB(1:NORDER) and off-diagonal
!!$!  coefficients B(1:NORDER-1) and MUZERO.
!!$!
!!$  call class ( kind, norder, alpha, beta, b, xtab, muzero )
!!$!
!!$!  The matrix of coefficients is assumed to be symmetric.
!!$!  The array XTAB contains the diagonal elements, the array
!!$!  B the off-diagonal elements.
!!$!  Make appropriate changes in the lower right 2 by 2 submatrix.
!!$!
!!$!  If KPTS = 1, only XTAB(NORDER) must be changed.
!!$!
!!$  if ( kpts == 1 ) then
!!$ 
!!$    xtab(norder) = endpts(1) + solve(endpts(1),norder,xtab,b) * b(norder-1)**2
!!$!
!!$!  If KPTS = 2, XTAB(NORDER) and B(NORDER-1) must be recomputed.
!!$!
!!$  else if ( kpts == 2 ) then
!!$ 
!!$    gam = solve ( endpts(1), norder, xtab, b )
!!$    t1 = ((endpts(1) - endpts(2)) / (solve(endpts(2),norder,xtab,b) - gam ) )
!!$    b(norder-1) = sqrt(t1)
!!$    xtab(norder) = endpts(1) + gam * t1
!!$ 
!!$  end if
!!$!
!!$!  The indices of the elements of B run from 1 to NORDER-1.
!!$!  The value of B(NORDER) is of no importance.
!!$!
!!$!  Now compute the eigenvalues of the symmetric tridiagonal
!!$!  matrix, which has been modified as necessary.
!!$!
!!$!  The method used is a QL-type method with origin shifting.
!!$!
!!$  weight(1) = 1.0E+00
!!$  weight(2:norder) = 0.0E+00
!!$ 
!!$  call gausq2 ( norder, xtab, b, weight, ierr )
!!$ 
!!$  do i = 1, norder
!!$    weight(i) = muzero * weight(i)**2
!!$  end do
!!$ 
!!$  return
!!$end
!!$subroutine hiordq ( n, y, delt, work, result )
!!$!
!!$!***********************************************************************
!!$!
!!$!! HIORDQ approximates the integral of a function using equally spaced data.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    The method applies the trapezoidal rule to various subsets of the
!!$!    data, and then applies Richardson extrapolation.
!!$!
!!$!  Author:
!!$!
!!$!    Alan Kaylor Cline,
!!$!    Department of Computer Science,
!!$!    University of Texas at Austin.
!!$!
!!$!  Modified:
!!$!
!!$!    19 February 2002
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, integer N, number of data points.
!!$!
!!$!    Input, real Y(N), the Y values of the data.
!!$!
!!$!    Input, real DELT, the spacing between the X values of the
!!$!    data.  The actual X values are not needed!
!!$!
!!$!    Work array, real WORK(2*(N-1)).  The actual minimum amount
!!$!    of workspace required is two times the number of integer
!!$!    divisors of N-1.
!!$!
!!$!    Output, real RESULT, the approximation to the integral.
!!$!
!!$  implicit none
!!$!
!!$  integer n
!!$!
!!$  real delt
!!$  real fac
!!$  integer i
!!$  integer j
!!$  integer jbak
!!$  integer jj
!!$  integer k
!!$  real result
!!$  real sum2
!!$  real sum1
!!$  real work(2*(n-1))
!!$  real y(n)
!!$!
!!$!  Determine initial trapezoidal rule
!!$!
!!$  sum1 = ( y(1) + y(n) ) / 2.0E+00
!!$  j = -1
!!$ 
!!$  do k = 1, n-1
!!$!
!!$!  Check if K divides N-1
!!$!
!!$    if ( ((n-1)/k)*k == n-1 ) then
!!$!
!!$!  Determine the K-point trapezoidal rule.
!!$!
!!$      sum2 = -sum1
!!$      do i = 1, n, (n-1)/k
!!$        sum2 = sum2 + y(i)
!!$      end do
!!$ 
!!$      j = j + 2
!!$      work(j) = delt * sum2 * real ( ( n - 1 ) / k )
!!$      work(j+1) = real ( ((n-1)/k)**2 )
!!$!
!!$!  Apply Richardson extrapolation.
!!$!
!!$      if ( k /= 1 ) then
!!$ 
!!$        do jj = 3, j, 2
!!$          jbak = j+1-jj
!!$          fac = work(j+1) / ( work(j+1) - work(jbak+1) )
!!$          work(jbak) = work(jbak+2) + fac * ( work(jbak) - work(jbak+2) )
!!$        end do
!!$ 
!!$      end if
!!$ 
!!$    end if
!!$ 
!!$  end do
!!$ 
!!$  result = work(1)
!!$ 
!!$  return
!!$end
!!$function i1mach ( i )
!!$!
!!$!*******************************************************************************
!!$!
!!$!! I1MACH returns integer machine constants.
!!$!
!!$!
!!$!  I/O unit numbers.
!!$!
!!$!    I1MACH(1) = the standard input unit.
!!$!    I1MACH(2) = the standard output unit.
!!$!    I1MACH(3) = the standard punch unit.
!!$!    I1MACH(4) = the standard error message unit.
!!$!
!!$!  Words.
!!$!
!!$!    I1MACH(5) = the number of bits per integer storage unit.
!!$!    I1MACH(6) = the number of characters per integer storage unit.
!!$!
!!$!  Integers.
!!$!
!!$!  Assume integers are represented in the S digit base A form:
!!$!
!!$!  Sign * (X(S-1)*A**(S-1) + ... + X(1)*A + X(0))
!!$!  where 0 <= X(I)<A for I=0 to S-1.
!!$!
!!$!    I1MACH(7) = A, the base.
!!$!    I1MACH(8) = S, the number of base A digits.
!!$!    I1MACH(9) = A**S-1, the largest integer.
!!$!
!!$!  Floating point numbers
!!$!
!!$!  Assume floating point numbers are represented in the T digit base B form:
!!$!
!!$!    Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B**T) )
!!$!
!!$!  where 0 <= X(I)<B for I=1 to T, 0<X(1) and EMIN <= E <= EMAX
!!$!
!!$!    I1MACH(10) = B, the base.
!!$!
!!$!  Single precision
!!$!
!!$!    I1MACH(11) = T, the number of base B digits.
!!$!    I1MACH(12) = EMIN, the smallest exponent E.
!!$!    I1MACH(13) = EMAX, the largest exponent E.
!!$!
!!$!  Double precision
!!$!
!!$!    I1MACH(14) = T, the number of base B digits.
!!$!    I1MACH(15) = EMIN, the smallest exponent E.
!!$!    I1MACH(16) = EMAX, the largest exponent E.
!!$!
!!$!  To alter this function for a particular environment, the desired set of DATA
!!$!  statements should be activated by removing the C from column 1.  On rare
!!$!  machines, a STATIC statement may need to be added, but probably more systems
!!$!  prohibit than require it.
!!$!
!!$!  Also, the values of I1MACH(1) through I1MACH(4) should be checked for
!!$!  consistency with the local operating system.  For FORTRAN 77, you may wish
!!$!  to adjust the data statement so imach(6) is set to 1, and then to comment
!!$!  out the executable test on I.EQ.6 below.
!!$!
!!$!  For IEEE-arithmetic machines (binary standard), the first set of constants
!!$!  below should be appropriate, except perhaps for IMACH(1) - IMACH(4).
!!$!
!!$  integer i
!!$  integer i1mach
!!$  integer imach(16)
!!$  integer output
!!$!
!!$  equivalence (imach(4),output)
!!$!
!!$!  IEEE arithmetic machines, such as the ATT 3B series, Motorola
!!$!  68000 based machines such as the SUN 3 and ATT PC 7300, and
!!$!  8087 based micros such asthe IBM PC and ATT 6300.
!!$!
!!$   data imach( 1) /    5 /
!!$   data imach( 2) /    6 /
!!$   data imach( 3) /    7 /
!!$   data imach( 4) /    6 /
!!$   data imach( 5) /   32 /
!!$   data imach( 6) /    4 /
!!$   data imach( 7) /    2 /
!!$   data imach( 8) /   31 /
!!$   data imach( 9) / 2147483647 /
!!$   data imach(10) /    2 /
!!$   data imach(11) /   24 /
!!$   data imach(12) / -125 /
!!$   data imach(13) /  128 /
!!$   data imach(14) /   53 /
!!$   data imach(15) / -1021 /
!!$   data imach(16) /  1024 /
!!$!
!!$!  ALLIANT FX/8 UNIX FORTRAN compiler.
!!$!
!!$!      data imach( 1) /     5 /
!!$!      data imach( 2) /     6 /
!!$!      data imach( 3) /     6 /
!!$!      data imach( 4) /     0 /
!!$!      data imach( 5) /    32 /
!!$!      data imach( 6) /     4 /
!!$!      data imach( 7) /     2 /
!!$!      data imach( 8) /    32 /
!!$!      data imach( 9) /2147483647/
!!$!      data imach(10) /     2 /
!!$!      data imach(11) /    24 /
!!$!      data imach(12) /  -126 /
!!$!      data imach(13) /   128 /
!!$!      data imach(14) /    53 /
!!$!      data imach(15) / -1022 /
!!$!      data imach(16) /  1024 /
!!$!
!!$!  AMDAHL machines.
!!$!
!!$!      data imach( 1) /   5 /
!!$!      data imach( 2) /   6 /
!!$!      data imach( 3) /   7 /
!!$!      data imach( 4) /   6 /
!!$!      data imach( 5) /  32 /
!!$!      data imach( 6) /   4 /
!!$!      data imach( 7) /   2 /
!!$!      data imach( 8) /  31 /
!!$!      data imach( 9) / 2147483647 /
!!$!      data imach(10) /  16 /
!!$!      data imach(11) /   6 /
!!$!      data imach(12) / -64 /
!!$!      data imach(13) /  63 /
!!$!      data imach(14) /  14 /
!!$!      data imach(15) / -64 /
!!$!      data imach(16) /  63 /
!!$!
!!$!  BURROUGHS 1700 system.
!!$!
!!$!      data imach( 1) /    7 /
!!$!      data imach( 2) /    2 /
!!$!      data imach( 3) /    2 /
!!$!      data imach( 4) /    2 /
!!$!      data imach( 5) /   36 /
!!$!      data imach( 6) /    4 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   33 /
!!$!      data imach( 9) / Z1FFFFFFFF /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   24 /
!!$!      data imach(12) / -256 /
!!$!      data imach(13) /  255 /
!!$!      data imach(14) /   60 /
!!$!      data imach(15) / -256 /
!!$!      data imach(16) /  255 /
!!$!
!!$!  BURROUGHS 5700 system.
!!$!
!!$!      data imach( 1) /   5 /
!!$!      data imach( 2) /   6 /
!!$!      data imach( 3) /   7 /
!!$!      data imach( 4) /   6 /
!!$!      data imach( 5) /  48 /
!!$!      data imach( 6) /   6 /
!!$!      data imach( 7) /   2 /
!!$!      data imach( 8) /  39 /
!!$!      data imach( 9) / O0007777777777777 /
!!$!      data imach(10) /   8 /
!!$!      data imach(11) /  13 /
!!$!      data imach(12) / -50 /
!!$!      data imach(13) /  76 /
!!$!      data imach(14) /  26 /
!!$!      data imach(15) / -50 /
!!$!      data imach(16) /  76 /
!!$!
!!$!  BURROUGHS 6700/7700 systems.
!!$!
!!$!      data imach( 1) /   5 /
!!$!      data imach( 2) /   6 /
!!$!      data imach( 3) /   7 /
!!$!      data imach( 4) /   6 /
!!$!      data imach( 5) /  48 /
!!$!      data imach( 6) /   6 /
!!$!      data imach( 7) /   2 /
!!$!      data imach( 8) /  39 /
!!$!      data imach( 9) / O0007777777777777 /
!!$!      data imach(10) /   8 /
!!$!      data imach(11) /  13 /
!!$!      data imach(12) / -50 /
!!$!      data imach(13) /  76 /
!!$!      data imach(14) /  26 /
!!$!      data imach(15) / -32754 /
!!$!      data imach(16) /  32780 /
!!$!
!!$!  CDC CYBER 170/180 series using NOS
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    7 /
!!$!      data imach( 4) /    6 /
!!$!      data imach( 5) /   60 /
!!$!      data imach( 6) /   10 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   48 /
!!$!      data imach( 9) / O"00007777777777777777" /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   48 /
!!$!      data imach(12) / -974 /
!!$!      data imach(13) / 1070 /
!!$!      data imach(14) /   96 /
!!$!      data imach(15) / -927 /
!!$!      data imach(16) / 1070 /
!!$!
!!$!  CDC CYBER 170/180 series using NOS/VE
!!$!
!!$!      data imach( 1) /     5 /
!!$!      data imach( 2) /     6 /
!!$!      data imach( 3) /     7 /
!!$!      data imach( 4) /     6 /
!!$!      data imach( 5) /    64 /
!!$!      data imach( 6) /     8 /
!!$!      data imach( 7) /     2 /
!!$!      data imach( 8) /    63 /
!!$!      data imach( 9) / 9223372036854775807 /
!!$!      data imach(10) /     2 /
!!$!      data imach(11) /    47 /
!!$!      data imach(12) / -4095 /
!!$!      data imach(13) /  4094 /
!!$!      data imach(14) /    94 /
!!$!      data imach(15) / -4095 /
!!$!      data imach(16) /  4094 /
!!$!
!!$!  CDC CYBER 200 series
!!$!
!!$!      data imach( 1) /      5 /
!!$!      data imach( 2) /      6 /
!!$!      data imach( 3) /      7 /
!!$!      data imach( 4) /      6 /
!!$!      data imach( 5) /     64 /
!!$!      data imach( 6) /      8 /
!!$!      data imach( 7) /      2 /
!!$!      data imach( 8) /     47 /
!!$!      data imach( 9) / X'00007FFFFFFFFFFF' /
!!$!      data imach(10) /      2 /
!!$!      data imach(11) /     47 /
!!$!      data imach(12) / -28625 /
!!$!      data imach(13) /  28718 /
!!$!      data imach(14) /     94 /
!!$!      data imach(15) / -28625 /
!!$!      data imach(16) /  28718 /
!!$!
!!$!  CDC 6000/7000 series using FTN4.
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    7 /
!!$!      data imach( 4) /    6 /
!!$!      data imach( 5) /   60 /
!!$!      data imach( 6) /   10 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   48 /
!!$!      data imach( 9) / 00007777777777777777B /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   47 /
!!$!      data imach(12) / -929 /
!!$!      data imach(13) / 1070 /
!!$!      data imach(14) /   94 /
!!$!      data imach(15) / -929 /
!!$!      data imach(16) / 1069 /
!!$!
!!$!  CDC 6000/7000 series using FTN5.
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    7 /
!!$!      data imach( 4) /    6 /
!!$!      data imach( 5) /   60 /
!!$!      data imach( 6) /   10 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   48 /
!!$!      data imach( 9) / O"00007777777777777777" /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   47 /
!!$!      data imach(12) / -929 /
!!$!      data imach(13) / 1070 /
!!$!      data imach(14) /   94 /
!!$!      data imach(15) / -929 /
!!$!      data imach(16) / 1069 /
!!$!
!!$!  CONVEX C-1.
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    7 /
!!$!      data imach( 4) /    6 /
!!$!      data imach( 5) /   32 /
!!$!      data imach( 6) /    4 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   31 /
!!$!      data imach( 9) / 2147483647 /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   24 /
!!$!      data imach(12) / -128 /
!!$!      data imach(13) /  127 /
!!$!      data imach(14) /   53 /
!!$!      data imach(15) /-1024 /
!!$!      data imach(16) / 1023 /
!!$!
!!$!  CONVEX C-120 (native mode) without -R8 option
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    0 /
!!$!      data imach( 4) /    6 /
!!$!      data imach( 5) /   32 /
!!$!      data imach( 6) /    4 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   31 /
!!$!      data imach( 9) / 2147483647 /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   24 /
!!$!      data imach(12) / -127 /
!!$!      data imach(13) /  127 /
!!$!      data imach(14) /   53 /
!!$!      data imach(15) / -1023 /
!!$!      data imach(16) /  1023 /
!!$!
!!$!  CONVEX C-120 (native mode) with -R8 option
!!$!
!!$!      data imach( 1) /     5 /
!!$!      data imach( 2) /     6 /
!!$!      data imach( 3) /     0 /
!!$!      data imach( 4) /     6 /
!!$!      data imach( 5) /    32 /
!!$!      data imach( 6) /     4 /
!!$!      data imach( 7) /     2 /
!!$!      data imach( 8) /    31 /
!!$!      data imach( 9) / 2147483647 /
!!$!      data imach(10) /     2 /
!!$!      data imach(11) /    53 /
!!$!      data imach(12) / -1023 /
!!$!      data imach(13) /  1023 /
!!$!      data imach(14) /    53 /
!!$!      data imach(15) / -1023 /
!!$!      data imach(16) /  1023 /
!!$!
!!$!  CONVEX C-120 (IEEE mode) without -R8 option
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    0 /
!!$!      data imach( 4) /    6 /
!!$!      data imach( 5) /   32 /
!!$!      data imach( 6) /    4 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   31 /
!!$!      data imach( 9) / 2147483647 /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   24 /
!!$!      data imach(12) / -125 /
!!$!      data imach(13) /  128 /
!!$!      data imach(14) /   53 /
!!$!      data imach(15) / -1021 /
!!$!      data imach(16) /  1024 /
!!$!
!!$!  CONVEX C-120 (IEEE mode) with -R8 option
!!$!
!!$!      data imach( 1) /     5 /
!!$!      data imach( 2) /     6 /
!!$!      data imach( 3) /     0 /
!!$!      data imach( 4) /     6 /
!!$!      data imach( 5) /    32 /
!!$!      data imach( 6) /     4 /
!!$!      data imach( 7) /     2 /
!!$!      data imach( 8) /    31 /
!!$!      data imach( 9) / 2147483647 /
!!$!      data imach(10) /     2 /
!!$!      data imach(11) /    53 /
!!$!      data imach(12) / -1021 /
!!$!      data imach(13) /  1024 /
!!$!      data imach(14) /    53 /
!!$!      data imach(15) / -1021 /
!!$!      data imach(16) /  1024 /
!!$!
!!$!  CRAY 1, 2, XMP and YMP.
!!$!
!!$!      data imach( 1) /     5 /
!!$!      data imach( 2) /     6 /
!!$!      data imach( 3) /   102 /
!!$!      data imach( 4) /     6 /
!!$!      data imach( 5) /    64 /
!!$!      data imach( 6) /     8 /
!!$!      data imach( 7) /     2 /
!!$!      data imach( 8) /    63 /
!!$!      data imach( 9) /  777777777777777777777B /
!!$!      data imach(10) /     2 /
!!$!      data imach(11) /    47 /
!!$!      data imach(12) / -8189 /
!!$!      data imach(13) /  8190 /
!!$!      data imach(14) /    94 /
!!$!      data imach(15) / -8099 /
!!$!      data imach(16) /  8190 /
!!$!
!!$!  DATA GENERAL ECLIPSE S/200.
!!$!
!!$!      data imach( 1) /   11 /
!!$!      data imach( 2) /   12 /
!!$!      data imach( 3) /    8 /
!!$!      data imach( 4) /   10 /
!!$!      data imach( 5) /   16 /
!!$!      data imach( 6) /    2 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   15 /
!!$!      data imach( 9) /32767 /
!!$!      data imach(10) /   16 /
!!$!      data imach(11) /    6 /
!!$!      data imach(12) /  -64 /
!!$!      data imach(13) /   63 /
!!$!      data imach(14) /   14 /
!!$!      data imach(15) /  -64 /
!!$!      data imach(16) /   63 /
!!$!
!!$!  ELXSI 6400
!!$!
!!$!      data imach( 1) /     5 /
!!$!      data imach( 2) /     6 /
!!$!      data imach( 3) /     6 /
!!$!      data imach( 4) /     6 /
!!$!      data imach( 5) /    32 /
!!$!      data imach( 6) /     4 /
!!$!      data imach( 7) /     2 /
!!$!      data imach( 8) /    32 /
!!$!      data imach( 9) / 2147483647 /
!!$!      data imach(10) /     2 /
!!$!      data imach(11) /    24 /
!!$!      data imach(12) /  -126 /
!!$!      data imach(13) /   127 /
!!$!      data imach(14) /    53 /
!!$!      data imach(15) / -1022 /
!!$!      data imach(16) /  1023 /
!!$!
!!$!  HARRIS 220
!!$!
!!$!      data imach( 1) /       5 /
!!$!      data imach( 2) /       6 /
!!$!      data imach( 3) /       0 /
!!$!      data imach( 4) /       6 /
!!$!      data imach( 5) /      24 /
!!$!      data imach( 6) /       3 /
!!$!      data imach( 7) /       2 /
!!$!      data imach( 8) /      23 /
!!$!      data imach( 9) / 8388607 /
!!$!      data imach(10) /       2 /
!!$!      data imach(11) /      23 /
!!$!      data imach(12) /    -127 /
!!$!      data imach(13) /     127 /
!!$!      data imach(14) /      38 /
!!$!      data imach(15) /    -127 /
!!$!      data imach(16) /     127 /
!!$!
!!$!  HARRIS SLASH 6 and SLASH 7.
!!$!
!!$!      data imach( 1) /       5 /
!!$!      data imach( 2) /       6 /
!!$!      data imach( 3) /       0 /
!!$!      data imach( 4) /       6 /
!!$!      data imach( 5) /      24 /
!!$!      data imach( 6) /       3 /
!!$!      data imach( 7) /       2 /
!!$!      data imach( 8) /      23 /
!!$!      data imach( 9) / 8388607 /
!!$!      data imach(10) /       2 /
!!$!      data imach(11) /      23 /
!!$!      data imach(12) /    -127 /
!!$!      data imach(13) /     127 /
!!$!      data imach(14) /      38 /
!!$!      data imach(15) /    -127 /
!!$!      data imach(16) /     127 /
!!$!
!!$!  HONEYWELL DPS 8/70 and 600/6000 series.
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /   43 /
!!$!      data imach( 4) /    6 /
!!$!      data imach( 5) /   36 /
!!$!      data imach( 6) /    4 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   35 /
!!$!      data imach( 9) / O377777777777 /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   27 /
!!$!      data imach(12) / -127 /
!!$!      data imach(13) /  127 /
!!$!      data imach(14) /   63 /
!!$!      data imach(15) / -127 /
!!$!      data imach(16) /  127 /
!!$!
!!$!  HP 2100, 3 word double precision option with FTN4
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    4 /
!!$!      data imach( 4) /    1 /
!!$!      data imach( 5) /   16 /
!!$!      data imach( 6) /    2 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   15 /
!!$!      data imach( 9) / 32767 /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   23 /
!!$!      data imach(12) / -128 /
!!$!      data imach(13) /  127 /
!!$!      data imach(14) /   39 /
!!$!      data imach(15) / -128 /
!!$!      data imach(16) /  127 /
!!$!
!!$!  HP 2100, 4 word double precision option with FTN4
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    4 /
!!$!      data imach( 4) /    1 /
!!$!      data imach( 5) /   16 /
!!$!      data imach( 6) /    2 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   15 /
!!$!      data imach( 9) / 32767 /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   23 /
!!$!      data imach(12) / -128 /
!!$!      data imach(13) /  127 /
!!$!      data imach(14) /   55 /
!!$!      data imach(15) / -128 /
!!$!      data imach(16) /  127 /
!!$!
!!$!  HP 9000
!!$!
!!$!      data imach( 1) /     5 /
!!$!      data imach( 2) /     6 /
!!$!      data imach( 3) /     6 /
!!$!      data imach( 4) /     7 /
!!$!      data imach( 5) /    32 /
!!$!      data imach( 6) /     4 /
!!$!      data imach( 7) /     2 /
!!$!      data imach( 8) /    32 /
!!$!      data imach( 9) / 2147483647 /
!!$!      data imach(10) /     2 /
!!$!      data imach(11) /    24 /
!!$!      data imach(12) /  -126 /
!!$!      data imach(13) /   127 /
!!$!      data imach(14) /    53 /
!!$!      data imach(15) / -1015 /
!!$!      data imach(16) /  1017 /
!!$!
!!$!  IBM 360/370 series, XEROX SIGMA 5/7/9, SEL systems 85/86, PERKIN ELMER 3230,
!!$!  and PERKIN ELMER (INTERDATA) 3230.
!!$!
!!$!      data imach( 1) /   5 /
!!$!      data imach( 2) /   6 /
!!$!      data imach( 3) /   7 /
!!$!      data imach( 4) /   6 /
!!$!      data imach( 5) /  32 /
!!$!      data imach( 6) /   4 /
!!$!      data imach( 7) /   2 /
!!$!      data imach( 8) /  31 /
!!$!      data imach( 9) / Z7FFFFFFF /
!!$!      data imach(10) /  16 /
!!$!      data imach(11) /   6 /
!!$!      data imach(12) / -64 /
!!$!      data imach(13) /  63 /
!!$!      data imach(14) /  14 /
!!$!      data imach(15) / -64 /
!!$!      data imach(16) /  63 /
!!$!
!!$!  IBM PC - Microsoft FORTRAN
!!$!
!!$!      data imach( 1) /     5 /
!!$!      data imach( 2) /     6 /
!!$!      data imach( 3) /     6 /
!!$!      data imach( 4) /     0 /
!!$!      data imach( 5) /    32 /
!!$!      data imach( 6) /     4 /
!!$!      data imach( 7) /     2 /
!!$!      data imach( 8) /    31 /
!!$!      data imach( 9) / 2147483647 /
!!$!      data imach(10) /     2 /
!!$!      data imach(11) /    24 /
!!$!      data imach(12) /  -126 /
!!$!      data imach(13) /   127 /
!!$!      data imach(14) /    53 /
!!$!      data imach(15) / -1022 /
!!$!      data imach(16) /  1023 /
!!$!
!!$!  IBM PC - Professional FORTRAN and Lahey FORTRAN
!!$!
!!$!      data imach( 1) /     4 /
!!$!      data imach( 2) /     7 /
!!$!      data imach( 3) /     7 /
!!$!      data imach( 4) /     0 /
!!$!      data imach( 5) /    32 /
!!$!      data imach( 6) /     4 /
!!$!      data imach( 7) /     2 /
!!$!      data imach( 8) /    31 /
!!$!      data imach( 9) / 2147483647 /
!!$!      data imach(10) /     2 /
!!$!      data imach(11) /    24 /
!!$!      data imach(12) /  -126 /
!!$!      data imach(13) /   127 /
!!$!      data imach(14) /    53 /
!!$!      data imach(15) / -1022 /
!!$!      data imach(16) /  1023 /
!!$!
!!$!  INTERDATA 8/32 with the UNIX system FORTRAN 77 compiler.
!!$!  For the INTERDATA FORTRAN VII compiler, replace the Z's specifying hex
!!$!  constants with Y's.
!!$!
!!$!      data imach( 1) /   5 /
!!$!      data imach( 2) /   6 /
!!$!      data imach( 3) /   6 /
!!$!      data imach( 4) /   6 /
!!$!      data imach( 5) /  32 /
!!$!      data imach( 6) /   4 /
!!$!      data imach( 7) /   2 /
!!$!      data imach( 8) /  31 /
!!$!      data imach( 9) / Z'7FFFFFFF' /
!!$!      data imach(10) /  16 /
!!$!      data imach(11) /   6 /
!!$!      data imach(12) / -64 /
!!$!      data imach(13) /  62 /
!!$!      data imach(14) /  14 /
!!$!      data imach(15) / -64 /
!!$!      data imach(16) /  62 /
!!$!
!!$!  PDP-10 (KA processor).
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    7 /
!!$!      data imach( 4) /    6 /
!!$!      data imach( 5) /   36 /
!!$!      data imach( 6) /    5 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   35 /
!!$!      data imach( 9) / "377777777777 /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   27 /
!!$!      data imach(12) / -128 /
!!$!      data imach(13) /  127 /
!!$!      data imach(14) /   54 /
!!$!      data imach(15) / -101 /
!!$!      data imach(16) /  127 /
!!$!
!!$!  PDP-10 (KI processor).
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    7 /
!!$!      data imach( 4) /    6 /
!!$!      data imach( 5) /   36 /
!!$!      data imach( 6) /    5 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   35 /
!!$!      data imach( 9) / "377777777777 /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   27 /
!!$!      data imach(12) / -128 /
!!$!      data imach(13) /  127 /
!!$!      data imach(14) /   62 /
!!$!      data imach(15) / -128 /
!!$!      data imach(16) /  127 /
!!$!
!!$!  PDP-11 FORTRANS supporting 32-bit integer arithmetic.
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    7 /
!!$!      data imach( 4) /    6 /
!!$!      data imach( 5) /   32 /
!!$!      data imach( 6) /    4 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   31 /
!!$!      data imach( 9) / 2147483647 /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   24 /
!!$!      data imach(12) / -127 /
!!$!      data imach(13) /  127 /
!!$!      data imach(14) /   56 /
!!$!      data imach(15) / -127 /
!!$!      data imach(16) /  127 /
!!$!
!!$!  PDP-11 FORTRANS supporting 16-bit integer arithmetic.
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    7 /
!!$!      data imach( 4) /    6 /
!!$!      data imach( 5) /   16 /
!!$!      data imach( 6) /    2 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   15 /
!!$!      data imach( 9) / 32767 /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   24 /
!!$!      data imach(12) / -127 /
!!$!      data imach(13) /  127 /
!!$!      data imach(14) /   56 /
!!$!      data imach(15) / -127 /
!!$!      data imach(16) /  127 /
!!$!
!!$!  PRIME 50 series systems with 32-bit integers and 64V MODE instructions,
!!$!  supplied by Igor Bray.
!!$!
!!$!      data imach( 1) /            1 /
!!$!      data imach( 2) /            1 /
!!$!      data imach( 3) /            2 /
!!$!      data imach( 4) /            1 /
!!$!      data imach( 5) /           32 /
!!$!      data imach( 6) /            4 /
!!$!      data imach( 7) /            2 /
!!$!      data imach( 8) /           31 /
!!$!      data imach( 9) / :17777777777 /
!!$!      data imach(10) /            2 /
!!$!      data imach(11) /           23 /
!!$!      data imach(12) /         -127 /
!!$!      data imach(13) /         +127 /
!!$!      data imach(14) /           47 /
!!$!      data imach(15) /       -32895 /
!!$!      data imach(16) /       +32637 /
!!$!
!!$!  SEQUENT BALANCE 8000.
!!$!
!!$!      data imach( 1) /     0 /
!!$!      data imach( 2) /     0 /
!!$!      data imach( 3) /     7 /
!!$!      data imach( 4) /     0 /
!!$!      data imach( 5) /    32 /
!!$!      data imach( 6) /     1 /
!!$!      data imach( 7) /     2 /
!!$!      data imach( 8) /    31 /
!!$!      data imach( 9) /  2147483647 /
!!$!      data imach(10) /     2 /
!!$!      data imach(11) /    24 /
!!$!      data imach(12) /  -125 /
!!$!      data imach(13) /   128 /
!!$!      data imach(14) /    53 /
!!$!      data imach(15) / -1021 /
!!$!      data imach(16) /  1024 /
!!$!
!!$!  SUN Microsystems UNIX F77 compiler.
!!$!
!!$!      data imach( 1) /     5 /
!!$!      data imach( 2) /     6 /
!!$!      data imach( 3) /     6 /
!!$!      data imach( 4) /     0 /
!!$!      data imach( 5) /    32 /
!!$!      data imach( 6) /     4 /
!!$!      data imach( 7) /     2 /
!!$!      data imach( 8) /    32 /
!!$!      data imach( 9) /2147483647/
!!$!      data imach(10) /     2 /
!!$!      data imach(11) /    24 /
!!$!      data imach(12) /  -126 /
!!$!      data imach(13) /   128 /
!!$!      data imach(14) /    53 /
!!$!      data imach(15) / -1022 /
!!$!      data imach(16) /  1024 /
!!$!
!!$!  SUN 3 (68881 or FPA)
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    6 /
!!$!      data imach( 4) /    0 /
!!$!      data imach( 5) /   32 /
!!$!      data imach( 6) /    4 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   31 /
!!$!      data imach( 9) / 2147483647 /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   24 /
!!$!      data imach(12) / -125 /
!!$!      data imach(13) /  128 /
!!$!      data imach(14) /   53 /
!!$!      data imach(15) / -1021 /
!!$!      data imach(16) /  1024 /
!!$!
!!$!  UNIVAC 1100 series.
!!$!  Note that the punch unit, I1MACH(3), has been set to 7, which is appropriate
!!$!  for the UNIVAC-FOR system.  If you have the UNIVAC-FTN system, set it to 1
!!$!  instead.
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    7 /
!!$!      data imach( 4) /    6 /
!!$!      data imach( 5) /   36 /
!!$!      data imach( 6) /    6 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   35 /
!!$!      data imach( 9) / O377777777777 /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   27 /
!!$!      data imach(12) / -128 /
!!$!      data imach(13) /  127 /
!!$!      data imach(14) /   60 /
!!$!      data imach(15) /-1024 /
!!$!      data imach(16) / 1023 /
!!$!
!!$!  VAX.
!!$!
!!$!      data imach( 1) /    5 /
!!$!      data imach( 2) /    6 /
!!$!      data imach( 3) /    7 /
!!$!      data imach( 4) /    6 /
!!$!      data imach( 5) /   32 /
!!$!      data imach( 6) /    4 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   31 /
!!$!      data imach( 9) / 2147483647 /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   24 /
!!$!      data imach(12) / -127 /
!!$!      data imach(13) /  127 /
!!$!      data imach(14) /   56 /
!!$!      data imach(15) / -127 /
!!$!      data imach(16) /  127 /
!!$!
!!$!  Z80 microprocessor.
!!$!
!!$!      data imach( 1) /    1 /
!!$!      data imach( 2) /    1 /
!!$!      data imach( 3) /    0 /
!!$!      data imach( 4) /    1 /
!!$!      data imach( 5) /   16 /
!!$!      data imach( 6) /    2 /
!!$!      data imach( 7) /    2 /
!!$!      data imach( 8) /   15 /
!!$!      data imach( 9) / 32767 /
!!$!      data imach(10) /    2 /
!!$!      data imach(11) /   24 /
!!$!      data imach(12) / -127 /
!!$!      data imach(13) /  127 /
!!$!      data imach(14) /   56 /
!!$!      data imach(15) / -127 /
!!$!      data imach(16) /  127 /
!!$!
!!$  if ( i < 1 .or. i > 16 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'I1MACH - Fatal error!'
!!$    write ( *, '(a,i6)' ) 'I is out of bounds:',i
!!$    i1mach = 0
!!$    stop
!!$  else
!!$    i1mach = imach(i)
!!$  end if
!!$
!!$  return
!!$end
!!$subroutine iratex ( func, a, b, epsin, epsout, result, ind )
!!$!
!!$!***********************************************************************
!!$!
!!$!! IRATEX estimates the integral of a function.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    IRATEX estimates the integral from A to B of F(X), using the
!!$!    trapezoidal rule for several stepsizes H.
!!$!
!!$!    Then a rational function of H*H is fitted to these results, and 
!!$!    an estimate of the integral is computed by extrapolating the 
!!$!    results to a stepsize of zero.
!!$!
!!$!  Reference:
!!$!
!!$!    R Bulirsch and J Stoer,
!!$!    Fehlerabschaetzungen und Extrapolation mit rationaled Funktionen
!!$!      bei Verfahren vom Richardson-Typus,
!!$!    (Error estimates and extrapolation with rational functions
!!$!      in processes of Richardson type),
!!$!    Numerische Mathematik,
!!$!    Volume 6 (1964), pages 413-427.
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real, external FUNC.
!!$!    FUNC is the name of the function to be
!!$!    integrated.  The user must declare the name an external
!!$!    parameter in the calling program, pass the name of the
!!$!    function in FUNC, and write a function of the form
!!$!
!!$!      FUNCTION FUNC(X)
!!$!
!!$!    which evaluates the function at the point X.
!!$!
!!$!    Input, real A, the lower limit of integration.
!!$!
!!$!    Input, real B, the upper limit of integration.
!!$!
!!$!    Input, real EPSIN, the requested relative error tolerance.
!!$!
!!$!    Output, real EPSOUT, an estimate of the error in the integration.
!!$!
!!$!    Output, real RESULT, the approximate value of the integral.
!!$!
!!$!    Output, integer IND, error return flag.
!!$!    IND = 0 if the requested accuracy was not achieved,
!!$!    IND = 1 if the accuracy was achieved.
!!$!
!!$  implicit none
!!$!
!!$  real a
!!$  real arg
!!$  real b
!!$  real ba
!!$  logical bo
!!$  logical bu
!!$  real c
!!$  real d(6)
!!$  real d1
!!$  real ddt
!!$  real den
!!$  real dt(7)
!!$  real e
!!$  real ent
!!$  real epsin
!!$  real epsout
!!$  real, external :: func
!!$  real gr
!!$  real hm
!!$  integer i
!!$  integer ind
!!$  integer m
!!$  integer mr
!!$  integer n
!!$  integer np1
!!$  logical odd
!!$  real rnderr
!!$  real result
!!$  real sm
!!$  real t
!!$  real t1
!!$  real t2
!!$  real t2a
!!$  real ta
!!$  real tab
!!$  real tb
!!$  real tnt
!!$  real v
!!$  real w
!!$!
!!$  if ( a == b ) then
!!$    result = 0.0E+00
!!$    return
!!$  end if
!!$ 
!!$  rnderr = epsilon ( 1.0E+00 )
!!$  epsin = max ( epsin, 8.0 * rnderr )
!!$  ind = 0
!!$  n = 2
!!$  np1 = 3
!!$  ba = b-a
!!$  t1 = 0.0E+00
!!$  gr = 0.0E+00
!!$  sm = 0.0E+00
!!$  t2a = 0.5E+00 * ( func ( a ) + func ( b ) )
!!$  t2 = t2a
!!$  tb = abs ( t2a )
!!$  c = t2 * ba
!!$
!!$  dt(1) = c
!!$  dt(2:7) = 0.0E+00
!!$ 
!!$  odd = .true.
!!$  bu = .false.
!!$ 
!!$  do m = 1, 15
!!$ 
!!$    bo = m>=7
!!$    hm = ba / real(n)
!!$!
!!$!  N+1 is odd.
!!$!
!!$    if ( odd ) then
!!$ 
!!$      do i = 1, n, 2
!!$        arg = a + real(i) * hm
!!$        w = func(arg)
!!$        t2 = t2 + w
!!$        tb = abs ( w ) + tb
!!$      end do
!!$ 
!!$      ent = t2
!!$      tab = tb * abs ( hm )
!!$      d(1) = 16.0E+00 / 9.0E+00
!!$      d(3) = 4.0E+00 * d(1)
!!$      d(5) = 4.0E+00 * d(3)
!!$!
!!$!  N+1 is even.
!!$!
!!$    else
!!$ 
!!$      do i = 1, n, 6
!!$        w = real(i) * hm
!!$        t1 = t1 + func(a+w)+func(b-w)
!!$      end do
!!$ 
!!$      ent = t1+t2a
!!$      t2a = t2
!!$      d(1) = 2.25E+00
!!$      d(3) = 9.0E+00
!!$      d(5) = 36.0E+00
!!$ 
!!$    end if
!!$ 
!!$    ddt = dt(1)
!!$    t = ent*hm
!!$    dt(1) = t
!!$    ent = t
!!$ 
!!$    if ( m < 7 ) then
!!$      mr = m
!!$      w = real ( n * n )
!!$      d(m) = w
!!$    else
!!$      mr = 6
!!$      d(6) = 64.0E+00
!!$      w = 144.0E+00
!!$    end if
!!$ 
!!$    do i = 1, mr
!!$ 
!!$      d1 = d(i) * ddt
!!$      den = d1 - ent
!!$      e = ent - ddt
!!$      tnt = ent
!!$      v = 0.0E+00
!!$      ent = 0.0E+00
!!$ 
!!$      if ( abs ( den ) >= epsin ) then
!!$        e = e / den
!!$        v = tnt * e
!!$        ent = d1 * e
!!$        t = v + t
!!$        ddt = dt(i+1)
!!$      end if
!!$ 
!!$      dt(i+1) = v
!!$ 
!!$    end do
!!$ 
!!$    ta = c
!!$    c = t
!!$    result = c
!!$    if ( .not. bo ) then
!!$      t = t-v
!!$    end if
!!$
!!$    v = t-ta
!!$    t = v+t
!!$    epsout = abs ( v )
!!$ 
!!$    if ( ta < t ) then
!!$      d1 = ta
!!$      ta = t
!!$      t = d1
!!$    end if
!!$ 
!!$    bo = bo .or. ( ta < gr .and. t > sm )
!!$ 
!!$    if ( bu .and. bo .and. epsout < epsin * w * tab ) then
!!$      go to 140
!!$    end if
!!$ 
!!$    gr = ta + epsin
!!$    sm = t - epsin
!!$    odd = .not. odd
!!$    n = np1
!!$    np1 = n+1
!!$    bu = bo
!!$    d(2) = 4.0E+00
!!$    d(4) = 16.0E+00
!!$ 
!!$  end do
!!$ 
!!$  bo = .false.
!!$ 
!!$  140 continue
!!$ 
!!$  v = rnderr*tab
!!$ 
!!$  epsout = max ( epsout, v )
!!$
!!$  if (bo) ind = 1
!!$
!!$  return
!!$end
!!$function pi ( )
!!$!
!!$!*******************************************************************************
!!$!
!!$!! PI returns the value of pi.
!!$!
!!$!
!!$!  Modified:
!!$!
!!$!    04 December 1998
!!$!
!!$!  Author:
!!$!
!!$!    John Burkardt
!!$!
!!$!  Parameters:
!!$!
!!$!    Output, real PI, the value of pi.
!!$!
!!$  implicit none
!!$!
!!$  real pi
!!$!
!!$  pi = 3.14159265358979323846264338327950288419716939937510E+00
!!$
!!$  return
!!$end
!!$subroutine plint ( ftab, xtab, ntab, a, b, result )
!!$!
!!$!***********************************************************************
!!$!
!!$!! PLINT approximates the integral of unequally spaced data.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    The method uses piecewise linear interpolation.
!!$!
!!$!  Reference:
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real FTAB(NTAB), the function values, FTAB(I) = F(XTAB(I)).
!!$!
!!$!    Input, real XTAB(NTAB), the abscissas at which the
!!$!    function values are given.  The XTAB's must be distinct
!!$!    and in ascending order.
!!$!
!!$!    Input, integer NTAB, the number of entries in FTAB and
!!$!    XTAB.  NTAB must be at least 2.
!!$!
!!$!    Input, real A, the lower limit of integration.  A should
!!$!    be, but need not be, near one endpoint of the interval
!!$!    (X(1), X(NTAB)).
!!$!
!!$!    Input, real B, the upper limit of integration.  B should
!!$!    be, but need not be, near one endpoint of the interval
!!$!    (X(1), X(NTAB)).
!!$!
!!$!    Output, real RESULT, the approximate value of the integral.
!!$!
!!$  implicit none
!!$!
!!$  integer ntab
!!$!
!!$  real a
!!$  real b
!!$  real fa
!!$  real fb
!!$  real ftab(ntab)
!!$  integer i
!!$  integer ihi
!!$  integer ilo
!!$  integer ind
!!$  real result
!!$  real slope
!!$  real syl
!!$  real xtab(ntab)
!!$!
!!$  if ( a == b ) then
!!$    result = 0.0E+00
!!$    return
!!$  end if
!!$ 
!!$  if ( ntab < 2 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'PLINT - Fatal error!'
!!$    write ( *, '(a,i6)' ) '  NTAB < 2, NTAB = ',ntab
!!$    stop
!!$  end if
!!$ 
!!$  do i = 2, ntab
!!$    if ( xtab(i) <= xtab(i-1) ) then
!!$      write ( *, '(a)' ) ' '
!!$      write ( *, '(a)' ) 'PLINT - Fatal error!'
!!$      write ( *, '(a)' ) '  Nodes not in strict increasing order.'
!!$      write ( *, '(a,i6)' ) '  XTAB(I) <= XTAB(I-1) for I=',i
!!$      write ( *, '(a,g14.6)' ) '  XTAB(I) = ',xtab(i)
!!$      write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ',xtab(i-1)
!!$      stop
!!$    end if
!!$  end do
!!$!
!!$!  If A > B, temporarily switch A and B, and store sign.
!!$!
!!$  if ( a > b ) then
!!$    syl = b
!!$    b = a
!!$    a = syl
!!$    ind = -1
!!$  else
!!$    syl = a
!!$    ind = 1
!!$  end if
!!$!
!!$!  Find ILO and IHI so that A <= XTAB(ILO) <= XTAB(IHI) <= B
!!$!  with the possible exception that A and B may be in the same
!!$!  interval, or completely to the right or left of the XTAB's.
!!$!
!!$  ilo = 1
!!$  ihi = ntab
!!$  do i = 1, ntab
!!$    if ( a <= xtab(i) ) then
!!$      exit
!!$    end if
!!$    ilo = ilo+1
!!$  end do
!!$  
!!$  do i = 1, ntab
!!$    if ( b >= xtab(i) ) then
!!$      exit
!!$    end if
!!$    ihi = ihi-1
!!$  end do
!!$!
!!$!  Treat special cases where A, B lie both to left or both to right
!!$!  of XTAB interval, or inbetween same pair of XTAB's.
!!$!
!!$  if ( ihi == 0 ) then
!!$    slope = (ftab(2)-ftab(1))/(xtab(2)-xtab(1))
!!$    fa = ftab(1) + slope*(a-xtab(1))
!!$    fb = ftab(1) + slope*(b-xtab(1))
!!$    result = 0.5 * (b-a) * (fa+fb)
!!$    go to 110
!!$  else if ( ilo == ntab+1 ) then
!!$    slope = (ftab(ntab)-ftab(ntab-1))/(xtab(ntab)-xtab(ntab-1))
!!$    fa = ftab(ntab-1)+slope*(a-xtab(ntab-1))
!!$    fb = ftab(ntab-1)+slope*(b-xtab(ntab-1))
!!$    result = 0.5 * (b-a) * (fa+fb)
!!$    go to 110
!!$  else if ( ihi+1 == ilo ) then
!!$    slope = (ftab(ilo)-ftab(ihi))/(xtab(ilo)-xtab(ihi))
!!$    fa = ftab(ihi)+slope*(a-xtab(ihi))
!!$    fb = ftab(ihi)+slope*(b-xtab(ihi))
!!$    result = 0.5 * (b-a) * (fa+fb)
!!$    go to 110
!!$  end if
!!$!
!!$!  Carry out approximate integration.  We know that ILO is no greater
!!$!  than IHI-1, but equality is possible; A and B may be on either side
!!$!  of a single XTAB(I).  That's OK, then the loop below won't be executed
!!$!  at all.
!!$!
!!$  result = 0.0E+00
!!$  do i = ilo, ihi-1
!!$    result = result + 0.5 * (xtab(i+1)-xtab(i))*(ftab(i)+ftab(i+1))
!!$  end do
!!$!
!!$!  Add contribution from A-ILO and IHI-B.
!!$!  Still have to watch out if ILO = 1 or IHI=NTAB...
!!$!
!!$  if ( ilo == 1 ) then
!!$    slope = (ftab(2)-ftab(1)) / (xtab(2)-xtab(1))
!!$    fa = ftab(1) + slope*(a-xtab(1))
!!$    result = result + 0.5 * (xtab(ilo)-a)*(fa+ftab(ilo))
!!$  else
!!$    slope = (ftab(ilo)-ftab(ilo-1)) / (xtab(ilo)-xtab(ilo-1))
!!$    fa = ftab(ilo-1) + slope*(a-xtab(ilo-1))
!!$    result = result + 0.5 * (xtab(ilo)-a)*(fa+ftab(ilo))
!!$  end if
!!$ 
!!$  if ( ihi == ntab ) then
!!$    slope = (ftab(ntab)-ftab(ntab-1)) / (xtab(ntab)-xtab(ntab-1))
!!$    fb = ftab(ntab-1) + slope*(b-xtab(ntab-1))
!!$    result = result + 0.5*(b-xtab(ntab))*(fb+ftab(ntab))
!!$  else
!!$    slope = (ftab(ihi+1)-ftab(ihi)) / (xtab(ihi+1)-xtab(ihi))
!!$    fb = ftab(ihi) + slope*(b-xtab(ihi))
!!$    result = result + 0.5*(b-xtab(ihi))*(fb+ftab(ihi))
!!$  end if
!!$!
!!$!  Restore original values of A and B, reverse sign of integral
!!$!  because of earlier switch.
!!$!
!!$110   continue
!!$ 
!!$  if ( ind /= 1 ) then
!!$    ind = 1
!!$    syl = b
!!$    b = a
!!$    a = syl
!!$    result = -result
!!$  end if
!!$ 
!!$  return
!!$end
!!$subroutine qnc79 ( func, a, b, err, result, ierr, k )
!!$!
!!$!***********************************************************************
!!$!
!!$!! QNC79 approximates the integral of F(X) using Newton-Cotes quadrature.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    QNC79 is a general purpose program for evaluation of one
!!$!    dimensional integrals  of user defined functions.  QNC79 will
!!$!    pick its own points for evaluation of the integrand and these
!!$!    will vary from problem to problem.
!!$!
!!$!    Thus QNC79 is not designed to integrate over data sets.
!!$!
!!$!    Moderately smooth integrands will be integrated efficiently
!!$!    and reliably.  For problems with strong singularities,
!!$!    oscillations etc., the user may wish to use more sophisticated
!!$!    routines such as those in QUADPACK.
!!$!
!!$!    One measure of the reliability of QNC79 is the output parameter
!!$!    K, giving the number of integrand evaluations that were needed.
!!$!
!!$!  Reference:
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real, external FUNC, the name of the function to be
!!$!    integrated.  The user must declare the name an external
!!$!    parameter in the calling program, pass the name of the
!!$!    function in FUNC, and write a function of the form
!!$!
!!$!      FUNCTION FUNC(X)
!!$!
!!$!    which evaluates the function at the point X.
!!$!
!!$!    Input, real A, lower limit of integral.
!!$!
!!$!    Input, real B, upper limit of integral.
!!$!
!!$!    Input, real ERR, is a requested error tolerance.
!!$!    Normally pick a value, 0 .LT. ERR .LT. 1.E-3.
!!$!
!!$!    Output, real RESULT, computed value of the integral.
!!$!    Hopefully RESULT is accurate to within ERR times the
!!$!    integral of ABS ( FUNC(X) ).
!!$!
!!$!    Output, integer IERR, a status code
!!$!     1 RESULT most likely meets requested error tolerance.
!!$!    -1 A and B are too nearly equal to allow normal integration.
!!$!     2 RESULT probably does not meet requested error tolerance.
!!$!
!!$!    Output, integer K, the number of function evaluations
!!$!    actually used to do the integration.
!!$!    A value of K .GT. 1000 indicates a difficult problem.
!!$!    Other programs may be more efficient.
!!$!    QNC79 will gracefully give up if K exceeds 2000.
!!$!
!!$  implicit none
!!$!
!!$  integer, parameter :: kmx = 2000
!!$!
!!$  real a
!!$  real aa(40)
!!$  real ae
!!$  real area
!!$  real b
!!$  real bank
!!$  real blocal
!!$  real c
!!$  real ce
!!$  real ee
!!$  real ef
!!$  real eps
!!$  real err
!!$  real f(13)
!!$  real f1(40)
!!$  real f2(40)
!!$  real f3(40)
!!$  real f4(40)
!!$  real f5(40)
!!$  real f6(40)
!!$  real f7(40)
!!$  real, external :: func
!!$  real hh(40)
!!$  integer i
!!$  integer i1mach
!!$  integer, save :: icall = 0
!!$  integer ierr
!!$  integer k
!!$  integer, save :: kml = 7
!!$  integer l
!!$  integer lmn
!!$  integer lmx
!!$  integer lr(40)
!!$  integer nbits
!!$  integer nib
!!$  integer, save :: nlmn = 2
!!$  integer nlmx
!!$  real q13
!!$  real q7
!!$  real q7l
!!$  real q7r(40)
!!$  real r1mach
!!$  real result
!!$  real test
!!$  real tol
!!$  real vl(40)
!!$  real vr
!!$  real w1
!!$  real w2
!!$  real w3
!!$  real w4
!!$!
!!$  if ( a == b ) then
!!$    result = 0.0E+00
!!$    return
!!$  end if
!!$ 
!!$  if ( icall /= 0 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'QNC79 - Fatal error!'
!!$    write ( *, '(a)' ) '  QNC79 was called recursively!'
!!$    stop
!!$  end if
!!$ 
!!$  icall = 1
!!$  w1 = 41.0E+00 / 140.0E+00
!!$  w2  = 216.0 / 140.0E+00
!!$  w3 = 27.0 / 140.0E+00
!!$  w4  = 272.0 / 140.0E+00
!!$  nbits = int ( r1mach(5) * real(i1mach(11)) / 0.30102000 )
!!$  nlmx = min ( 40, (nbits*4)/5 )
!!$  result = 0.0E+00
!!$  ierr = 1
!!$  ce = 0.0E+00
!!$  lmx = nlmx
!!$  lmn = nlmn
!!$  if ( b == 0.0E+00 ) go to 3
!!$  if ( sign ( 1.0E+00, b ) * a <= 0.0E+00 ) go to 3
!!$  c = abs ( 1.0E+00 - a / b )
!!$
!!$  if ( c > 0.1E+00 ) then
!!$    go to 3
!!$  end if
!!$ 
!!$  if ( c <= 0.0E+00 ) then
!!$    ierr  = -1
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'QNC79 - Fatal error!'
!!$    write ( *, '(a)' ) '  A and B are too close.'
!!$    stop
!!$  end if
!!$ 
!!$  nib = int ( 0.5E+00 - log(c) / log(2.0E+00) )
!!$  lmx = min ( nlmx, nbits-nib-4)
!!$
!!$  if ( lmx < 2 ) then
!!$    go to 32
!!$  end if
!!$
!!$  lmn = min(lmn,lmx)
!!$ 
!!$3 continue
!!$ 
!!$  tol = max ( abs ( err ), 2.0E+00**(5-nbits) )
!!$  if ( err == 0.0E+00 ) tol = sqrt ( epsilon ( tol ) )
!!$  eps = tol
!!$  hh(1) = (b-a) / 12.0E+00
!!$  aa(1) = a
!!$  lr(1) = 1
!!$ 
!!$  do i = 1, 11, 2
!!$    f(i) = func(a+real(i-1)*hh(1))
!!$  end do
!!$ 
!!$  blocal = b
!!$  f(13) = func(blocal)
!!$  k = 7
!!$  l = 1
!!$  area = 0.0E+00
!!$  q7 = 0.0E+00
!!$  ef = 256.0E+00 / 255.0E+00
!!$  bank = 0.0E+00
!!$!
!!$!  Compute refined estimates, estimate the error, etc.
!!$!
!!$5 continue
!!$ 
!!$  do i = 2, 12, 2
!!$    f(i) = func ( aa(l)+real(i-1)*hh(l) )
!!$  end do
!!$ 
!!$  k = k+6
!!$!
!!$!  Compute left and right half estimates.
!!$!
!!$  q7l = hh(l)*( ( w1*(f(1)+f(7)) + w2*(f(2)+f(6)) ) &
!!$              + ( w3*(f(3)+f(5)) + w4*f(4) ) )
!!$
!!$  q7r(l) = hh(l)*( ( w1*(f(7)+f(13)) + w2*(f(8)+f(12)) ) &
!!$                + ( w3*(f(9)+f(11)) + w4*f(10) ) )
!!$!
!!$!  Update estimate of integral of absolute value.
!!$!
!!$  area = area + ( abs ( q7l ) + abs ( q7r(l) ) - abs ( q7 ) )
!!$!
!!$!  Do not bother to test convergence before minimum refinement level.
!!$!
!!$  if (l<lmn) go to 11
!!$!
!!$!  Estimate the error in new value for whole interval, Q13.
!!$!
!!$  q13 = q7l+q7r(l)
!!$  ee = abs ( q7 - q13 ) * ef
!!$!
!!$!  Compute nominal allowed error.
!!$!
!!$  ae = eps*area
!!$!
!!$!  Borrow from bank account, but not too much.
!!$!
!!$  test = min(ae+0.8E+00*bank, 10.0E+00*ae)
!!$!
!!$!  Don't ask for excessive accuracy.
!!$!
!!$  test = max ( test, tol * abs ( q13 ), 0.00003E+00 * tol * area )
!!$!
!!$!  Now, did this interval pass or not?
!!$!
!!$  if ( ee <= test )go to 8
!!$  go to 10
!!$!
!!$!  Have hit max refinement level - penalize the cumulative error.
!!$!
!!$    7 continue
!!$ 
!!$  ce = ce+(q7-q13)
!!$  go to 9
!!$!
!!$!  On good intervals accumulate the theoretical estimate.
!!$!
!!$    8 ce = ce+(q7-q13) / 255.0E+00
!!$!
!!$!  Update the bank account.  Don't go into debt.
!!$!
!!$9 continue
!!$
!!$  bank = bank + (ae-ee)
!!$  if ( bank < 0.0E+00 ) bank = 0.0E+00
!!$!
!!$!  Did we just finish a left half or a right half?
!!$!
!!$  if ( lr(l) <= 0 ) go to 15
!!$  go to 20
!!$!
!!$!  Consider the left half of the next deeper level.
!!$!
!!$10 continue
!!$ 
!!$  if ( k > kmx ) lmx = min(kml,lmx)
!!$  if ( l >= lmx ) go to 7
!!$
!!$11 continue
!!$
!!$  l = l+1
!!$  eps = eps * 0.5E+00
!!$  if ( l <= 17 ) ef = ef/sqrt(2.0E+00)
!!$  hh(l) = hh(l-1)*0.5E+00
!!$  lr(l) = -1
!!$  aa(l) = aa(l-1)
!!$  q7 = q7l
!!$  f1(l) = f(7)
!!$  f2(l) = f(8)
!!$  f3(l) = f(9)
!!$  f4(l) = f(10)
!!$  f5(l) = f(11)
!!$  f6(l) = f(12)
!!$  f7(l) = f(13)
!!$  f(13) = f(7)
!!$  f(11) = f(6)
!!$  f(9)  = f(5)
!!$  f(7)  = f(4)
!!$  f(5)  = f(3)
!!$  f(3)  = f(2)
!!$  go to 5
!!$!
!!$!  Proceed to right half at this level
!!$!
!!$   15 vl(l) = q13
!!$   16 q7 = q7r(l-1)
!!$  lr(l) = 1
!!$  aa(l) = aa(l)+12.0E+00 * hh(l)
!!$  f(1)  = f1(l)
!!$  f(3)  = f2(l)
!!$  f(5)  = f3(l)
!!$  f(7)  = f4(l)
!!$  f(9)  = f5(l)
!!$  f(11) = f6(l)
!!$  f(13) = f7(l)
!!$  go to 5
!!$!
!!$!  Left and right halves are done, so go back up a level
!!$!
!!$   20 vr = q13
!!$   22 if ( l <= 1 ) go to 30
!!$
!!$  if ( l <= 17 ) then
!!$    ef = ef * sqrt ( 2.0E+00 )
!!$  end if
!!$
!!$  eps = eps * 2.0E+00
!!$  l = l-1
!!$ 
!!$  if ( lr(l) <= 0 ) then
!!$    vl(l) = vl(l+1)+vr
!!$    go to 16
!!$  else
!!$    vr = vl(l+1)+vr
!!$    go to 22
!!$  end if
!!$ 
!!$   30 continue
!!$ 
!!$  if ( abs ( ce ) > 2.0E+00 * tol * area ) then
!!$    ierr = 2
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'QNC79 - Warning!'
!!$    write ( *, '(a)' ) '  RESULT is probably insufficiently accurate.'
!!$  end if
!!$ 
!!$32    continue
!!$ 
!!$  result = vr
!!$ 
!!$  icall = 0
!!$ 
!!$  return
!!$end
!!$subroutine quad ( func, a, b, abserr, relerr, nleast, nmost, work, &
!!$  result )
!!$!
!!$!***********************************************************************
!!$!
!!$!! QUAD approximates the integral of F(X) by Romberg integration.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    The integration is repeated until convergence of the results, 
!!$!    or the maximum number of steps is taken.  The Romberg method 
!!$!    uses the trapezoidal rule, subinterval halving, and Richardson 
!!$!    extrapolation.
!!$!
!!$!    Convergence is declared if either of the following occurs:
!!$!
!!$!      ABS ( RESULT - INTEGRAL ) < ABSERR
!!$!
!!$!    or
!!$!
!!$!      RESULT = integral from A to B (1+Y(X))*FUNC(X)DX for some
!!$!      function Y with ABS ( Y(X) ) < RELERR+RNDERR  with RNDERR the
!!$!      machine rounding error.
!!$!
!!$!  Reference:
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!    C F Dunkl,
!!$!    Romberg quadrature to prescribed accuracy,
!!$!    SHARE file number 7090-1481 TYQUAD
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real, external FUNC, the name of the function to be integrated.
!!$!    The user must declare the name an external parameter in the
!!$!    calling program, pass the name of the function in FUNC,
!!$!    and write a function of the form FUNCTION FUNC(X) which
!!$!    evaluates the function at the point X.
!!$!
!!$!    Input, real A, the lower limit of integration.
!!$!
!!$!    Input, real B, the upper limit of integration.
!!$!
!!$!    Input, real ABSERR, the absolute error tolerance.
!!$!
!!$!    Input, real RELERR, the relative error tolerance.
!!$!
!!$!    Input, integer NLEAST, the least number of times the integration
!!$!    is to be carried out before the convergence test is made.
!!$!    A value 3 <= NLEAST <= 15 is suggested.
!!$!
!!$!    Input, integer NMOST, the most number of times the
!!$!    integration is to be carried out.
!!$!
!!$!    Workspace, real WORK(NMOST+1).
!!$!
!!$!    Output, real RESULT, the approximate value of the integral.
!!$!
!!$  implicit none
!!$!
!!$  integer nmost
!!$!
!!$  real a
!!$  real abserr
!!$  real b
!!$  real f
!!$  real fcna
!!$  real fcnb
!!$  real fcnxi
!!$  real, external :: func
!!$  real h
!!$  integer i
!!$  integer j
!!$  integer k
!!$  integer nleast
!!$  integer nx
!!$  real qx1
!!$  real qx2
!!$  real relerr
!!$  real result
!!$  real rnderr
!!$  real sum1
!!$  real sumabs
!!$  real t
!!$  real tabs
!!$  real work(nmost+1)
!!$  real x
!!$!
!!$  if ( a == b ) then
!!$    result = 0.0E+00
!!$    return
!!$  end if
!!$ 
!!$  rnderr = epsilon ( 1.0E+00 )
!!$  qx1 = 0
!!$  qx2 = 0
!!$  h = b-a
!!$  fcna = func(a)
!!$  fcnb = func(b)
!!$  tabs = abs ( h ) * ( abs ( fcna ) + abs ( fcnb ) ) / 2.0E+00
!!$  t = h*(fcna+fcnb) / 2.0E+00
!!$  nx = 1
!!$ 
!!$  do i = 1, nmost
!!$ 
!!$    h = 0.5*h
!!$ 
!!$    sum1 = 0.0E+00
!!$    sumabs = 0.0E+00
!!$    do j = 1, nx
!!$      fcnxi = func(a+h*(2*j-1))
!!$      sumabs = sumabs + abs ( fcnxi )
!!$      sum1 = sum1+fcnxi
!!$    end do
!!$ 
!!$    t = 0.5E+00 * t + h * sum1
!!$    tabs = tabs/2.0E+00 + abs ( h ) * sumabs
!!$    work(i) = 2.0E+00*(t+h*sum1)/3.0E+00
!!$
!!$    if ( i > 1 ) then
!!$!
!!$!  Construct difference table for Richardson extrapolation.
!!$!
!!$      f = 4.0E+00
!!$      do j = 2, i
!!$        k = i+1-j
!!$        f = f*4.0E+00
!!$        work(k) = work(k+1)+(work(k+1)-work(k) ) / ( f - 1.0E+00 )
!!$      end do
!!$!
!!$!  Perform acceptance check.
!!$!
!!$      if ( i >= nleast ) then
!!$        x = abs ( work(1)-qx2 ) + abs ( qx2 - qx1 )
!!$        if ( x <= 3.0E+00 * tabs * ( abs ( relerr ) + rnderr ) ) go to 10
!!$        if ( x <= 3.0E+00 * abs ( abserr ) ) go to 10
!!$      end if
!!$!
!!$!  Save old result, perform bisection, repeat.
!!$!
!!$      qx1 = qx2
!!$    end if
!!$ 
!!$    qx2 = work(1)
!!$    nx = nx*2
!!$ 
!!$  end do
!!$!
!!$!  Accept result.
!!$!
!!$10    continue
!!$  result = work(1)
!!$  return
!!$end
!!$subroutine r_swap ( x, y )
!!$!
!!$!*******************************************************************************
!!$!
!!$!! R_SWAP swaps two real values.
!!$!
!!$!
!!$!  Modified:
!!$!
!!$!    01 May 2000
!!$!
!!$!  Author:
!!$!
!!$!    John Burkardt
!!$!
!!$!  Parameters:
!!$!
!!$!    Input/output, real X, Y.  On output, the values of X and
!!$!    Y have been interchanged.
!!$!
!!$  implicit none
!!$!
!!$  real x
!!$  real y
!!$  real z
!!$!
!!$  z = x
!!$  x = y
!!$  y = z
!!$
!!$  return
!!$end
!!$function r1mach(i)
!!$!
!!$!*******************************************************************************
!!$!
!!$!! R1MACH returns single precision machine constants.
!!$!
!!$!
!!$!  Assume that single precision numbers are stored with a mantissa of T digits
!!$!  in base B, with an exponent whose value must lie between EMIN and EMAX.  Then
!!$!  for values of I between 1 and 5, R1MACH will return the following values:
!!$!
!!$!    R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!!$!    R1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
!!$!    R1MACH(3) = B**(-T), the smallest relative spacing.
!!$!    R1MACH(4) = B**(1-T), the largest relative spacing.
!!$!    R1MACH(5) = log10(B)
!!$!
!!$!  To alter this function for a particular environment, the desired set of data
!!$!  statements should be activated by removing the C from column 1.
!!$!
!!$!  On rare machines a STATIC statement may need to be added.  But probably more
!!$!  systems prohibit it that require it.
!!$!
!!$!  For IEEE-arithmetic machines (binary standard), the first set of constants
!!$!  below should be appropriate.
!!$!
!!$!  Where possible, octal or hexadecimal constants have been used to specify the
!!$!  constants exactly which has in some cases required the use of EQUIVALENCED
!!$!  integer arrays.  If your compiler uses half-word integers by default
!!$!  (sometimes called integer*2), you may need to change INTEGER to INTEGER*4 or
!!$!  otherwise instruct your compiler to use full-word integers in the next 5
!!$!  declarations.
!!$!
!!$  integer diver(2)
!!$  integer i
!!$  integer large(2)
!!$  integer log10(2)
!!$  real r1mach
!!$  integer right(2)
!!$  real rmach(5)
!!$  integer small(2)
!!$!
!!$  equivalence (rmach(1),small(1))
!!$  equivalence (rmach(2),large(1))
!!$  equivalence (rmach(3),right(1))
!!$  equivalence (rmach(4),diver(1))
!!$  equivalence (rmach(5),log10(1))
!!$!
!!$!  IEEE arithmetic machines, such as the ATT 3B series, Motorola 68000 based
!!$!  machines such as the SUN 3 and ATT PC 7300, and 8087 based micros such as
!!$!  the IBM PC and ATT 6300.
!!$!
!!$   data small(1) /     8388608 /
!!$   data large(1) /  2139095039 /
!!$   data right(1) /   864026624 /
!!$   data diver(1) /   872415232 /
!!$   data log10(1) /  1050288283 /
!!$!
!!$!  ALLIANT FX/8 UNIX Fortran compiler with the -r8 command line option.  This
!!$!  option causes all variables declared with 'real' to be of type 'REAL*8' or
!!$!  DOUBLE PRECISION.  This option does not override the 'real*4' declarations.
!!$!  These R1MACH numbers below and the coresponding I1MACH are simply the DOUBLE
!!$!  PRECISION or 'real*8' numbers.  If you use the -r8 your whole code (and the
!!$!  user libraries you link with, the system libraries are taken care of
!!$!  automagicly) must be compiled with this option.
!!$!
!!$!      data rmach(1) / 2.22507385850721D-308 /
!!$!      data rmach(2) / 1.79769313486231D+308 /
!!$!      data rmach(3) / 1.1101827117665D-16 /
!!$!      data rmach(4) / 2.2203654423533D-16 /
!!$!      data rmach(5) / 3.01029995663981E-1 /
!!$!
!!$!  AMDAHL machines.
!!$!
!!$!      data small(1) /    1048576 /
!!$!      data large(1) / 2147483647 /
!!$!      data right(1) /  990904320 /
!!$!      data diver(1) / 1007681536 /
!!$!      data log10(1) / 1091781651 /
!!$!
!!$!  BURROUGHS 1700 system.
!!$!
!!$!      data rmach(1) / Z400800000 /
!!$!      data rmach(2) / Z5FFFFFFFF /
!!$!      data rmach(3) / Z4E9800000 /
!!$!      data rmach(4) / Z4EA800000 /
!!$!      data rmach(5) / Z500E730E8 /
!!$!
!!$!  BURROUGHS 5700/6700/7700 systems.
!!$!
!!$!      data rmach(1) / O1771000000000000 /
!!$!      data rmach(2) / O0777777777777777 /
!!$!      data rmach(3) / O1311000000000000 /
!!$!      data rmach(4) / O1301000000000000 /
!!$!      data rmach(5) / O1157163034761675 /
!!$!
!!$!  CDC CYBER 170/180 series using NOS
!!$!
!!$!      data rmach(1) / O"00014000000000000000" /
!!$!      data rmach(2) / O"37767777777777777777" /
!!$!      data rmach(3) / O"16404000000000000000" /
!!$!      data rmach(4) / O"16414000000000000000" /
!!$!      data rmach(5) / O"17164642023241175720" /
!!$!
!!$!  CDC CYBER 170/180 series using NOS/VE
!!$!
!!$!      data rmach(1) / Z"3001800000000000" /
!!$!      data rmach(2) / Z"4FFEFFFFFFFFFFFE" /
!!$!      data rmach(3) / Z"3FD2800000000000" /
!!$!      data rmach(4) / Z"3FD3800000000000" /
!!$!      data rmach(5) / Z"3FFF9A209A84FBCF" /
!!$!
!!$!  CDC CYBER 200 series
!!$!
!!$!      data rmach(1) / X'9000400000000000' /
!!$!      data rmach(2) / X'6FFF7FFFFFFFFFFF' /
!!$!      data rmach(3) / X'FFA3400000000000' /
!!$!      data rmach(4) / X'FFA4400000000000' /
!!$!      data rmach(5) / X'FFD04D104D427DE8' /
!!$!
!!$!  CDC 6000/7000 series using FTN4.
!!$!
!!$!      data rmach(1) / 00564000000000000000B /
!!$!      data rmach(2) / 37767777777777777776B /
!!$!      data rmach(3) / 16414000000000000000B /
!!$!      data rmach(4) / 16424000000000000000B /
!!$!      data rmach(5) / 17164642023241175720B /
!!$!
!!$!  CDC 6000/7000 series using FTN5.
!!$!
!!$!      data rmach(1) / O"00564000000000000000" /
!!$!      data rmach(2) / O"37767777777777777776" /
!!$!      data rmach(3) / O"16414000000000000000" /
!!$!      data rmach(4) / O"16424000000000000000" /
!!$!      data rmach(5) / O"17164642023241175720" /
!!$!
!!$!  CONVEX C-1.
!!$!
!!$!      data rmach(1) / '00800000'X /
!!$!      data rmach(2) / '7FFFFFFF'X /
!!$!      data rmach(3) / '34800000'X /
!!$!      data rmach(4) / '35000000'X /
!!$!      data rmach(5) / '3F9A209B'X /
!!$!
!!$!  CONVEX C-120 (native mode) without -R8 option
!!$!
!!$!      data rmach(1) / 2.9387360E-39 /
!!$!      data rmach(2) / 1.7014117E+38 /
!!$!      data rmach(3) / 5.9604645E-08 /
!!$!      data rmach(4) / 1.1920929E-07 /
!!$!      data rmach(5) / 3.0102999E-01 /
!!$!
!!$!  CONVEX C-120 (native mode) with -R8 option
!!$!
!!$!      data rmach(1) / 5.562684646268007D-309 /
!!$!      data rmach(2) / 8.988465674311577D+307 /
!!$!      data rmach(3) / 1.110223024625157D-016 /
!!$!      data rmach(4) / 2.220446049250313D-016 /
!!$!      data rmach(5) / 3.010299956639812D-001 /
!!$!
!!$!  CONVEX C-120 (IEEE mode) without -R8 option
!!$!
!!$!      data rmach(1) / 1.1754945E-38 /
!!$!      data rmach(2) / 3.4028234E+38 /
!!$!      data rmach(3) / 5.9604645E-08 /
!!$!      data rmach(4) / 1.1920929E-07 /
!!$!      data rmach(5) / 3.0102999E-01 /
!!$!
!!$!  CONVEX C-120 (IEEE mode) with -R8 option
!!$!
!!$!      data rmach(1) / 2.225073858507202D-308 /
!!$!      data rmach(2) / 1.797693134862315D+308 /
!!$!      data rmach(3) / 1.110223024625157D-016 /
!!$!      data rmach(4) / 2.220446049250313D-016 /
!!$!      data rmach(5) / 3.010299956639812D-001 /
!!$!
!!$!  CRAY 1, 2, XMP and YMP.
!!$!
!!$!      data rmach(1) / 200034000000000000000B /
!!$!      data rmach(2) / 577767777777777777776B /
!!$!      data rmach(3) / 377224000000000000000B /
!!$!      data rmach(4) / 377234000000000000000B /
!!$!      data rmach(5) / 377774642023241175720B /
!!$!
!!$!  DATA GENERAL ECLIPSE S/200.
!!$!  Note - It may be appropriate to include the line: STATIC RMACH(5)
!!$!
!!$!      data small /20K,0/
!!$!      data large /77777K,177777K/
!!$!      data right /35420K,0/
!!$!      data diver /36020K,0/
!!$!      data log10 /40423K,42023K/
!!$!
!!$!  ELXSI 6400, assuming real*4 is the default real type.
!!$!
!!$!      data small(1) / '00800000'X /
!!$!      data large(1) / '7F7FFFFF'X /
!!$!      data right(1) / '33800000'X /
!!$!      data diver(1) / '34000000'X /
!!$!      data log10(1) / '3E9A209B'X /
!!$!
!!$!  HARRIS 220
!!$!
!!$!      data small(1),small(2) / '20000000, '00000201 /
!!$!      data large(1),large(2) / '37777777, '00000177 /
!!$!      data right(1),right(2) / '20000000, '00000352 /
!!$!      data diver(1),diver(2) / '20000000, '00000353 /
!!$!      data log10(1),log10(2) / '23210115, '00000377 /
!!$!
!!$!  HARRIS SLASH 6 and SLASH 7.
!!$!
!!$!      data small(1),small(2) / '20000000, '00000201 /
!!$!      data large(1),large(2) / '37777777, '00000177 /
!!$!      data right(1),right(2) / '20000000, '00000352 /
!!$!      data diver(1),diver(2) / '20000000, '00000353 /
!!$!      data log10(1),log10(2) / '23210115, '00000377 /
!!$!
!!$!  HONEYWELL DPS 8/70 and 600/6000 series.
!!$!
!!$!      data rmach(1) / O402400000000 /
!!$!      data rmach(2) / O376777777777 /
!!$!      data rmach(3) / O714400000000 /
!!$!      data rmach(4) / O716400000000 /
!!$!      data rmach(5) / O776464202324 /
!!$!
!!$!  HP 2100, 3 word double precision with FTN4
!!$!
!!$!      data small(1), small(2) / 40000B,       1 /
!!$!      data large(1), large(2) / 77777B, 177776B /
!!$!      data right(1), right(2) / 40000B,    325B /
!!$!      data diver(1), diver(2) / 40000B,    327B /
!!$!      data log10(1), log10(2) / 46420B,  46777B /
!!$!
!!$!  HP 2100, 4 word double precision with FTN4
!!$!
!!$!      data small(1), small(2) / 40000B,       1 /
!!$!      data large91), large(2) / 77777B, 177776B /
!!$!      data right(1), right(2) / 40000B,    325B /
!!$!      data diver(1), diver(2) / 40000B,    327B /
!!$!      data log10(1), log10(2) / 46420B,  46777B /
!!$!
!!$!  HP 9000
!!$!
!!$!      r1mach(1) = 1.17549435E-38
!!$!      r1mach(2) = 1.70141163E+38
!!$!      r1mach(3) = 5.960464478E-8
!!$!      r1mach(4) = 1.119209290E-7
!!$!      r1mach(5) = 3.01030010E-1
!!$!
!!$!      data small(1) / 00040000000B /
!!$!      data large(1) / 17677777777B /
!!$!      data right(1) / 06340000000B /
!!$!      data diver(1) / 06400000000B /
!!$!      data log10(1) / 07646420233B /
!!$!
!!$!  IBM 360/370 series, XEROX SIGMA 5/7/9, SEL systems 85/86, PERKIN ELMER 3230,
!!$!  and PERKIN ELMER (INTERDATA) 3230.
!!$!
!!$!      data rmach(1) / Z00100000 /
!!$!      data rmach(2) / Z7FFFFFFF /
!!$!      data rmach(3) / Z3B100000 /
!!$!      data rmach(4) / Z3C100000 /
!!$!      data rmach(5) / Z41134413 /
!!$!
!!$!  IBM PC - Microsoft FORTRAN
!!$!
!!$!      data small(1) / #00800000 /
!!$!      data large(1) / #7F7FFFFF /
!!$!      data right(1) / #33800000 /
!!$!      data diver(1) / #34000000 /
!!$!      data log10(1) / #3E9A209A /
!!$!
!!$!  IBM PC - Professional FORTRAN and Lahey FORTRAN
!!$!
!!$!      data small(1)/ Z'00800000' /
!!$!      data large(1)/ Z'7F7FFFFF' /
!!$!      data right(1)/ Z'33800000' /
!!$!      data diver(1)/ Z'34000000' /
!!$!      data log10(1)/ Z'3E9A209A' /
!!$!
!!$!  INTERDATA 8/32 with the UNIX system FORTRAN 77 compiler.
!!$!  For the INTERDATA FORTRAN VII compiler replace the Z'S specifying HEX
!!$!  constants with Y'S.
!!$!
!!$!      data rmach(1) / Z'00100000' /
!!$!      data rmach(2) / Z'7EFFFFFF' /
!!$!      data rmach(3) / Z'3B100000' /
!!$!      data rmach(4) / Z'3C100000' /
!!$!      data rmach(5) / Z'41134413' /
!!$!
!!$!  PDP-10 (KA or KI processor).
!!$!
!!$!      data rmach(1) / "000400000000 /
!!$!      data rmach(2) / "377777777777 /
!!$!      data rmach(3) / "146400000000 /
!!$!      data rmach(4) / "147400000000 /
!!$!      data rmach(5) / "177464202324 /
!!$!
!!$!  PDP-11 FORTRANS supporting 32-bit integers (integer version).
!!$!
!!$!      data small(1) /    8388608 /
!!$!      data large(1) / 2147483647 /
!!$!      data right(1) /  880803840 /
!!$!      data diver(1) /  889192448 /
!!$!      data log10(1) / 1067065499 /
!!$!
!!$!  PDP-11 FORTRANS supporting 32-bit integers (octal version).
!!$!
!!$!      data rmach(1) / O00040000000 /
!!$!      data rmach(2) / O17777777777 /
!!$!      data rmach(3) / O06440000000 /
!!$!      data rmach(4) / O06500000000 /
!!$!      data rmach(5) / O07746420233 /
!!$!
!!$!  PDP-11 FORTRANS supporting 16-bit integers (integer version).
!!$!
!!$!      data small(1),small(2) /   128,     0 /
!!$!      data large(1),large(2) / 32767,    -1 /
!!$!      data right(1),right(2) / 13440,     0 /
!!$!      data diver(1),diver(2) / 13568,     0 /
!!$!      data log10(1),log10(2) / 16282,  8347 /
!!$!
!!$!  PDP-11 FORTRANS supporting 16-bit integers (octal version).
!!$!
!!$!      data small(1),small(2) / O000200, O000000 /
!!$!      data large(1),large(2) / O077777, O177777 /
!!$!      data right(1),right(2) / O032200, O000000 /
!!$!      data diver(1),diver(2) / O032400, O000000 /
!!$!      data log10(1),log10(2) / O037632, O020233 /
!!$!
!!$!  SEQUENT BALANCE 8000.
!!$!
!!$!      data small(1) / $00800000 /
!!$!      data large(1) / $7F7FFFFF /
!!$!      data right(1) / $33800000 /
!!$!      data diver(1) / $34000000 /
!!$!      data log10(1) / $3E9A209B /
!!$!
!!$!  SUN Microsystems UNIX F77 compiler.
!!$!
!!$!      data rmach(1) / 1.17549435E-38 /
!!$!      data rmach(2) / 3.40282347E+38 /
!!$!      data rmach(3) / 5.96016605E-08 /
!!$!      data rmach(4) / 1.19203321E-07 /
!!$!      data rmach(5) / 3.01030010E-01 /
!!$!
!!$!  SUN 3 (68881 or FPA)
!!$!
!!$!      data small(1) / X'00800000' /
!!$!      data large(1) / X'7F7FFFFF' /
!!$!      data right(1) / X'33800000' /
!!$!      data diver(1) / X'34000000' /
!!$!      data log10(1) / X'3E9A209B' /
!!$!
!!$!  UNIVAC 1100 series.
!!$!
!!$!      data rmach(1) / O000400000000 /
!!$!      data rmach(2) / O377777777777 /
!!$!      data rmach(3) / O146400000000 /
!!$!      data rmach(4) / O147400000000 /
!!$!      data rmach(5) / O177464202324 /
!!$!
!!$!  VAX/ULTRIX F77 compiler.
!!$!
!!$!      data small(1) /       128 /
!!$!      data large(1) /    -32769 /
!!$!      data right(1) /     13440 /
!!$!      data diver(1) /     13568 /
!!$!      data log10(1) / 547045274 /
!!$!
!!$!  VAX-11 with FORTRAN IV-PLUS compiler.
!!$!
!!$!      data rmach(1) / Z00000080 /
!!$!      data rmach(2) / ZFFFF7FFF /
!!$!      data rmach(3) / Z00003480 /
!!$!      data rmach(4) / Z00003500 /
!!$!      data rmach(5) / Z209B3F9A /
!!$!
!!$!  VAX/VMS version 2.2.
!!$!
!!$!      data rmach(1) /       '80'X /
!!$!      data rmach(2) / 'FFFF7FFF'X /
!!$!      data rmach(3) /     '3480'X /
!!$!      data rmach(4) /     '3500'X /
!!$!      data rmach(5) / '209B3F9A'X /
!!$!
!!$!  VAX/VMS 11/780
!!$!
!!$!      data small(1) / Z00000080 /
!!$!      data large(1) / ZFFFF7FFF /
!!$!      data right(1) / Z00003480 /
!!$!      data diver(1) / Z00003500 /
!!$!      data log10(1) / Z209B3F9A /
!!$!
!!$!  Z80 microprocessor.
!!$!
!!$!      data small(1), small(2) /     0,    256 /
!!$!      data large(1), large(2) /    -1,   -129 /
!!$!      data right(1), right(2) /     0,  26880 /
!!$!      data diver(1), diver(2) /     0,  27136 /
!!$!      data log10(1), log10(2) /  8347,  32538 /
!!$!
!!$  if ( i < 1 .or. i > 5 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'R1MACH - Fatal error!'
!!$    write ( *, '(a,i6)' ) 'I is out of bounds = ',i
!!$    r1mach = 0.0E+00
!!$    stop
!!$  else
!!$    r1mach = rmach(i)
!!$  end if
!!$
!!$  return
!!$end
!!$subroutine rminsp ( func, a, b, epsin, epsout, iop, result )
!!$!
!!$!***********************************************************************
!!$!
!!$!! RMINSP approximates the integral of a function using Romberg integration.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    Both the midpoint and trapezoidal rule are used,
!!$!    the intervals are repeatedly bisected, and Richardson
!!$!    extrapolation is carried out to achieve a high accuracy.
!!$!
!!$!    RMINSP can carry out a cosine-transform of the integral.  The
!!$!    only effect this has is to handle a function F(X) which has
!!$!    singularities near the endpoints.  This transform is based on
!!$!    the fact that
!!$!
!!$!      Integral from -1 to 1  ( F(    X))       DX
!!$!
!!$!    equals
!!$!
!!$!      Integral from  0 to PI ( F(COS(Z))*SIN(Z) )  DZ
!!$!
!!$!    If suitable accuracy is not achieved, the internal variable
!!$!    NUPPER might be increased.  Its current value of 9 corresponds
!!$!    to a maximum of 1024 subintervals and 1025 function evaluations.
!!$!
!!$!  Reference:
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!    T Havie,
!!$!    Algorithm 257,
!!$!    Communications of the Association for Computing Machinery,
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real, external FUNC, the name of the function to be
!!$!    integrated.  The user must declare the name an external
!!$!    parameter in the calling program, pass the name of the
!!$!    function in FUNC, and write a function of the form
!!$!
!!$!       FUNCTION FUNC(X)
!!$!
!!$!    which evaluates the function at the point X.
!!$!
!!$!    Input, real A, lower limit of integration.
!!$!
!!$!    Input, real B, upper limit of integration.
!!$!
!!$!    Input, real EPSIN, requested relative error tolerance.
!!$!
!!$!    Output, real EPSOUT, estimated achieved relative error.
!!$!
!!$!    Input, integer IOP, method switch:
!!$!    1, Use ordinary algorithm.
!!$!    2, Use cosine transformation to decrease effect of
!!$!       singularities near the endpoints.
!!$!
!!$!    Output, real RESULT, the approximate value of the integral.
!!$!
!!$  implicit none
!!$!
!!$  integer, parameter :: nupper = 9
!!$!
!!$  real a
!!$  real acof(11)
!!$  real alf
!!$  real alfnj
!!$  real alfno
!!$  real ar
!!$  real b
!!$  real bcof(nupper+1)
!!$  real bet
!!$  real betnj
!!$  real betno
!!$  real const1
!!$  real const2
!!$  real deltan
!!$  real endpts
!!$  real epsin
!!$  real epsout
!!$  real error
!!$  real, parameter :: fac1 = 0.411233516712057E+00
!!$  real, parameter :: fac2 = 0.822467033441132E+00
!!$  real factor
!!$  real, external :: func
!!$  real gamman
!!$  real hnstep
!!$  integer i
!!$  integer index
!!$  integer iop
!!$  integer iout
!!$  integer j
!!$  integer n
!!$  integer nhalf
!!$  real pi
!!$  real r1
!!$  real r2
!!$  real rn
!!$  real rnderr
!!$  real result
!!$  real rounde
!!$  real tend
!!$  real triarg
!!$  real umid
!!$  real xmin
!!$  real xplus
!!$!
!!$  result = 0.0E+00
!!$ 
!!$  if ( a == b ) then
!!$    return
!!$  end if
!!$!
!!$!  Set coefficients in formula for accumulated roundoff error,
!!$!  rounde = rnderr*(r1+r2*n), where r1, r2 are two empirical constants
!!$!  and n is the current number of function values used.
!!$!
!!$  rnderr = epsilon ( 1.0E+00 )
!!$ 
!!$  r1 = 1.0E+00
!!$  r2 = 2.0E+00
!!$  if ( iop==2 ) r1 = 50.0E+00
!!$  if ( iop==1 ) r2 = 0.01E+00 * r2
!!$  error = epsin
!!$!
!!$!  Initial calculations.
!!$!
!!$  alf = 0.5E+00 * (b-a)
!!$  bet = 0.5E+00 * (b+a)
!!$  acof(1) = func(a)+func(b)
!!$  bcof(1) = func(bet)
!!$!
!!$!  Modified Romberg algorithm, ordinary case.
!!$!
!!$  if ( iop /= 2 ) then
!!$
!!$    hnstep = 2.0E+00
!!$    bcof(1) = hnstep*bcof(1)
!!$    factor = 1.0E+00
!!$!
!!$!  Modified Romberg, cosine transformed case.
!!$!
!!$  else
!!$    hnstep = pi()
!!$    ar = fac1
!!$    endpts = acof(1)
!!$    acof(1) = fac2*acof(1)
!!$    bcof(1) = hnstep*bcof(1)-ar*endpts
!!$    factor = 4.0E+00
!!$    ar = ar/4.0E+00
!!$    triarg = pi() / 4.0E+00
!!$    alfno = -1.0E+00
!!$  end if
!!$ 
!!$  hnstep = 0.5E+00 * hnstep
!!$  nhalf = 1
!!$  n = 2
!!$  rn = 2.0E+00
!!$  acof(1) = 0.5E+00 * (acof(1)+bcof(1))
!!$  acof(2) = acof(1)-(acof(1)-bcof(1))/(4.0E+00*factor-1.0E+00)
!!$!
!!$!  End of initial calculation.
!!$!
!!$!  Start actual calculations.
!!$!
!!$  do i = 1, nupper
!!$ 
!!$    umid = 0.0E+00
!!$!
!!$!  Modified Romberg algorithm, ordinary case.
!!$!  compute first element in mid-point formula for ordinary case
!!$!
!!$    if ( iop == 1 ) then
!!$ 
!!$      alfnj = 0.5E+00*hnstep
!!$ 
!!$      do j = 1, nhalf
!!$        xplus = alf*alfnj+bet
!!$        xmin = -alf*alfnj+bet
!!$        umid = umid + func(xplus)+func(xmin)
!!$        alfnj = alfnj+hnstep
!!$      end do
!!$ 
!!$      umid = hnstep*umid
!!$!
!!$!  Modified Romberg algorithm, cosine transformed case
!!$!  compute first element in mid-point formula for cosine transformed
!!$!  Romberg algorithm
!!$!
!!$    else if ( iop == 2 ) then
!!$ 
!!$      const1 = -sin(triarg)
!!$      const2 = 0.5E+00 * alfno/const1
!!$ 
!!$      alfno = const1
!!$      betno = const2
!!$      gamman = 1.0E+00 - 2.0E+00 * alfno**2
!!$      deltan = -2.0E+00 * alfno*betno
!!$ 
!!$      do j = 1, nhalf
!!$        alfnj = gamman*const1+deltan*const2
!!$        betnj = gamman*const2-deltan*const1
!!$        xplus = alf*alfnj+bet
!!$        xmin = -alf*alfnj+bet
!!$        umid = umid+betnj*(func(xplus)+func(xmin))
!!$        const1 = alfnj
!!$        const2 = betnj
!!$      end do
!!$ 
!!$      umid = hnstep*umid-ar*endpts
!!$      ar = ar / 4.0E+00
!!$ 
!!$    end if
!!$!
!!$!  Modified Romberg algorithm, calculate (i+1)-th row in the U table
!!$!
!!$    const1 = 4.0E+00 * factor
!!$    index = i+1
!!$ 
!!$    do j = 2, i+1
!!$      tend = umid + ( umid - bcof(j-1) ) / ( const1 - 1.0E+00 )
!!$      bcof(j-1) = umid
!!$      umid = tend
!!$      const1 = 4.0E+00 * const1
!!$    end do
!!$ 
!!$    bcof(i+1) = tend
!!$    xplus = const1
!!$!
!!$!  Calculation of (i+1)-th row in the U table is finished
!!$!
!!$!  Test to see if the required accuracy has been obtained.
!!$!
!!$    epsout = 1.0E+00
!!$    iout = 1
!!$ 
!!$    do j = 1, index
!!$ 
!!$      const1 = 0.5E+00 * ( acof(j) + bcof(j) )
!!$      const2 = 0.5E+00 * abs ( ( acof(j) - bcof(j) ) / const1 )
!!$ 
!!$      if ( const2 <= epsout ) then
!!$        epsout = const2
!!$        iout = j
!!$      end if
!!$ 
!!$      acof(j) = const1
!!$ 
!!$    end do
!!$!
!!$!  Testing on accuracy finished
!!$!
!!$    if (iout==index) iout=iout+1
!!$    acof(index+1) = acof(index)-(acof(index)-bcof(index))/(xplus-1.0)
!!$    rounde = rnderr*(r1+r2*rn)
!!$
!!$    epsout = max ( epsout, rounde )
!!$    error = max ( error, rounde )
!!$ 
!!$    if ( epsout <= error ) go to 10
!!$ 
!!$    nhalf = n
!!$    n = 2 * n
!!$    rn = 2.0E+00 * rn
!!$    hnstep = 0.5E+00 * hnstep
!!$
!!$    if ( iop > 1 ) then
!!$      triarg = 0.5 * triarg
!!$    end if
!!$ 
!!$  end do
!!$!
!!$!  Accuracy not reached with maximum number of subdivisions
!!$!
!!$  n = nhalf
!!$!
!!$!  Calculation for modified Romberg algorithm finished
!!$!
!!$  10  continue
!!$ 
!!$  n = 2*n
!!$  index = index-1
!!$  n = n+1
!!$  j = iout
!!$  if ((j-1)>=index) j = index
!!$  tend = alf * (2.0E+00*acof(j)-bcof(j))
!!$  umid = alf * bcof(j)
!!$  result = alf * acof(iout)
!!$
!!$  return
!!$end
!!$subroutine rvec_even ( alo, ahi, n, a )
!!$!
!!$!*******************************************************************************
!!$!
!!$!! RVEC_EVEN returns N real values, evenly spaced between ALO and AHI.
!!$!
!!$!
!!$!  Modified:
!!$!
!!$!    31 October 2000
!!$!
!!$!  Author:
!!$!
!!$!    John Burkardt
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real ALO, AHI, the low and high values.
!!$!
!!$!    Input, integer N, the number of values.
!!$!
!!$!    Output, real A(N), N evenly spaced values.
!!$!    Normally, A(1) = ALO and A(N) = AHI.
!!$!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!!$!
!!$!
!!$  implicit none
!!$!
!!$  integer n
!!$!
!!$  real a(n)
!!$  real ahi
!!$  real alo
!!$  integer i
!!$!
!!$  if ( n == 1 ) then
!!$
!!$    a(1) = 0.5E+00 * ( alo + ahi )
!!$
!!$  else
!!$
!!$    do i = 1, n
!!$      a(i) = ( real ( n - i ) * alo + real ( i - 1 ) * ahi ) / real ( n - 1 )
!!$    end do
!!$
!!$  end if
!!$
!!$  return
!!$end
!!$subroutine simp ( func, a, b, eps, result )
!!$!
!!$!***********************************************************************
!!$!
!!$!! SIMP approximates the integral of a function using an adaptive Simpson's rule.
!!$!
!!$!
!!$!  Reference:
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!    J N Lyness,
!!$!    Algorithm 379,
!!$!    SQUANK - Simpson quadrature used adaptively, noise killed,
!!$!    Communications of the Association for Computing Machinery,
!!$!    Volume 13 (1970), pages 260-263.
!!$!
!!$!    W M McKeeman and L Tesler,
!!$!    Algorithm 182,
!!$!    Nonrecursive adaptive integration,
!!$!    Communications of the Association for Computing Machinery,
!!$!    Volume 6 (1963), page 315.
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real, external FUNC, the name of the function to be
!!$!    integrated.  The user must declare the name an external
!!$!    parameter in the calling program, pass the name of the
!!$!    function in FUNC, and write a function of the form
!!$!
!!$!    FUNCTION FUNC(X)
!!$!
!!$!    which evaluates the function at the point X.
!!$!
!!$!    Input, real A, the lower limit of integration.
!!$!
!!$!    Input, real B, the upper limit of integration.
!!$!
!!$!    Input, real EPS, the requested error tolerance.
!!$!
!!$!    Output, real RESULT, the approximation to the integral.
!!$!
!!$  implicit none
!!$!
!!$  integer, parameter :: maxlev = 30
!!$!
!!$  real a
!!$  real a1
!!$  real absar
!!$  real b
!!$  real da
!!$  real dx(maxlev)
!!$  real ep
!!$  real eps
!!$  real epsp(maxlev)
!!$  real est
!!$  real est1
!!$  real est2(maxlev)
!!$  real est3(maxlev)
!!$  real f1
!!$  real f2(maxlev)
!!$  real f3(maxlev)
!!$  real f4(maxlev)
!!$  real fa
!!$  real fb
!!$  real fbp(maxlev)
!!$  real fm
!!$  real fmp(maxlev)
!!$  real, external :: func
!!$  integer i
!!$  integer j
!!$  integer l
!!$  integer lvl
!!$  integer nrtr(maxlev)
!!$  real pval(maxlev,3)
!!$  real result
!!$  real sum1
!!$  real sx
!!$  real x2(maxlev)
!!$  real x3(maxlev)
!!$!
!!$  result = 0.0E+00
!!$ 
!!$  if ( a == b ) then
!!$    return
!!$  end if
!!$ 
!!$  ep = eps
!!$  a1 = a
!!$  nrtr(1:maxlev) = 0
!!$  pval(1:maxlev,1:3) = 0.0E+00
!!$ 
!!$  lvl = 0
!!$  absar = 0.0E+00
!!$  est = 0.0E+00
!!$  da = b - a1
!!$
!!$  fa = func ( a1 )
!!$  fm = 4.0E+00 * func ( (a1+b) * 0.5E+00 )
!!$  fb = func ( b )
!!$!
!!$!  1 = RECUR
!!$!
!!$   30 continue
!!$ 
!!$  lvl = lvl + 1
!!$  dx(lvl) = da / 3.0E+00
!!$  sx = dx(lvl)/6.0E+00
!!$  f1 = 4.0E+00 * func(0.5*dx(lvl)+a1)
!!$  x2(lvl) = a1+dx(lvl)
!!$  f2(lvl) = func(x2(lvl))
!!$  x3(lvl) = x2(lvl)+dx(lvl)
!!$  f3(lvl) = func(x3(lvl))
!!$  epsp(lvl) = ep
!!$  f4(lvl) = 4.0E+00 * func(dx(lvl)*0.5E+00+x3(lvl))
!!$  fmp(lvl) = fm
!!$  est1 = sx*(fa+f1+f2(lvl))
!!$  fbp(lvl) = fb
!!$  est2(lvl) = sx * (f2(lvl)+f3(lvl)+fm)
!!$  est3(lvl) = sx * (f3(lvl)+f4(lvl)+fb)
!!$  sum1 = est1+est2(lvl)+est3(lvl)
!!$  absar = absar - abs ( est ) + abs ( est1 ) + abs ( est2(lvl) ) &
!!$    + abs ( est3(lvl) )
!!$  if ( abs ( est - sum1 ) <= epsp(lvl) * absar ) go to 40
!!$  if ( lvl >= maxlev ) go to 50
!!$!
!!$!  2 = UP
!!$!
!!$40 continue
!!$ 
!!$  if ( lvl > 1 ) then
!!$    lvl = lvl-1
!!$  end if
!!$
!!$  l = nrtr(lvl)
!!$
!!$  if ( l == 0 ) then
!!$    go to 50
!!$  end if
!!$
!!$  pval(lvl,l) = sum1
!!$
!!$  if ( l == 1 ) go to 60
!!$  if ( l == 2 ) go to 70
!!$  if ( l == 3 ) go to 80
!!$ 
!!$50 continue
!!$ 
!!$  nrtr(lvl) = 1
!!$  est = est1
!!$  fm = f1
!!$  fb = f2(lvl)
!!$  ep = epsp(lvl) / 1.7E+00
!!$  da = dx(lvl)
!!$  go to 30
!!$ 
!!$60 continue
!!$ 
!!$  nrtr(lvl) = 2
!!$  fa = f2(lvl)
!!$  fm = fmp(lvl)
!!$  fb = f3(lvl)
!!$  est = est2(lvl)
!!$  a1 = x2(lvl)
!!$  ep = epsp(lvl) / 1.7E+00
!!$  da = dx(lvl)
!!$  go to 30
!!$ 
!!$70 continue
!!$ 
!!$  nrtr(lvl) = 3
!!$  fa = f3(lvl)
!!$  fm = f4(lvl)
!!$  fb = fbp(lvl)
!!$  est = est3(lvl)
!!$  a1 = x3(lvl)
!!$  ep = epsp(lvl) / 1.7E+00
!!$  da = dx(lvl)
!!$  go to 30
!!$ 
!!$80 continue
!!$ 
!!$  sum1 = pval(lvl,1)+pval(lvl,2)+pval(lvl,3)
!!$
!!$  if ( lvl > 1 ) then
!!$    go to 40
!!$  end if
!!$ 
!!$90 continue
!!$ 
!!$  result = sum1
!!$ 
!!$  return
!!$end
!!$subroutine simpne ( x, y, num, result )
!!$!
!!$!***********************************************************************
!!$!
!!$!! SIMPNE approximates the integral of unevenly spaced data.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    The routine repeatedly interpolates a 3-point Lagrangian polynomial 
!!$!    to the data and integrates that exactly.
!!$!
!!$!  Reference:
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real X(NUM), contains the X values of the data, in order.
!!$!
!!$!    Input, real Y(NUM), contains the Y values of the data.
!!$!
!!$!    Input, integer NUM, number of data points.  NUM must be at least 3.
!!$!
!!$!    Output, real RESULT.
!!$!    RESULT is the approximate value of the integral.
!!$!
!!$  implicit none
!!$!
!!$  integer num
!!$!
!!$  real del(3)
!!$  real e
!!$  real f
!!$  real feints
!!$  real g(3)
!!$  integer i
!!$  integer n
!!$  real pi(3)
!!$  real result
!!$  real sum1
!!$  real x(num)
!!$  real x1
!!$  real x2
!!$  real x3
!!$  real y(num)
!!$!
!!$  result = 0.0E+00
!!$ 
!!$  if ( num <= 2 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'SIMPNE - Fatal error!'
!!$    write ( *, '(a)' ) '  NUM <= 2.'
!!$    stop
!!$  end if
!!$ 
!!$  n = 1
!!$ 
!!$  do
!!$ 
!!$    x1 = x(n)
!!$    x2 = x(n+1)
!!$    x3 = x(n+2)
!!$    e = x3*x3-x1*x1
!!$    f = x3*x3*x3-x1*x1*x1
!!$    feints = x3-x1
!!$    del(1) = x3-x2
!!$    del(2) = x1-x3
!!$    del(3) = x2-x1
!!$    g(1) = x2+x3
!!$    g(2) = x1+x3
!!$    g(3) = x1+x2
!!$    pi(1) = x2*x3
!!$    pi(2) = x1*x3
!!$    pi(3) = x1*x2
!!$ 
!!$    sum1 = 0.0E+00
!!$    do i = 1, 3
!!$      sum1 = sum1 + y(n-1+i)*del(i)*(f/3.0E+00-g(i)*0.5E+00*e+pi(i)*feints)
!!$    end do
!!$    result = result - sum1 / ( del(1) * del(2) * del(3) )
!!$ 
!!$    n = n+2
!!$
!!$    if ( n + 1 >= num ) then
!!$      exit
!!$    end if
!!$
!!$  end do
!!$ 
!!$  if ( mod(num,2) /= 0 ) then
!!$    return
!!$  end if
!!$
!!$  n = num-2
!!$  x3 = x(num)
!!$  x2 = x(num-1)
!!$  x1 = x(num-2)
!!$  e = x3*x3-x2*x2
!!$  f = x3*x3*x3-x2*x2*x2
!!$  feints = x3-x2
!!$  del(1) = x3-x2
!!$  del(2) = x1-x3
!!$  del(3) = x2-x1
!!$  g(1) = x2+x3
!!$  g(2) = x1+x3
!!$  g(3) = x1+x2
!!$  pi(1) = x2*x3
!!$  pi(2) = x1*x3
!!$  pi(3) = x1*x2
!!$ 
!!$  sum1 = 0.0E+00
!!$  do i = 1, 3
!!$    sum1 = sum1 + y(n-1+i) * del(i) * &
!!$      ( f / 3.0E+00 - g(i) * 0.5E+00 * e + pi(i) * feints )
!!$  end do
!!$ 
!!$  result = result - sum1 / ( del(1) * del(2) * del(3) )
!!$ 
!!$  return
!!$end
!!$subroutine simpsn ( h, y, num, result )
!!$!
!!$!***********************************************************************
!!$!
!!$!! SIMPSN approximates the integral of evenly spaced data.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    Simpson's rule is used.
!!$!
!!$!  Reference:
!!$!
!!$!    Philip Davis and Philip Rabinowitz,
!!$!    Methods of Numerical Integration,
!!$!    Blaisdell Publishing, 1967.
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real H, specifies the increment between the
!!$!    X values.  Note that the actual X values are not needed,
!!$!    just the constant spacing!
!!$!
!!$!    Input, real Y(NUM), the data.
!!$!
!!$!    Input, integer NUM, the number of data points.  NUM must be at least 3.
!!$!
!!$!    Output, real RESULT, the value of the integral
!!$!    from the first to the last point.
!!$!
!!$  implicit none
!!$!
!!$  integer num
!!$!
!!$  real del(3)
!!$  real f
!!$  real g(3)
!!$  real h
!!$  integer i
!!$  integer n
!!$  real pii(3)
!!$  real result
!!$  real sum1
!!$  real y(num)
!!$!
!!$  result = 0.0E+00
!!$ 
!!$  if ( num <= 2 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'SIMPSN - Fatal error!'
!!$    write ( *, '(a,i6)' ) '  NUM < 2, NUM = ', num
!!$    stop
!!$  end if
!!$ 
!!$  if ( mod ( num, 2 ) == 0 ) then
!!$    n = num-1
!!$  else
!!$    n = num
!!$  end if
!!$ 
!!$  result = y(1) + y(n) + 4.0E+00 * y(n-1)
!!$  do i = 2, n-2, 2
!!$    result = result + 4.0E+00 * y(i) + 2.0E+00 * y(i+1)
!!$  end do
!!$  result = h * result / 3.0E+00
!!$ 
!!$  if ( mod(num,2) == 1 ) then
!!$    return
!!$  end if
!!$ 
!!$  f = h*h*h
!!$  del(1) = h
!!$  del(2) = -2.0E+00 * h
!!$  del(3) = h
!!$  g(1) = h
!!$  g(2) = 0.0E+00
!!$  g(3) = -h
!!$  pii(1) = 0.0E+00
!!$  pii(2) = -h*h
!!$  pii(3) = 0.0E+00
!!$  n = n-1
!!$ 
!!$  sum1 = 0.0E+00
!!$  do i = 1, 3
!!$    sum1 = sum1 + y(n-1+i) * del(i) * &
!!$      ( f / 3.0E+00 - g(i) * 0.5E+00 * h * h + pii(i) * h )
!!$  end do
!!$ 
!!$  result = result + 0.5E+00 * sum1 / h**3
!!$ 
!!$  return
!!$end
!!$function solve ( shift, n, a, b )
!!$!
!!$!***********************************************************************
!!$!
!!$!! SOLVE solves a special linear system.
!!$!
!!$!
!!$!  Discussion:
!!$!
!!$!    SOLVE solves for the N-th component of the solution DELTA to the equation
!!$!
!!$!      (Jn - shift*Identity) * DELTA  = En,
!!$!
!!$!    En is the vector of all zeroes except for 1 in the N-th position.
!!$!
!!$!    The matrix Jn is symmetric tridiagonal, with diagonal
!!$!    elements A(I), off-diagonal elements B(I).  This equation
!!$!    must be solved to obtain the appropriate changes in the lower
!!$!    2 by 2 submatrix of coefficients for orthogonal polynomials.
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real SHIFT, the value of the factor that multiplies
!!$!    the identity matrix, in the definition of the system matrix.
!!$!
!!$!    Input, integer N, the index of the desired component.
!!$!
!!$  implicit none
!!$!
!!$  integer n
!!$!
!!$  real a(n)
!!$  real alpha
!!$  real b(n)
!!$  integer i
!!$  real shift
!!$  real solve
!!$!
!!$  alpha = a(1) - shift
!!$  do i = 2, n-1
!!$    alpha = a(i) - shift - b(i-1)**2 / alpha
!!$  end do
!!$ 
!!$  solve = 1.0E+00 / alpha
!!$ 
!!$  return
!!$end
!!$subroutine timestamp ( )
!!$!
!!$!*******************************************************************************
!!$!
!!$!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!!$!
!!$!
!!$!  Example:
!!$!
!!$!    May 31 2001   9:45:54.872 AM
!!$!
!!$!  Modified:
!!$!
!!$!    31 May 2001
!!$!
!!$!  Author:
!!$!
!!$!    John Burkardt
!!$!
!!$!  Parameters:
!!$!
!!$!    None
!!$!
!!$  implicit none
!!$!
!!$  character ( len = 8 ) ampm
!!$  integer d
!!$  character ( len = 8 ) date
!!$  integer h
!!$  integer m
!!$  integer mm
!!$  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
!!$    'January  ', 'February ', 'March    ', 'April    ', &
!!$    'May      ', 'June     ', 'July     ', 'August   ', &
!!$    'September', 'October  ', 'November ', 'December ' /)
!!$  integer n
!!$  integer s
!!$  character ( len = 10 )  time
!!$  integer values(8)
!!$  integer y
!!$  character ( len = 5 ) zone
!!$!
!!$  call date_and_time ( date, time, zone, values )
!!$
!!$  y = values(1)
!!$  m = values(2)
!!$  d = values(3)
!!$  h = values(5)
!!$  n = values(6)
!!$  s = values(7)
!!$  mm = values(8)
!!$
!!$  if ( h < 12 ) then
!!$    ampm = 'AM'
!!$  else if ( h == 12 ) then
!!$    if ( n == 0 .and. s == 0 ) then
!!$      ampm = 'Noon'
!!$    else
!!$      ampm = 'PM'
!!$    end if
!!$  else
!!$    h = h - 12
!!$    if ( h < 12 ) then
!!$      ampm = 'PM'
!!$    else if ( h == 12 ) then
!!$      if ( n == 0 .and. s == 0 ) then
!!$        ampm = 'Midnight'
!!$      else
!!$        ampm = 'AM'
!!$      end if
!!$    end if
!!$  end if
!!$
!!$  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
!!$    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )
!!$
!!$  return
!!$end
!!$subroutine wedint ( ftab, h, ntab, result )
!!$!
!!$!***********************************************************************
!!$!
!!$!! WEDINT uses Weddle's rule to integrate data at equally spaced points.
!!$!
!!$!
!!$!  Modified:
!!$!
!!$!    30 October 2000
!!$!
!!$!  Author:
!!$!
!!$!    John Burkardt
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real FTAB(NTAB), contains the tabulated data values.
!!$!
!!$!    Input, real H, is the spacing between the points at which the data
!!$!    was evaluated.
!!$!
!!$!    Input, integer NTAB, is the number of data points.  (NTAB-1) must be
!!$!    divisible by 6.
!!$!
!!$!    Output, real RESULT, is the approximation to the integral.
!!$!
!!$  implicit none
!!$!
!!$  integer ntab
!!$!
!!$  real ftab(ntab)
!!$  real h
!!$  integer i
!!$  real result
!!$!
!!$  result = 0.0E+00
!!$ 
!!$  if ( ntab <= 1 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'WEDINT - Fatal error!'
!!$    write ( *, '(a)' ) '  NTAB < 2'
!!$    write ( *, '(a,i6)' ) '  NTAB = ', ntab
!!$    stop
!!$  end if
!!$ 
!!$  if ( mod ( ntab, 6 ) /= 1 ) then
!!$    write ( *, '(a)' ) ' '
!!$    write ( *, '(a)' ) 'WEDINT - Fatal error!'
!!$    write ( *, '(a)' ) '  NTAB must equal 6*N+1 for some N!'
!!$    stop
!!$  end if
!!$ 
!!$  do i = 1, ntab-6, 6
!!$    result = result & 
!!$      +           ftab(i) &
!!$      + 5.0E+00 * ftab(i+1) &
!!$      +           ftab(i+2) &
!!$      + 6.0E+00 * ftab(i+3) &
!!$      +           ftab(i+4) &
!!$      + 5.0E+00 * ftab(i+5) &
!!$      +           ftab(i+6)
!!$  end do
!!$ 
!!$  result = 3.0E+00 * h * result / 10.0E+00
!!$ 
!!$  return
!!$end

END MODULE INTLIB
