!********************************************************************************
! Automatic Differentiation (AD) module, Version 1 (2012),
!
!        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!
! the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!
!
!
! This module contains the derivative data type (ddt) and associated operators,
! specifically written for storing 5 derivatives (thus, ddt5).
!
! To use the data type and operators, declare "use derivative_data".
! Then, declare a scalar derivative_data_type variable, say 'a', by
!
!  type(derivative_data_type_df5) :: a
!
! Then, 'a' has components, a%f to store a function value, and the array of
! dimension 5, a%df(1:5) to store derivatives.
!
! Operators are overloaded as defined in this module. That is, the operators
! such as + or * perform not only + or * on %f, but also differentiations in
! %df components automatically. Currently available operators are
!
! [ assignment, +, -, *, /, **, ddt_sqrt, ddt_abs, ddt_sign, ddt_max, ddt_min,
!   <, <=, >, >=, == ]
!
! More operators may be added as needed.
!
! 
! Example:
!
! It may be easier to understand if the variable is declared as array of
! dimension of 5 as follows:
!
!  type(derivative_data_type_df5), dimension(5) :: a
!
! Then, the variable 'a' contains function values in %f:
!
!    [ a(1)%f, a(2)%f, a(3)%f, a(4)%f, a(5)%f ]
!
! and 5x5 derivatives in %df
!
!    | a(1)%df(1), a(1)%df(2), a(1)%df(3), a(1)%df(4), a(1)%df(5) |
!    | a(2)%df(1), a(2)%df(2), a(2)%df(3), a(2)%df(4), a(2)%df(5) |
!    | a(3)%df(1), a(3)%df(2), a(3)%df(3), a(3)%df(4), a(3)%df(5) |
!    | a(4)%df(1), a(4)%df(2), a(4)%df(3), a(4)%df(4), a(4)%df(5) |
!    | a(5)%df(1), a(5)%df(2), a(5)%df(3), a(5)%df(4), a(5)%df(5) |
!
! Suppose we have a vector of primitive variables in a real array, w(:),
!
!    real, dimension(5) :: w
!    w = [ 1.1, 2.0, 3.0, 4.0, 5.0 ] 
!
! and, we give it to the ddt array, w_ddt, by '=',
!
!    type(derivative_data_type_df5), dimension(5) :: w_ddt
!    w_ddt = w
!
! Then, we have
!
!    w_ddt(1:5)%f       = [ 1.0, 2.1, 0.1, 0.0, 1.2 ] 
!    w_ddt(1:5)%df(1:5) = [ 5x5 zero matrix ]
!
! To compute the derivatives of some quantities with respect to the
! primitive variables, du/dw, we first set
!
!    call ddt_seed(w_ddt)
!
! This will set
!
!    w_ddt(1:5)%df(1:5) = [ 5x5 identity matrix ]
!
! indicating that %df is now defined as the derivatives of w with respect to w:
!
!    w_ddt(1:5)%df(1:5)
!
!     = | dw1/dw1 dw1/dw2 dw1/dw3 dw1/dw4 dw1/dw5 |
!       | dw2/dw1 dw2/dw2 dw2/dw3 dw2/dw4 dw2/dw5 |
!       | dw3/dw1 dw3/dw2 dw3/dw3 dw3/dw4 dw3/dw5 |
!       | dw4/dw1 dw4/dw2 dw4/dw3 dw4/dw4 dw4/dw5 |
!       | dw5/dw1 dw5/dw2 dw5/dw3 dw5/dw4 dw5/dw5 |
!
!     = | 1 0 0 0 0 |
!       | 0 1 0 0 0 |
!       | 0 0 1 0 0 |
!       | 0 0 0 1 0 |
!       | 0 0 0 0 1 |
!
! Any subsequent operations will be based on this definition, and at the end,
! the resulting derivatives will be with respect to w. For example, if we compute
! the conservative variables, which is declared by
!
!   type(derivative_data_type_df5), dimension(5) :: u_ddt
!
! by the following usual operations:
!
!  u_ddt(1) = w_ddt(1)
!  u_ddt(2) = w_ddt(1)*w_ddt(2)
!  u_ddt(3) = w_ddt(1)*w_ddt(3)
!  u_ddt(4) = w_ddt(1)*w_ddt(4)
!  u_ddt(5) = w_ddt(5)/(gamma-one)+half*w_ddt(1)*(w_ddt(2)*w_ddt(2)+w_ddt(3)*w_ddt(3)+w_ddt(4)*w_ddt(4))
!
! Then, we obtain
!
!   u_ddt(1:5)%f       = [ 1.0, 2.2, 3.3, 4.4, 5/(gamma-1)+29/2 ] !conservative variables
!
!   u_ddt(1:5)%df(1:5) = | du1/dw1 du1/dw2 du1/dw3 du1/dw4 du1/dw5 |
!                        | du2/dw1 du2/dw2 du2/dw3 du2/dw4 du2/dw5 |
!                        | du3/dw1 du3/dw2 du3/dw3 du3/dw4 du3/dw5 |
!                        | du4/dw1 du4/dw2 du4/dw3 du4/dw4 du4/dw5 |
!                        | du5/dw1 du5/dw2 du5/dw3 du5/dw4 du5/dw5 |   !derivatives w.r.t. w
!
! So, the derivatives are computed automatically.
!
!
! Note: This module can be easily extended to accomodate more derivatives.
!       Replace "df5" by "df20", for example, to define the derivative data type
!       having %df of dimension 20 ( also modify dimension(5) -> dimension(20) ).
!
! Katate Masatsuka, November 2012. http://www.cfdbooks.com
!********************************************************************************
 module derivative_data_df5

  implicit none

! Select the precision:
!  selected_real_kind( 5) = single precision
!  selected_real_kind(15) = double precision

  integer           , parameter :: my_precision = selected_real_kind(15)
  real(my_precision), parameter :: zero = 0.0_my_precision
  real(my_precision), parameter :: half = 0.5_my_precision
  real(my_precision), parameter ::  one = 1.0_my_precision

  private

!--------------------------------------------------------------------------------
! Derivative data type definition

  public :: derivative_data_type_df5

  type derivative_data_type_df5
   real(my_precision)               ::  f ! function value
   real(my_precision), dimension(5) :: df ! df/du
  end type derivative_data_type_df5

!--------------------------------------------------------------------------------
! Set derivatives w.r.t. which derivatives of a functional are computed.
! Basically, this sets variable(i)%df(i)=1.0
! E.g.,  call ddt_seed(ddt variable)

  public :: ddt_seed
  interface ddt_seed
   module procedure ddt_seed_scalar
   module procedure ddt_seed_vector
  end interface

!--------------------------------------------------------------------------------
! Assignment operator
! E.g.,  ddt_variable1 = ddt_variable2
! E.g.,  ddt_variable  = real value
! E.g.,    real value  = ddt_variable

  public :: assignment(=)
  interface assignment(=)
   module procedure d_assign_d
   module procedure d_assign_r
   module procedure r_assign_d
  end interface

!--------------------------------------------------------------------------------
! Addition operator
! E.g.,  ddt_variable1 + ddt_variable2
! E.g.,  ddt_variable  + real value
! E.g.,    real value  + ddt_variable

  public :: operator(+)
  interface operator(+)
   module procedure d_plus_d
   module procedure d_plus_r
   module procedure r_plus_d
  end interface

!--------------------------------------------------------------------------------
! Subtraction operator
! E.g.,    ddt_variable1 - ddt_variable2
! E.g.,    ddt_variable  - real value
! E.g.,      real value  - ddt_variable
! E.g.,  - ddt_variable

  public :: operator(-)
  interface operator(-)
   module procedure d_minus_d
   module procedure d_minus_r
   module procedure r_minus_d
   module procedure   minus_d
  end interface

!--------------------------------------------------------------------------------
! Multiplication operator
! E.g.,    ddt_variable1 * ddt_variable2
! E.g.,    ddt_variable  * real value
! E.g.,      real value  * ddt_variable

  public :: operator(*)
  interface operator(*)
   module procedure d_times_d
   module procedure d_times_r
   module procedure r_times_d
  end interface

!--------------------------------------------------------------------------------
! Division operator
! E.g.,    ddt_variable1 / ddt_variable2
! E.g.,    ddt_variable  / real value
! E.g.,      real value  / ddt_variable

  public :: operator(/)
  interface operator(/)
   module procedure d_divided_by_d
   module procedure d_divided_by_r
   module procedure r_divided_by_d
  end interface

!--------------------------------------------------------------------------------
! Exponent operator
! E.g.,    ddt_variable1 ** ddt_variable2
! E.g.,    ddt_variable  ** (real value)
! E.g.,     (real value) ** ddt_variable
! E.g.,    ddt_variable  ** (integer)

  public :: operator(**)
  interface operator(**)
   module procedure d_power_d
   module procedure d_power_r
   module procedure r_power_d
   module procedure d_power_i
  end interface

!--------------------------------------------------------------------------------
! Square root operator
! E.g.,   ddt_sqrt(ddt_variable)

  public :: ddt_sqrt
  interface ddt_sqrt
   module procedure sqrt_d
  end interface

!--------------------------------------------------------------------------------
! Square root operator
! E.g.,   ddt_tanh(ddt_variable)

  public :: ddt_tanh
  interface ddt_tanh
   module procedure tanh_d
  end interface

!--------------------------------------------------------------------------------
! Absolute value operator
! E.g.,  ddt_abs(ddt_variable)

  public :: ddt_abs
  interface ddt_abs
   module procedure abs_d
  end interface

!--------------------------------------------------------------------------------
! Sign operator
! E.g.,  ddt_sign(    one, ddt_variable)
! E.g.,  ddt_sign(rael value, ddt_variable)

  public :: ddt_sign
  interface ddt_sign
   module procedure sign_rd
  end interface

!--------------------------------------------------------------------------------
! Max operator
! E.g.,  ddt_max(ddt_variable1, ddt_variable2)
! E.g.,  ddt_max(ddt_variable1,    real value)
! E.g.,  ddt_max(   real value, ddt_variable2)

  public :: ddt_max
  interface ddt_max
   module procedure max_of_d_and_d
   module procedure max_of_d_and_r
   module procedure max_of_r_and_d
  end interface

!--------------------------------------------------------------------------------
! Min operator
! E.g.,  ddt_min(ddt_variable1, ddt_variable2)
! E.g.,  ddt_min(ddt_variable1,    real value)
! E.g.,  ddt_min(   real value, ddt_variable2)

  public :: ddt_min
  interface ddt_min
   module procedure min_of_d_and_d
   module procedure min_of_d_and_r
   module procedure min_of_r_and_d
  end interface

!--------------------------------------------------------------------------------
! Less-than (logical). It compares the function values.
! E.g.,  ddt_variable1 < ddt_variable2
! E.g.,  ddt_variable  < real value
! E.g.,    real value  < ddt_variable

  public :: operator(<)
  interface operator(<)
   module procedure d_less_than_d
   module procedure d_less_than_r
   module procedure r_less_than_d
  end interface

!--------------------------------------------------------------------------------
! Less-than-equal (logical). It compares the function values.
! E.g.,  ddt_variable1 <= ddt_variable2
! E.g.,  ddt_variable  <= real value
! E.g.,    real value  <= ddt_variable

  public :: operator(<=)
  interface operator(<=)
   module procedure d_less_than_equal_d
   module procedure d_less_than_equal_r
   module procedure r_less_than_equal_d
  end interface

!--------------------------------------------------------------------------------
! Greater-than (logical). It compares the function values.
! E.g.,  ddt_variable1 > ddt_variable2
! E.g.,  ddt_variable  > real value
! E.g.,    real value  > ddt_variable

  public :: operator(>)
  interface operator(>)
   module procedure d_greater_than_d
   module procedure d_greater_than_r
   module procedure r_greater_than_d
  end interface

!--------------------------------------------------------------------------------
! Greater-than-equal (logical). It compares the function values.
! E.g.,  ddt_variable1 >= ddt_variable2
! E.g.,  ddt_variable  >= real value
! E.g.,    real value  >= ddt_variable

  public :: operator(>=)
  interface operator(>=)
   module procedure d_greater_than_equal_d
   module procedure d_greater_than_equal_r
   module procedure r_greater_than_equal_d
  end interface

!--------------------------------------------------------------------------------
! Equal (logical). It compares the function values.
! E.g.,  ddt_variable1 == ddt_variable2
! E.g.,  ddt_variable  == real value
! E.g.,    real value  == ddt_variable

  public :: operator(==)
  interface operator(==)
   module procedure d_equal_d
   module procedure d_equal_r
   module procedure r_equal_d
  end interface

 contains

!*******************************************************************************
! Seeding: Set the diagonal components to be 1.0, so that derivatives are
!          computed w.r.t. these variables. In a way, set the correct derivative.
!
!  Input: d = the derivative_data_type variable to be seeded.
! Output: d with d%df = 1.0 (only the diagonals, d(f_k)/d(f_k) = 1.0)
!
!*******************************************************************************
  subroutine ddt_seed_scalar(d)

   type(derivative_data_type_df5), intent(inout) :: d

    d%df = one

  end subroutine ddt_seed_scalar
!-------------------------------------------------------------------------------
  subroutine ddt_seed_vector(d)

   type(derivative_data_type_df5), dimension(:), intent(inout) :: d
   integer                                   :: n, i

    n = size(d(1)%df)

    do i = 1, n
     d(i)%df(i) = one
    end do

  end subroutine ddt_seed_vector

!*******************************************************************************
! Assignment
!*******************************************************************************
  pure elemental subroutine d_assign_d(d1,d2)

   type(derivative_data_type_df5),  intent(out) :: d1
   type(derivative_data_type_df5),  intent( in) :: d2

    d1%f  = d2%f
    d1%df = d2%df

  end subroutine d_assign_d
!-------------------------------------------------------------------------------
  pure elemental subroutine d_assign_r(d,r)

   type(derivative_data_type_df5),  intent(out) :: d
   real(my_precision)        ,  intent( in) :: r

    d%f  = r
    d%df = zero

  end subroutine d_assign_r
!-------------------------------------------------------------------------------
  pure elemental subroutine r_assign_d(r,d)

   real(my_precision)        ,  intent(out) :: r
   type(derivative_data_type_df5),  intent( in) :: d

    r = d%f

  end subroutine r_assign_d

!*******************************************************************************
! Addition
!*******************************************************************************
  pure elemental function d_plus_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   type(derivative_data_type_df5)             :: d_plus_d

    d_plus_d%f  = d1%f  + d2%f
    d_plus_d%df = d1%df + d2%df

  end function d_plus_d
!-------------------------------------------------------------------------------
  pure elemental function d_plus_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   type(derivative_data_type_df5)             :: d_plus_r

    d_plus_r%f  = d%f  + r
    d_plus_r%df = d%df

  end function d_plus_r
!-------------------------------------------------------------------------------
  pure elemental function r_plus_d(r, d)

   type(derivative_data_type_df5),  intent(in) :: d
   real(my_precision)        ,  intent(in) :: r
   type(derivative_data_type_df5)              :: r_plus_d

    r_plus_d%f  = d%f  + r
    r_plus_d%df = d%df

  end function r_plus_d

!*******************************************************************************
! Subtraction
!*******************************************************************************
  pure elemental function d_minus_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   type(derivative_data_type_df5)             :: d_minus_d

    d_minus_d%f  = d1%f  - d2%f
    d_minus_d%df = d1%df - d2%df

  end function d_minus_d
!-------------------------------------------------------------------------------
  pure elemental function d_minus_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   type(derivative_data_type_df5)             :: d_minus_r

    d_minus_r%f  = d%f - r
    d_minus_r%df = d%df

  end function d_minus_r
!-------------------------------------------------------------------------------
  pure elemental function r_minus_d(r, d)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   type(derivative_data_type_df5)             :: r_minus_d

    r_minus_d%f  = r - d%f
    r_minus_d%df =   - d%df

  end function r_minus_d
!-------------------------------------------------------------------------------
  pure elemental function minus_d(d)

   type(derivative_data_type_df5), intent(in) :: d
   type(derivative_data_type_df5)             :: minus_d

    minus_d%f  = -d%f
    minus_d%df = -d%df

  end function minus_d

!*******************************************************************************
! Multiplication
!*******************************************************************************
  pure elemental function d_times_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   type(derivative_data_type_df5)             :: d_times_d

    d_times_d%f  = d1%f*d2%f
    d_times_d%df = d1%f*d2%df + d2%f*d1%df

  end function d_times_d
!-------------------------------------------------------------------------------
  pure elemental function d_times_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   type(derivative_data_type_df5)             :: d_times_r

    d_times_r%f  = r*d%f
    d_times_r%df = r*d%df

  end function d_times_r
!-------------------------------------------------------------------------------
  pure elemental function r_times_d(r, d)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   type(derivative_data_type_df5)             :: r_times_d

    r_times_d%f  = r*d%f
    r_times_d%df = r*d%df

  end function r_times_d

!*******************************************************************************
! Division (NOTE: derivative of r is zero: r'=0)
!*******************************************************************************
  pure elemental function d_divided_by_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   type(derivative_data_type_df5)             :: d_divided_by_d

    d_divided_by_d%f  = d1%f / d2%f
    d_divided_by_d%df = ( d1%df*d2%f - d2%df*d1%f ) / (d2%f*d2%f)

  end function d_divided_by_d
!-------------------------------------------------------------------------------
  pure elemental function d_divided_by_r(d, r)

   type(derivative_data_type_df5),  intent(in) :: d
   real(my_precision)        ,  intent(in) :: r
   type(derivative_data_type_df5)              :: d_divided_by_r

    d_divided_by_r%f  = d%f  / r
    d_divided_by_r%df = d%df / r

  end function d_divided_by_r
!-------------------------------------------------------------------------------
  pure elemental function r_divided_by_d(r, d)

   real(my_precision),    intent(in) :: r
   type(derivative_data_type_df5),  intent(in) :: d
   type(derivative_data_type_df5)              :: r_divided_by_d

    r_divided_by_d%f  = r / d%f
    r_divided_by_d%df = -d%df*r / (d%f*d%f) ! (r/f)'=(r'*f-f'*r)/f^2=-f'*r/f^2

  end function r_divided_by_d

!*******************************************************************************
! Square root
!*******************************************************************************
  pure elemental function sqrt_d(d)
 
   type(derivative_data_type_df5), intent(in) :: d
   type(derivative_data_type_df5)             :: sqrt_d

    sqrt_d%f = sqrt(d%f)

    avoid_zero_denominator : if ( sqrt_d%f > epsilon(one) ) then
      sqrt_d%df = half * d%df / sqrt_d%f
    else
      sqrt_d%df = half * d%df / epsilon(one)
    end if avoid_zero_denominator

  end function sqrt_d

!*******************************************************************************
! Square root
!*******************************************************************************
  pure elemental function tanh_d(d)
 
   type(derivative_data_type_df5), intent(in) :: d
   type(derivative_data_type_df5)             :: tanh_d

    tanh_d%f = tanh(d%f)

    tanh_d%df = (one - tanh(d%f)**2)*d%df

  end function tanh_d

!*******************************************************************************
! Power: Let F=f1^f2. We want F'.
!  1. Set Z=log(f1^f2)=f2*log(f1) -> Differentiate -> Z'=f2'*log(f1)+f2/f1*f1'
!  2. Since Z=log(F), we have Z'=F'/F, so, F'=Z'*F.
!  1. Hence, F' = f1^f2*( f2'*log(f1) + f2/f1*f1' )
!*******************************************************************************
  pure elemental function d_power_d(d1, d2)

  type(derivative_data_type_df5), intent(in) :: d1, d2
  type(derivative_data_type_df5)             :: d_power_d

   d_power_d%f  = d1%f**d2%f
   d_power_d%df = d1%f**d2%f * ( d2%df*log(d1%f) + d2%f*d1%df/d1%f  )
 
  end function d_power_d
!-------------------------------------------------------------------------------
  pure elemental function d_power_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   type(derivative_data_type_df5)             :: d_power_r

    d_power_r%f  = d%f**r
    d_power_r%df = r*d%f**(r-one)*d%df

  end function d_power_r
!-------------------------------------------------------------------------------
  pure elemental function r_power_d(r, d)

  real(my_precision)        , intent(in) :: r
  type(derivative_data_type_df5), intent(in) :: d
  type(derivative_data_type_df5)             :: r_power_d

   r_power_d%f  = r**d%f
   r_power_d%df = r**d%f * d%df*log(r)
 
  end function r_power_d
!-------------------------------------------------------------------------------
  pure elemental function d_power_i(d, i)

   type(derivative_data_type_df5), intent(in) :: d
   integer                   , intent(in) :: i
   type(derivative_data_type_df5)             :: d_power_i

    d_power_i%f  = d%f**i
    d_power_i%df = real(i,my_precision)*d%f**(i-1)*d%df

  end function d_power_i

!*******************************************************************************
! Absolute value fucntion
!*******************************************************************************
  pure elemental function abs_d(d)

   type(derivative_data_type_df5), intent(in) :: d
   type(derivative_data_type_df5)             :: abs_d

    if ( d%f < zero ) then
      abs_d%f  = -d%f
      abs_d%df = -d%df
    else
      abs_d%f  =  d%f
      abs_d%df =  d%df
    endif

  end function abs_d

!*******************************************************************************
! Sign function
!*******************************************************************************
  pure elemental function sign_rd(r,d)

   real(my_precision)        , intent(in) :: r
   type(derivative_data_type_df5), intent(in) :: d
   type(derivative_data_type_df5)             :: sign_rd

    if ( d%f >= zero) then
      sign_rd%f  =  abs(r)
      sign_rd%df =  zero
    else
      sign_rd%f  = -abs(r)
      sign_rd%df =  zero
    end if

  end function sign_rd

!*******************************************************************************
! Max
!*******************************************************************************
  pure elemental function max_of_d_and_d(d1,d2)
 
   type(derivative_data_type_df5), intent(in) :: d1, d2
   type(derivative_data_type_df5)             :: max_of_d_and_d

    if ( d1%f > d2%f ) then
      max_of_d_and_d = d1
    else
      max_of_d_and_d = d2
    endif

  end function max_of_d_and_d
!-------------------------------------------------------------------------------
  pure elemental function max_of_d_and_r(d,r)

   real(my_precision)        , intent(in) :: r
   type(derivative_data_type_df5), intent(in) :: d
   type(derivative_data_type_df5)             :: max_of_d_and_r

    if ( r > d%f ) then
     max_of_d_and_r%f  = r
     max_of_d_and_r%df = zero
    else
     max_of_d_and_r = d
    endif

  end function max_of_d_and_r
!-------------------------------------------------------------------------------
  pure elemental function max_of_r_and_d(r,d2)

   real(my_precision)        , intent(in) :: r
   type(derivative_data_type_df5), intent(in) :: d2
   type(derivative_data_type_df5)             :: max_of_r_and_d

    if ( r > d2%f ) then
     max_of_r_and_d%f  = r
     max_of_r_and_d%df = zero
    else
      max_of_r_and_d = d2
    endif

  end function max_of_r_and_d

!*******************************************************************************
! Min
!*******************************************************************************
  pure elemental function min_of_d_and_d(d1,d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   type(derivative_data_type_df5)             :: min_of_d_and_d

    if ( d1%f < d2%f ) then
      min_of_d_and_d = d1
    else
      min_of_d_and_d = d2
    endif

  end function min_of_d_and_d
!-------------------------------------------------------------------------------
  pure elemental function min_of_d_and_r(d,r)

   real(my_precision)        , intent(in) :: r
   type(derivative_data_type_df5), intent(in) :: d
   type(derivative_data_type_df5)             :: min_of_d_and_r

    if ( r < d%f ) then
     min_of_d_and_r%f  = r
     min_of_d_and_r%df = zero
    else
     min_of_d_and_r = d
    endif

  end function min_of_d_and_r
!-------------------------------------------------------------------------------
  pure elemental function min_of_r_and_d(r,d2)

   real(my_precision)        , intent(in) :: r
   type(derivative_data_type_df5), intent(in) :: d2
   type(derivative_data_type_df5)             :: min_of_r_and_d

    if ( r < d2%f ) then
     min_of_r_and_d%f  = r
     min_of_r_and_d%df = zero
    else
      min_of_r_and_d = d2
    endif

  end function min_of_r_and_d

!*******************************************************************************
! Less than (Logical)
!*******************************************************************************
  pure elemental function d_less_than_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   logical                                :: d_less_than_d

    d_less_than_d = ( d1%f < d2%f )

  end function d_less_than_d
!-------------------------------------------------------------------------------
  pure elemental function d_less_than_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   logical                                :: d_less_than_r

    d_less_than_r = ( d%f < r )

  end function d_less_than_r
!-------------------------------------------------------------------------------
  pure elemental function r_less_than_d(r, d)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)                  , intent(in) :: r
   logical                                :: r_less_than_d

    r_less_than_d = ( r < d%f )

  end function r_less_than_d

!*******************************************************************************
! Less than or equal (Logical)
!*******************************************************************************
  pure elemental function d_less_than_equal_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   logical                                :: d_less_than_equal_d

    d_less_than_equal_d = ( d1%f <= d2%f )

  end function d_less_than_equal_d
!-------------------------------------------------------------------------------
  pure elemental function d_less_than_equal_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)                  , intent(in) :: r
   logical                                :: d_less_than_equal_r

    d_less_than_equal_r = ( d%f <= r )

  end function d_less_than_equal_r
!-------------------------------------------------------------------------------
  pure elemental function r_less_than_equal_d(r, d)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   logical                                :: r_less_than_equal_d

    r_less_than_equal_d = ( r <= d%f )

  end function r_less_than_equal_d

!*******************************************************************************
! Greater than (Logical)
!*******************************************************************************
  pure elemental function d_greater_than_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   logical                                :: d_greater_than_d

    d_greater_than_d = ( d1%f > d2%f )

  end function d_greater_than_d
!-------------------------------------------------------------------------------
  pure elemental function d_greater_than_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   logical                                :: d_greater_than_r

    d_greater_than_r = ( d%f > r )

  end function d_greater_than_r
!-------------------------------------------------------------------------------
  pure elemental function r_greater_than_d(r, d)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   logical                                :: r_greater_than_d

    r_greater_than_d = ( r > d%f )

  end function r_greater_than_d

!*******************************************************************************
! Greater than or equal (Logical)
!*******************************************************************************
  pure elemental function d_greater_than_equal_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   logical                                :: d_greater_than_equal_d

    d_greater_than_equal_d = ( d1%f >= d2%f )

  end function d_greater_than_equal_d
!-------------------------------------------------------------------------------
  pure elemental function d_greater_than_equal_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   logical                                :: d_greater_than_equal_r

    d_greater_than_equal_r = ( d%f >= r )

  end function d_greater_than_equal_r
!-------------------------------------------------------------------------------
  pure elemental function r_greater_than_equal_d(r, d)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   logical                                :: r_greater_than_equal_d

    r_greater_than_equal_d = ( r >= d%f )

  end function r_greater_than_equal_d

!*******************************************************************************
! Equal (Logical)
!*******************************************************************************
  pure elemental function d_equal_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   logical                                :: d_equal_d

    d_equal_d = ( abs(d1%f-d2%f) < epsilon(d1%f) )

  end function d_equal_d
!-------------------------------------------------------------------------------
  pure elemental function d_equal_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   logical                                :: d_equal_r

    d_equal_r = ( abs(d%f-r) < epsilon(d%f) )

  end function d_equal_r
!-------------------------------------------------------------------------------
  pure elemental function r_equal_d(r, d2)

   type(derivative_data_type_df5), intent(in) :: d2
   real(my_precision)        , intent(in) :: r
   logical                                :: r_equal_d

    r_equal_d = ( abs(r-d2%f) < epsilon(r) )

  end function r_equal_d

 end module derivative_data_df5

