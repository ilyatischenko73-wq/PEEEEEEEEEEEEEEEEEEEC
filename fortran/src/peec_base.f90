module peec_base
  use iso_fortran_env, only: dp => real64, error_unit
  use, intrinsic :: ieee_arithmetic
  implicit none
  real(dp), parameter :: pi=acos(-1.0_dp), c0=299792458.0_dp
  real(dp), parameter :: eps0=8.8541878128e-12_dp, mu0=1.25663706212e-6_dp
  complex(dp), parameter :: iu=(0.0_dp,1.0_dp)
contains
  subroutine require(ok,message)
    logical, intent(in) :: ok
    character(*), intent(in) :: message
    if (.not.ok) then
      write(error_unit,'(a)') 'ERROR: '//message
      error stop 1
    end if
  end subroutine
  pure function cross(a,b) result(c)
    real(dp), intent(in) :: a(3),b(3)
    real(dp) :: c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function
  pure real(dp) function norm(a)
    real(dp), intent(in) :: a(:)
    norm=sqrt(sum(a*a))
  end function
  pure real(dp) function tri_area(p)
    real(dp), intent(in) :: p(3,3)
    tri_area=0.5_dp*norm(cross(p(:,2)-p(:,1),p(:,3)-p(:,1)))
  end function
  pure real(dp) function patch_area(p)
    real(dp), intent(in) :: p(3,4)
    patch_area=tri_area(p(:,[1,2,3]))+tri_area(p(:,[1,3,4]))
  end function
  pure real(dp) function sinc(x)
    real(dp), intent(in) :: x
    if(abs(x)<1e-4_dp) then
      sinc=1-x*x/6+x**4/120
    else
      sinc=sin(x)/x
    end if
  end function
  pure real(dp) function pulse(t,k,alpha,beta)
    real(dp), intent(in) :: t,k,alpha,beta
    pulse=0
    if(t>=0) pulse=k*(exp(-alpha*t)-exp(-beta*t))
  end function
  function integer_text(i) result(s)
    integer, intent(in) :: i
    character(:), allocatable :: s
    character(40) :: b
    write(b,'(i0)') i
    s=trim(b)
  end function
end module
