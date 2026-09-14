module peec_linalg
  use peec_base
  implicit none
  private
  public :: real_lu, complex_lu, factor_real, factor_complex, solve_real, solve_complex
  type real_lu
    real(dp), allocatable :: a(:,:), row(:), col(:)
    integer, allocatable :: piv(:)
  end type
  type complex_lu
    complex(dp), allocatable :: a(:,:)
    real(dp), allocatable :: row(:), col(:)
    integer, allocatable :: piv(:)
  end type
  interface
    subroutine dgetrf(m,n,a,lda,ipiv,info)
      import dp
      integer :: m,n,lda,ipiv(*),info
      real(dp) :: a(lda,*)
    end subroutine
    subroutine dgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
      import dp
      character :: trans
      integer :: n,nrhs,lda,ipiv(*),ldb,info
      real(dp) :: a(lda,*),b(ldb,*)
    end subroutine
    subroutine zgetrf(m,n,a,lda,ipiv,info)
      import dp
      integer :: m,n,lda,ipiv(*),info
      complex(dp) :: a(lda,*)
    end subroutine
    subroutine zgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
      import dp
      character :: trans
      integer :: n,nrhs,lda,ipiv(*),ldb,info
      complex(dp) :: a(lda,*),b(ldb,*)
    end subroutine
  end interface
contains
  subroutine factor_real(a,lu,info)
    real(dp), intent(in) :: a(:,:)
    type(real_lu), intent(out) :: lu
    integer, intent(out) :: info
    integer :: n,j
    n=size(a,1); info=-1
    if(n==0.or.size(a,2)/=n) return
    if(.not.all(ieee_is_finite(a))) return
    lu%a=a
    lu%row=maxval(abs(a),dim=2)
    if(any(lu%row<=0)) return
    do j=1,n
      lu%a(:,j)=lu%a(:,j)/lu%row
    end do
    lu%col=maxval(abs(lu%a),dim=1)
    if(any(lu%col<=0)) return
    do j=1,n
      lu%a(:,j)=lu%a(:,j)/lu%col(j)
    end do
    allocate(lu%piv(n))
    call dgetrf(n,n,lu%a,n,lu%piv,info)
    if(info/=0) return
    do j=1,n
      if(abs(lu%a(j,j))<=1e-14_dp) info=j
    end do
  end subroutine
  subroutine factor_complex(a,lu,info)
    complex(dp), intent(in) :: a(:,:)
    type(complex_lu), intent(out) :: lu
    integer, intent(out) :: info
    integer :: n,j
    n=size(a,1); info=-1
    if(n==0.or.size(a,2)/=n) return
    if(.not.all(ieee_is_finite(abs(a)))) return
    lu%a=a
    lu%row=maxval(abs(a),dim=2)
    if(any(lu%row<=0)) return
    do j=1,n
      lu%a(:,j)=lu%a(:,j)/lu%row
    end do
    lu%col=maxval(abs(lu%a),dim=1)
    if(any(lu%col<=0)) return
    do j=1,n
      lu%a(:,j)=lu%a(:,j)/lu%col(j)
    end do
    allocate(lu%piv(n))
    call zgetrf(n,n,lu%a,n,lu%piv,info)
    if(info/=0) return
    do j=1,n
      if(abs(lu%a(j,j))<=1e-14_dp) info=j
    end do
  end subroutine
  subroutine solve_real(lu,b,x)
    type(real_lu), intent(in) :: lu
    real(dp), intent(in) :: b(:)
    real(dp), intent(out) :: x(:)
    real(dp) :: work(size(b),1)
    integer :: n,info
    n=size(b)
    call require(allocated(lu%piv),'LU is not factorized')
    call require(size(x)==n.and.size(lu%a,1)==n,'LU solve dimensions')
    work(:,1)=b/lu%row
    call dgetrs('N',n,1,lu%a,n,lu%piv,work,n,info)
    call require(info==0,'LAPACK dgetrs failed')
    x=work(:,1)/lu%col
    call require(all(ieee_is_finite(x)),'nonfinite real solution')
  end subroutine
  subroutine solve_complex(lu,b,x)
    type(complex_lu), intent(in) :: lu
    complex(dp), intent(in) :: b(:)
    complex(dp), intent(out) :: x(:)
    complex(dp) :: work(size(b),1)
    integer :: n,info
    n=size(b)
    call require(allocated(lu%piv),'LU is not factorized')
    call require(size(x)==n.and.size(lu%a,1)==n,'LU solve dimensions')
    work(:,1)=b/lu%row
    call zgetrs('N',n,1,lu%a,n,lu%piv,work,n,info)
    call require(info==0,'LAPACK zgetrs failed')
    x=work(:,1)/lu%col
    call require(all(ieee_is_finite(abs(x))),'nonfinite complex solution')
  end subroutine
end module
