module peec_model
  use peec_base
  use peec_mesh
  use peec_quadrature
  use peec_linalg
  implicit none
  type model_type
    type(mesh_type) :: mesh
    real(dp), allocatable :: l(:,:),p(:,:),a(:,:),r(:),pat(:,:),f(:),s(:,:),as(:,:)
  end type
contains
  real(dp) function region_kernel(a,b) result(v)
    type(patch_list), intent(in) :: a,b
    integer :: i,j
    v=0
    do i=1,size(a%p,3)
      do j=1,size(b%p,3)
        v=v+kernel(a%p(:,:,i),b%p(:,:,j))
      end do
    end do
  end function
  subroutine assemble(path,m,threads)
    character(*), intent(in) :: path
    type(model_type), intent(out) :: m
    integer, intent(in) :: threads
    integer :: ne,nv,i,j
    real(dp) :: d,v,wa,wb
    call read_mesh(path,m%mesh)
    ne=size(m%mesh%edges,2); nv=size(m%mesh%xyz,2)
    allocate(m%l(ne,ne),m%p(nv,nv),m%a(ne,nv),m%r(ne))
    m%a=0; m%r=0
    do i=1,ne
      m%a(i,m%mesh%edges(1,i))=-1; m%a(i,m%mesh%edges(2,i))=1
    end do
    print '(a)', 'Assembling L ...'
    !$omp parallel do schedule(dynamic) num_threads(threads) private(j,d,v,wa,wb)
    do i=1,ne
      do j=i,ne
        d=dot_product(m%mesh%direction(:,i),m%mesh%direction(:,j))
        v=0
        if(abs(d)>=1e-10_dp) then
          wa=m%mesh%branches(i)%area/m%mesh%length(i)
          wb=m%mesh%branches(j)%area/m%mesh%length(j)
          v=1e-7_dp*d*region_kernel(m%mesh%branches(i),m%mesh%branches(j))/(wa*wb)
        end if
        m%l(i,j)=v; m%l(j,i)=v
      end do
    end do
    !$omp end parallel do
    call require(all(ieee_is_finite(m%l)),'inductance quadrature did not converge')
    print '(a)', 'Assembling P ...'
    !$omp parallel do schedule(dynamic) num_threads(threads) private(j,v)
    do i=1,nv
      do j=i,nv
        v=8.987551792261171e9_dp*region_kernel(m%mesh%nodes(i),m%mesh%nodes(j)) &
          /(m%mesh%nodes(i)%area*m%mesh%nodes(j)%area)
        m%p(i,j)=v; m%p(j,i)=v
      end do
    end do
    !$omp end parallel do
    call require(all(ieee_is_finite(m%p)),'potential quadrature did not converge')
    call operators(m)
  end subroutine
  subroutine operators(m)
    type(model_type), intent(inout) :: m
    integer :: nv,j
    nv=size(m%p,1)
    call require(size(m%p,2)==nv.and.size(m%a,2)==nv,'PEEC matrix dimensions')
    allocate(m%f(nv),m%s(nv,nv))
    do j=1,nv
      call require(m%p(j,j)>0,'P diagonal must be positive')
      m%f(j)=1/m%p(j,j)
      m%s(:,j)=m%p(:,j)*m%f(j)
    end do
    m%pat=matmul(m%p,transpose(m%a))
    m%as=matmul(m%a,m%s)
  end subroutine
  subroutine harmonic(m,frequency,u,current,phi,charge,backward_error)
    type(model_type), intent(in) :: m
    real(dp), intent(in) :: frequency
    complex(dp), intent(in) :: u(:)
    complex(dp), allocatable, intent(out) :: current(:),phi(:),charge(:)
    real(dp), intent(out) :: backward_error
    type(complex_lu) :: lu
    type(real_lu) :: plu
    complex(dp), allocatable :: h(:,:),b(:),x(:)
    real(dp), allocatable :: qr(:),qi(:)
    real(dp) :: omega,den
    integer :: ne,nv,n,i,info
    ne=size(m%l,1); nv=size(m%p,1); n=ne+nv
    call require(frequency>0.and.size(u)==ne,'invalid harmonic arguments')
    omega=2*pi*frequency
    allocate(h(n,n),b(n),x(n),qr(nv),qi(nv)); h=0; b=0
    h(:ne,:ne)=iu*omega*m%l
    do i=1,ne
      h(i,i)=h(i,i)+m%r(i)
    end do
    h(:ne,ne+1:)=cmplx(m%a,0.0_dp,dp)
    h(ne+1:,:ne)=cmplx(-m%pat,0.0_dp,dp)
    do i=1,nv
      h(ne+i,ne+i)=iu*omega
    end do
    b(:ne)=-u
    call factor_complex(h,lu,info)
    call require(info==0,'singular harmonic matrix')
    call solve_complex(lu,b,x)
    den=maxval(sum(abs(h),dim=2))*maxval(abs(x))+maxval(abs(b))
    backward_error=maxval(abs(matmul(h,x)-b))/max(den,tiny(1.0_dp))
    call require(ieee_is_finite(backward_error).and.backward_error<=1e-10_dp,'harmonic residual too large')
    current=x(:ne); phi=x(ne+1:)
    call factor_real(m%p,plu,info)
    call require(info==0,'singular potential matrix')
    call solve_real(plu,real(phi,dp),qr)
    call solve_real(plu,aimag(phi),qi)
    charge=cmplx(qr,qi,dp)
  end subroutine
end module
