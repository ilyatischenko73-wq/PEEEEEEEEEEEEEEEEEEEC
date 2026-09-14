program test_core
  use peec_io
  implicit none
  integer :: passed=0
  call test_quadrature()
  call test_lu()
  call test_harmonic()
  call test_time(1)
  call test_time(2)
  call test_shunt()
  call test_slots()
  call test_aperture()
  write(*,'(a,i0,a)') 'PASS: ',passed,' numerical checks'
contains
  subroutine check(ok,name)
    logical, intent(in) :: ok
    character(*), intent(in) :: name
    call require(ok,'test failed: '//name)
    passed=passed+1
    print '(a)', 'PASS '//name
  end subroutine
  real(dp) function rectangle(a,b) result(v)
    real(dp), intent(in) :: a,b
    real(dp) :: h,i00,ix,iy,ixy
    h=sqrt(a*a+b*b)
    i00=a*asinh(b/a)+b*asinh(a/b)
    ix=(b*h+a*a*asinh(b/a)-b*b)/2
    iy=(a*h+b*b*asinh(a/b)-a*a)/2
    ixy=(h**3-a**3-b**3)/3
    v=4*(a*b*i00-b*ix-a*iy+ixy)
  end function
  subroutine test_quadrature()
    real(dp) :: a(3,4),b(3,4),ref,ab,ba,total,s,e(3),f(3),x(16),w(16)
    integer :: k,j,n
    call check(.not.order_supported(7),'reject unsupported quadrature order')
    do n=2,16
      if(.not.order_supported(n)) cycle
      call gauss_rule(n,x(:n),w(:n))
      call check(abs(sum(w(:n))-2)<1e-14_dp,'Gauss weights '//integer_text(n))
    end do
    do k=0,2
      s=10.0_dp**k
      a=reshape([0.0_dp,0.0_dp,0.0_dp,s,0.0_dp,0.0_dp,s,1.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp],[3,4])
      ref=rectangle(s,1.0_dp)
      call check(abs(kernel(a,a)/ref-1)<=1e-3_dp,'rectangle aspect '//integer_text(10**k))
    end do
    a=reshape([0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp],[3,4])
    b=a; b(1,:)=b(1,:)+1
    ab=kernel(a,b); ba=kernel(b,a)
    total=kernel(a,a)+kernel(b,b)+ab+ba
    call check(abs(total/rectangle(2.0_dp,1.0_dp)-1)<1e-3_dp,'touching partition')
    call check(abs(ab-ba)<1e-12_dp,'quadrature reciprocity')
    b=a; b(3,:)=1e-7_dp
    call check(abs(kernel(a,b)/rectangle(1.0_dp,1.0_dp)-1)<1e-3_dp,'near-to-self limit')
    b(3,:)=100
    s=kernel(a,b)
    call check(s>=1/sqrt(10002.0_dp).and.s<=0.01_dp,'far-field bounds')
    e=[0.6_dp,0.8_dp,0.0_dp]; f=[-0.48_dp,0.36_dp,0.8_dp]
    do j=1,4
      b(:,j)=3.7_dp*(a(1,j)*e+a(2,j)*f)+[1.0_dp,2.0_dp,3.0_dp]
    end do
    call check(abs(kernel(b,b)/(rectangle(1.0_dp,1.0_dp)*3.7_dp**3)-1)<1e-3_dp,'rigid transform and scale')
  end subroutine
  subroutine test_lu()
    type(real_lu) :: rlu
    type(complex_lu) :: clu
    real(dp) :: a(3,3),exact(3),x(3),b(3)
    complex(dp) :: z(2,2),zx(2),ze(2),zb(2)
    integer :: info
    a=transpose(reshape([0.0_dp,2.0_dp,1.0_dp,1.0_dp,-2.0_dp,-3.0_dp,2.0_dp,3.0_dp,1.0_dp],[3,3]))
    exact=[1.0_dp,-2.0_dp,0.5_dp]; b=matmul(a,exact)
    call factor_real(a,rlu,info); call check(info==0,'real LU pivoting factorization')
    call solve_real(rlu,b,x)
    call check(maxval(abs(x-exact))<1e-11_dp,'real LU solution')
    z=transpose(reshape([(0.0_dp,0.0_dp),(2.0_dp,1.0_dp),(1.0_dp,-1.0_dp),(3.0_dp,0.5_dp)],[2,2]))
    ze=[(1.0_dp,2.0_dp),(-0.5_dp,0.75_dp)]; zb=matmul(z,ze)
    call factor_complex(z,clu,info); call check(info==0,'complex LU pivoting factorization')
    call solve_complex(clu,zb,zx)
    call check(maxval(abs(zx-ze))<1e-11_dp,'complex LU solution')
    z=reshape([(1.0_dp,1.0_dp),(2.0_dp,2.0_dp),(2.0_dp,2.0_dp),(4.0_dp,4.0_dp)],[2,2])
    call factor_complex(z,clu,info); call check(info/=0,'singular complex LU rejected')
  end subroutine
  subroutine oscillator_model(m)
    type(model_type), intent(out) :: m
    m%l=reshape([1.0_dp],[1,1]); m%r=[0.0_dp]
    m%a=reshape([-1.0_dp,1.0_dp],[1,2])
    m%p=reshape([1.0_dp,0.25_dp,0.25_dp,1.0_dp],[2,2])
    call operators(m)
  end subroutine
  real(dp) function energy(m,t,cells) result(v)
    type(model_type), intent(in) :: m
    type(transient_type), intent(in) :: t
    type(slot_cell), intent(in) :: cells(:)
    real(dp) :: q(t%nv),phi(t%nv),current(t%ne)
    integer :: j
    current=t%x(:t%ne)
    call recover(m,t,phi,q)
    v=0.5_dp*(dot_product(current,matmul(m%l,current))+dot_product(q,phi))
    do j=1,size(cells)
      v=v+0.5_dp*cells(j)%l*t%x(t%ne+t%nv+j)**2 &
        +0.5_dp*cells(j)%c*(phi(cells(j)%a)-phi(cells(j)%b))**2
    end do
  end function
  subroutine test_harmonic()
    type(model_type) :: m
    complex(dp), allocatable :: current(:),phi(:),q(:)
    complex(dp) :: u(3)
    real(dp) :: w,den,f,err,worst,v(3),e(3)
    integer :: i,k
    allocate(m%l(3,3),m%p(3,3),m%r(3)); m%l=0; m%p=0; m%r=0
    do i=1,3
      m%l(i,i)=1e-7_dp; m%p(i,i)=1e11_dp
    end do
    m%a=transpose(reshape([-1.0_dp,1.0_dp,0.0_dp,0.0_dp,-1.0_dp,1.0_dp,1.0_dp,0.0_dp,-1.0_dp],[3,3]))
    call operators(m); worst=0; e=[1.0_dp,-1.0_dp,0.0_dp]
    do k=0,2
      f=10.0_dp**(3*k); w=2*pi*f; den=3e11_dp-w*w*1e-7_dp
      u=cmplx(e,-w*1e-7_dp*0.003_dp,dp)
      call harmonic(m,f,u,current,phi,q,err)
      v=[1e11_dp,-2e11_dp,1e11_dp]/den
      worst=max(worst,maxval(abs(current-cmplx(0.003_dp,-w*e/den,dp)))/0.003_dp,maxval(abs(phi-v)))
    end do
    call check(worst<1e-6_dp,'analytic harmonic loop over 1 Hz to 1 MHz')
  end subroutine
  subroutine test_time(order)
    integer, intent(in) :: order
    type(model_type) :: m
    type(transient_type) :: t
    type(shunt_type) :: shunt
    type(slot_cell) :: empty(0)
    real(dp) :: errors(3),dt,omega,ref_v,previous,en,ratio,expected
    real(dp) :: zn(2)=0,ze(1)=0
    integer :: k,steps,j
    call oscillator_model(m); omega=sqrt(1.5_dp); ref_v=sin(omega)/omega
    do k=1,3
      steps=10*2**(k-1); dt=1.0_dp/steps
      call build_transient(m,empty,shunt,dt,order,t); t%x(1)=1
      previous=0.5_dp
      do j=1,steps
        call step_transient(t,zn,zn,ze,ze)
        en=energy(m,t,empty)
        call require(abs(t%x(2)+t%x(3))<1e-11_dp,'oscillator charge conservation')
        if(order==2) call require(abs(en-0.5_dp)<1e-11_dp,'trapezoidal energy conservation')
        if(order==1) call require(en<=previous+1e-11_dp,'BE energy monotonicity')
        previous=en
      end do
      errors(k)=max(abs(t%x(1)-cos(omega)),abs(t%x(2)+ref_v),abs(t%x(3)-ref_v))
    end do
    expected=merge(2.0_dp,4.0_dp,order==1)
    do k=2,3
      ratio=errors(k-1)/errors(k)
      call check(abs(ratio-expected)<0.1_dp*expected,'time convergence order '//integer_text(order))
    end do
  end subroutine
  subroutine test_shunt()
    type(model_type) :: m
    type(transient_type) :: t
    type(shunt_type) :: shunt
    type(slot_cell) :: empty(0)
    real(dp) :: zn(2)=0,ze(1)=0,en,prev,old_i,balance,ohm,dt
    integer :: j
    call oscillator_model(m)
    shunt=shunt_type(1,2,2.0_dp,.true.); dt=0.025_dp
    call build_transient(m,empty,shunt,dt,2,t)
    t%x=[0.0_dp,1.0_dp,-1.0_dp,0.75_dp]; balance=0; ohm=0
    do j=1,40
      prev=energy(m,t,empty); old_i=t%x(4)
      call step_transient(t,zn,zn,ze,ze); en=energy(m,t,empty)
      balance=max(balance,abs(en-prev+dt*2*(0.5_dp*(old_i+t%x(4)))**2))
      ohm=max(ohm,abs(0.75_dp*(t%x(2)-t%x(3))-2*t%x(4)))
    end do
    call check(balance<1e-11_dp.and.ohm<1e-11_dp,'shunt energy balance and Ohm law')
  end subroutine
  subroutine test_slots()
    type(model_type) :: m
    type(transient_type) :: t
    type(shunt_type) :: shunt
    type(slot_cell) :: cells(1)
    real(dp) :: zn(2)=0,ze(1)=0,lp,cp
    integer :: j
    call oscillator_model(m)
    cells(1)=slot_cell(1,2,0.5_dp,0.4_dp,0.1_dp)
    call build_transient(m,cells,shunt,0.025_dp,2,t); t%x(1)=1
    do j=1,40
      call step_transient(t,zn,zn,ze,ze)
      call require(abs(energy(m,t,cells)-0.5_dp)<1e-11_dp,'slot energy conservation')
    end do
    call check(.true.,'LC slot transient energy')
    call slot_line(0.1_dp,2.0_dp,0.002_dp,1.0_dp,1.0_dp,lp,cp)
    call check(abs(lp*cp/(4e-7_dp*pi*eps0)-1)<1e-12_dp,'slot transmission-line identity')
  end subroutine
  subroutine test_aperture()
    type(mesh_type) :: closed,opened
    type(coupling_type) :: coupling
    integer :: info
    real(dp) :: source(1)
    closed%xyz=reshape([0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp],[3,2])
    opened%xyz=reshape([0.0_dp,0.0_dp,0.0_dp],[3,1])
    closed%edges=reshape([1,2],[2,1])
    call build_coupling(closed,opened,reshape([1,1],[2,1]),[1],coupling,info)
    call check(info==0,'valid aperture mapping')
    call apply_coupling(coupling,[2.0_dp],source)
    call check(abs(source(1)+2)<1e-11_dp,'aperture current orientation')
    opened%xyz(1,1)=100
    call build_coupling(closed,opened,reshape([1,1],[2,1]),[1],coupling,info)
    call check(info/=0.and..not.allocated(coupling%edge),'invalid aperture geometry rejected cleanly')
  end subroutine
end program
