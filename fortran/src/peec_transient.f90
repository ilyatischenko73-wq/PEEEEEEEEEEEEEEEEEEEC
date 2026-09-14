module peec_transient
  use peec_physics
  implicit none
  type transient_type
    integer :: ne=0,nv=0,ns=0,n=0
    real(dp) :: theta=1
    real(dp), allocatable :: matrix(:,:),history(:,:),x(:)
    type(real_lu) :: lu
  end type
contains
  subroutine build_transient(m,cells,shunt,dt,order,t)
    type(model_type), intent(in) :: m
    type(slot_cell), intent(in) :: cells(:)
    type(shunt_type), intent(in) :: shunt
    real(dp), intent(in) :: dt
    integer, intent(in) :: order
    type(transient_type), intent(out) :: t
    real(dp), allocatable :: mass(:,:),stiff(:,:),op(:),b(:)
    integer :: ne,nv,n,ns,i,j,s,a1,b1,info
    call require(dt>0.and.(order==1.or.order==2),'invalid time integration parameters')
    ne=size(m%l,1); nv=size(m%p,1); ns=size(cells); n=ne+nv+ns
    if(shunt%enabled) n=n+1
    t%ne=ne; t%nv=nv; t%ns=ns; t%n=n
    t%theta=merge(1.0_dp,0.5_dp,order==1)
    allocate(mass(n,n),stiff(n,n),op(nv),b(nv)); mass=0; stiff=0
    mass(:ne,:ne)=m%l
    do i=1,ne
      stiff(i,i)=m%r(i)
    end do
    stiff(:ne,ne+1:ne+nv)=m%as
    stiff(ne+1:ne+nv,:ne)=-transpose(m%a)
    do i=1,nv
      mass(ne+i,ne+i)=m%f(i)
    end do
    do j=1,ns
      a1=cells(j)%a; b1=cells(j)%b; s=ne+nv+j
      call require(a1>=1.and.b1>=1.and.a1<=nv.and.b1<=nv.and.a1/=b1,'invalid slot terminals')
      call require(cells(j)%l>0.and.cells(j)%c>0,'invalid slot L/C')
      op=(m%p(b1,:)-m%p(a1,:))*m%f
      mass(ne+a1,ne+1:ne+nv)=mass(ne+a1,ne+1:ne+nv)-cells(j)%c*op
      mass(ne+b1,ne+1:ne+nv)=mass(ne+b1,ne+1:ne+nv)+cells(j)%c*op
      stiff(ne+a1,s)=1; stiff(ne+b1,s)=-1
      stiff(s,ne+1:ne+nv)=op
      mass(s,s)=cells(j)%l
    end do
    if(shunt%enabled) then
      a1=shunt%a; b1=shunt%b
      call require(a1>=1.and.b1>=1.and.a1<=nv.and.b1<=nv.and.a1/=b1,'invalid shunt terminals')
      b=0; b(a1)=-1; b(b1)=1
      stiff(ne+1:ne+nv,n)=-b
      stiff(n,ne+1:ne+nv)=(m%p(b1,:)-m%p(a1,:))*m%f
      stiff(n,n)=shunt%r
    end if
    t%matrix=mass/dt+t%theta*stiff
    t%history=mass/dt-(1-t%theta)*stiff
    ! The shunt is an algebraic constraint imposed at n+1, not time averaged.
    if(shunt%enabled) then
      t%matrix(n,:)=stiff(n,:); t%history(n,:)=0
    end if
    allocate(t%x(n)); t%x=0
    call factor_real(t%matrix,t%lu,info)
    call require(info==0,'singular transient matrix')
  end subroutine
  subroutine step_transient(t,node_n,node_next,edge_n,edge_next)
    type(transient_type), intent(inout) :: t
    real(dp), intent(in) :: node_n(:),node_next(:),edge_n(:),edge_next(:)
    real(dp) :: rhs(t%n),next(t%n),den,residual
    integer :: ne,nv
    ne=t%ne; nv=t%nv
    rhs=matmul(t%history,t%x)
    rhs(:ne)=rhs(:ne)-(1-t%theta)*edge_n-t%theta*edge_next
    rhs(ne+1:ne+nv)=rhs(ne+1:ne+nv)+(1-t%theta)*node_n+t%theta*node_next
    call solve_real(t%lu,rhs,next)
    den=maxval(sum(abs(t%matrix),dim=2))*maxval(abs(next))+maxval(abs(rhs))
    residual=maxval(abs(matmul(t%matrix,next)-rhs))/max(den,tiny(1.0_dp))
    call require(ieee_is_finite(residual).and.residual<=1e-10_dp,'transient residual too large')
    t%x=next
  end subroutine
  subroutine recover(m,t,phi,charge)
    type(model_type), intent(in) :: m
    type(transient_type), intent(in) :: t
    real(dp), intent(out) :: phi(:),charge(:)
    charge=m%f*t%x(t%ne+1:t%ne+t%nv)
    phi=matmul(m%p,charge)
  end subroutine
end module
