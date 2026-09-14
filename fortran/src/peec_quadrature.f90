module peec_quadrature
  use peec_base
  implicit none
  private
  public :: kernel, gauss_rule, order_supported, patch_phase
  type triangle_potential
    real(dp) :: p(3,3), unit(3,3), outward(3,3), length(3), normal(3), twice_area
  end type
contains
  pure logical function order_supported(n)
    integer, intent(in) :: n
    order_supported=any(n==[2,3,4,5,6,8,10,12,16])
  end function
  subroutine gauss_rule(n,x,w)
    integer, intent(in) :: n
    real(dp), intent(out) :: x(n),w(n)
    call require(order_supported(n),'unsupported Gauss order')
    select case(n)
    case(2)
      x=[ &
        -0.57735026918962576451_dp, &
        0.57735026918962576451_dp]
      w=[ &
        1.0_dp, &
        1.0_dp]
    case(3)
      x=[ &
        -0.77459666924148337704_dp, &
        0.0_dp, &
        0.77459666924148337704_dp]
      w=[ &
        0.55555555555555555556_dp, &
        0.88888888888888888889_dp, &
        0.55555555555555555556_dp]
    case(4)
      x=[ &
        -0.86113631159405257522_dp, &
        -0.33998104358485626480_dp, &
        0.33998104358485626480_dp, &
        0.86113631159405257522_dp]
      w=[ &
        0.34785484513745385737_dp, &
        0.65214515486254614263_dp, &
        0.65214515486254614263_dp, &
        0.34785484513745385737_dp]
    case(5)
      x=[ &
        -0.90617984593866399280_dp, &
        -0.53846931010568309104_dp, &
        0.0_dp, &
        0.53846931010568309104_dp, &
        0.90617984593866399280_dp]
      w=[ &
        0.23692688505618908751_dp, &
        0.47862867049936646804_dp, &
        0.56888888888888888889_dp, &
        0.47862867049936646804_dp, &
        0.23692688505618908751_dp]
    case(6)
      x=[ &
        -0.93246951420315202781_dp, &
        -0.66120938646626451366_dp, &
        -0.23861918608319690863_dp, &
        0.23861918608319690863_dp, &
        0.66120938646626451366_dp, &
        0.93246951420315202781_dp]
      w=[ &
        0.17132449237917034504_dp, &
        0.36076157304813860757_dp, &
        0.46791393457269104739_dp, &
        0.46791393457269104739_dp, &
        0.36076157304813860757_dp, &
        0.17132449237917034504_dp]
    case(8)
      x=[ &
        -0.96028985649753623168_dp, &
        -0.79666647741362673959_dp, &
        -0.52553240991632898582_dp, &
        -0.18343464249564980494_dp, &
        0.18343464249564980494_dp, &
        0.52553240991632898582_dp, &
        0.79666647741362673959_dp, &
        0.96028985649753623168_dp]
      w=[ &
        0.10122853629037625915_dp, &
        0.22238103445337447054_dp, &
        0.31370664587788728734_dp, &
        0.36268378337836198297_dp, &
        0.36268378337836198297_dp, &
        0.31370664587788728734_dp, &
        0.22238103445337447054_dp, &
        0.10122853629037625915_dp]
    case(10)
      x=[ &
        -0.97390652851717172008_dp, &
        -0.86506336668898451073_dp, &
        -0.67940956829902440623_dp, &
        -0.43339539412924719080_dp, &
        -0.14887433898163121088_dp, &
        0.14887433898163121088_dp, &
        0.43339539412924719080_dp, &
        0.67940956829902440623_dp, &
        0.86506336668898451073_dp, &
        0.97390652851717172008_dp]
      w=[ &
        0.06667134430868813759_dp, &
        0.14945134915058059315_dp, &
        0.21908636251598204399_dp, &
        0.26926671930999635509_dp, &
        0.29552422471475287017_dp, &
        0.29552422471475287017_dp, &
        0.26926671930999635509_dp, &
        0.21908636251598204399_dp, &
        0.14945134915058059315_dp, &
        0.06667134430868813759_dp]
    case(12)
      x=[ &
        -0.98156063424671925069_dp, &
        -0.90411725637047485668_dp, &
        -0.76990267419430468704_dp, &
        -0.58731795428661744730_dp, &
        -0.36783149899818019375_dp, &
        -0.12523340851146891547_dp, &
        0.12523340851146891547_dp, &
        0.36783149899818019375_dp, &
        0.58731795428661744730_dp, &
        0.76990267419430468704_dp, &
        0.90411725637047485668_dp, &
        0.98156063424671925069_dp]
      w=[ &
        0.04717533638651182719_dp, &
        0.10693932599531843096_dp, &
        0.16007832854334622633_dp, &
        0.20316742672306592175_dp, &
        0.23349253653835480876_dp, &
        0.24914704581340278500_dp, &
        0.24914704581340278500_dp, &
        0.23349253653835480876_dp, &
        0.20316742672306592175_dp, &
        0.16007832854334622633_dp, &
        0.10693932599531843096_dp, &
        0.04717533638651182719_dp]
    case(16)
      x=[ &
        -0.98940093499164993260_dp, &
        -0.94457502307323257608_dp, &
        -0.86563120238783174388_dp, &
        -0.75540440835500303390_dp, &
        -0.61787624440264374845_dp, &
        -0.45801677765722738634_dp, &
        -0.28160355077925891323_dp, &
        -0.09501250983763744019_dp, &
        0.09501250983763744019_dp, &
        0.28160355077925891323_dp, &
        0.45801677765722738634_dp, &
        0.61787624440264374845_dp, &
        0.75540440835500303390_dp, &
        0.86563120238783174388_dp, &
        0.94457502307323257608_dp, &
        0.98940093499164993260_dp]
      w=[ &
        0.02715245941175409485_dp, &
        0.06225352393864789286_dp, &
        0.09515851168249278481_dp, &
        0.12462897125553387205_dp, &
        0.14959598881657673208_dp, &
        0.16915651939500253819_dp, &
        0.18260341504492358887_dp, &
        0.18945061045506849629_dp, &
        0.18945061045506849629_dp, &
        0.18260341504492358887_dp, &
        0.16915651939500253819_dp, &
        0.14959598881657673208_dp, &
        0.12462897125553387205_dp, &
        0.09515851168249278481_dp, &
        0.06225352393864789286_dp, &
        0.02715245941175409485_dp]
    end select
  end subroutine

  subroutine tri_map(a,u,v,r,jac)
    real(dp), intent(in) :: a(3,3),u,v
    real(dp), intent(out) :: r(3),jac
    r=a(:,1)+u*(a(:,2)-a(:,1))+(1-u)*v*(a(:,3)-a(:,1))
    jac=2*tri_area(a)*(1-u)
  end subroutine
  subroutine prepare(p,t)
    real(dp), intent(in) :: p(3,3)
    type(triangle_potential), intent(out) :: t
    integer :: i
    t%p=p
    do i=1,3
      t%unit(:,i)=p(:,mod(i,3)+1)-p(:,i)
    end do
    t%normal=cross(t%unit(:,1),t%unit(:,2))
    t%twice_area=norm(t%normal)
    if(t%twice_area==0) return
    t%normal=t%normal/t%twice_area
    do i=1,3
      t%length(i)=norm(t%unit(:,i))
      t%unit(:,i)=t%unit(:,i)/t%length(i)
      t%outward(:,i)=cross(t%unit(:,i),t%normal)
    end do
  end subroutine
  real(dp) function potential_at(t,r) result(value)
    type(triangle_potential), intent(in) :: t
    real(dp), intent(in) :: r(3)
    real(dp) :: v(3,3),radius(3),height,d,along,perp,den,omega
    integer :: i
    do i=1,3
      v(:,i)=t%p(:,i)-r
      radius(i)=norm(v(:,i))
    end do
    height=dot_product(v(:,1),t%normal)
    value=0
    do i=1,3
      d=dot_product(v(:,i),t%outward(:,i))
      if(d==0) cycle
      along=dot_product(v(:,i),t%unit(:,i))
      perp=sqrt(d*d+height*height)
      value=value+d*(asinh((along+t%length(i))/perp)-asinh(along/perp))
    end do
    if(height/=0) then
      den=product(radius)+dot_product(v(:,1),v(:,2))*radius(3) &
        +dot_product(v(:,2),v(:,3))*radius(1)+dot_product(v(:,3),v(:,1))*radius(2)
      omega=2*atan2(dot_product(v(:,1),cross(v(:,2),v(:,3))),den)
      value=value-abs(height*omega)
    end if
  end function
  real(dp) function potential_rule(a,b,order) result(value)
    real(dp), intent(in) :: a(3,3)
    type(triangle_potential), intent(in) :: b
    integer, intent(in) :: order
    real(dp) :: x(order),w(order),r(3),jac
    integer :: i,j
    call gauss_rule(order,x,w)
    value=0
    do i=1,order
      do j=1,order
        call tri_map(a,(x(i)+1)/2,(x(j)+1)/2,r,jac)
        value=value+0.25_dp*w(i)*w(j)*jac*potential_at(b,r)
      end do
    end do
  end function
  recursive real(dp) function adaptive(a,b,order,tol,depth) result(value)
    real(dp), intent(in) :: a(3,3),tol
    type(triangle_potential), intent(in) :: b
    integer, intent(in) :: order,depth
    real(dp) :: low,high,lengths(3),left(3,3),right(3,3)
    integer :: e,j,k,i
    low=potential_rule(a,b,merge(3,4,order<=6))
    high=potential_rule(a,b,order)
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    if(.not.ieee_is_finite(low).or..not.ieee_is_finite(high)) return
    if(abs(high-low)<=tol+2e-5_dp*abs(high)) then
      value=high
      return
    end if
    if(depth>=24) return
    do i=1,3
      lengths(i)=norm(a(:,i)-a(:,mod(i,3)+1))
    end do
    e=maxloc(lengths,dim=1); j=mod(e,3)+1; k=mod(e+1,3)+1
    left(:,1)=a(:,e); left(:,2)=(a(:,e)+a(:,j))/2; left(:,3)=a(:,k)
    right(:,1)=left(:,2); right(:,2)=a(:,j); right(:,3)=a(:,k)
    value=adaptive(left,b,order,tol/2,depth+1)+adaptive(right,b,order,tol/2,depth+1)
  end function
  real(dp) function regular(a,b) result(value)
    real(dp), intent(in) :: a(3,3),b(3,3)
    real(dp) :: x(4),w(4),ra(3),rb(3),ja,jb,d
    integer :: i,j,k,l
    call gauss_rule(4,x,w)
    value=0
    do i=1,4
      do j=1,4
        call tri_map(a,(x(i)+1)/2,(x(j)+1)/2,ra,ja)
        do k=1,4
          do l=1,4
            call tri_map(b,(x(k)+1)/2,(x(l)+1)/2,rb,jb)
            d=norm(ra-rb)
            if(d>0) value=value+w(i)*w(j)*w(k)*w(l)*ja*jb/(16*d)
          end do
        end do
      end do
    end do
  end function
  real(dp) function triangle_pair(a,b,kind) result(value)
    real(dp), intent(in) :: a(3,3),b(3,3)
    integer, intent(in) :: kind
    type(triangle_potential) :: ta,tb
    integer :: i,j,n
    logical :: equal,found
    real(dp) :: scale
    if(kind==0) then
      value=regular(a,b)
      return
    end if
    call prepare(a,ta); call prepare(b,tb)
    value=0
    if(ta%twice_area==0.or.tb%twice_area==0) return
    equal=.true.
    do i=1,3
      found=.false.
      do j=1,3
        if(norm(a(:,i)-b(:,j))<=1e-12_dp) found=.true.
      end do
      equal=equal.and.found
    end do
    if(equal) then
      do i=1,3
        value=value+potential_at(tb,a(:,i))
      end do
      value=value*ta%twice_area/3
      return
    end if
    n=10
    if(kind==2) n=12
    if(kind==3) n=8
    scale=ta%twice_area*tb%twice_area/sqrt(max(ta%twice_area,tb%twice_area))
    value=0.5_dp*(adaptive(a,tb,n,2e-5_dp*scale,0)+adaptive(b,ta,n,2e-5_dp*scale,0))
  end function
  real(dp) function kernel(a,b) result(value)
    real(dp), intent(in) :: a(3,4),b(3,4)
    real(dp) :: gap(3),h
    integer :: i,j,kind,common
    integer, parameter :: ix(3,2)=reshape([1,2,3,1,3,4],[3,2])
    logical :: used(4)
    used=.false.; common=0; h=0
    do i=1,4
      h=max(h,norm(a(:,i)-a(:,mod(i,4)+1)),norm(b(:,i)-b(:,mod(i,4)+1)))
      do j=1,4
        if(used(j)) cycle
        if(norm(a(:,i)-b(:,j))<=1e-12_dp) then
          used(j)=.true.; common=common+1
          exit
        end if
      end do
    end do
    gap=max(0.0_dp,max(minval(a,dim=2)-maxval(b,dim=2),minval(b,dim=2)-maxval(a,dim=2)))
    kind=0
    if(norm(gap)<0.5_dp*h) kind=1
    if(common>0) kind=2
    if(common==4) kind=3
    value=0
    do i=1,2
      do j=1,2
        value=value+triangle_pair(a(:,ix(:,i)),b(:,ix(:,j)),kind)
      end do
    end do
  end function
  complex(dp) function patch_phase(p,k,direction,order) result(value)
    real(dp), intent(in) :: p(3,4),k,direction(3)
    integer, intent(in) :: order
    real(dp) :: x(order),w(order),r(3),jac
    integer :: t,i,j
    integer, parameter :: ix(3,2)=reshape([1,2,3,1,3,4],[3,2])
    call gauss_rule(order,x,w)
    value=0
    do t=1,2
      do i=1,order
        do j=1,order
          call tri_map(p(:,ix(:,t)),(x(i)+1)/2,(x(j)+1)/2,r,jac)
          value=value+0.25_dp*w(i)*w(j)*jac*exp(iu*k*dot_product(direction,r))
        end do
      end do
    end do
  end function
end module
