module peec_physics
  use peec_base
  use peec_model
  use peec_config
  implicit none
  type slot_cell
    integer :: a=0,b=0
    real(dp) :: length=0,l=0,c=0
  end type
  type shunt_type
    integer :: a=0,b=0
    real(dp) :: r=0
    logical :: enabled=.false.
  end type
  type coupling_type
    integer, allocatable :: edge(:),node(:)
    real(dp), allocatable :: sign(:)
  end type
contains
  subroutine basis(theta,phi,polarization,khat,ehat)
    real(dp), intent(in) :: theta,phi
    character(*), intent(in) :: polarization
    real(dp), intent(out) :: khat(3),ehat(3)
    real(dp) :: t,p
    t=theta*pi/180; p=phi*pi/180
    khat=-[sin(t)*cos(p),sin(t)*sin(p),cos(t)]
    if(polarization=='horizontal') then
      ehat=[-sin(p),cos(p),0.0_dp]
    else
      ehat=[cos(t)*cos(p),cos(t)*sin(p),-sin(t)]
    end if
  end subroutine
  function harmonic_excitation(m,c,frequency) result(u)
    type(mesh_type), intent(in) :: m
    type(config_type), intent(in) :: c
    real(dp), intent(in) :: frequency
    complex(dp), allocatable :: u(:)
    real(dp) :: khat(3),ehat(3),k
    integer :: e
    call basis(c%theta_deg,c%phi_deg,c%polarization,khat,ehat)
    k=2*pi*frequency/c0
    allocate(u(size(m%edges,2)))
    do e=1,size(u)
      u(e)=c%field_amplitude*dot_product(ehat,m%vector(:,e))* &
        exp(iu*(c%phase_rad-k*dot_product(khat,m%center(:,e))))* &
        sinc(0.5_dp*k*dot_product(khat,m%vector(:,e)))
    end do
  end function
  subroutine time_excitation(m,c,t,strike,return_node,node_source,edge_voltage)
    type(mesh_type), intent(in) :: m
    type(config_type), intent(in) :: c
    real(dp), intent(in) :: t
    integer, intent(in) :: strike,return_node
    real(dp), intent(out) :: node_source(:),edge_voltage(:)
    real(dp) :: khat(3),ehat(3),tau,val
    integer :: e
    node_source=0; edge_voltage=0
    if(c%excitation=='lightning-current') then
      val=pulse(t-c%lightning_delay,c%lightning_K,c%lightning_alpha,c%lightning_beta)
      node_source(strike)=val; node_source(return_node)=-val
    else
      call basis(c%theta_deg,c%phi_deg,c%polarization,khat,ehat)
      do e=1,size(edge_voltage)
        tau=t-dot_product(khat,m%center(:,e))/c0
        edge_voltage(e)=dot_product(ehat,m%vector(:,e))*pulse(tau,c%pulse_K,c%pulse_alpha,c%pulse_beta)
      end do
    end if
  end subroutine
  real(dp) function rcs(m,current,k,amplitude,theta,phi,dual,order) result(sigma)
    type(mesh_type), intent(in) :: m
    complex(dp), intent(in) :: current(:)
    real(dp), intent(in) :: k,amplitude,theta,phi
    logical, intent(in) :: dual
    integer, intent(in) :: order
    real(dp) :: d(3),t,p
    complex(dp) :: f(3),transverse(3),g
    integer :: e,j
    call require(order_supported(order),'unsupported RCS order')
    call require(k>0.and.amplitude>0,'invalid RCS k/amplitude')
    t=theta*pi/180; p=phi*pi/180
    d=[sin(t)*cos(p),sin(t)*sin(p),cos(t)]
    f=0
    do e=1,size(current)
      if(dual) then
        g=0
        do j=1,size(m%branches(e)%p,3)
          g=g+patch_phase(m%branches(e)%p(:,:,j),k,d,order)
        end do
        g=g*current(e)*m%length(e)/m%branches(e)%area
      else
        g=current(e)*m%length(e)*exp(iu*k*dot_product(d,m%center(:,e)))
      end if
      f=f+g*m%direction(:,e)
    end do
    transverse=f-d*sum(d*f)
    sigma=k*k*(mu0/eps0)*sum(abs(transverse)**2)/(4*pi*amplitude**2)
    call require(ieee_is_finite(sigma).and.sigma>=0,'nonfinite RCS')
  end function
  subroutine configure_shunt(m,c,shunt)
    type(mesh_type), intent(in) :: m
    type(config_type), intent(in) :: c
    type(shunt_type), intent(out) :: shunt
    shunt%enabled=c%use_internal_shunt
    if(.not.shunt%enabled) return
    shunt%a=resolve_node(m,c%shunt_a_node,c%shunt_a_xyz,c%shunt_a_xyz_set)
    shunt%b=resolve_node(m,c%shunt_b_node,c%shunt_b_xyz,c%shunt_b_xyz_set)
    call require(shunt%a/=shunt%b,'shunt terminals coincide')
    shunt%r=c%shunt_resistance
  end subroutine
  subroutine slot_line(width,span,thickness,epsr,mur,lp,cp)
    real(dp), intent(in) :: width,span,thickness,epsr,mur
    real(dp), intent(out) :: lp,cp
    real(dp) :: we,ratio,q,den,z,vel,mu
    call require(width>0.and.span>0.and.thickness>=0.and.epsr>0.and.mur>0,'invalid slot parameters')
    we=width
    if(thickness>0) we=width-5*thickness/(4*pi)*(1+log(4*pi*width/thickness))
    call require(we>0.and.we<span/sqrt(2.0_dp),'slot outside narrow-slot model range')
    ratio=we/span; q=sqrt(sqrt(1-ratio**2))
    call require(q<1,'slot ratio underflow')
    den=log(2*(1+q)/(1-q))
    ! Match the original slot module's permeability constant exactly.
    mu=4e-7_dp*pi
    z=120*pi*pi/den*sqrt(mur/epsr)
    vel=1/sqrt(mu*mur*eps0*epsr)
    lp=z/vel; cp=1/(z*vel)
  end subroutine
  subroutine read_slots(c,nv,cells)
    type(config_type), intent(in) :: c
    integer, intent(in) :: nv
    type(slot_cell), allocatable, intent(out) :: cells(:)
    type(slot_cell) :: cell
    type(token), allocatable :: args(:)
    character(4096) :: line
    real(dp) :: lp,cp
    integer :: u,ios
    call slot_line(c%slot_width,c%slot_wall_span,c%slot_wall_thickness,c%slot_epsilon_r,c%slot_mu_r,lp,cp)
    allocate(cells(0))
    open(newunit=u,file=c%slot_cells_file,status='old',action='read',iostat=ios)
    call require(ios==0,'cannot read slot-cell file')
    do
      read(u,'(a)',iostat=ios) line
      if(ios<0) exit
      call require(ios==0,'slot-cell read failed')
      allocate(args(0)); call tokenize(trim(line),args)
      if(size(args)>0) then
        call require(size(args)==3,'slot-cell record requires node_a node_b support_length')
        cell%a=parse_int(args(1)%s)+1; cell%b=parse_int(args(2)%s)+1
        cell%length=parse_real(args(3)%s)
        call require(cell%a>=1.and.cell%a<=nv.and.cell%b>=1.and.cell%b<=nv.and. &
          cell%a/=cell%b.and.cell%length>0,'invalid slot cell')
        cell%l=lp*cell%length; cell%c=cp*cell%length
        cells=[cells,cell]
      end if
      deallocate(args)
    end do
    close(u)
    call require(size(cells)>0,'empty slot-cell file')
  end subroutine
  subroutine build_coupling(closed,opened,mapping,cover,coupling,info)
    type(mesh_type), intent(in) :: closed,opened
    integer, intent(in) :: mapping(:,:),cover(:)
    type(coupling_type), intent(out) :: coupling
    integer, intent(out) :: info
    type(coupling_type) :: work
    integer, allocatable :: map(:)
    real(dp) :: scale,tol
    integer :: k,a,b,e,j
    info=1
    if(size(mapping,1)/=2.or.size(mapping,2)==0.or.size(cover)==0) return
    allocate(map(size(closed%xyz,2))); map=0
    scale=max(maxval(maxval(closed%xyz,dim=2)-minval(closed%xyz,dim=2)), &
      maxval(maxval(opened%xyz,dim=2)-minval(opened%xyz,dim=2)))
    if(scale<=0) scale=1
    tol=max(1e-12_dp,1e-10_dp*scale)
    do k=1,size(mapping,2)
      a=mapping(1,k); b=mapping(2,k)
      if(a<1.or.a>size(map).or.b<1.or.b>size(opened%xyz,2)) return
      if(map(a)/=0.or.any(map==b)) return
      if(norm(closed%xyz(:,a)-opened%xyz(:,b))>tol) return
      map(a)=b
    end do
    allocate(work%edge(0),work%node(0),work%sign(0))
    do k=1,size(cover)
      e=cover(k)
      if(e<1.or.e>size(closed%edges,2)) return
      if(count(cover==e)/=1) return
      a=closed%edges(1,e); b=closed%edges(2,e)
      if(a<1.or.a>size(map).or.b<1.or.b>size(map)) return
      if((map(a)>0).eqv.(map(b)>0)) cycle
      work%edge=[work%edge,e]
      if(map(a)>0) then
        work%node=[work%node,map(a)]; work%sign=[work%sign,-1.0_dp]
      else
        work%node=[work%node,map(b)]; work%sign=[work%sign,1.0_dp]
      end if
    end do
    coupling=work; info=0
  end subroutine
  subroutine apply_coupling(coupling,current,source)
    type(coupling_type), intent(in) :: coupling
    real(dp), intent(in) :: current(:)
    real(dp), intent(out) :: source(:)
    integer :: k
    source=0
    do k=1,size(coupling%edge)
      source(coupling%node(k))=source(coupling%node(k))+coupling%sign(k)*current(coupling%edge(k))
    end do
  end subroutine
end module
