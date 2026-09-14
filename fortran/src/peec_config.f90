module peec_config
  use peec_base
  use peec_quadrature, only: order_supported
  implicit none
  type token
    character(:), allocatable :: s
  end type
  type config_type
    real(dp) :: frequency=300e6_dp
    real(dp) :: field_amplitude=1.0_dp
    real(dp) :: theta_deg=90.0_dp
    real(dp) :: phi_deg=45.0_dp
    real(dp) :: phase_rad=0.0_dp
    real(dp) :: ka=1.0_dp
    real(dp) :: characteristic_length=1.0_dp
    real(dp) :: rcs_observation_theta=90.0_dp
    real(dp) :: rcs_phi_start=0.0_dp
    real(dp) :: rcs_phi_end=360.0_dp
    real(dp) :: rcs_phi_step=1.0_dp
    real(dp) :: slot_width=0.0_dp
    real(dp) :: slot_wall_span=0.0_dp
    real(dp) :: slot_wall_thickness=0.0_dp
    real(dp) :: slot_epsilon_r=1.0_dp
    real(dp) :: slot_mu_r=1.0_dp
    real(dp) :: dt=0.0_dp
    real(dp) :: t_end=0.0_dp
    real(dp) :: lightning_K=102900.0_dp
    real(dp) :: lightning_alpha=1500.0_dp
    real(dp) :: lightning_beta=1e6_dp
    real(dp) :: lightning_delay=0.0_dp
    real(dp) :: pulse_K=1.0_dp
    real(dp) :: pulse_alpha=1500.0_dp
    real(dp) :: pulse_beta=1e6_dp
    real(dp) :: shunt_resistance=1000.0_dp
    integer :: parallel_threads=1
    integer :: rcs_order=6
    integer :: vtk_every=1
    integer :: time_order=2
    integer :: strike_node=-1
    integer :: return_node=-1
    integer :: shunt_a_node=-1
    integer :: shunt_b_node=-1
    character(:), allocatable :: task
    character(:), allocatable :: physics
    character(:), allocatable :: mesh_file
    character(:), allocatable :: polarization
    character(:), allocatable :: output_file
    character(:), allocatable :: matrix_directory
    character(:), allocatable :: vtk_directory
    character(:), allocatable :: lightning_case
    character(:), allocatable :: excitation
    character(:), allocatable :: closed_mesh_file
    character(:), allocatable :: open_mesh_file
    character(:), allocatable :: aperture_map_file
    character(:), allocatable :: slot_cells_file
    logical :: rcs_use_dual=.true.,save_matrices=.false.,write_vtk=.false.,use_internal_shunt=.false.
    real(dp) :: strike_xyz(3)=0,return_xyz(3)=0,shunt_a_xyz(3)=0,shunt_b_xyz(3)=0
    logical :: strike_xyz_set=.false.,return_xyz_set=.false.,shunt_a_xyz_set=.false.,shunt_b_xyz_set=.false.
  end type
contains
  subroutine append_token(args,word)
    type(token), allocatable, intent(inout) :: args(:)
    character(*), intent(in) :: word
    type(token), allocatable :: work(:)
    integer :: n
    n=size(args)
    allocate(work(n+1))
    work(:n)=args
    work(n+1)%s=word
    call move_alloc(work,args)
  end subroutine
  subroutine tokenize(line,args)
    character(*), intent(in) :: line
    type(token), allocatable, intent(inout) :: args(:)
    integer :: i,n
    character :: quote,ch
    character(:), allocatable :: word
    n=len_trim(line); i=1
    do while(i<=n)
      if(line(i:i)=='#') exit
      if(iachar(line(i:i))<=32) then
        i=i+1
        cycle
      end if
      word=''; quote=achar(0)
      do while(i<=n)
        ch=line(i:i)
        if(quote==achar(0)) then
          if(ch=='#'.or.iachar(ch)<=32) exit
          if(ch=='"'.or.ch=="'") then
            quote=ch
          else
            word=word//ch
          end if
        else if(ch==quote) then
          quote=achar(0)
        else
          word=word//ch
        end if
        i=i+1
      end do
      call require(quote==achar(0),'unclosed quote in configuration')
      call append_token(args,word)
      if(i<=n) then
        if(line(i:i)=='#') exit
      end if
    end do
  end subroutine
  function next_arg(args,i) result(s)
    type(token), intent(in) :: args(:)
    integer, intent(inout) :: i
    character(:), allocatable :: s
    call require(i<size(args),'missing value for '//args(i)%s)
    i=i+1; s=args(i)%s
  end function
  real(dp) function parse_real(s) result(v)
    character(*), intent(in) :: s
    integer :: ios
    call require(len(s)>0,'empty numeric value')
    call require(verify(s,'0123456789+-.eEdD')==0,'invalid number: '//s)
    read(s,*,iostat=ios) v
    call require(ios==0,'invalid number: '//s)
    call require(ieee_is_finite(v),'nonfinite number: '//s)
  end function
  integer function parse_int(s) result(v)
    character(*), intent(in) :: s
    integer :: ios
    call require(len(s)>0,'empty integer')
    call require(verify(s,'0123456789+-')==0,'invalid integer: '//s)
    read(s,*,iostat=ios) v
    call require(ios==0,'invalid integer: '//s)
  end function
  subroutine help()
    print '(a)', 'PeecSolverFortran --config cases/run_scattering.cfg'
    print '(a)', 'or: --task mesh-info|scattering|rcs|lightning [original PEEC options]'
    print '(a)', 'Fortran solver; ASCII Gmsh 4.1; all external node/edge indices are zero-based.'
    print '(a)', 'Run from the repository root. See fortran/README_RU.md.'
  end subroutine
  subroutine read_config(cfg)
    type(config_type), intent(out) :: cfg
    type(token), allocatable :: args(:)
    character(:), allocatable :: arg,path
    character(16384) :: line
    integer :: n,i,j,u,ios,length
    cfg%task='mesh-info'
    cfg%physics='quasistatic'
    cfg%mesh_file=''
    cfg%polarization='horizontal'
    cfg%output_file='results/rcs.csv'
    cfg%matrix_directory='cache'
    cfg%vtk_directory='results/vtk'
    cfg%lightning_case='closed'
    cfg%excitation='lightning-current'
    cfg%closed_mesh_file=''
    cfg%open_mesh_file=''
    cfg%aperture_map_file=''
    cfg%slot_cells_file=''
    allocate(args(0))
    n=command_argument_count()
    if(n==0) then
      call help()
      stop
    end if
    do i=1,n
      call get_command_argument(i,length=length)
      allocate(character(length) :: arg)
      call get_command_argument(i,arg)
      call append_token(args,arg)
      deallocate(arg)
    end do
    if(args(1)%s=='--config') then
      call require(n==2,'--config FILE must be used without CLI overrides')
      path=args(2)%s
      open(newunit=u,file=path,status='old',action='read',iostat=ios)
      call require(ios==0,'cannot open config: '//path)
      deallocate(args); allocate(args(0))
      do
        read(u,'(a)',iostat=ios) line
        if(ios<0) exit
        call require(ios==0,'config read failed')
        call tokenize(trim(line),args)
      end do
      close(u)
    end if
    i=1
    do while(i<=size(args))
      arg=args(i)%s
      select case(arg)
      case('--task')
        cfg%task=next_arg(args,i)
      case('--physics')
        cfg%physics=next_arg(args,i)
      case('--mesh')
        cfg%mesh_file=next_arg(args,i)
      case('--parallel')
        cfg%parallel_threads=parse_int(next_arg(args,i))
      case('--frequency')
        cfg%frequency=parse_real(next_arg(args,i))
      case('--E','--field-amplitude')
        cfg%field_amplitude=parse_real(next_arg(args,i))
      case('--theta')
        cfg%theta_deg=parse_real(next_arg(args,i))
      case('--phi')
        cfg%phi_deg=parse_real(next_arg(args,i))
      case('--phase')
        cfg%phase_rad=parse_real(next_arg(args,i))
      case('--polarization')
        cfg%polarization=next_arg(args,i)
      case('--ka')
        cfg%ka=parse_real(next_arg(args,i))
      case('--a','--characteristic-length')
        cfg%characteristic_length=parse_real(next_arg(args,i))
      case('--rcs-order')
        cfg%rcs_order=parse_int(next_arg(args,i))
      case('--output')
        cfg%output_file=next_arg(args,i)
      case('--obs-theta','--rcs-observation-theta')
        cfg%rcs_observation_theta=parse_real(next_arg(args,i))
      case('--phi-start','--rcs-phi-start')
        cfg%rcs_phi_start=parse_real(next_arg(args,i))
      case('--phi-end','--rcs-phi-end')
        cfg%rcs_phi_end=parse_real(next_arg(args,i))
      case('--phi-step','--rcs-phi-step')
        cfg%rcs_phi_step=parse_real(next_arg(args,i))
      case('--matrix-dir')
        cfg%matrix_directory=next_arg(args,i)
      case('--vtk-dir')
        cfg%vtk_directory=next_arg(args,i)
      case('--vtk-every')
        cfg%vtk_every=parse_int(next_arg(args,i))
      case('--excitation')
        cfg%excitation=next_arg(args,i)
      case('--lightning-case')
        cfg%lightning_case=next_arg(args,i)
      case('--closed-mesh')
        cfg%closed_mesh_file=next_arg(args,i)
      case('--open-mesh')
        cfg%open_mesh_file=next_arg(args,i)
      case('--aperture-map')
        cfg%aperture_map_file=next_arg(args,i)
      case('--slot-cells')
        cfg%slot_cells_file=next_arg(args,i)
      case('--slot-width')
        cfg%slot_width=parse_real(next_arg(args,i))
      case('--slot-wall-span')
        cfg%slot_wall_span=parse_real(next_arg(args,i))
      case('--slot-wall-thickness')
        cfg%slot_wall_thickness=parse_real(next_arg(args,i))
      case('--slot-epsilon-r')
        cfg%slot_epsilon_r=parse_real(next_arg(args,i))
      case('--slot-mu-r')
        cfg%slot_mu_r=parse_real(next_arg(args,i))
      case('--dt')
        cfg%dt=parse_real(next_arg(args,i))
      case('--t-end')
        cfg%t_end=parse_real(next_arg(args,i))
      case('--time-order')
        cfg%time_order=parse_int(next_arg(args,i))
      case('--lightning-K')
        cfg%lightning_K=parse_real(next_arg(args,i))
      case('--lightning-alpha')
        cfg%lightning_alpha=parse_real(next_arg(args,i))
      case('--lightning-beta')
        cfg%lightning_beta=parse_real(next_arg(args,i))
      case('--lightning-delay')
        cfg%lightning_delay=parse_real(next_arg(args,i))
      case('--pulse-K')
        cfg%pulse_K=parse_real(next_arg(args,i))
      case('--pulse-alpha')
        cfg%pulse_alpha=parse_real(next_arg(args,i))
      case('--pulse-beta')
        cfg%pulse_beta=parse_real(next_arg(args,i))
      case('--strike-node')
        cfg%strike_node=parse_int(next_arg(args,i))
      case('--return-node')
        cfg%return_node=parse_int(next_arg(args,i))
      case('--shunt-R')
        cfg%shunt_resistance=parse_real(next_arg(args,i))
        cfg%use_internal_shunt=.true.
      case('--shunt-a-node')
        cfg%shunt_a_node=parse_int(next_arg(args,i))
      case('--shunt-b-node')
        cfg%shunt_b_node=parse_int(next_arg(args,i))
      case('--strike')
        do j=1,3
          cfg%strike_xyz(j)=parse_real(next_arg(args,i))
        end do
        cfg%strike_xyz_set=.true.
      case('--return')
        do j=1,3
          cfg%return_xyz(j)=parse_real(next_arg(args,i))
        end do
        cfg%return_xyz_set=.true.
      case('--shunt-a')
        do j=1,3
          cfg%shunt_a_xyz(j)=parse_real(next_arg(args,i))
        end do
        cfg%shunt_a_xyz_set=.true.
      case('--shunt-b')
        do j=1,3
          cfg%shunt_b_xyz(j)=parse_real(next_arg(args,i))
        end do
        cfg%shunt_b_xyz_set=.true.
      case('--rcs-dual','--rcs-use-dual')
        cfg%rcs_use_dual=.true.
      case('--rcs-no-dual')
        cfg%rcs_use_dual=.false.
      case('--save-matrices')
        cfg%save_matrices=.true.
      case('--vtk')
        cfg%write_vtk=.true.
      case('-h','--help')
        call help()
        stop
      case default
        call require(.false.,'unknown option: '//arg)
      end select
      i=i+1
    end do
    call validate_config(cfg)
  end subroutine
  subroutine validate_pair(a,b,xa,xb)
    integer, intent(in) :: a,b
    logical, intent(in) :: xa,xb
    call require(((a>=0.and.b>=0).and..not.(xa.or.xb)).or. &
      ((a<0.and.b<0).and.xa.and.xb),'provide either two node indices or two coordinate triples')
  end subroutine
  subroutine validate_config(c)
    type(config_type), intent(inout) :: c
    call require(c%physics=='quasistatic','retarded mode is not implemented in the original solver or this port')
    call require(c%parallel_threads>=1.and.c%vtk_every>=1,'parallel/vtk-every must be >= 1')
    call require(c%polarization=='horizontal'.or.c%polarization=='vertical','invalid polarization')
    call require(c%task=='mesh-info'.or.c%task=='scattering'.or.c%task=='rcs'.or.c%task=='lightning','invalid task')
    if(c%task/='lightning') call require(len(c%mesh_file)>0,'--mesh is required')
    if(c%task=='scattering'.or.c%task=='rcs') then
      call require(c%frequency>0.and.c%field_amplitude>=0,'invalid frequency or field amplitude')
    end if
    if(c%task=='rcs') then
      call require(c%ka>0.and.c%characteristic_length>0.and.c%field_amplitude>0,'invalid RCS parameters')
      call require(order_supported(c%rcs_order),'unsupported RCS quadrature order')
      call require(c%rcs_phi_step>0.and.c%rcs_phi_end>=c%rcs_phi_start,'invalid RCS sweep')
      call require(c%rcs_phi_end-c%rcs_phi_start<=1e6_dp*c%rcs_phi_step,'RCS sweep is too large')
      call require(c%rcs_phi_start+c%rcs_phi_step>c%rcs_phi_start,'RCS step is too small')
    end if
    if(c%task/='lightning') return
    call require(c%dt>0.and.c%t_end>=c%dt,'invalid dt or t-end')
    call require(c%time_order==1.or.c%time_order==2,'time-order must be 1 or 2')
    call require(c%t_end/c%dt<real(huge(1)-1,dp),'too many time steps')
    call require(abs(c%t_end/c%dt-real(nint(c%t_end/c%dt),dp))<= &
      1e-10_dp*max(1.0_dp,c%t_end/c%dt),'t-end must be an integer multiple of dt')
    call require(c%lightning_case=='closed'.or.c%lightning_case=='two-stage'.or. &
      c%lightning_case=='slot-cells','invalid lightning-case')
    if(c%lightning_case=='slot-cells') then
      if(len(c%open_mesh_file)==0) c%open_mesh_file=c%mesh_file
      call require(len(c%open_mesh_file)>0.and.len(c%slot_cells_file)>0,'open-mesh and slot-cells are required')
      call require(c%slot_width>0.and.c%slot_wall_span>0.and.c%slot_wall_thickness>=0.and. &
        c%slot_epsilon_r>0.and.c%slot_mu_r>0,'invalid slot geometry/material')
    else
      if(len(c%closed_mesh_file)==0) c%closed_mesh_file=c%mesh_file
      call require(len(c%closed_mesh_file)>0,'closed-mesh is required')
      if(c%lightning_case=='two-stage') then
        call require(len(c%open_mesh_file)>0.and.len(c%aperture_map_file)>0,'open-mesh and aperture-map are required')
      end if
    end if
    if(c%excitation=='lightning-current') then
      call require(c%lightning_alpha>0.and.c%lightning_beta>0.and.c%lightning_delay>=0,'invalid lightning pulse')
      call validate_pair(c%strike_node,c%return_node,c%strike_xyz_set,c%return_xyz_set)
    else
      call require(c%excitation=='incident-pulse','invalid excitation')
      call require(c%pulse_alpha>0.and.c%pulse_beta>0,'invalid incident pulse')
    end if
    if(c%lightning_case=='closed') c%use_internal_shunt=.false.
    if(c%use_internal_shunt) then
      call require(c%shunt_resistance>=0,'shunt resistance must be >= 0')
      call validate_pair(c%shunt_a_node,c%shunt_b_node,c%shunt_a_xyz_set,c%shunt_b_xyz_set)
    end if
  end subroutine
end module
