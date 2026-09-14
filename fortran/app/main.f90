program peec_solver_fortran
  use peec_io
  use peec_aperture_json
  implicit none
  ! 0: menu; 1: mesh; 2: scattering; 3: RCS; 4: closed lightning;
  ! 5: two-stage lightning; 6: slot cells; 7: custom_config.
  ! Command-line arguments take precedence over selected_task.
  integer, parameter :: selected_task=0
  character(*), parameter :: custom_config='cases/run_scattering.cfg'
  type(config_type) :: cfg
  integer :: choice,ios
  character(4096) :: input,config_path
  if(command_argument_count()>0) then
    call read_config(cfg)
  else
    choice=selected_task
    config_path=custom_config
    if(choice==0) then
      do
        print '(a)', 'PEEC Fortran: select task'
        print '(a)', '1 Mesh information (plate)'
        print '(a)', '2 Harmonic scattering'
        print '(a)', '3 RCS'
        print '(a)', '4 Lightning: closed body'
        print '(a)', '5 Lightning: two stages'
        print '(a)', '6 Lightning: slot cells'
        print '(a)', '7 Custom configuration file'
        print '(a)', '0 Exit'
        read(*,'(a)',iostat=ios) input
        if(ios<0) stop
        if(ios/=0) error stop 'Cannot read menu input'
        input=trim(adjustl(input))
        if(len_trim(input)/=1) cycle
        if(verify(trim(input),'01234567')/=0) cycle
        read(input,*) choice
        exit
      end do
      if(choice==0) stop
      if(choice==7) then
        print '(a)', 'Configuration path (without surrounding quotes):'
        read(*,'(a)',iostat=ios) config_path
        if(ios/=0) error stop 'Cannot read configuration path'
      end if
    end if
    select case(choice)
    case(1)
      config_path='fortran/cases/mesh_info.cfg'
    case(2)
      config_path='cases/run_scattering.cfg'
    case(3)
      config_path='cases/run_rcs.cfg'
    case(4)
      config_path='fortran/cases/lightning_closed.cfg'
    case(5)
      config_path='cases/run_lightning.cfg'
    case(6)
      config_path='cases/run_lightning_slot_cells.cfg'
    case(7)
      ! Use custom_config or the path entered in the menu.
    case default
      error stop 'selected_task must be between 0 and 7'
    end select
    print '(a)', 'Configuration: '//trim(config_path)
    call read_config(cfg,trim(config_path))
  end if
  select case(cfg%task)
  case('mesh-info')
    block
      type(mesh_type) :: mesh
      call read_mesh(cfg%mesh_file,mesh)
    end block
  case('scattering','rcs')
    call run_harmonic(cfg)
  case('lightning')
    call run_lightning(cfg)
  end select
contains
  subroutine run_harmonic(c)
    type(config_type), intent(in) :: c
    type(model_type) :: model
    complex(dp), allocatable :: u(:),current(:),phi(:),charge(:)
    real(dp) :: f,err
    character(:), allocatable :: prefix
    call assemble(c%mesh_file,model,c%parallel_threads)
    if(c%save_matrices) call export_matrices(c%matrix_directory,model)
    f=c%frequency
    if(c%task=='rcs') f=c0*c%ka/(2*pi*c%characteristic_length)
    u=harmonic_excitation(model%mesh,c,f)
    call harmonic(model,f,u,current,phi,charge,err)
    write(*,'(a,es18.9e3)') 'Frequency          : ',f
    write(*,'(a,es18.9e3)') 'Maximum |I_e|      : ',maxval(abs(current))
    write(*,'(a,es18.9e3)') 'Maximum |V_j|      : ',maxval(abs(phi))
    write(*,'(a,es18.9e3)') 'Backward error     : ',err
    prefix='scattering'
    if(c%task=='rcs') prefix='rcs_solution'
    if(c%write_vtk) call vtk_state(c%vtk_directory,prefix,-1,0.0_dp,model%mesh,current,phi,charge,.true.)
    if(c%task=='rcs') call write_rcs(c,model%mesh,current,f)
  end subroutine
  subroutine summary_row(u,time,source,m,t,shunt)
    integer, intent(in) :: u
    real(dp), intent(in) :: time,source
    type(model_type), intent(in) :: m
    type(transient_type), intent(in) :: t
    type(shunt_type), intent(in) :: shunt
    real(dp) :: shunt_current,slot_current,total_charge
    shunt_current=0; slot_current=0
    if(shunt%enabled) shunt_current=t%x(t%n)
    if(t%ns>0) slot_current=maxval(abs(t%x(t%ne+t%nv+1:t%ne+t%nv+t%ns)))
    total_charge=sum(m%f*t%x(t%ne+1:t%ne+t%nv))
    write(u,"(*(es26.17e3,:,','))") time,source,maxval(abs(t%x(:t%ne))),total_charge,shunt_current,slot_current
  end subroutine
  subroutine run_lightning(c)
    type(config_type), intent(in) :: c
    type(model_type) :: m1,m2
    type(transient_type) :: t1,t2
    type(shunt_type) :: shunt1,shunt2
    type(slot_cell), allocatable :: cells(:),empty(:)
    type(coupling_type) :: coupling
    real(dp), allocatable :: jn(:),jp(:),un(:),up(:),gn(:),gp(:),zero(:)
    character(:), allocatable :: prefix,path
    logical :: two_stage
    integer :: strike,ret,step,steps,frame,ne,nv,out1,out2
    real(dp) :: time,source
    allocate(empty(0),cells(0))
    two_stage=c%lightning_case=='two-stage'
    if(c%lightning_case=='slot-cells') then
      path=c%open_mesh_file; prefix='lightning_slot_cells'
    else
      path=c%closed_mesh_file; prefix='lightning_closed'
      if(two_stage) prefix='lightning_stage1'
    end if
    call assemble(path,m1,c%parallel_threads)
    ne=size(m1%l,1); nv=size(m1%p,1)
    if(c%lightning_case=='slot-cells') then
      call read_slots(c,nv,cells)
      call configure_shunt(m1%mesh,c,shunt1)
    end if
    call build_transient(m1,cells,shunt1,c%dt,c%time_order,t1)
    if(c%save_matrices) then
      path=c%matrix_directory
      if(two_stage) path=path//'/stage1'
      call export_matrices(path,m1)
    end if
    strike=0; ret=0
    if(c%excitation=='lightning-current') then
      strike=resolve_node(m1%mesh,c%strike_node,c%strike_xyz,c%strike_xyz_set)
      ret=resolve_node(m1%mesh,c%return_node,c%return_xyz,c%return_xyz_set)
      call require(strike/=ret,'strike and return nodes coincide')
    end if
    allocate(jn(nv),jp(nv),un(ne),up(ne))
    call time_excitation(m1%mesh,c,0.0_dp,strike,ret,jn,un)
    if(two_stage) then
      call assemble(c%open_mesh_file,m2,c%parallel_threads)
      if(c%save_matrices) call export_matrices(c%matrix_directory//'/stage2',m2)
      call read_coupling(c%aperture_map_file,m1%mesh,m2%mesh,coupling)
      call configure_shunt(m2%mesh,c,shunt2)
      call build_transient(m2,empty,shunt2,c%dt,c%time_order,t2)
      allocate(gn(t2%nv),gp(t2%nv),zero(t2%ne)); zero=0
      call apply_coupling(coupling,t1%x(:ne),gn)
    end if
    call open_output(c%vtk_directory//'/'//prefix//'.csv',out1)
    write(out1,'(a)') 'time_s,source_A,max_edge_current_A,total_charge_C,shunt_current_A,max_slot_current_A'
    call summary_row(out1,0.0_dp,0.0_dp,m1,t1,shunt1)
    call write_frame(c,prefix,0,0.0_dp,m1,t1,cells,shunt1)
    if(two_stage) then
      call open_output(c%vtk_directory//'/lightning_stage2.csv',out2)
      write(out2,'(a)') 'time_s,source_A,max_edge_current_A,total_charge_C,shunt_current_A,max_slot_current_A'
      call summary_row(out2,0.0_dp,0.0_dp,m2,t2,shunt2)
      call write_frame(c,'lightning_stage2',0,0.0_dp,m2,t2,empty,shunt2)
    end if
    steps=nint(c%t_end/c%dt); frame=0
    do step=1,steps
      time=step*c%dt
      call time_excitation(m1%mesh,c,time,strike,ret,jp,up)
      call step_transient(t1,jn,jp,un,up)
      source=0
      if(strike>0) source=jp(strike)
      call summary_row(out1,time,source,m1,t1,shunt1)
      if(two_stage) then
        call apply_coupling(coupling,t1%x(:ne),gp)
        call step_transient(t2,gn,gp,zero,zero)
        call summary_row(out2,time,0.5_dp*sum(abs(gp)),m2,t2,shunt2)
        gn=gp
      end if
      if(mod(step,c%vtk_every)==0.or.step==steps) then
        frame=frame+1
        call write_frame(c,prefix,frame,time,m1,t1,cells,shunt1)
        if(two_stage) call write_frame(c,'lightning_stage2',frame,time,m2,t2,empty,shunt2)
      end if
      if(mod(step,max(1,steps/10))==0.or.step==steps) then
        write(*,'(a,i0,a,i0,a,es14.6e3)') 'Step ',step,' / ',steps,'; t=',time
      end if
      jn=jp; un=up
    end do
    close(out1)
    if(two_stage) close(out2)
    print '(a)', 'Transient results: '//c%vtk_directory
  end subroutine
end program
