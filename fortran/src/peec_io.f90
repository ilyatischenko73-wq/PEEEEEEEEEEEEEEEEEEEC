module peec_io
  use peec_transient
  implicit none
contains
  function shell_quote(path) result(q)
    character(*), intent(in) :: path
    character(:), allocatable :: q
    integer :: i
    q=achar(39)
    do i=1,len(path)
      if(path(i:i)==achar(39)) then
        q=q//achar(39)//achar(34)//achar(39)//achar(34)//achar(39)
      else
        q=q//path(i:i)
      end if
    end do
    q=q//achar(39)
  end function
  subroutine mkdir(path)
    character(*), intent(in) :: path
    integer :: status,cmdstatus,i
    character(32) :: os
    character(:), allocatable :: native
    if(len_trim(path)==0.or.path=='.') return
    call get_environment_variable('OS',os,status=cmdstatus)
    if(cmdstatus==0.and.trim(os)=='Windows_NT') then
      native=path
      do i=1,len(native)
        if(native(i:i)=='/') native(i:i)=achar(92)
      end do
      ! Prevent cmd.exe expansion; spaces and ordinary quoted paths are accepted.
      call require(scan(native,'"%!&|<>^')==0,'unsupported character in Windows output path')
      call execute_command_line('if not exist "'//native//'\\." mkdir "'//native//'"', &
        exitstat=status,cmdstat=cmdstatus)
    else
      call execute_command_line('mkdir -p -- '//shell_quote(path),exitstat=status,cmdstat=cmdstatus)
    end if
    call require(cmdstatus==0,'cannot execute mkdir')
    call require(status==0,'cannot create output directory: '//path)
  end subroutine
  subroutine open_output(path,u)
    character(*), intent(in) :: path
    integer, intent(out) :: u
    integer :: k,ios
    k=max(index(path,'/',back=.true.),index(path,achar(92),back=.true.))
    if(k>1) call mkdir(path(:k-1))
    open(newunit=u,file=path,status='replace',action='write',iostat=ios)
    call require(ios==0,'cannot write '//path)
  end subroutine
  subroutine write_matrix(path,a)
    character(*), intent(in) :: path
    real(dp), intent(in) :: a(:,:)
    integer :: u,j
    call open_output(path,u)
    write(u,'(a)') '%%MatrixMarket matrix array real general'
    write(u,'(i0,1x,i0)') size(a,1),size(a,2)
    do j=1,size(a,2)
      write(u,'(es26.17e3)') a(:,j)
    end do
    close(u)
  end subroutine
  subroutine export_matrices(path,m)
    character(*), intent(in) :: path
    type(model_type), intent(in) :: m
    integer :: u,e,ne,nv
    ne=size(m%a,1); nv=size(m%a,2)
    call write_matrix(path//'/L.mtx',m%l)
    call write_matrix(path//'/P.mtx',m%p)
    call write_matrix(path//'/PAT.mtx',m%pat)
    call open_output(path//'/A.mtx',u)
    write(u,'(a)') '%%MatrixMarket matrix coordinate real general'
    write(u,'(3(i0,1x))') ne,nv,2*ne
    do e=1,ne
      write(u,'(2(i0,1x),a)') e,m%mesh%edges(1,e),'-1.0'
      write(u,'(2(i0,1x),a)') e,m%mesh%edges(2,e),'1.0'
    end do
    close(u)
    call open_output(path//'/R.mtx',u)
    write(u,'(a)') '%%MatrixMarket matrix coordinate real general'
    write(u,'(3(i0,1x))') ne,ne,ne
    do e=1,ne
      write(u,'(2(i0,1x),es26.17e3)') e,e,m%r(e)
    end do
    close(u)
  end subroutine
  subroutine scalar(u,name,value)
    integer, intent(in) :: u
    character(*), intent(in) :: name
    real(dp), intent(in) :: value(:)
    write(u,'(a)') 'SCALARS '//name//' double 1','LOOKUP_TABLE default'
    write(u,'(es26.17e3)') value
  end subroutine
  subroutine complex_scalar(u,stem,unit,value)
    integer, intent(in) :: u
    character(*), intent(in) :: stem,unit
    complex(dp), intent(in) :: value(:)
    call scalar(u,stem//'_re_'//unit,real(value,dp))
    call scalar(u,stem//'_im_'//unit,aimag(value))
    call scalar(u,stem//'_abs_'//unit,abs(value))
  end subroutine
  function frame_name(prefix,kind,frame) result(name)
    character(*), intent(in) :: prefix,kind
    integer, intent(in) :: frame
    character(:), allocatable :: name
    character(20) :: number
    if(frame<0) then
      name=prefix//'_'//kind
    else
      write(number,'(i6.6)') frame
      name=prefix//'_'//kind//'_'//trim(number)
    end if
  end function
  subroutine vtk_state(directory,prefix,frame,time,m,current,phi,charge,harmonic_state)
    character(*), intent(in) :: directory,prefix
    integer, intent(in) :: frame
    real(dp), intent(in) :: time
    type(mesh_type), intent(in) :: m
    complex(dp), intent(in) :: current(:),phi(:),charge(:)
    logical, intent(in) :: harmonic_state
    integer :: kind,u,j,part,nv,ne,nq
    character(:), allocatable :: label,name
    real(dp) :: area(size(phi)),jvec(3),v
    nv=size(phi); ne=size(current); nq=size(m%quads,2)
    area=m%nodes%area
    do kind=1,2
      if(kind==1) then
        label='surface'
      else
        label='edges'
      end if
      name=frame_name(prefix,label,frame)//'.vtk'
      call open_output(directory//'/'//name,u)
      write(u,'(a)') '# vtk DataFile Version 3.0','PEEC Fortran','ASCII','DATASET POLYDATA'
      write(u,'(a,i0,a)') 'POINTS ',nv,' double'
      do j=1,nv
        write(u,'(3(es26.17e3,1x))') m%xyz(:,j)
      end do
      if(kind==1) then
        write(u,'(a,i0,1x,i0)') 'POLYGONS ',nq,5*nq
        do j=1,nq
          write(u,'(5(i0,1x))') 4,m%quads(:,j)-1
        end do
        write(u,'(a,i0)') 'POINT_DATA ',nv
        if(harmonic_state) then
          call complex_scalar(u,'phi','V',phi)
          call complex_scalar(u,'charge','C',charge)
          call complex_scalar(u,'sigma','C_m2',charge/area)
        else
          call scalar(u,'time',[(time,j=1,nv)])
          call scalar(u,'phi_V',real(phi,dp))
          call scalar(u,'charge_C',real(charge,dp))
          call scalar(u,'sigma_C_m2',real(charge,dp)/area)
        end if
      else
        write(u,'(a,i0,1x,i0)') 'LINES ',ne,3*ne
        do j=1,ne
          write(u,'(3(i0,1x))') 2,m%edges(:,j)-1
        end do
        write(u,'(a,i0)') 'CELL_DATA ',ne
        if(harmonic_state) then
          call complex_scalar(u,'current','A',current)
        else
          call scalar(u,'time',[(time,j=1,ne)])
          call scalar(u,'current_A',real(current,dp))
        end if
        do part=1,merge(2,1,harmonic_state)
          label='J_A_m'
          if(harmonic_state) then
            if(part==1) label='J_re_A_m'
            if(part==2) label='J_im_A_m'
          end if
          write(u,'(a)') 'VECTORS '//label//' double'
          do j=1,ne
            v=real(current(j),dp)
            if(part==2) v=aimag(current(j))
            jvec=v*m%length(j)/m%branches(j)%area*m%direction(:,j)
            write(u,'(3(es26.17e3,1x))') jvec
          end do
        end do
      end if
      close(u)
    end do
  end subroutine
  subroutine write_frame(c,prefix,frame,time,m,t,cells,shunt)
    type(config_type), intent(in) :: c
    character(*), intent(in) :: prefix
    integer, intent(in) :: frame
    real(dp), intent(in) :: time
    type(model_type), intent(in) :: m
    type(transient_type), intent(in) :: t
    type(slot_cell), intent(in) :: cells(:)
    type(shunt_type), intent(in) :: shunt
    real(dp) :: phi(t%nv),charge(t%nv)
    integer :: u,j,s
    if(.not.c%write_vtk) return
    call recover(m,t,phi,charge)
    call vtk_state(c%vtk_directory,prefix,frame,time,m%mesh,cmplx(t%x(:t%ne),0.0_dp,dp), &
      cmplx(phi,0.0_dp,dp),cmplx(charge,0.0_dp,dp),.false.)
    if(shunt%enabled) then
      call open_output(c%vtk_directory//'/'//frame_name(prefix,'shunt',frame)//'.dat',u)
      write(u,'(a,es26.17e3)') 'time ',time
      write(u,'(a,i0)') 'node_a ',shunt%a-1
      write(u,'(a,3(es26.17e3,1x))') 'xyz_a ',m%mesh%xyz(:,shunt%a)
      write(u,'(a,i0)') 'node_b ',shunt%b-1
      write(u,'(a,3(es26.17e3,1x))') 'xyz_b ',m%mesh%xyz(:,shunt%b)
      write(u,'(a,es26.17e3)') 'shunt_current_A ',t%x(t%n)
      write(u,'(a,es26.17e3)') 'shunt_resistance_ohm ',shunt%r
      close(u)
    end if
    if(size(cells)>0) then
      call open_output(c%vtk_directory//'/'//frame_name(prefix,'slot',frame)//'.dat',u)
      write(u,'(a,es26.17e3)') 'time ',time
      write(u,'(a,i0)') 'n_cells ',size(cells)
      do j=1,size(cells)
        s=t%ne+t%nv+j
        write(u,'(3(a,i0),4(a,es26.17e3))') 'cell ',j-1,' node_a ',cells(j)%a-1, &
          ' node_b ',cells(j)%b-1,' I_L_A ',t%x(s),' L_H ',cells(j)%l,' C_F ',cells(j)%c,' dl_m ',cells(j)%length
      end do
      close(u)
    end if
  end subroutine
  subroutine write_rcs(c,m,current,frequency)
    type(config_type), intent(in) :: c
    type(mesh_type), intent(in) :: m
    complex(dp), intent(in) :: current(:)
    real(dp), intent(in) :: frequency
    real(dp) :: phi,sigma,db,k
    integer :: u,j,n
    k=2*pi*frequency/c0
    n=floor((c%rcs_phi_end-c%rcs_phi_start)/c%rcs_phi_step+0.5_dp)
    call open_output(c%output_file,u)
    write(u,'(a)') 'theta_deg,phi_deg,sigma_m2,sigma_dbsm,frequency_hz,k_rad_m,ka,a_m'
    do j=0,n
      phi=c%rcs_phi_start+j*c%rcs_phi_step
      sigma=rcs(m,current,k,c%field_amplitude,c%rcs_observation_theta,phi,c%rcs_use_dual,c%rcs_order)
      db=ieee_value(0.0_dp,ieee_negative_inf)
      if(sigma>tiny(1.0_dp)) db=10*log10(sigma)
      write(u,"(*(es26.17e3,:,','))") c%rcs_observation_theta,phi,sigma,db,frequency,k,c%ka,c%characteristic_length
    end do
    close(u)
    print '(a)', 'RCS saved: '//c%output_file
  end subroutine
end module
