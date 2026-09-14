module peec_aperture_json
  use peec_physics
  implicit none
contains
  subroutine lex_json(text,args)
    character(*), intent(in) :: text
    type(token), allocatable, intent(out) :: args(:)
    integer :: i,j,n
    character :: c
    allocate(args(0)); i=1; n=len(text)
    do while(i<=n)
      c=text(i:i)
      if(iachar(c)<=32) then
        i=i+1
      else if(index('{}[],:',c)>0) then
        call append_token(args,c); i=i+1
      else if(c=='"') then
        j=i+1
        do while(j<=n)
          if(text(j:j)=='"') exit
          call require(text(j:j)/=achar(92),'escaped JSON keys are unsupported')
          j=j+1
        end do
        call require(j<=n,'unterminated JSON string')
        call append_token(args,text(i:j)); i=j+1
      else
        j=i
        do while(j<=n)
          if(index('0123456789-',text(j:j))==0) exit
          j=j+1
        end do
        call require(j>i,'invalid JSON token')
        call append_token(args,text(i:j-1)); i=j
      end if
    end do
  end subroutine
  subroutine expect(args,i,s)
    type(token), intent(in) :: args(:)
    integer, intent(inout) :: i
    character(*), intent(in) :: s
    call require(i<=size(args),'unexpected end of aperture JSON')
    call require(args(i)%s==s,'expected '//s//' in aperture JSON')
    i=i+1
  end subroutine
  logical function at(args,i,s)
    type(token), intent(in) :: args(:)
    integer, intent(in) :: i
    character(*), intent(in) :: s
    at=.false.
    if(i<=size(args)) at=args(i)%s==s
  end function
  integer function json_index(args,i) result(n)
    type(token), intent(in) :: args(:)
    integer, intent(inout) :: i
    call require(i<=size(args),'missing JSON index')
    n=parse_int(args(i)%s)
    call require(n>=0.and.n<huge(1),'JSON indices must be nonnegative')
    n=n+1; i=i+1
  end function
  subroutine read_coupling(path,closed,opened,coupling)
    character(*), intent(in) :: path
    type(mesh_type), intent(in) :: closed,opened
    type(coupling_type), intent(out) :: coupling
    type(token), allocatable :: args(:)
    character(:), allocatable :: text,key,field
    character(16384) :: line
    integer, allocatable :: mapping(:,:),cover(:)
    integer :: u,ios,i,a,b,info,value
    logical :: got_map,got_cover
    open(newunit=u,file=path,status='old',action='read',iostat=ios)
    call require(ios==0,'cannot read aperture map: '//path)
    text=''
    do
      read(u,'(a)',iostat=ios) line
      if(ios<0) exit
      call require(ios==0,'aperture JSON read failed')
      text=text//trim(line)//' '
    end do
    close(u)
    call lex_json(text,args); i=1
    allocate(mapping(2,0),cover(0)); got_map=.false.; got_cover=.false.
    call expect(args,i,'{')
    do
      call require(i<=size(args),'incomplete aperture JSON')
      key=args(i)%s; i=i+1
      call expect(args,i,':')
      select case(key)
      case('"node_mapping"')
        call require(.not.got_map,'duplicate node_mapping')
        got_map=.true.; call expect(args,i,'[')
        if(.not.at(args,i,']')) then
          do
            call expect(args,i,'{'); a=0; b=0
            do
              call require(i<=size(args),'incomplete node mapping')
              field=args(i)%s; i=i+1
              call expect(args,i,':'); value=json_index(args,i)
              select case(field)
              case('"closed_node"')
                call require(a==0,'duplicate closed_node'); a=value
              case('"open_node"')
                call require(b==0,'duplicate open_node'); b=value
              case default
                call require(.false.,'unknown node_mapping field')
              end select
              if(at(args,i,'}')) exit
              call expect(args,i,',')
            end do
            call expect(args,i,'}')
            call require(a>0.and.b>0,'node mapping needs both indices')
            mapping=reshape([reshape(mapping,[size(mapping)]),a,b],[2,size(mapping,2)+1])
            if(at(args,i,']')) exit
            call expect(args,i,',')
          end do
        end if
        call expect(args,i,']')
      case('"cover_edges"')
        call require(.not.got_cover,'duplicate cover_edges')
        got_cover=.true.; call expect(args,i,'[')
        if(.not.at(args,i,']')) then
          do
            value=json_index(args,i); cover=[cover,value]
            if(at(args,i,']')) exit
            call expect(args,i,',')
          end do
        end if
        call expect(args,i,']')
      case default
        call require(.false.,'unknown aperture JSON key: '//key)
      end select
      if(at(args,i,'}')) exit
      call expect(args,i,',')
    end do
    call expect(args,i,'}')
    call require(i>size(args).and.got_map.and.got_cover,'invalid aperture JSON structure')
    call build_coupling(closed,opened,mapping,cover,coupling,info)
    call require(info==0,'invalid aperture mapping geometry or indices')
  end subroutine
end module
