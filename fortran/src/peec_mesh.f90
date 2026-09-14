module peec_mesh
  use peec_base
  implicit none
  type patch_list
    real(dp), allocatable :: p(:,:,:)
    real(dp) :: area=0
  end type
  type mesh_type
    real(dp), allocatable :: xyz(:,:), vector(:,:), center(:,:), length(:), direction(:,:)
    integer, allocatable :: quads(:,:), edges(:,:), cell_edges(:,:), cell_sign(:,:)
    type(patch_list), allocatable :: nodes(:), branches(:)
  end type
contains
  subroutine read_mesh(path,m)
    character(*), intent(in) :: path
    type(mesh_type), intent(out) :: m
    character(4096) :: line
    integer :: u,ios,blocks,total,lo,hi,block,dim,entity,kind,n,j,k,cursor,id,e,idx
    integer, allocatable :: tags(:),alltags(:),q(:,:)
    real(dp) :: version
    logical :: have_format,have_nodes,have_elements
    have_format=.false.; have_nodes=.false.; have_elements=.false.
    open(newunit=u,file=path,status='old',action='read',iostat=ios)
    call require(ios==0,'cannot open mesh: '//path)
    do
      read(u,'(a)',iostat=ios) line
      if(ios<0) exit
      call require(ios==0,'mesh read failed')
      select case(trim(line))
      case('$MeshFormat')
        read(u,*,iostat=ios) version,kind,n
        call require(ios==0,'invalid Gmsh header')
        call require(version>=4.1_dp.and.version<5.and.kind==0,'expected ASCII Gmsh 4.1')
        have_format=.true.
      case('$Nodes')
        call require(have_format.and..not.have_nodes,'invalid or duplicate Nodes section')
        read(u,*,iostat=ios) blocks,total,lo,hi
        call require(ios==0.and.total>0.and.blocks>0,'invalid Nodes header')
        allocate(m%xyz(3,total),alltags(total)); cursor=0
        do block=1,blocks
          read(u,*,iostat=ios) dim,entity,kind,n
          call require(ios==0.and.kind==0.and.n>=0,'invalid or parametric node block')
          call require(n<=total-cursor,'node count exceeds header')
          allocate(tags(n))
          if(n>0) then
            read(u,*,iostat=ios) tags
            call require(ios==0.and.all(tags>0),'invalid node tags')
          end if
          do j=1,n
            read(u,*,iostat=ios) m%xyz(:,cursor+j)
            call require(ios==0,'invalid node coordinates')
            call require(all(ieee_is_finite(m%xyz(:,cursor+j))),'nonfinite coordinates')
            alltags(cursor+j)=tags(j)
          end do
          cursor=cursor+n; deallocate(tags)
        end do
        call require(cursor==total,'node count mismatch')
        do j=1,total
          call require(count(alltags==alltags(j))==1,'duplicate node tag')
        end do
        have_nodes=.true.
      case('$Elements')
        call require(have_nodes.and..not.have_elements,'invalid Elements section')
        read(u,*,iostat=ios) blocks,total,lo,hi
        call require(ios==0.and.total>0.and.blocks>0,'invalid Elements header')
        allocate(q(4,total),tags(4)); cursor=0; e=0
        do block=1,blocks
          read(u,*,iostat=ios) dim,entity,kind,n
          call require(ios==0.and.n>=0.and.n<=total-e,'invalid element block')
          do j=1,n
            read(u,'(a)',iostat=ios) line
            call require(ios==0,'missing element record')
            if(kind/=3) cycle
            read(line,*,iostat=ios) id,tags
            call require(ios==0,'invalid quadrangle')
            cursor=cursor+1
            do k=1,4
              idx=findloc(alltags,tags(k),dim=1)
              call require(idx>0,'unknown node tag in quadrangle')
              q(k,cursor)=idx
            end do
            do k=1,4
              call require(count(q(:,cursor)==q(k,cursor))==1,'repeated quadrangle vertex')
            end do
          end do
          e=e+n
        end do
        call require(e==total.and.cursor>0,'element count mismatch or no quadrangles')
        m%quads=q(:,:cursor)
        deallocate(q,tags); have_elements=.true.
      end select
    end do
    close(u)
    call require(have_nodes.and.have_elements,'incomplete mesh')
    call build_geometry(m)
  end subroutine
  subroutine build_geometry(m)
    type(mesh_type), intent(inout) :: m
    integer :: nq,nv,ne,c,k,j,a,b,e,nxt,prv,opp
    integer, allocatable :: edges(:,:),nn(:),nb(:)
    real(dp) :: p(3,4),mid(3,4),center(3),area
    nq=size(m%quads,2); nv=size(m%xyz,2); ne=0
    allocate(edges(2,4*nq),m%cell_edges(4,nq),m%cell_sign(4,nq))
    do c=1,nq
      do k=1,4
        a=m%quads(k,c); b=m%quads(mod(k,4)+1,c); e=0
        do j=1,ne
          if((edges(1,j)==a.and.edges(2,j)==b).or.(edges(1,j)==b.and.edges(2,j)==a)) then
            e=j
            exit
          end if
        end do
        if(e==0) then
          ne=ne+1; e=ne; edges(:,e)=[a,b]
        end if
        m%cell_edges(k,c)=e
        m%cell_sign(k,c)=merge(1,-1,edges(1,e)==a)
      end do
    end do
    m%edges=edges(:,:ne)
    allocate(m%vector(3,ne),m%center(3,ne),m%length(ne),m%direction(3,ne))
    do e=1,ne
      a=m%edges(1,e); b=m%edges(2,e)
      m%vector(:,e)=m%xyz(:,b)-m%xyz(:,a)
      m%length(e)=norm(m%vector(:,e))
      call require(m%length(e)>0,'zero-length edge')
      m%direction(:,e)=m%vector(:,e)/m%length(e)
      m%center(:,e)=(m%xyz(:,b)+m%xyz(:,a))/2
    end do
    allocate(m%nodes(nv),m%branches(ne),nn(nv),nb(ne))
    nn=0; nb=0
    do c=1,nq
      do k=1,4
        a=m%quads(k,c); e=m%cell_edges(k,c)
        nn(a)=nn(a)+1; nb(e)=nb(e)+1
      end do
    end do
    call require(all(nn>0),'mesh contains nodes not belonging to quadrangles')
    do a=1,nv
      allocate(m%nodes(a)%p(3,4,nn(a)))
    end do
    do e=1,ne
      allocate(m%branches(e)%p(3,4,nb(e)))
    end do
    nn=0; nb=0
    do c=1,nq
      p=m%xyz(:,m%quads(:,c)); center=sum(p,dim=2)/4
      area=patch_area(p)
      call require(area>0,'zero-area quadrangle')
      do k=1,4
        mid(:,k)=(p(:,k)+p(:,mod(k,4)+1))/2
      end do
      do k=1,4
        nxt=mod(k,4)+1; prv=mod(k+2,4)+1; opp=mod(k+1,4)+1
        a=m%quads(k,c); e=m%cell_edges(k,c)
        nn(a)=nn(a)+1; nb(e)=nb(e)+1
        m%nodes(a)%p(:,:,nn(a))=reshape([p(:,k),mid(:,k),center,mid(:,prv)],[3,4])
        m%branches(e)%p(:,:,nb(e))=reshape([p(:,k),p(:,nxt),mid(:,nxt),mid(:,prv)],[3,4])
      end do
    end do
    do a=1,nv
      do j=1,size(m%nodes(a)%p,3)
        m%nodes(a)%area=m%nodes(a)%area+patch_area(m%nodes(a)%p(:,:,j))
      end do
      call require(m%nodes(a)%area>0,'invalid node support')
    end do
    do e=1,ne
      do j=1,size(m%branches(e)%p,3)
        m%branches(e)%area=m%branches(e)%area+patch_area(m%branches(e)%p(:,:,j))
      end do
      call require(m%branches(e)%area>0,'invalid edge support')
    end do
    write(*,'(a,i0,a,i0,a,i0)') 'Mesh: nodes=',nv,', quads=',nq,', edges=',ne
  end subroutine
  integer function resolve_node(m,id,xyz,has_xyz) result(node)
    type(mesh_type), intent(in) :: m
    integer, intent(in) :: id
    real(dp), intent(in) :: xyz(3)
    logical, intent(in) :: has_xyz
    real(dp) :: best,d
    integer :: j
    if(id>=0) then
      node=id+1
    else
      call require(has_xyz,'specify a node index or coordinates')
      best=huge(0.0_dp); node=0
      do j=1,size(m%xyz,2)
        d=norm(m%xyz(:,j)-xyz)
        if(d<best) then
          best=d; node=j
        end if
      end do
    end if
    call require(node>=1.and.node<=size(m%xyz,2),'node index outside mesh')
  end function
end module
