module li_com

  implicit none

! li.in primary input variables
  integer :: nx,ny,nt
  real(8) :: tol,lx,ly,dt,theta,amp
  integer :: ni,initphi,ninit
  real(8) :: kapn,kapt,tets,memi
  integer :: bounded,enlin,wnlin,openmp,isolate,zflow
  real(8) :: xshape,yshape
  integer :: nrec,nprint,nmode
  integer,dimension(:,:),allocatable :: modeindices

! additional parameters
  complex(8) :: IU
  real(8) :: dx,dy
  real(8) :: sdt,cdt,sth,cth
  real(8) :: totvol,n0
  integer :: timestep
  integer :: iseed
  real(8) :: pi,pi2
  real(8) :: res

! timing variables
  integer :: wall_start,wall_finish
  real(8) :: wall_total

! fft variables
  complex(8),dimension(:),allocatable :: tmpx
  complex(8),dimension(:),allocatable :: tmpy

! particle array declarations
  real(8),dimension(:),allocatable :: x0,y0,vx0,vy0,vz0
  real(8),dimension(:),allocatable :: x1,y1,vx1,vy1,vz1
  real(8),dimension(:),allocatable :: w0,w1

! grid array declarations
  real(8),dimension(:,:),allocatable :: den0,den1,dent
  real(8),dimension(:,:),allocatable :: phi0,phi1,coeff
  real(8),dimension(:,:),allocatable :: ex0,ey0
  real(8),dimension(:,:),allocatable :: ex1,ey1
  real(8),dimension(:,:),allocatable :: tempxy

! ky=0 quantities
  real(8),dimension(:),allocatable :: uy,uz,px,py,dp,er

! mode history allocation
  complex(8),dimension(:),allocatable :: phihist,denhist,temphist


  save

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine init_com

  implicit none

! fft allocation
  allocate(tmpx(0:nx-1),tmpy(0:ny-1))

! particle allocation
  allocate(x0(1:ni),y0(1:ni))
  allocate(vx0(1:ni),vy0(1:ni),vz0(1:ni))
  allocate(x1(1:ni),y1(1:ni))
  allocate(vx1(1:ni),vy1(1:ni),vz1(1:ni))
  allocate(w0(1:ni),w1(1:ni))

! grid allocation
  allocate(den0(0:nx,0:ny),den1(0:nx,0:ny),dent(0:nx,0:ny))
  allocate(phi0(0:nx,0:ny),phi1(0:nx,0:ny),coeff(0:nx-1,0:ny-1))
  allocate(ex0(0:nx,0:ny),ey0(0:nx,0:ny))
  allocate(ex1(0:nx,0:ny),ey1(0:nx,0:ny))
  allocate(tempxy(0:nx,0:ny))

! ky=0 quantities
  allocate(uy(0:nx),uz(0:nx),px(0:nx),py(0:nx),dp(0:nx),er(0:nx))

! history allocation
  allocate(phihist(1:nmode),denhist(1:nmode),temphist(1:nmode))
  allocate(modeindices(2,nmode))

end

!-----------------------------------------------------------------------

subroutine finalize_com

  implicit none

! fft deallocation
  deallocate(tmpx,tmpy)

! particle deallocation
  deallocate(x0,y0)
  deallocate(vx0,vy0,vz0)
  deallocate(x1,y1)
  deallocate(vx1,vy1,vz1)
  deallocate(w0,w1)

! grid deallocation
  deallocate(den0,den1,dent)
  deallocate(phi0,phi1,coeff)
  deallocate(ex0,ey0)
  deallocate(ex1,ey1)
  deallocate(tempxy)

! ky=0 quantities
  deallocate(uy,uz,px,py,dp,er)

! history deallocation
  deallocate(phihist,denhist,temphist)
  deallocate(modeindices)


  end

!-----------------------------------------------------------------------

end
