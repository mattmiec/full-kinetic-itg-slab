!-----------------------------------------------------------------------
!---------2D Lorentz Ion ITG Delta-f Code-------------------------------
!-----------------------------------------------------------------------

program li
  use li_com
  use fft_wrapper
  use fcnt

  implicit none
  integer :: doprint

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,nproc,ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)

  !set start time in sec
  wall_start = mpi_wtime()

  !initialization
  call initialize
  call load
  call update
  call accumulate
  call field

  !main loop
  do timestep=1,nt

    doprint=0
    if (myid==0.and.mod(timestep,nprint).eq.0) doprint=1
    if (doprint.eq.1) print *
    if (doprint.eq.1) print *, 'timestep', timestep

    call epush

    !iterate over ipush
    res=1
    do while (res.gt.tol)
      call ipush
      call accumulate
      call field
      call residual
      if (doprint.eq.1) print *,'residual =',res
    end do

    call update

    call temperature

    !file output
    call modeout(phihist,'phist',11)
    call modeout(denhist,'dhist',12)
    call modeout(temphist,'thist',13)
    call diagnostics
    !call zdiagnostics
    if (mod(timestep,nrec).eq.0) then
      call gridout(phi,'phixy',14)
      call gridout(den,'denxy',15)
      call gridout(deni,'denii',17)
      call gridout(dene,'denee',18)
      call gridout(tempxy,'temxy',16)
    end if

  end do

  call finalize_com

  !set end time in sec
  wall_finish = mpi_wtime()

  call mpi_finalize(ierr)

  !time elapsed in sec
  wall_total=wall_finish-wall_start
  if (myid==0) then  
    print *
    print *, 'Wall Clock Seconds Elapsed = ',wall_total
  endif

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!---------initialization subroutines------------------------------------
!-----------------------------------------------------------------------

subroutine initialize

  implicit none
  character*(72) dumchar
  integer :: idum,i,j,ki,kj
  real(8) :: kx,ky,kp2,filt

  do i=0,nproc-1
    if (myid==i) then
      !read parameters from li.in
      open(115,file='li.in')
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) ni,ne
      read(115,*) dumchar
      read(115,*) nx,ny,nt,tol
      read(115,*) dumchar
      read(115,*) lx,ly,dt,theta
      read(115,*) dumchar
      read(115,*) amp,initphi,ninit
      read(115,*) dumchar
      read(115,*) kapni,kapti,kapne,kapte
      read(115,*) dumchar
      read(115,*) teti,memif,memip
      read(115,*) dumchar
      read(115,*) bounded,enlin,wnlin,odd
      read(115,*) dumchar
      read(115,*) isolate,zflow,xshape,yshape
      read(115,*) dumchar
      read(115,*) nrec,nprint,nmode
      tni=ni*nproc !total ions
      tne=ne*nproc !total electrons
      call init_com
      read(115,*) dumchar
      read(115,*) modeindices
      close(115)
    endif
    call mpi_barrier(mpi_comm_world,ierr)
  enddo

! compute remaining parameters
  IU=cmplx(0.,1.)
  pi=4.0*datan(dble(1.0))
  pi2=pi*2.0
  dx=lx/float(nx)
  dy=ly/float(ny)
  sdt=dsin(dt)
  cdt=dcos(dt)
  sth=dsin(theta)
  cth=dcos(theta)

  iseed=-(1777)
  idum=ran2(iseed)

! intialize fourier transform subroutines
  call ccfft('x',0,nx,tmpx)
  call ccfft('y',0,ny,tmpy)

! calculate coefficients for potential
  do i=0,nx-1
    do j=0,ny-1
      coeff(i,j)=0.
      ki = min(i,nx-i)
      kj = min(j,ny-j)
      kx = 2*pi*ki/lx
      ky = 2*pi*kj/ly
      kp2 = kx*kx + ky*ky
      filt = exp(-1*(xshape**2*kx**2+yshape**2*ky**2)**2)
      coeff(i,j) = filt/(memif*teti*kp2)
      ! zonal flow excluded if zflow != 1
      if ((zflow /= 1) .and. kj==0) coeff(i,j)=0.
      ! isolate 1,1 and 2,0 if isolate == 1
      if ((isolate == 1) .and. (.not.(((ki == 1).and.(kj == 1)) .or. ((ki == 2).and.(kj == 0))))) coeff(i,j)=0.
      !if (myid==0) then
      !  print*,'i = ',i,', j = ',j
      !  print*,'coeff =',coeff(i,j)
      !end if
    end do
  end do


end

!-----------------------------------------------------------------------

subroutine load

  implicit none
  integer :: m

  wi1 = 0.
  we1 = 0.

  ! ions
  do m=1,ni
!   load particle positions
    xi(m)=lx*revers(myid*ni+m,2)
    yi(m)=ly*(dble(myid*ni+m)-0.5)/dble(tni)
!   load maxwellian velocities
    vxi(m)=dinvnorm(revers(myid*ni+m,3))
    vyi(m)=dinvnorm(revers(myid*ni+m,5))
    vzi(m)=dinvnorm(revers(myid*ni+m,7))
!   initialize weights
    if (initphi /= 1) wi1(m)=amp*dsin(pi2*xi(m)/lx)*dsin(pi2*yi(m)/ly)
  end do

  ! electrons
  do m=1,ne
!   load particle positions
    xe1(m)=lx*revers(myid*ne+m,2)
    ye1(m)=ly*(dble(myid*ne+m)-0.5)/dble(tne)
!   load maxwellian velocities
    vpe(m)=dinvnorm(revers(myid*ne+m,3))/sqrt(memip)
!   initialize weights
    if (initphi /= 1) we1(m)=amp*dsin(pi2*xe1(m)/lx)*dsin(pi2*ye1(m)/ly)
  end do

end

!-----------------------------------------------------------------------
!---------deposit subroutines-------------------------------------------
!-----------------------------------------------------------------------

subroutine accumulate

  implicit none
  real(8) :: xpdx,ypdy,wx,wy
  real(8) :: mydeni(0:nx,0:ny),mydene(0:nx,0:ny)
  integer :: i,j,m

  denlast=den
  den=0
  mydeni=0
  mydene=0

  ! ions
  do m=1,ni
    xpdx=xi(m)/dx
    ypdy=yi(m)/dy
    i=int(xpdx)
    j=int(ypdy)
    wx=dble(i+1)-xpdx
    wy=dble(j+1)-ypdy
    mydeni(i,j)=mydeni(i,j)+wi1(m)*wx*wy
    mydeni(i+1,j)=mydeni(i+1,j)+wi1(m)*(1.0-wx)*wy
    mydeni(i,j+1)=mydeni(i,j+1)+wi1(m)*wx*(1.0-wy)
    mydeni(i+1,j+1)=mydeni(i+1,j+1)+wi1(m)*(1.0-wx)*(1.0-wy)
  end do

  ! electrons
  do m=1,ne
    xpdx=xe1(m)/dx
    ypdy=ye1(m)/dy
    i=int(xpdx)
    j=int(ypdy)
    wx=dble(i+1)-xpdx
    wy=dble(j+1)-ypdy
    mydene(i,j)=mydene(i,j)+we1(m)*wx*wy
    mydene(i+1,j)=mydene(i+1,j)+we1(m)*(1.0-wx)*wy
    mydene(i,j+1)=mydene(i,j+1)+we1(m)*wx*(1.0-wy)
    mydene(i+1,j+1)=mydene(i+1,j+1)+we1(m)*(1.0-wx)*(1.0-wy)
  end do

  call mpi_allreduce(mydeni,deni,(nx+1)*(ny+1),mpi_real8,mpi_sum,mpi_comm_world,ierr)
  call mpi_allreduce(mydene,dene,(nx+1)*(ny+1),mpi_real8,mpi_sum,mpi_comm_world,ierr)

  !divide by particles per cell
  deni=deni*dble(nx)*dble(ny)/dble(tni)
  dene=dene*dble(nx)*dble(ny)/dble(tne)

  do i=0,nx
    deni(i,0)=deni(i,0)+deni(i,ny)
    deni(i,ny)=deni(i,0)
    dene(i,0)=dene(i,0)+dene(i,ny)
    dene(i,ny)=dene(i,0)
  end do

  do j=0,ny
    deni(0,j)=deni(0,j)+deni(nx,j)
    deni(nx,j)=deni(0,j)
    dene(0,j)=dene(0,j)+dene(nx,j)
    dene(nx,j)=dene(0,j)
  end do

  den = deni - dene

end

!-----------------------------------------------------------------------
!--------field subroutines----------------------------------------------
!-----------------------------------------------------------------------

subroutine field

  implicit none
  integer :: i,j,ki,kj
  real(8) :: kx,ky
  complex(8) :: phit(0:nx-1,0:ny-1)
  complex(8) :: ext(0:nx-1,0:ny-1),eyt(0:nx-1,0:ny-1)

  !set potential equal to density and transform to k-space
  phit = den(0:nx-1,0:ny-1)

  do j=0,ny-1
    call ccfft('x',-1,nx,phit(:,j))
  end do

  do i=0,nx-1
    call ccfft('y',-1,ny,phit(i,:))
  end do

  !normalize
  phit=phit/nx/ny

  !record selected density modes
  do i=1,nmode
    denhist(i) = phit(modeindices(1,i),modeindices(2,i))
  end do

  !(optional) enforce phi=0 at x=0 and x=lx/2
  if (bounded==1) then
    do i=1,nx/2
      do j=0,ny-1
        phit(i,j)=(phit(i,j)-phit(nx-i,j))/2
        phit(nx-i,j)=-1*phit(i,j)
      end do
    end do
  end if

  !calculate phi with coefficients calculated during initialization
  do i=1,nx-1
    do j=0,ny-1
      phit(i,j) = coeff(i,j)*phit(i,j)
    end do
  end do

  !initialize if initphi
  if ((timestep.le.ninit).and.(initphi.eq.1)) then
      phit = 0.
      phit(1,1)=amp
      phit(nx-1,1)=-amp
      phit(1,ny-1)=-amp
      phit(nx-1,ny-1)=amp
  end if

  !record selected modes
  do i=1,nmode
    phihist(i) = phit(modeindices(1,i),modeindices(2,i))
  end do

  !calculate e-field
  do i=0,nx-1
    do j=0,ny-1
      if (i<=nx/2) kx = 2*pi*i/lx
      if (i>nx/2) kx = -2*pi*(nx-i)/lx
      if (j<=ny/2) ky = 2*pi*j/ly
      if (j>ny/2) ky = -2*pi*(ny-j)/ly
      ext(i,j) = -IU*kx*teti*phit(i,j)
      eyt(i,j) = -IU*ky*teti*phit(i,j)
    end do
  end do

  do i=0,nx-1
    call ccfft('y',1,ny,phit(i,:))
    call ccfft('y',1,ny,ext(i,:))
    call ccfft('y',1,ny,eyt(i,:))

  end do

  do j=0,ny-1
    call ccfft('x',1,nx,phit(:,j))
    call ccfft('x',1,nx,ext(:,j))
    call ccfft('x',1,nx,eyt(:,j))

  end do

  !store final phi,e-field
  phi(0:nx-1,0:ny-1)=real(phit)
  ex(0:nx-1,0:ny-1)=real(ext)
  ey(0:nx-1,0:ny-1)=real(eyt)

  !periodic boundaries
  phi(:,ny)=phi(:,0)
  phi(nx,:)=phi(0,:)
  ex(:,ny)=ex(:,0)
  ex(nx,:)=ex(0,:)
  ey(:,ny)=ey(:,0)
  ey(nx,:)=ey(0,:)

  return

end

!-----------------------------------------------------------------------

subroutine residual

  implicit none
  integer :: i,j
  real(8) :: norm

  res=0
  norm=0

  do i=0,nx
    do j=0,ny
      res = res + (den(i,j)-denlast(i,j))**2
      norm = norm + den(i,j)**2
    end do
  end do
  res = (res/norm)**.5

end

!-----------------------------------------------------------------------
!--------particle push subroutines--------------------------------------
!-----------------------------------------------------------------------

subroutine epush

  implicit none
  integer :: m,i,j
  real(8) :: vdv,kap,edv
  real(8) :: ax,ay
  real(8) :: wx,wy
  real(8) :: xpdx,ypdy
  real(8) :: vxt,vyt !temp velocity storage

  ! ions
  do m=1,ni
    ! interpolation weights
    xpdx=xi(m)/dx
    ypdy=yi(m)/dy
    i=int(xpdx)
    j=int(ypdy)
    wx=dble(i+1)-xpdx
    wy=dble(j+1)-ypdy
    ! interpolate e-field
    ax=ex(i,j)*wx*wy+ex(i+1,j)*(1.0-wx)*wy+&
      ex(i,j+1)*wx*(1.0-wy)+ex(i+1,j+1)*(1.0-wx)*(1.0-wy)
    ay=ey(i,j)*wx*wy+ey(i+1,j)*(1.0-wx)*wy+&
      ey(i,j+1)*wx*(1.0-wy)+ey(i+1,j+1)*(1.0-wx)*(1.0-wy)
    ! 1/2 velocity push
    vxi(m)=vxi(m)+.5*dt*ax*enlin
    vyi(m)=vyi(m)+.5*dt*ay*enlin
    ! full velocity rotation (theta,dt,-theta)
    vxt=cdt*vxi(m)+sdt*cth*vyi(m)-sdt*sth*vzi(m)
    vyt=-1.0*sdt*cth*vxi(m)+(cdt*cth**2+sth**2)*vyi(m)&
      +(-1.0*cdt*sth*cth+sth*cth)*vzi(m)
    vzi(m)=sdt*sth*vxi(m)+(-1.0*cdt*sth*cth+sth*cth)*vyi(m)&
      +(cdt*sth**2+cth**2)*vzi(m)
    vxi(m)=vxt
    vyi(m)=vyt
    ! 1/2 velocity push
    vxi(m)=vxi(m)+.5*dt*ax*enlin
    vyi(m)=vyi(m)+.5*dt*ay*enlin
    ! weight equation terms
    vdv=vxi(m)**2+vyi(m)**2+vzi(m)**2
    edv=vxi(m)*ax+vyi(m)*ay
    kap=kapni+kapti*(.5*vdv-1.5)
    ! explicit 1/2 weight advance
    wi0(m)=wi0(m)+.5*dt*(1-wi0(m)*wnlin)*(edv+cth*ay*kap)
    ! full position advance
    xi(m)=xi(m)+dt*vxi(m)
    yi(m)=yi(m)+dt*vyi(m)
    ! periodic boundaries
    xi(m)=xi(m)-lx*dble(floor(xi(m)/lx))
    yi(m)=yi(m)-ly*dble(floor(yi(m)/ly))
  end do

  ! electrons
  do m=1,ne
    ! interpolation weights
    xpdx=xe0(m)/dx
    ypdy=ye0(m)/dy
    i=int(xpdx)
    j=int(ypdy)
    wx=dble(i+1)-xpdx
    wy=dble(j+1)-ypdy
    ! interpolate e-field
    ax=ex(i,j)*wx*wy+ex(i+1,j)*(1.0-wx)*wy+&
      ex(i,j+1)*wx*(1.0-wy)+ex(i+1,j+1)*(1.0-wx)*(1.0-wy)
    ay=ey(i,j)*wx*wy+ey(i+1,j)*(1.0-wx)*wy+&
      ey(i,j+1)*wx*(1.0-wy)+ey(i+1,j+1)*(1.0-wx)*(1.0-wy)
    ! full parallel velocity push
    vpe(m)=vpe(m)-dt*ay*sth*enlin/memip
    ! weight equation terms
    vdv=vpe(m)**2
    kap=kapne+kapte*(.5*memip*vdv-1.5)
    ! explicit 1/2 weight advance
    we0(m)=we0(m)-.5*dt*(1-we0(m)*wnlin)*(sth*ay*vpe(m)-cth*ay*kap)
    ! explicit part of position advance
    xe0(m)=xe0(m)+.5*dt*ay*cth
    ye0(m)=ye0(m)-.5*dt*ax*cth+dt*sth*vpe(m)
    ! periodic boundaries
    xe0(m)=xe0(m)-lx*dble(floor(xe0(m)/lx))
    ye0(m)=ye0(m)-ly*dble(floor(ye0(m)/ly))
  end do

end

!-----------------------------------------------------------------------

subroutine ipush

  implicit none
  integer :: m,i,j
  real(8) :: vdv,kap,edv
  real(8) :: ax,ay
  real(8) :: wx,wy
  real(8) :: xpdx,ypdy

  ! ions
  do m=1,ni
    xpdx=xi(m)/dx
    ypdy=yi(m)/dy
    i=int(xpdx)
    j=int(ypdy)
    wx=dble(i+1)-xpdx
    wy=dble(j+1)-ypdy
    ! interpolate e-field
    ax=ex(i,j)*wx*wy+ex(i+1,j)*(1.0-wx)*wy+&
      ex(i,j+1)*wx*(1.0-wy)+ex(i+1,j+1)*(1.0-wx)*(1.0-wy)
    ay=ey(i,j)*wx*wy+ey(i+1,j)*(1.0-wx)*wy+&
      ey(i,j+1)*wx*(1.0-wy)+ey(i+1,j+1)*(1.0-wx)*(1.0-wy)
    ! weight equation terms
    vdv=vxi(m)**2+vyi(m)**2+vzi(m)**2
    edv=vxi(m)*ax+vyi(m)*ay
    kap=kapni+kapti*(.5*vdv-1.5)
    ! implicit weight advance
    wi1(m)=wi0(m)+.5*dt*(1-wi1(m)*wnlin)*(edv+cth*ay*kap)
  end do

  ! electrons
  do m=1,ne
    ! interpolation weights
    xpdx=xe1(m)/dx
    ypdy=ye1(m)/dy
    i=int(xpdx)
    j=int(ypdy)
    wx=dble(i+1)-xpdx
    wy=dble(j+1)-ypdy
    ! interpolate e-field
    ax=ex(i,j)*wx*wy+ex(i+1,j)*(1.0-wx)*wy+&
      ex(i,j+1)*wx*(1.0-wy)+ex(i+1,j+1)*(1.0-wx)*(1.0-wy)
    ay=ey(i,j)*wx*wy+ey(i+1,j)*(1.0-wx)*wy+&
      ey(i,j+1)*wx*(1.0-wy)+ey(i+1,j+1)*(1.0-wx)*(1.0-wy)
    ! weight equation terms
    vdv=vpe(m)**2
    kap=kapne+kapte*(.5*memip*vdv-1.5)
    ! implicit part of weight advance
    we1(m)=we0(m)-.5*dt*(1-we1(m)*wnlin)*(sth*ay*vpe(m)-cth*ay*kap)
    ! implicit part of position advance
    xe1(m)=xe0(m)+.5*dt*ay*cth
    ye1(m)=ye0(m)-.5*dt*ax*cth
    ! periodic boundaries
    xe1(m)=xe1(m)-lx*dble(floor(xe1(m)/lx))
    ye1(m)=ye1(m)-ly*dble(floor(ye1(m)/ly))
  end do

end

!-----------------------------------------------------------------------

subroutine update

  implicit none

  wi0 = wi1
  xe0 = xe1
  ye0 = ye1
  we0 = we1

end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine temperature

  !calculate tperp on grid

  implicit none
  real(8) :: xpdx,ypdy,wx,wy
  real(8) :: mytempxy(0:nx,0:ny)
  integer :: i,j,m
  complex(8) :: ktemp(0:nx-1,0:ny-1)

  tempxy=0

  do m=1,ni
    xpdx=xi(m)/dx
    ypdy=yi(m)/dy
    i=int(xpdx)
    j=int(ypdy)
    wx=dble(i+1)-xpdx
    wy=dble(j+1)-ypdy
    mytempxy(i,j)=mytempxy(i,j)+wi1(m)*wx*wy*(vxi(m)**2+(cth*vyi(m)-sth*vzi(m))**2)
    mytempxy(i+1,j)=mytempxy(i+1,j)+wi1(m)*(1.0-wx)*wy*(vxi(m)**2+(cth*vyi(m)-sth*vzi(m))**2)
    mytempxy(i,j+1)=mytempxy(i,j+1)+wi1(m)*wx*(1.0-wy)*(vxi(m)**2+(cth*vyi(m)-sth*vzi(m))**2)
    mytempxy(i+1,j+1)=mytempxy(i+1,j+1)+wi1(m)*(1.0-wx)*(1.0-wy)*(vxi(m)**2+(cth*vyi(m)-sth*vzi(m))**2)
  end do
 
  call mpi_allreduce(mytempxy,tempxy,(nx+1)*(ny+1),mpi_real8,mpi_sum,mpi_comm_world,ierr)

  !divide by particles per cell
  tempxy=tempxy*dble(nx)*dble(ny)/dble(tni)

  do i=0,nx
    tempxy(i,0)=tempxy(i,0)+tempxy(i,ny)
    tempxy(i,ny)=tempxy(i,0)
  end do

  do j=0,ny
    tempxy(0,j)=tempxy(0,j)+tempxy(nx,j)
    tempxy(nx,j)=tempxy(0,j)
  end do

  !transform to record fourier space components
  ktemp = tempxy(0:nx-1,0:ny-1)

  do j=0,ny-1
    call ccfft('x',-1,nx,ktemp(:,j))
  end do

  do i=0,nx-1
    call ccfft('y',-1,ny,ktemp(i,:))
  end do

  ktemp = ktemp/nx/ny

  !record selected modes
  do i=1,nmode
    temphist(i) = ktemp(modeindices(1,i),modeindices(2,i))
  end do


end

!-----------------------------------------------------------------------
!--------history subroutines--------------------------------------------
!-----------------------------------------------------------------------

subroutine modeout(hist,fl,id)

  !record components of phi given by modehist

  implicit none
  integer :: i,id
  character*5 :: fl
  character*70 :: flnm
  complex(8) :: hist(1:nmode)

  if (myid==0) then
    flnm=fl//'.out'
    open(id,file=flnm,form='formatted',status='unknown',&
      position='append')
    write(id,'(f8.2)',advance="no") dt*timestep
    do i=1,nmode
      write(id,'(a2,e13.6,a2,e13.6)',advance="no") '  ',real(hist(i)),&
        '  ',imag(hist(i))
    end do
    write(id,*)
    endfile id
    close(id)
  endif

end

!-----------------------------------------------------------------------

subroutine gridout(u,fl,id)

  !record all values of grid quantity

  implicit none
  integer :: i,j,id
  real(8) :: u(0:nx,0:ny)
  character*5 :: fl
  character*70 :: flnm

  if (myid==0) then
    flnm=fl//'.out'
    open(id,file=flnm,form='formatted',status='unknown',&
      position='append')
    do j=0,ny
      do i=0,nx
        write(id,'(e11.4)') u(i,j)
      end do
    end do
    endfile id
    close(id)
  endif
  
end

!-----------------------------------------------------------------------

subroutine diagnostics

  implicit none
  integer :: id,m,i,j
  real(8) ::xpdx,ypdy,wx,wy,qx,myqx,w2i,myw2i,w2e,myw2e,mykei,kei,mykee,kee
  character*5 :: fl
  character*70 :: flnm

  id=89

  qx=0
  w2i=0
  w2e=0
  kei=0
  kee=0

  myqx=0
  myw2i=0
  myw2e=0
  mykei=0
  mykee=0
  do m=1,ni
    !net ion heat flux in x-direction
    myqx = myqx + wi1(m)*vxi(m)*(vxi(m)**2+vyi(m)**2+vzi(m)**2)
    !weight squared sum
    myw2i = myw2i + wi1(m)**2
    myw2e = myw2e + we1(m)**2
    !kinetic energies
    mykei = mykei + 0.5*wi1(m)*(vxi(m)**2+vyi(m)**2+vzi(m)**2)
    mykee = mykee + 0.5*we1(m)*memip*vpe(m)**2
  end do

  call mpi_allreduce(myqx,qx,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
  call mpi_allreduce(myw2i,w2i,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
  call mpi_allreduce(myw2e,w2e,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
  call mpi_allreduce(mykei,kei,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
  call mpi_allreduce(mykee,kee,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

  qx=qx/dble(tni)
  w2i=w2i/dble(tni)
  w2e=w2e/dble(tne)
  kei=kei/dble(tni)
  kee=kee/dble(tne)

  if (myid==0) then
    flnm='diagn.out'
    open(id,file=flnm,form='formatted',status='unknown',&
      position='append')
    if (timestep==1) write(id) 't  qx  w2i  w2e  kei  kee'
    write(id,'(f8.2)',advance="no") dt*timestep
    write(id,'(a2,e13.6,a2,e13.6,a2,e13.6,a2,e13.6,a2,e13.6)') '  ',qx,'  ',w2i,'  ',w2e,'  ',kei,'  ',kee
    endfile id
    close(id)
  endif

end

!-----------------------------------------------------------------------

subroutine zdiagnostics

  implicit none
  integer :: i,m,il,ir,uyid,uzid,pxid,pyid
  character*70 :: uyflnm,uzflnm,pxflnm,pyflnm
  real(8) :: xpdx,wx,kx,p(0:nx)
  real(8) :: myuy(0:nx),myuz(0:nx),mypx(0:nx),mypy(0:nx)

  uy=0
  myuy=0
  uz=0
  myuz=0
  px=0
  mypx=0
  py=0
  mypy=0

  do m=1,ni
    xpdx=xi(m)/dx
    i=int(xpdx)
    wx=dble(i+1)-xpdx
    myuy(i)=myuy(i)+wi1(m)*wx*vyi(m)
    myuy(i+1)=myuy(i+1)+wi1(m)*(1.0-wx)*vyi(m)
    myuz(i)=myuz(i)+wi1(m)*wx*vzi(m)
    myuz(i+1)=myuz(i+1)+wi1(m)*(1.0-wx)*vzi(m)
    mypx(i)=mypx(i)+wi1(m)*wx*vxi(m)**2
    mypx(i+1)=mypx(i+1)+wi1(m)*(1.0-wx)*vxi(m)**2
    mypy(i)=mypy(i)+wi1(m)*wx*(cth*vyi(m)-sth*vzi(m))**2
    mypy(i+1)=mypy(i+1)+wi1(m)*(1.0-wx)*(cth*vyi(m)-sth*vzi(m))**2
  enddo

  call mpi_allreduce(myuy,uy,nx+1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
  call mpi_allreduce(myuz,uz,nx+1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
  call mpi_allreduce(mypx,px,nx+1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
  call mpi_allreduce(mypy,py,nx+1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

  !divide by particles per column
  uy=uy*dble(nx)/dble(tni)
  uz=uz*dble(nx)/dble(tni)
  px=px*dble(nx)/dble(tni)
  py=py*dble(nx)/dble(tni)

  !periodic boundaries
  uy(0)=uy(0)+uy(nx)
  uz(0)=uz(0)+uz(nx)
  px(0)=px(0)+px(nx)
  py(0)=py(0)+py(nx)
  uy(nx)=uy(0)
  uz(nx)=uz(0)
  px(nx)=px(0)
  py(nx)=py(0)

  uyid=99
  uzid=98
  pxid=97
  pyid=96

  uyflnm='uy.out'
  uzflnm='uz.out'
  pxflnm='px.out'
  pyflnm='py.out'

  if (myid==0) then
    open(uyid,file=uyflnm,form='formatted',status='unknown',&
      position='append')
    write(uyid,'(f8.2)',advance='no') dt*timestep
    do i=0,nx
      write(uyid,'(a2,e13.6)',advance='no') '  ',uy(i)
    end do
    write(uyid,*)
    endfile uyid
    close(uyid)
  
    open(uzid,file=uzflnm,form='formatted',status='unknown',&
      position='append')
    write(uzid,'(f8.2)',advance='no') dt*timestep
    do i=0,nx
      write(uzid,'(a2,e13.6)',advance='no') '  ',uz(i)
    end do
    write(uzid,*)
    endfile uzid
    close(uzid)
  
    open(pxid,file=pxflnm,form='formatted',status='unknown',&
      position='append')
    write(pxid,'(f8.2)',advance='no') dt*timestep
    do i=0,nx
      write(pxid,'(a2,e13.6)',advance='no') '  ',px(i)
    end do
    write(pxid,*)
    endfile pxid
    close(pxid)
  
    open(pyid,file=pyflnm,form='formatted',status='unknown',&
      position='append')
    write(pyid,'(f8.2)',advance='no') dt*timestep
    do i=0,nx
      write(pyid,'(a2,e13.6)',advance='no') '  ',py(i)
    end do
    write(pyid,*)
    endfile pyid
    close(pyid)
  endif
 
  end

!-----------------------------------------------------------------------

end
