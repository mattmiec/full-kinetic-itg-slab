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

  tstep = 0
  call output

  !main loop
  do tstep=1,nt

    doprint=0
    if (myid==0.and.mod(tstep,nprint).eq.0) doprint=1
    if (doprint.eq.1) print *
    if (doprint.eq.1) print *, 'tstep', tstep

    !explicit part of push
    call epush

    !iterate over implicit part of push
    res=1
    do while (res.gt.tol)
      call ipush
      call accumulate
      call field
      call residual
      if (doprint.eq.1) print *,'residual =',res
    end do

    !set new t0 arrays to old t1 arrays
    call update

    !output

    ! modeout, gridout, diagnostics
    call output

!   call modeout(phihist,'phist',11)
!   call modeout(denhist,'dhist',12)
!   call modeout(temphist,'thist',13)
!   if (mod(tstep,nrec).eq.0) then
!     call gridout(phi,'phixy',14)
!     if (dke /= 1) then
!       call gridout(den,'denxy',15)
!     else
!       call gridout(deni,'denii',17)
!       call gridout(dene,'denee',18)
!     end if
!     call gridout(tempxy,'temxy',16)
!     call gendiagnostics
!     !call ionvelocity
!   end if

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
      read(115,*) dumchar
      read(115,*) ni,ne,nt,dt,tol
      read(115,*) dumchar
      read(115,*) nx,ny,lx,ly,theta,order
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) amp,initphi,ninit,rand
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) kapni,kapti,kapne,kapte
      read(115,*) dumchar
      read(115,*) teti,memi,fki,gke
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) eperpi,epari,weighti,eperpe,epare,weighte
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) reflect,isolate,zflow,oddmodes,xshape,yshape
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) nrec,nprint,nmode
      tni=ni*nproc !total ions
      tne=ne*nproc !total electrons
      call init_com
      read(115,*) dumchar
      read(115,*) modeindices
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) bamp
      read(115,*) dumchar
      read(115,*) bmode
      close(115)
    endif
    call mpi_barrier(mpi_comm_world,ierr)
  enddo

! compute remaining parameters
  IU=cmplx(0.,1.)
  pi=4.0*datan(dble(1.0))
  pi2=pi*2.0
  dx=lx/dble(nx)
  dy=ly/dble(ny)
  sdt=dsin(memi*dt)
  cdt=dcos(memi*dt)
  sth=dsin(theta)
  cth=dcos(theta)

  iseed=-(1777)
  idum=ran2(iseed)

! intialize fourier transform subroutines
  call ccfft('x',0,nx)
  call ccfft('y',0,ny)

  call dsinf(1,nx-1)
  call dcosf(1,nx+1)

! calculate coefficients for potential
  do i=0,nx-1
    do j=0,ny-1
      coeff(i,j)=0.
      if (reflect == 1) then
        ki = i
        kx = pi*dble(ki)/lx
      else
        ki = min(i,nx-i)
        kx = pi2*dble(ki)/lx
      end if
      kj = min(j,ny-j)
      ky = pi2*dble(kj)/ly
      kp2 = kx*kx + ky*ky
      filt = exp(-1*(xshape**2*kx**2+yshape**2*ky**2)**2)
      !default solution to Poisson equation
      if (ki /= 0) coeff(i,j) = filt/(teti*kp2)
      ! use adiabatic response for k_par /= 0
      if ((fki /= 1) .and. ((kj /= 0) .or. (ki /= 0))) coeff(i,j) = filt/(1.0 + teti*kp2)
      ! zonal flow excluded if zflow != 1
      if ((zflow /= 1) .and. kj==0) coeff(i,j)=0.
      ! isolate 1,1 and 2,0 if isolate == 1
      if ((isolate == 1) .and. (.not.(((ki == 1).and.(kj == 1)) .or. ((ki == 2).and.(kj == 0))))) coeff(i,j)=0.
      ! remove odd kx modes for ky=0
      if ((oddmodes /= 1) .and. (kj == 0) .and. (mod(ki,2) == 1)) coeff(i,j) = 0
      if (myid==0) then
        print*,'coeff(',i,',',j,') = ',coeff(i,j)
      end if
    end do
  end do


end

!-----------------------------------------------------------------------

subroutine load

  implicit none
  integer :: m

  ! ions
  if (fki == 1) then
    wi1 = 0.
    wpi1 = 1.
    do m=1,ni
  !   load particle positions
      xi(m)=lx*revers(myid*ni+m,2)
      yi(m)=ly*(dble(myid*ni+m)-0.5)/dble(tni)
  !   load maxwellian velocities
      vxi(m)=dinvnorm(revers(myid*ni+m,3))*sqrt(memi)
      vyi(m)=dinvnorm(revers(myid*ni+m,5))*sqrt(memi)
      vpari(m)=dinvnorm(revers(myid*ni+m,7))*sqrt(memi)
  !   initialize weights
      if (initphi /= 1) then
        if (rand /= 0) then
         wi1(m) = 2. * amp * (ran2(rand) - 0.5)
        else
          if (reflect /= 1) wi1(m)=amp*dsin(pi2*xi(m)/lx)*dsin(pi2*yi(m)/ly)
          if (reflect == 1) wi1(m)=amp*dsin(pi*xi(m)/lx)*dsin(pi2*yi(m)/ly)
        end if
      end if
    end do
  end if

! electrons
  we1 = 0.
  wpe1 = 1.
  do m=1,ne
!   load particle positions
    xe1(m)=lx*revers(myid*ne+m,11)
    ye1(m)=ly*(dble(myid*ne+m)-0.5)/dble(tne)
!   load maxwellian velocities
    vpare(m)=dinvnorm(revers(myid*ne+m,13))
    if (gke == 1) then
      mue(m) = 2.*(dinvnorm(revers(myid*ne+m,17))**2)
    else
      mue(m) = 1.0
    end if
!   initialize weights
    if (initphi /= 1) then
      if (rand /= 0) then
       we1(m) = 2. * amp * (ran2(rand) - 0.5)
      else
        if (reflect /=1) we1(m)=amp*dsin(pi2*xe1(m)/lx)*dsin(pi2*ye1(m)/ly)
        if (reflect ==1) we1(m)=amp*dsin(pi*xe1(m)/lx)*dsin(pi2*ye1(m)/ly)
      end if
    end if
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
  if (fki == 1) then
    do m=1,ni
      call spline(1,xi(m),yi(m),wi1(m),mydeni)
    end do
  end if

  ! electrons
  if (gke == 1) then
    do m=1,ne
      call gyrospline(1,xe1(m),ye1(m),sqrt(mue(m)),we1(m),mydene)
    end do
  else
    do m=1,ne
      call spline(1,xe1(m),ye1(m),we1(m),mydene)
    end do
  end if

  if (fki == 1) call mpi_allreduce(mydeni,deni,(nx+1)*(ny+1),mpi_real8,mpi_sum,mpi_comm_world,ierr)
  call mpi_allreduce(mydene,dene,(nx+1)*(ny+1),mpi_real8,mpi_sum,mpi_comm_world,ierr)

  !divide by particles per cell
  if (fki == 1) deni=deni*dble(nx)*dble(ny)/dble(tni)
  dene=dene*dble(nx)*dble(ny)/dble(tne)

  do i=0,nx
    if (fki == 1) then
      deni(i,0)=deni(i,0)+deni(i,ny)
      deni(i,ny)=deni(i,0)
    end if
    dene(i,0)=dene(i,0)+dene(i,ny)
    dene(i,ny)=dene(i,0)
  end do

  if (reflect == 1) then
    do j=0,ny
      if (fki == 1) then
        deni(0,j)=0.
        deni(nx,j)=0.
      end if
      dene(0,j)=0.
      dene(nx,j)=0.
    end do
  else
    do j=0,ny
      if (fki == 1) then
        deni(0,j)=deni(0,j)+deni(nx,j)
        deni(nx,j)=deni(0,j)
      end if
      dene(0,j)=dene(0,j)+dene(nx,j)
      dene(nx,j)=dene(0,j)
    end do
  end if

  if (fki == 1) then
    den = deni - dene
  else
    den = -dene
  end if

end

!-----------------------------------------------------------------------
!--------field subroutines----------------------------------------------
!-----------------------------------------------------------------------

subroutine field

  implicit none
  integer :: i,j,ki,kj
  real(8) :: kx,ky
  real(8) :: phitr(1:nx-1,0:ny-1),extr(0:nx,0:ny-1),eytr(1:nx-1,0:ny-1)
  complex(8) :: phit(0:nx-1,0:ny-1),ext(0:nx-1,0:ny-1),eyt(0:nx-1,0:ny-1)

  !set potential equal to density and transform to k-space
  ! transform in x
  if (reflect == 1) then
    phitr = den(1:nx-1,0:ny-1)
    do j=0,ny-1
      call dsinf(0,nx-1,phitr(:,j))
    end do
    ! pad i=0 slice with zeros
    phit = 0.
    phit(1:nx-1,:) = phitr
  else
    phit = den(0:nx-1,0:ny-1)
    do j=0,ny-1
      call ccfft('x',-1,nx,phit(:,j))
    end do
  end if
  ! transform in y
  do i=0,nx-1
    call ccfft('y',-1,ny,phit(i,:))
  end do

  !normalize
  phit=phit/dble(nx)/dble(ny)
  if (reflect == 1) phit = 2.*phit

  !record selected density modes
  do i=1,nmode
    denhist(i) = phit(modeindices(1,i),modeindices(2,i))
  end do

  !enforce phi=0 at x=0 and x=lx/2
  if (reflect /= 1) then
    do j=0,ny-1
      phit(0,j) = 0.
      do i=1,nx/2
        phit(i,j)=(phit(i,j)-phit(nx-i,j))/2.
        phit(nx-i,j)=-1*phit(i,j)
      end do
    end do
  end if

  !calculate phi with coefficients calculated during initialization
  do i=0,nx-1
    do j=0,ny-1
      phit(i,j) = coeff(i,j)*phit(i,j)
    end do
  end do

  !initialize if initphi
  if ((tstep.le.ninit).and.(initphi.eq.1)) then
      phit = 0.
      phit(1,1)=amp
      if (reflect /= 1) phit(nx-1,1)=-amp
      phit(1,ny-1)=-amp
      if (reflect /= 1) phit(nx-1,ny-1)=amp
  end if

  !record selected modes
  do i=1,nmode
    phihist(i) = phit(modeindices(1,i),modeindices(2,i))
  end do

  !add background
  phit(bmode(1), bmode(2)) = phit(bmode(1), bmode(2)) + bamp
  phit(bmode(1), ny - bmode(2)) = phit(bmode(1), ny - bmode(2)) - bamp
  if (reflect /= 1) then
      phit(nx - bmode(1), bmode(2)) = phit(nx - bmode(1), bmode(2)) - bamp
      phit(nx - bmode(1), ny - bmode(2)) = phit(nx - bmode(1), ny - bmode(2)) + bamp
  end if

!  print *, 'bamp = ', bamp
!  print *, 'bmode(1) = ', bmode(1)
!  print *, 'bmode(2) = ', bmode(2)


  !calculate e-field, field energy
  fe = 0.0
  do i=0,nx-1
    do j=0,ny-1
      !calculate kx, ky
      if (reflect == 1) then
        kx = pi*dble(i)/lx
      else
        if (i<=nx/2) kx = pi2*dble(i)/lx
        if (i>nx/2) kx = -pi2*dble(nx-i)/lx
      end if
      if (j<=ny/2) ky = pi2*dble(j)/ly
      if (j>ny/2) ky = -pi2*dble(ny-j)/ly
      !calculate ex, ey
      ext(i,j) = -kx*teti*phit(i,j)
      if (reflect /= 1) ext(i,j) = IU*ext(i,j)
      eyt(i,j) = -IU*ky*teti*phit(i,j)
      !calculate fe
      fe = fe + .5 * (kx**2 + ky**2) * abs(phit(i,j))**2
    end do
  end do

  ! transform back
  ! transform in y
  do i=0,nx-1
    call ccfft('y',1,ny,phit(i,:))
    call ccfft('y',1,ny,ext(i,:))
    call ccfft('y',1,ny,eyt(i,:))
  end do

  if (reflect == 1) then

    ! ignore padded i=0 zeros, except extr add zeros at i=0, nx
    phitr          = real(phit(1:nx-1,:))
    extr(0,     :) = 0.
    extr(1:nx-1,:) = real(ext(1:nx-1,:))
    extr(nx,    :) = 0.
    eytr           = real(eyt(1:nx-1,:))

    ! transform in x
    do j=0,ny-1
      call dsinf(0, nx-1, phitr(:,j))
      call dcosf(0, nx+1, extr(:,j))
      call dsinf(0, nx-1, eytr(:,j))
    end do
       
    !store final phi,e-field
    phi(0,      0:ny-1) = 0.
    phi(1:nx-1, 0:ny-1) = phitr
    phi(nx,     0:ny-1) = 0.

    ex(:,       0:ny-1) = extr

    ey(0,       0:ny-1) = 0.
    ey(1:nx-1,  0:ny-1) = eytr
    ey(nx,      0:ny-1) = 0.

  else

    ! transform in x
    do j=0,ny-1
      call ccfft('x', 1,nx, phit(:,j))
      call ccfft('x', 1,nx, ext(:,j))
      call ccfft('x', 1,nx, eyt(:,j))
    end do
       
    !store final phi,e-field
    phi(0:nx-1, 0:ny-1) = real(phit)
    ex(0:nx-1,  0:ny-1) = real(ext)
    ey(0:nx-1,  0:ny-1) = real(eyt)
    phi(nx,     0:ny-1) = phi(0, 0:ny-1)
    ex(nx,      0:ny-1) = ex(0,  0:ny-1)
    ey(nx,      0:ny-1) = ey(0,  0:ny-1)

  end if

  ! enforce periodic boundary in y
  phi(:, ny) = phi(:, 0)
  ex(:,  ny) = ex(:,  0)
  ey(:,  ny) = ey(:,  0)

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

subroutine spline(io, x, y, w, grid)

  ! io = 0 for retrieve
  ! io = 1 for deposit
  ! x, y position
  ! w is weight or output field depending on io
  ! grid is field quantitiy to interpolate

  implicit none

  ! arguments
  integer :: io
  real(8) :: x,y,w
  real(8),dimension(0:nx,0:ny) :: grid

  ! temporary vars
  real(8) :: xpdx,ypdy,wx0,wx1,wx2,wy0,wy1,wy2
  integer :: i0,i1,i2,j0,j1,j2

  ! position in grid units
  xpdx = x / dx
  ypdy = y / dy

  ! interpolate e-field
  if (order == 2) then
    !quadratic see Birdsall Langdon p. 169

    i1 = nint(xpdx) ! central point of spline in x
    j1 = nint(ypdy) ! central point of spline in y

    ! calc weights
    wx0 = .5 * (.5 - (xpdx - dble(i1))) ** 2.
    wx1 = .75 - (xpdx - dble(i1)) ** 2.
    wx2 = .5 * (.5 + (xpdx - dble(i1))) ** 2.

    wy0 = .5 * (.5 - (ypdy - dble(j1))) ** 2.
    wy1 = .75 - (ypdy - dble(j1)) ** 2.
    wy2 = .5 * (.5 + (ypdy - dble(j1))) ** 2.

    ! calc indices
    if (reflect == 1) then
        if (i1 <= 0) then
            i0 = -1.0 * (i1 - 1)
            i2 = -1.0 * (i1 + 1)
            i1 = -i1
        else if (i1 >= nx) then
            i0 = nx - ((i1 - 1) - nx)
            i2 = nx - ((i1 + 1) - nx)
            i1 = nx - (i1 - nx)
        else
            i0 = i1 - 1
            i2 = i1 + 1
        end if
    else
        i0 = modulo(i1-1,nx)
        i1 = modulo(i1,nx)
        i2 = modulo(i1+1,nx)
    end if
    j0 = modulo(j1-1,ny)
    j1 = modulo(j1,ny)
    j2 = modulo(j1+1,ny)

    if (i1 <= 0) then
      if (reflect == 1) then
        i0 = modulo(i1+1,nx)
      else
        i0 = modulo(i1-1,nx)
      end if
      i2 = 1
    else if (i1 == nx) then
      if (reflect == 1) then
        i2 = nx-1
      else
        i2 = 1
      end if
      i0 = nx-1
    else
      i0 = i1 - 1
      i2 = i1 + 1
    end if

    if (j1 == 0) then
      j0 = ny-1
      j2 = 1
    else if (j1 == ny) then
      j2 = 1
      j0 = ny-1
    else
      j0 = j1 - 1
      j2 = j1 + 1
    end if

    if (io == 0) then
      w = grid(i0,j0) * wx0 * wy0 &
        + grid(i0,j1) * wx0 * wy1 &
        + grid(i0,j2) * wx0 * wy2 &
        + grid(i1,j0) * wx1 * wy0 &
        + grid(i1,j1) * wx1 * wy1 &
        + grid(i1,j2) * wx1 * wy2 &
        + grid(i2,j0) * wx2 * wy0 &
        + grid(i2,j1) * wx2 * wy1 &
        + grid(i2,j2) * wx2 * wy2
    else
      grid(i0,j0) = grid(i0,j0) + wx0 * wy0 * w
      grid(i0,j1) = grid(i0,j1) + wx0 * wy1 * w
      grid(i0,j2) = grid(i0,j2) + wx0 * wy2 * w
      grid(i1,j0) = grid(i1,j0) + wx1 * wy0 * w
      grid(i1,j1) = grid(i1,j1) + wx1 * wy1 * w
      grid(i1,j2) = grid(i1,j2) + wx1 * wy2 * w
      grid(i2,j0) = grid(i2,j0) + wx2 * wy0 * w
      grid(i2,j1) = grid(i2,j1) + wx2 * wy1 * w
      grid(i2,j2) = grid(i2,j2) + wx2 * wy2 * w
    end if

  else !linear

    i0 = int(xpdx)
    j0 = int(ypdy)

    !interpolation weights
    wx0 = dble(i0 + 1) - xpdx
    wx1 = 1. - wx0
    wy0 = dble(j0 + 1) - ypdy
    wy1 = 1. - wy0

    !in-bound indices
    i0 = modulo(i0,nx)
    i1 = modulo(i0+1,nx)
    j0 = modulo(j0,ny)
    j1 = modulo(j0+1,ny)

    if (io == 0) then
      w  = grid(i0,j0) * wx0 * wy0 &
         + grid(i1,j0) * wx1 * wy0 &
         + grid(i0,j1) * wx0 * wy1 &
         + grid(i1,j1) * wx1 * wy1
    else
      grid(i0,j0) = grid(i0,j0) + wx0 * wy0 * w
      grid(i0,j1) = grid(i0,j1) + wx0 * wy1 * w
      grid(i1,j0) = grid(i1,j0) + wx1 * wy0 * w
      grid(i1,j1) = grid(i1,j1) + wx1 * wy1 * w
    end if

  end if

end

!-----------------------------------------------------------------------

subroutine gyrospline(io, x, y, rho, w, grid)

  ! io = 0 for retrieve
  ! io = 1 for deposit
  ! x, y position
  ! rho is gyroradius
  ! w is weight or output field depending on io
  ! grid is field quantitiy to interpolate

  implicit none

  ! arguments
  integer :: io
  real(8) :: x,y,w,rho
  real(8),dimension(0:nx,0:ny) :: grid

  ! unit vectors for gyroaverage
  !real(8), parameter, dimension(4) :: unitvecx = (/1.0,-1.0,0.0,0.0/)
  !real(8), parameter, dimension(4) :: unitvecy = (/0.0,0.0,1.0,-1.0/)
  real(8), parameter, dimension(4) :: unitvecx = (/0.0,0.0,0.0,0.0/)
  real(8), parameter, dimension(4) :: unitvecy = (/0.0,0.0,0.0,0.0/)

  ! local weight
  real(8) :: myw

  ! loop control
  integer :: i

  do i=1,4
    if (io == 0) then
      w = 0
      !call spline(io, x + rho*unitvecx(i), y + rho*unitvecy(i), myw, grid)
      call spline(io, x, y, myw, grid)
      w = w + 0.25*myw
    else
      !call spline(io, x + rho*unitvecx(i), y + rho*unitvecy(i), 0.25*w, grid)
      call spline(io, x, y, 0.25*w, grid)
    end if
  end do

end subroutine gyrospline

!-----------------------------------------------------------------------
!--------particle push subroutines--------------------------------------
!-----------------------------------------------------------------------

subroutine epush

  implicit none
  integer :: m
  real(8) :: vdv,kap,edv
  real(8) :: ax,ay
  real(8) :: vxt,vyt !temp velocity storage

  ! ions
  if (fki == 1) then
    do m=1,ni
      call spline(0,xi(m),yi(m),ax,ex)
      call spline(0,xi(m),yi(m),ay,ey)
      ! 1/2 perp velocity push (note that vy is in rotated frame)
      vxi(m)=vxi(m) + .5*dt*ax*eperpi*memi
      vyi(m)=vyi(m) + .5*dt*ay*cth*eperpi*memi
      ! full velocity rotation (note that vy is in rotated frame)
      vxt = cdt*vxi(m) + sdt*vyi(m)
      vyt = -1.0*sdt*vxi(m) + cdt*vyi(m)
      vxi(m) = vxt
      vyi(m) = vyt
      ! 1/2 perp velocity push (note that vy is in rotated frame)
      vxi(m) = vxi(m) + .5*dt*ax*eperpi*memi
      vyi(m) = vyi(m) + .5*dt*ay*cth*eperpi*memi
      ! parallel velocity push
      vpari(m) = vpari(m) + dt*ay*sth*epari*memi
      ! weight equation terms
      !!!!!!!!!!!!!!!!!!!! CHECK THIS !!!!!!!!!!!!!!!!!!!!!
      vdv = vxi(m)**2+vyi(m)**2+vpari(m)**2
      edv = vxi(m)*ax + vyi(m)*ay*cth + vpari(m)*ay*sth
      kap = kapni+kapti*(.5*vdv/memi-1.5)
      ! explicit 1/2 weight advance
      wi0(m)=wi0(m)+.5*dt*wpi0(m)*(edv+cth*ay*kap)
      wpi0(m)=wpi0(m)-.5*dt*weighti*edv*wpi0(m)
      ! full position advance
      xi(m) = xi(m) + dt*vxi(m)
      yi(m) = yi(m) + dt*cth*vyi(m) + dt*sth*vpari(m)
      ! boundaries
      call enforce_bounds(xi(m),yi(m),vxi(m),vyi(m))
    end do
  end if

  ! electrons
  do m=1,ne
    if (gke == 1) then
      call gyrospline(0,xe0(m),ye0(m),sqrt(mue(m)),ax,ex)
      call gyrospline(0,xe0(m),ye0(m),sqrt(mue(m)),ay,ey)
    else
      call spline(0,xe0(m),ye0(m),ax,ex)
      call spline(0,xe0(m),ye0(m),ay,ey)
    end if
    ! full parallel velocity push
    vpare(m)=vpare(m) - dt*ay*sth*epare
    ! weight equation terms
    !!!!!!!!!!!!!!!!!!!! CHECK THIS !!!!!!!!!!!!!!!!!!!!!
    vdv=vpare(m)**2
    !kap=kapne+kapte*(.5*vdv-1.5)
    kap=kapne+kapte*(.5*vdv+mue(m)-1.5)
    ! explicit 1/2 weight advance
    we0(m)=we0(m)-.5*dt*wpe0(m)*(sth*ay*vpare(m)-cth*ay*kap)
    wpe0(m)=wpe0(m)+.5*dt*weighte*sth*ay*vpare(m)*wpe0(m)
    ! explicit part of position advance
    xe0(m) = xe0(m) + .5*dt*ay*cth*eperpe
    ye0(m) = ye0(m) - .5*dt*ax*cth*eperpe + dt*sth*vpare(m)
    ! boundaries
    call enforce_bounds(xe0(m),ye0(m))
    ! initial guess for implicit position
    xe1(m) = xe0(m) + .5*dt*ay*cth*eperpe
    ye1(m) = ye0(m) - .5*dt*ax*cth*eperpe
  end do

end

!-----------------------------------------------------------------------

subroutine ipush

  implicit none
  integer :: m
  real(8) :: vdv,kap,edv
  real(8) :: ax,ay

  ! ions
  if (fki == 1) then
    do m=1,ni
      call spline(0,xi(m),yi(m),ax,ex)
      call spline(0,xi(m),yi(m),ay,ey)
      ! weight equation terms
      vdv = vxi(m)**2 + vyi(m)**2 + vpari(m)**2
      edv = vxi(m)*ax + vyi(m)*ay*cth + vpari(m)*ay*sth
      kap = kapni + kapti*(.5*vdv/memi-1.5)
      ! implicit weight advance
      wpi1(m) = wpi0(m) - .5*dt*weighti*edv*wpi1(m)
      wi1(m) = wi0(m) + .5*dt*wpi1(m)*(edv+cth*ay*kap)
    end do
  end if

  ! electrons
  do m=1,ne
    if (gke == 1) then
      call gyrospline(0,xe1(m),ye1(m),sqrt(mue(m)),ax,ex)
      call gyrospline(0,xe1(m),ye1(m),sqrt(mue(m)),ay,ey)
    else
      call spline(0,xe1(m),ye1(m),ax,ex)
      call spline(0,xe1(m),ye1(m),ay,ey)
    end if
    ! weight equation terms
    vdv=vpare(m)**2
    !kap=kapne+kapte*(.5*vdv-1.5)
    kap=kapne+kapte*(.5*vdv+mue(m)-1.5)
    ! implicit part of weight advance
    wpe1(m)=wpe0(m)+.5*dt*weighte*sth*ay*vpare(m)*wpe1(m)
    we1(m)=we0(m)-.5*dt*wpe1(m)*(sth*ay*vpare(m)-cth*ay*kap)
    ! implicit part of position advance
    xe1(m)=xe0(m)+.5*dt*ay*cth*eperpe
    ye1(m)=ye0(m)-.5*dt*ax*cth*eperpe
    ! boundaries
    call enforce_bounds(xe1(m),ye1(m))
  end do

end

!-----------------------------------------------------------------------

subroutine enforce_bounds(x,y,vx,vy)

  implicit none
  real(8) :: x,y
  real(8), optional :: vx,vy

  !reflecting x-boundaries if reflect==1 otherwise periodic
  if (reflect == 1) then
    if (x > lx) then
      x = x - 2.*(x-lx)
      if(present(vx)) vx = -vx
    elseif (x < 0.) then
      x = -x
      if(present(vx)) vx = -vx
    elseif (x == lx) then
      x = 0.9999*lx
      if(present(vx)) vx = -vx
    endif
  else
    x=x-lx*dble(floor(x/lx))
  end if
  y=y-ly*dble(floor(y/ly))

end

!-----------------------------------------------------------------------

subroutine update

  implicit none

  if (fki == 1) then
    wi0 = wi1
    wpi0 = wpi1
  end if
  xe0 = xe1
  ye0 = ye1
  we0 = we1
  wpe0 = wpe1

end

!-----------------------------------------------------------------------
!--------output subroutines---------------------------------------------
!-----------------------------------------------------------------------

  ! OUTPUT SUBROUTINES

  subroutine output

    implicit none

    ! selected modes of potential and density
    call mode_output(phihist, 'phi_hist.out')
    call mode_output(denhist, 'den_hist.out')
    call mode_output(denhist, 'temp_hist.out')

    ! various diagnostic quantities
    if (mod(tstep, nrec) == 0) call diagnostics

    ! grid quantities
    if (mod(tstep, nrec) == 0) then
      call grid_output(phi, 'phi_grid.out')
      !call grid_output(phi_kmod, 'phi_kmod_grid.out')
      call grid_output(den, 'den_grid.out')
      if (fki == 1) then
        call grid_output(deni, 'denii_grid.out')
        call grid_output(dene, 'denee_grid.out')
      end if
      call grid_output(tempxy, 'temp_grid.out')
    end if

  end subroutine output


  subroutine mode_output(hist, filename)
    ! record real, imaginary components given by hist

    implicit none
    integer :: i, id = 99
    complex(8) :: hist(1: nmode)
    character(*) :: filename

    if (myid == 0) then ! master writes to files

      open(id, file = filename, form = 'formatted', status = 'unknown', position = 'append')

      if (tstep == 0) then
          write(id, '(a12, a4)', advance = 'no') 't', '    '
          do i = 1, nmode
            write(id, '(a4, i3.3, a1, i3.3, a1, a4)', advance = 'no') 'Re[', modeindices(1,i), ',', modeindices(2,i), ']', '    '
            write(id, '(a4, i3.3, a1, i3.3, a1, a4)', advance = 'no') 'Im[', modeindices(1,i), ',', modeindices(2,i), ']', '    '
          end do
          write(id, *)
      end if

      write(id, '(f12.4)', advance = 'no') dt * tstep
      do i = 1, nmode
        write(id, '(a4, e12.5)', advance = 'no') '    ', real(hist(i))
        write(id, '(a4, e12.5)', advance = 'no') '    ', imag(hist(i))
      end do
      write(id, *)

      endfile id
      close(id)

    endif

  end subroutine mode_output


  subroutine grid_output(grid_array, filename)
    ! record all values of grid quantity

    implicit none
    integer :: i, j, id = 99
    real(8) :: grid_array(0: nx, 0: ny)
    character(*) :: filename

    if (myid == 0) then ! master writes to files

      open(id, file = filename, form = 'formatted', status = 'unknown', position = 'append')

      do j = 0, ny - 1
        do i = 0, nx - 1
          write(id, '(e12.5)') grid_array(i, j)
        end do
      end do

      write(id, *)

      endfile id
      close(id)

    endif

  end subroutine grid_output


  subroutine diagnostics
    ! calculate and record various diagnostic quantities

    implicit none
    integer :: i, j, m, id = 99
    real(8) :: qxi, my_qxi ! heat flux in x - direction
    real(8) :: qxe, my_qxe ! heat flux in x - direction
    real(8) :: wi, my_wi, w2i, my_w2i ! mean weight, mean squared weight
    real(8) :: we, my_we, w2e, my_w2e ! mean weight, mean squared weight
    real(8) :: my_kei, kei, fe ! kinetic and field energies
    real(8) :: my_kee, kee ! kinetic and field energies
    real(8) :: te
    real(8) :: myupari(0:nx,0:ny), myupare(0:nx,0:ny) ! parallel flows
    real(8) :: upari(0:nx,0:ny), upare(0:nx,0:ny) ! parallel flows
    real(8) :: upari_fsa(0:nx), upare_fsa(0:nx)
    character(72) :: filename

    qxi = 0.0
    my_qxi = 0.0
    qxe = 0.0
    my_qxe = 0.0
    wi  = 0.0
    my_wi  = 0.0
    w2i = 0.0
    my_w2i = 0.0
    we  = 0.0
    my_we  = 0.0
    w2e = 0.0
    my_w2e = 0.0
    my_kei = 0.0
    kei = 0.0
    my_kee = 0.0
    kee = 0.0
    myupari = 0.0
    myupare = 0.0
    upari = 0.0
    upare = 0.0
    upari_fsa = 0.0
    upare_fsa = 0.0

    !!$acc parallel loop reduction(+:my_qxi,my_wi,my_w2i,my_kei)
    if (fki == 1) then
      do m = 1, ni
        ! x - direction heat flux
        my_qxi = my_qxi + 0.25 * (wi0(m) + wi1(m)) * vxi(m) &
            * (vxi(m) ** 2 + vyi(m) ** 2 + vpari(m) ** 2)
        ! weight sum
        my_wi = my_wi + wi1(m)
        ! weight squared sum
        my_w2i = my_w2i + wi1(m) ** 2
        ! kinetic energy
        my_kei = my_kei + 0.25 * (wi0(m) + wi1(m)) * (vxi(m) ** 2 + vyi(m) ** 2 + vpari(m) ** 2)
        ! parallel flow
        call spline(1,xi(m),yi(m),0.5*(wi0(m)+wi1(m))*vpari(m),myupari)
      end do
    end if

    do m=1,ne
      ! weight sum
      my_we = my_we + we1(m)
      !weight squared sum
      my_w2e = my_w2e + we1(m)**2
      !kinetic energy
      my_kee = my_kee + 0.25 * (we0(m) + we1(m)) *memi*vpare(m)**2
      ! parallel flow
      call spline(1,0.5*(xe0(m)+xe1(m)),0.5*(ye0(m)+ye1(m)),0.5*(we0(m)+we1(m))*vpare(m),myupare)
    end do

    if (fki == 1) then
      call mpi_allreduce(my_qxi, qxi, 1, mpi_real8, mpi_sum, mpi_comm_world, ierr)
      call mpi_allreduce(my_wi,  wi,  1, mpi_real8, mpi_sum, mpi_comm_world, ierr)
      call mpi_allreduce(my_w2i, w2i, 1, mpi_real8, mpi_sum, mpi_comm_world, ierr)
      call mpi_allreduce(my_kei, kei, 1, mpi_real8, mpi_sum, mpi_comm_world, ierr)
      call mpi_allreduce(myupari, upari, nx*ny, mpi_real8, mpi_sum, mpi_comm_world, ierr)
    end if
    call mpi_allreduce(my_we,we,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(my_w2e,w2e,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(my_kee,kee,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(myupare, upare, nx*ny, mpi_real8, mpi_sum, mpi_comm_world, ierr)

    if (fki == 1) then
      qxi = qxi / dble(tni)
      wi  = wi  / dble(tni)
      w2i = w2i / dble(tni)
      kei = kei / dble(tni)
      upari = upari * nx * ny / dble(tni)
    end if

    we=we/dble(tne)
    w2e=w2e/dble(tne)
    kee=kee/dble(tne)
    upare = upare * nx * ny / dble(tne)

    ! enforce periodic boundaries on parallel flow
    if (reflect /= 1) then
      do j=0,ny
        upari(0,j)=upari(0,j)+upari(nx,j)
        upari(nx,j)=upari(0,j)
        upare(0,j)=upare(0,j)+upare(nx,j)
        upare(nx,j)=upare(0,j)
      end do
    end if

    ! take fsa of parallel flow
    upari_fsa = sum(upari, 2)/ny
    upare_fsa = sum(upare, 2)/ny

    ! field energy
    fe = 0.
    do j = 0, ny - 1
      do i = 0, nx - 1
        fe = fe + 0.5 * (ex(i,j) ** 2 + ey(i,j) ** 2)
      end do
    end do

    fe = fe/nx/ny ! factor of 8 smaller than Scott's old GK code

    ! total energy
    te = kei + kee + fe

    if (myid == 0) then ! master writes to files

      ! heat flux

      filename = 'qflux.out'
      open(id, file = filename, form = 'formatted', status = 'unknown', position = 'append')

      if (tstep == 0) then
          write(id, '(a12, a4)', advance = 'no') 't', '    '
          write(id, '(a12, a4)', advance = 'no') 'ion_qflux_x', '    '
          write(id, *)
      end if

      write(id,   '(f12.4)', advance = 'no') dt * tstep
      write(id, '(a4, e12.5)', advance = 'no') '    ', qxi
      write(id, *)

      endfile id
      close(id)

      ! weight

      filename = 'weight.out'
      open(id, file = filename, form = 'formatted', status = 'unknown', position = 'append')

      if (tstep == 0) then
          write(id, '(a12, a4)', advance = 'no') 't', '    '
          write(id, '(a12, a4)', advance = 'no') '<w_i>', '    '
          write(id, '(a12, a4)', advance = 'no') '<w_i^2>', '    '
          write(id, '(a12, a4)', advance = 'no') '<w_e>', '    '
          write(id, '(a12, a4)', advance = 'no') '<w_e^2>', '    '
          write(id, *)
      end if

      write(id,   '(f12.4)', advance = 'no') dt * tstep
      write(id, '(a4, e12.5)', advance = 'no') '    ', wi
      write(id, '(a4, e12.5)', advance = 'no') '    ', w2i
      write(id, '(a4, e12.5)', advance = 'no') '    ', we
      write(id, '(a4, e12.5)', advance = 'no') '    ', w2e
      write(id, *)

      endfile id
      close(id)

      ! energy

      filename = 'energy.out'
      open(id, file = filename, form = 'formatted', status = 'unknown', position = 'append')

      if (tstep == 0) then
          write(id, '(a12, a4)', advance = 'no') 't', '    '
          write(id, '(a12, a4)', advance = 'no') 'ikinetic', '    '
          write(id, '(a12, a4)', advance = 'no') 'ekinetic', '    '
          write(id, '(a12, a4)', advance = 'no') 'field', '    '
          write(id, '(a12, a4)', advance = 'no') 'total', '    '
          write(id, *)
      end if

      write(id,    '(f12.4)', advance = 'no') dt * tstep
      write(id, '(a4, e12.5)', advance = 'no') '    ', kei
      write(id, '(a4, e12.5)', advance = 'no') '    ', kee
      write(id, '(a4, e12.5)', advance = 'no') '    ', fe
      write(id, '(a4, e12.5)', advance = 'no') '    ', te
      write(id, *)

      endfile id
      close(id)

      ! parallel flow

      filename = 'upari_fsa.out'

      open(id, file = filename, form = 'formatted', status = 'unknown', position = 'append')

      write(id,   '(f12.4)', advance = 'no') dt * tstep
      do i = 0, nx
        write(id, '(a4, e12.5)', advance = 'no') '    ', upari_fsa(i)
      end do
      write(id,*)

      endfile id
      close(id)

      filename = 'upare_fsa.out'

      open(id, file = filename, form = 'formatted', status = 'unknown', position = 'append')

      write(id,   '(f12.4)', advance = 'no') dt * tstep
      do i = 0, nx
        write(id, '(a4, e12.5)', advance = 'no') '    ', upare_fsa(i)
      end do
      write(id,*)

      endfile id
      close(id)

    endif

  end subroutine diagnostics


  ! PRINT SUBROUTINES FOR DEBUGGING

  subroutine print_particle_array(particle_array, name)

    implicit none
    integer :: m
    real(8), intent(in) :: particle_array(1: ni)
    character(*) :: name

    do m = 1, ni
      print *, name, '(', m, ') = ', particle_array(m)
    end do

    print *

  end subroutine print_particle_array


  subroutine print_grid_array(grid_array, name)

    implicit none
    integer :: i, j
    real(8), intent(in) :: grid_array(0: nx - 1, 0: ny - 1)
    character(*) :: name

    do j = 0, ny - 1
      do i = 0, nx - 1
        print *, name, '(', i, ',', j, ') = ', grid_array(i, j)
      end do
    end do

    print *

  end subroutine print_grid_array


  subroutine print_grid_array_dft(grid_array_dft, name)

    implicit none
    integer :: i, j
    complex(8), intent(in) :: grid_array_dft(0: nx - 1, 0: ny - 1)
    character(*) :: name

    do j = 0, ny - 1
      do i = 0, nx - 1
        print *, name, '(', i, ',', j, ') = ', grid_array_dft(i, j)
      end do
    end do

    print *

  end subroutine print_grid_array_dft


!subroutine modeout(hist,fl,id)
!
!  !record components of phi given by modehist
!
!  implicit none
!  integer :: i,id
!  character*5 :: fl
!  character*70 :: flnm
!  complex(8) :: hist(1:nmode)
!
!  if (myid==0) then
!    flnm=fl//'.out'
!    open(id,file=flnm,form='formatted',status='unknown',&
!      position='append')
!    write(id,'(f8.2)',advance="no") dt*tstep
!    do i=1,nmode
!      write(id,'(a2,e13.6,a2,e13.6)',advance="no") '  ',real(hist(i)),&
!        '  ',imag(hist(i))
!    end do
!    write(id,*)
!    endfile id
!    close(id)
!  endif
!
!end
!
!!-----------------------------------------------------------------------
!
!subroutine gridout(u,fl,id)
!
!  !record all values of grid quantity
!
!  implicit none
!  integer :: i,j,id
!  real(8) :: u(0:nx,0:ny)
!  character*5 :: fl
!  character*70 :: flnm
!
!  if (myid==0) then
!    flnm=fl//'.out'
!    open(id,file=flnm,form='formatted',status='unknown',&
!      position='append')
!    do j=0,ny
!      do i=0,nx
!        write(id,'(e11.4)') u(i,j)
!      end do
!    end do
!    endfile id
!    close(id)
!  endif
!  
!end
!
!!-----------------------------------------------------------------------
!
!subroutine gendiagnostics
!
!  implicit none
!  integer :: id,m,i,j
!  real(8) :: qx,myqx,w2i,myw2i,w2e,myw2e,mykei,kei,mykee,kee,te
!  character*5 :: fl
!  character*70 :: flnm
!
!  id=89
!
!  qx=0
!  w2i=0
!  w2e=0
!  kei=0
!  kee=0
!
!  myqx=0
!  myw2i=0
!  myw2e=0
!  mykei=0
!  mykee=0
!
!  do m=1,ni
!    !net ion heat flux in x-direction
!    myqx = myqx + wi1(m)*vxi(m)*(vxi(m)**2+vyi(m)**2+vpari(m)**2)
!    !weight squared sum
!    myw2i = myw2i + wi1(m)**2
!    !kinetic energy
!    mykei = mykei + 0.5*wi1(m)*(vxi(m)**2+vyi(m)**2+vpari(m)**2)
!  end do
!
!  if (dke == 1) then
!    do m=1,ne
!      myw2e = myw2e + we1(m)**2 !weight squared sum
!      mykee = mykee + 0.5*we1(m)*memi*vpare(m)**2 !kinetic energy
!    end do
!  end if
!
!
!  call mpi_allreduce(myqx,qx,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
!  call mpi_allreduce(myw2i,w2i,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
!  if (dke == 1) call mpi_allreduce(myw2e,w2e,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
!  call mpi_allreduce(mykei,kei,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
!  if (dke == 1) call mpi_allreduce(mykee,kee,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
!
!  qx=qx/dble(tni)
!  w2i=w2i/dble(tni)
!  kei=kei/dble(tni)
!
!  if (dke == 1) then
!    w2e=w2e/dble(tne)
!    kee=kee/dble(tne)
!  else
!    w2e = 0.
!    kee = 0.
!  end if
!
!  ! total energy
!  te = kei + fe
!  if (dke == 1) te = te + kee
!
!  if (myid==0) then
!    flnm='diagn.out'
!    open(id,file=flnm,form='formatted',status='unknown', &
!      position='append')
!    if (tstep==1) write(id,'(a31)') 't  qx  w2i  w2e  kei  kee fe te'
!    write(id,'(f8.2)',advance="no") dt*tstep
!    write(id,'(a2,e13.6,a2,e13.6,a2,e13.6,a2,e13.6,a2,e13.6,a2,e13.6,a2,e13.6)') &
!      '  ',qx,'  ',w2i,'  ',w2e,'  ',kei,'  ',kee,'  ',fe,'  ',te
!    endfile id
!    close(id)
!  endif
!
!end

!-----------------------------------------------------------------------

end
