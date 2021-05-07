! ``Volume penalization for inhomogeneous Neumann boundary conditions modeling scalar flux in complicated geometry'' [J. Comput. Phys. 390 (2019) 452-469]
! and its Corrigendum.
! by T. Sakurai, K. Yoshimatsu, N. Okamoto and K. Schneider
! This is the fortran 90 code for Fig. 12 in 
!  "Corrigendum to `Volume penalization for inhomogeneous Neumann boundary conditions modeling scalar flux in complicated geometry'
!cc remake by KY.
!cc 2D periodic domain
!c L2norm
program main
implicit none
integer,parameter :: n=2**10, itmax=10**6
integer :: i, j, itr
double precision, dimension(0:n+1,0:n+1) :: vp, tmp, chi, theta
double precision, dimension(n,n) :: f, xbeta, ybeta, w
double precision, dimension(n,n) :: err_v
!, err_test
double precision :: rhs, tvp, sum_uni
double precision :: dchidx, dchidy, rhs0,left1x,left1y
double precision :: th_ip, th_im, th_jp, th_jm
double precision :: facv, sor_inc
double precision :: error, var, ern, exact_s, errL1, errL2, atmp
double precision :: x, y, dx, dy, r, rb, rb3, b, c, pi, xrange, yrange, gc, dx2, dy2
double precision, parameter :: alpha = 1.0d0, eta = 1.0d-8, eps = 1.0d-12
! for n=32 -> omg=1.2
!double precision, parameter :: omg=1.2d0, conv = 1.0d-9
double precision, parameter :: omg=1.9d0, conv = 1.0d-9
!double precision, parameter :: omg=1.9d0, conv = 1.0d-10
!*************************
 pi = 4.0d0*datan(1.0d0)
!write(*,*) pi
! 3.1415926535897931
 rb = pi/4.0d0
 rb3 = 3.0d0*rb
 xrange = 2.0d0*pi
 yrange = 2.0d0*pi
 dx = xrange/dble(n)
 dy = yrange/dble(n)
!*************************
! zero set
  xbeta = 0.0d0
  ybeta = 0.0d0
  w = 0.0d0
  f = 0.0d0
! for exact solution (non-pena) 
  b = 3.0d0*pi/4.0d0*alpha
  c =-3.0d0/32.0d0*pi*alpha* &
 &   (9.0d0*dlog(3.0d0*pi/4.0d0)-dlog(pi/4.0d0) - 4.0d0)
!== set conditions ==
 do j=1,n
  do i=1,n
   x=dx*dble(i-1)
   y=dy*dble(j-1)
   r=dsqrt((x-pi)**2+(y-pi)**2) 
!c
   if ((r>= 0.0d0 - eps).and.(r <= pi + eps)) then
    gc = alpha*((4.0d0*r/3.0d0/pi)**2)    &
 &             *(4.0d0*(1.0d0-r/pi))**3
    if (r>0.0d0) then
     xbeta(i,j) = gc*(x-pi)/r
     ybeta(i,j) = gc*(y-pi)/r
    end if
   endif
!   if (r > 0.0d0) then
   if ((r>rb -eps).and.(r < rb3 + eps)) then
!c  the exact solution
    w(i,j) = dcos(4.0d0*r) + b*dlog(r) + c
!c  source term
    f(i,j) = 16.0d0*dcos(4.d0*r)+4.0d0*dsin(4.d0*r)/r
!c  later, we multiply 1-chi -> (1-chi)*f
   end if
  end do
 end do
!=======  mask function =====================
 chi = 1.0d0
 do j=1,n
  do i=1,n
   x=dx*dble(i-1)
   y=dy*dble(j-1)
   r=dsqrt((x-pi)**2+(y-pi)**2) 
   if( (r > rb) .and. (r < rb3)) chi(i,j)=0.0d0
   if( (dabs(r-rb)<eps) .or. (dabs(r-rb3)<eps)) chi(i,j)=0.5d0
  end do
 end do
!c===========================================
! periodic bc (chi, theta)
 do i=1,n
!   theta(i,0)   = theta(i,n)
!   theta(i,n+1) = theta(i,1)
   chi(i,0)   = chi(i,n)
   chi(i,n+1) = chi(i,1)
 end do
 do j=1,n
!   theta(0,  j) = theta(n,j)
!   theta(n+1,j) = theta(1,j)
   chi(0,  j) = chi(n,j)
   chi(n+1,j) = chi(1,j)
 end do
!c  chi=1.0d0 -> theta = eta
! theta=eta
! do j=1,n
!  do i=1,n
 do j=0,n+1
  do i=0,n+1
   theta(i,j) = 1.0d0-chi(i,j) + eta*chi(i,j)
  end do
 end do
!c*************** SOR **************************************
!c* int_(no solid region) 1 dxdy, int_(no solid region) w dxdy 
!c   w(x,0)=w(x,2pi), w(0,y)=w(2pi,y) : 0.5*w(x,0)+0.5*w(x,2pi)=w(x,0)
  sum_uni = 0.0d0
  exact_s = 0.0d0
  do j=1,n
   do i=1,n
    sum_uni = sum_uni + (1.0d0-chi(i,j))*dx*dy
    exact_s = exact_s + (1.0d0-chi(i,j))*dx*dy*w(i,j)
   end do
  end do
!c**********************************************************
! vp=0.0d0
  do j=1,n
   do i=1,n
    vp(i,j) = w(i,j)
   end do
  end do
 tmp = 0.0d0
 dx2 = dx**2
 dy2 = dy**2
do itr=1,itmax
!c++  periodic bc. for vp: after update of tmp ++
   do j=0,n+1
    do i=0,n+1
     tmp(i,j)=vp(i,j)
    end do
   end do
!c=========================
  do j=1,n
   do i=1,n
!c nabla (theta nabla v) = -(1-chi) f - beta*(nabla chi)
!c  beta=(xbeta, ybeta)
    dchidx =(-chi(i-1,j) + chi(i+1,j))/(2.0d0*dx)
    dchidy =(-chi(i,j-1) + chi(i,j+1))/(2.0d0*dy)
!c
    rhs0 = - (1.0d0-chi(i,j))*f(i,j)                                &
 &         - xbeta(i,j)*dchidx - ybeta(i,j)*dchidy  
!c ip = i+1/2, im = i-1/2, jp = j+1/2, jm = j-1/2
    th_ip = 0.5d0*(theta(i,j)   + theta(i+1,j))
    th_im = 0.5d0*(theta(i-1,j) + theta(i,j)  )
    th_jp = 0.5d0*(theta(i,j)   + theta(i,j+1))
    th_jm = 0.5d0*(theta(i,j-1) + theta(i,j)  )
!c dx (theta dxv) no v_ij part
    left1x = th_ip*tmp(i+1,j)/dx2 + th_im*tmp(i-1,j)/dx2
!c dy (theta dyv) no v_ij part
    left1y = th_jp*tmp(i,j+1)/dy2 + th_jm*tmp(i,j-1)/dy2
!c  rhs    
     rhs = rhs0 - ( left1x + left1y )
!c factor of v_ij part 
     facv = - th_ip/dx2 - th_im/dx2 - th_jp/dy2 - th_jm/dy2
!c
     sor_inc = rhs/facv
!c--------------------
    tmp(i,j)=tmp(i,j)+omg*(sor_inc - tmp(i,j))
   end do
  end do
!c=========================
!c int_(no solid domain) tmp dx dy -> 0
  tvp=0.d0
  do j=1,n
   do i=1,n
    tvp = tvp + (1.0d0-chi(i,j))*tmp(i,j)*dx*dy
   end do
  end do
!c  write(*,*) tvp
!c
  do j=1,n
   do i=1,n
    tmp(i,j) = tmp(i,j) - tvp/sum_uni
   end do
  end do
!c
! tvp=0.d0
!  do j=1,n
!   do i=1,n
!    tvp = tvp + (1.0d0-chi(i,j))*tmp(i,j)*dx*dy
!   end do
!  end do
!  write(*,*) tvp 
!c tvp : from  e-14 to e-16    confirmed
!c----------------------
  var = 0.d0
  ern = 0.d0
  do j=1,n
   do i=1,n
    var = var + tmp(i,j) **2
    ern = ern + (vp(i,j) - tmp(i,j))**2
   end do
  end do
  error = dsqrt(ern)/dsqrt(var)
  if (mod(itr,10).eq.0) write(*,*) itr, real(error)
!c----------------------------
 if (error < conv) then
   exit
 endif
  if(itr==itmax) write(*,*) 'iteration max'
!c----------- update of v -----------------
  do j=1,n
   do i=1,n
    vp(i,j)=tmp(i,j)
   end do
  end do
!c---- periodic bc for vp ----
   do i=1,n
     vp(i,0)   = vp(i,n)
     vp(i,n+1) = vp(i,1)
   end do
   do j=1,n
     vp(0,j)   = vp(n,j)
     vp(n+1,j) = vp(1,j)
   end do
!c
end do
!+++++++++++++++++++++++++++++++++++++++++++
!    open(101,file='res.dat')
! do j=1,n
!  do i=1,n
!   x=dx*dble(i-1)
!   y=dy*dble(j-1)
!   if (dabs(chi(i,j)-1.0d0) > eps) then
!    write(101,*) x, y, vp(i,j), w(i,j)
!   else
!    write(101,*) x, y, vp(i,j), 0.0d0
!   endif
!  end do
!    write(101,*)
! end do
! close(101)
!c-----------------------
!c error in the fluid region
  err_v = 0.0d0
!  err_test = 0.0d0
  do j=1,n
   do i=1,n
    if (dabs(chi(i,j)) < eps) then
          err_v(i,j)=dabs( w(i,j) - vp(i,j) )
!       err_test(i,j)=dabs( w(i,j) - vp(i,j) - exact_s/sum_uni )
    end if
   end do
  end do
!+++++++++++++++++++++++++++++++++++++++++++
!    sum_uni = sum_uni - pi*(rb3**2 - rb**2)
!+++++++++++++++++++++++++++++++++++++++++++
  errL2=0.0d0
  errL1=0.0d0
!  atmp=0.0d0
  do j=1,n
   do i=1,n
!    if (err_v(i,j) > atmp ) atmp = err_v(i,j)
!    if (dabs(chi(i,j)) < eps) then
       errL2 = errL2 +  err_v(i,j)**2*dx*dy
       errL1 = errL1 +  err_v(i,j)*dx*dy
!       err_test(i,j)=dabs( w(i,j) - vp(i,j) - exact_s/sum_uni )
!    end if
   end do
  end do
!c
 open(111,file='./err.dat')
  write(111,*) n, dx, maxval(err_v), dsqrt(errL2), errL1
!c, atmp .. ok
 close(111)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 open(102,file='./res_x_pena.dat')
 open(202,file='./res_x_theo.dat')
!c
 i = n / 2
 do j=1,n
  x=dx*dble(i-1)
  y=dy*dble(j-1)
  if (dabs(chi(i,j)-1.0d0) > eps) then
   write(102,*) x, y, vp(i,j), w(i,j), dabs(w(i,j) -vp(i,j))
   write(202,*) x, y, vp(i,j), w(i,j), dabs(w(i,j) -vp(i,j))
  else
   write(102,*) x, y, vp(i,j), 0.0d0, dabs(w(i,j) -vp(i,j))
  endif
 end do
 close(102)
 close(202)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 open(103,file='./res_y_pena.dat')
 open(203,file='./res_y_theo.dat')
 j = n / 2
 do i=1,n
  x=dx*dble(i-1)
  y=dy*dble(j-1)
  if (dabs(chi(i,j)-1.0d0) > eps) then
   write(103,*) x, y, vp(i,j), w(i,j), dabs(w(i,j) -vp(i,j))
   write(203,*) x, y, vp(i,j), w(i,j), dabs(w(i,j) -vp(i,j))
  else
   write(103,*) x, y, vp(i,j), 0.0d0, dabs(w(i,j) -vp(i,j))
  endif
 end do
 close(103)
 close(203)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 stop
 end
!c///////////////////////////////////////////////
