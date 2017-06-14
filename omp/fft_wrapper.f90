module fft_wrapper
implicit none
real(8) :: coefxn(20000),coefyn(20000),coefzn(20000)
real(8) :: scoef(50000)

contains

subroutine ccfft(c,isign,n,x)
   character :: c
   integer :: isign,n      
   complex(8),dimension(0:) :: x

!initialize
   if(isign==0)then
      if(c=='x')then
 call zffti(n,coefxn)
 call zffti(n,coefxn)
      end if
      if(c=='y')then
 call zffti(n,coefyn)
 call zffti(n,coefyn)
      end if
      if(c=='z')then
 call zffti(n,coefzn)
 call zffti(n,coefzn)
      end if
   end if

!backward transform?
   if(isign==1)then
      if(c=='x')then
 call zfftb(n,x,coefxn)
      end if
      if(c=='y')then
 call zfftb(n,x,coefyn)
      end if
      if(c=='z')then
 call zfftb(n,x,coefzn)
      end if
   end if

!forward transform?
   if(isign==-1)then
      if(c=='x')then
 call zfftf(n,x,coefxn)
      end if
      if(c=='y')then
 call zfftf(n,x,coefyn)
      end if
      if(c=='z')then
 call zfftf(n,x,coefzn)
      end if
   end if

   return
end subroutine ccfft

subroutine dsinf(init,x,n)
   integer :: init,n,i
   real(8),dimension(:) :: x

   if(init/=0)then
 call dsinti(n,scoef)
   end if

   if(init==0)then
 call dsint(n,x,scoef)
 do i=1,n
    x(i)=0.5d0*x(i) !dsint outputs 2 * sine transform
 end do
   end if

   return
end subroutine dsinf

end module fft_wrapper

