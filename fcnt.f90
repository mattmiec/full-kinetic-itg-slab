!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!   Functions, etc.
!   (things that don't change)

module fcnt

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      real(8) function revers(num,n)

      implicit none
      integer :: num,n
      integer :: inum,iquot,irem
      real(8) :: rev,power

!     function to reverse the digits of num in base n.

      rev = 0.
      inum = num
      power = 1.

   11 continue
      iquot = int(inum/n)
      irem = inum - n*iquot
      power = power/n
      rev = rev + irem*power
      inum = iquot
      if(inum.gt.0) goto 11

      revers = rev
      return
      end


!--------------------------------------------------

      real(8) function ran2(idum)
      parameter( IM1=2147483563,  &
                IM2=2147483399, &
                AM=1.0/IM1,&
                IMM1=IM1-1,&
                IA1=40014,&
                IA2=40692,&
                IQ1=53668,&
                IQ2=52774,&
                IR1=12211,&
                IR2=3791,&
                NTAB=32,&
                NDIV=1+IMM1/NTAB,&
                EPS=1.2e-7,&
                RNMX=1.0-EPS &
               )
      integer :: j,k,idum2=123456789,iy=0,iv(0:NTAB-1)
      real(8) :: temp

      save idum2, iy,iv
!      write(*,*)'idum2,iy  ',idum2,iy
      if(idum.le.0)then
         if(-idum.lt.1)then
            idum=1
         else
            idum = -idum
         end if
         idum2 = idum
         do j = NTAB+7,0,-1
            k = idum/IQ1
            idum = IA1*(idum-k*IQ1)-k*IR1
            if(idum.lt.0)idum = idum+IM1
            if(j.lt.NTAB)iv(j) = idum
         end do
         iy = iv(0)
      end if

      k = idum/IQ1
      idum = IA1*(idum-k*IQ1)-k*IR1
      if(idum.lt.0)idum = idum+IM1
      k = idum2/IQ2
      idum2 = IA2*(idum2-k*IQ2)-k*IR2
      if(idum2.lt.0)idum2 = idum2+IM2
      j = iy/NDIV
      iy = iv(j)-idum2
      iv(j) = idum
      if(iy<1)iy = iy+IMM1
      temp = AM*iy
      if(temp>RNMX)then
         ran2 = RNMX
      else
         ran2 = temp
      end if
      return
      end

!-------------------------------------------------

      real(8) function dinvnorm(p)
      implicit none
      real(8) p,p_low,p_high
      real(8) a1,a2,a3,a4,a5,a6
      real(8) b1,b2,b3,b4,b5
      real(8) c1,c2,c3,c4,c5,c6
      real(8) d1,d2,d3,d4
      real(8) z,q,r
      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1-p_low
      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
201   q=dsqrt(-2*dlog(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/&
      &((((d1*q+d2)*q+d3)*q+d4)*q+1)
      goto 204
301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
202   q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/&
      &(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      goto 204
302   if((p.gt.p_high).and.(p.lt.1)) goto 203
203   q=dsqrt(-2*dlog(1-p))
      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/&
      &((((d1*q+d2)*q+d3)*q+d4)*q+1)
204   dinvnorm=z
      return
      end

!-------------------------------------------------
      real(8) function invteq(e,u)
      
      implicit none
      real(8) :: e,u
      real(8) :: x
      real(8) :: f,fp
      integer :: i

      ! function for solving equation:
      ! x + e*sin(x) = u
      ! via 6 newton iterations, assuming e < 1
      ! 0 < u < 1

      x=u
      do i=1,6
        f=x+e*dsin(x)-u
        fp=1.0+e*dcos(x)
        x=x-f/fp
      end do

      invteq=x

      return
      end

!-------------------------------------------------
end module fcnt
