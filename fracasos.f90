subroutine fracasos(x,frac,xsolv)
use system
use const
 implicit none
integer i
real*16 frac(8),x(2)
real*16 xphiA,xphiB,xmphiA,xmphiB,xsolv
real*16 eta,M ,etap,Kap,Kbp,Kop,xOHmin,xHplus
real*16 KNap,KClp
xphiA=x(1);
xphiB=x(2);
xmphiA=xphiA/(Ma*vp);
xmphiB=xphiB/(Mb*vp);

xHplus=xsolv*expmuHplus
xOHmin=xsolv*expmuOHmin
xNaplus=(xsolv**vsal)*expmupos
xClmin=(xsolv**vsal)*expmuneg


!print* ,x
eta=xphiB/xphiA
etap=xphiA/xphiB
Kap=KaA*xsolv/xHplus
Kbp=KaB*xsolv/xOHmin
Kop=Keo*xmphiA*Ma*vs

KNap=xNaplus*KaA*xsolv/(vsal*xHplus*KaNa*(xsolv**vsal))
KClp=xClmin*KaB*xsolv/(vsal*xOHmin*KaCl*(xsolv**vsal))

Knew=Kop*Kap*Kbp/((1.0+Kap+KNap)*(1.0+Kbp+KClp))
frac(2)=0.5*(1.0+etap+etap/Knew )-( (0.5*(1.0+etap+etap/Knew))**2-etap)**0.5;!fB_a
frac(1)=frac(2)*eta;!fA_a
!frac(3)=(1.-frac(1))/(1.+Kap) ;!fAnc
!frac(4)=(1.-frac(2))/(1.+Kbp) ;!fBnc
frac(5)=(1.-frac(1))*Kap/(1.+Kap+KNap);!fAc
frac(6)=(1.-frac(2))*Kbp/(1.+Kbp+KClp);!fBc
frac(7)=KNap*(1.-frac(1)-frac(5))/(1.+KNap) !fANa
frac(8)=KClp*(1.-frac(2)-frac(6))/(1.+KClp) !fBCl
frac(3)=(1.-frac(1)-frac(5)-frac(7)) !fAnc
frac(4)=(1.-frac(2)-frac(6)-frac(8)) !fBnc

!print* ,xsolv
!print* , 'f',frac
!stop
!do i = 1,6
!  if (frac(i) .lt. 10**(-15))then
!    frac(i)=10**(-15)
! endif
!enddo
                                  
end subroutine
