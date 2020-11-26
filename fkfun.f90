subroutine fkfun(x,f,ier)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! User provided routine for kinsol
! x is the input vector
! f is the output vector, kinsol will change x in order to get f = 0
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use newton, only: solve
use functions, only: falpha, fbeta, fp
use solver
use system
use const
use results
implicit none
 
integer ntot
real*8 x(3),f(3)
real*16 protemp, protemp1
integer i,j, k, ix, iy, iz, ii, ax, ay, az, temp, iiZ
!real*16 xpotA(dimz),fdisbc,fdisAC
!real*16 xpotB(dimz)
!real*16 psi2(0:dimz+1) ! psi plus boundaries at z=0 and dimz+1
!real*8 xtotal(-Xulimit:dimz+Xulimit) ! xtotal for poor solvent
real*16 m, eta, Penality,Penality2

! Kinsol
integer*4 ier
real*16 vectalpha(2),vectbeta(2)
real*16 mu2alpha,mu2beta,mu3alpha,mu3beta,fealpha,febeta
real*16 potquim2,elib,potquim3
real*16 xsolventalpha,xsolventbeta,xHplusalpha,xHplusbeta
real* 16 xOHminalpha,xOHminbeta
real*8 xsol, xsol0
logical debug
integer iters
! Recovers xh and psi from x
real*16 fracalpha(8), fracbeta(8)

debug = .false.

x2alpha=10**(-x(1))
x3alpha=x2alpha
x2beta=10**(-x(2))
x3beta=x2beta

vectalpha(1)=x2alpha
vectalpha(2)=x3alpha
vectbeta(1)=x2beta
vectbeta(2)=x3beta

iters = 0
xsol0=(1.0-x2alpha-x3alpha)/(1+expmuOHmin+expmuHplus+expmupos+expmuneg) ! initial guess vs = vsalt
call solve(falpha, fp, xsol0, xsol, iters, debug)
xsolventalpha = xsol

iters = 0
xsol0=(1.0-x2beta-x3beta)/(1+expmuOHmin+expmuHplus+expmupos+expmuneg) ! initial guess vs = vsalt
call solve(fbeta, fp, xsol0, xsol, iters, debug)
xsolventbeta = xsol


xHplusalpha=xsolventalpha*cHplusbulk
xOHminalpha=xsolventalpha*cOHminbulk
xHplusbeta=xsolventbeta*cHplusbulk
xOHminbeta=xsolventbeta*cOHminbulk


call fracasos(vectalpha,fracalpha,xsolventalpha)
call fracasos(vectbeta,fracbeta,xsolventbeta)
 
!print* , fracalpha
!print* , fracbeta
!stop 
! Pot quimico respecto de phiA 
Penality=(x2alpha-x2beta)**2/(x2alpha+x2beta)**2
Penality2=((x2alpha-x3alpha)**2+(x2beta-x3beta)**2)

call mu2(vectalpha,potquim2,xsolventalpha)
mu2alpha = potquim2

call mu2(vectbeta,potquim2,xsolventbeta)
mu2beta=potquim2

call mu3(vectalpha,potquim3,xsolventalpha)
mu3alpha = potquim3

call mu3(vectbeta,potquim3,xsolventbeta)
mu3beta=potquim3
 
call fe(vectalpha,elib,xsolventalpha)
fealpha=elib

call fe(vectbeta,elib,xsolventbeta)
febeta=elib

!print*, mu2beta,mu2alpha,mu3alpha,mu3beta
!stop
! ### EQUATIONS TO SOLVE

 f(1)= mu2alpha-mu2beta
 f(1)=f(1)/Penality

! Recta tangente

 f(2)= (fealpha-febeta)/Penality !&

iter = iter + 1
norma = 0.0

do i = 1, 2
norma = norma +(f(i))**2    
enddo

!print*,iter, x(1), x(2), f(1), f(2), Penality
ier = 0.0
return
end subroutine
