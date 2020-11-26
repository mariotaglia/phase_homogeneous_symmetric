  !###############################################################################!     
  !     Simple brush: Standard Molecular Theory Program 
  !    
  !     Calculates a weak polyelectrolyte brush in poor sv conditions 
  !     Calculates free-energy
  !     Solve BINODAL METHOD
  !     MARCH 2020
  !         
  !###############################################################################
!use pks
use system
use const

  implicit none
  integer i, flagcrash
  real*16 pKw,K0A,K0EO,K0B,K0ANa,K0BCl,kte
  real*16 xposbulk,xnegbulk,xHplusbulk,xOHminbulk,Kw,xsalt

  !print*, 'Program Simple Brush'
  !print*, 'GIT Version: ', _VERSION
yes=0 ! es para  chequear si encuentra o no xalpha, xbeta
flagcrash=1
  call readinput  ! reads input variables from file
!  call allocation ! allocates memory

  K0eo=10**(-pKeo)
  K0A=10**(-pKaA)
  K0B=10**(-pKaB)
  K0ANa=10**(-pKaNa)
  K0BCl=10**(-pKaCl)

  pKw=14.0
  kW=10**(-pKw)
  pOHbulk=pKw-pHbulk
  cOHminbulk=10**(-pOHbulk)
  cHplusbulk=10**(-pHbulk)
  !cNaplus=10**(-pNabulk)
  !cClmin=10**(-pClbulk)
  vp=vpol!
  
  xposbulk=phi_sal  !NUEVO
  xnegbulk=phi_sal  !NUEVO
  
  xHplusbulk = (cHplusbulk*Na/(1.0d24))*(vs)
  xOHminbulk = (cOHminbulk*Na/(1.0d24))*(vs) 

        
  do i = 1, ncsal ! loop in csal

  csal = csalini + (csalfin-csalini)/float(ncsal-1)*float(i-1) 


  xsalt=csal*vsal*vs*6.02/10.0

  if(pHbulk.le.7) then  ! pH<= 7
     xposbulk=xsalt
     xnegbulk= xsalt +(xHplusbulk -xOHminbulk)*vsal! NaCl+ HCl  
  else                  ! pH >7 
     xposbulk=xsalt +(xOHminbulk -xHplusbulk)*vsal ! NaCl+ NaOH   
     xnegbulk=xsalt
  endif

  xsolbulk=1.0 -xHplusbulk -xOHminbulk - xnegbulk -xposbulk

  KaA = (K0A*vs/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant 
  KaNa = (K0ANa*vs/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant 
  KaCl = (K0BCl*vs/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant 
  KaB = (Kw/K0B*vs/xsolbulk)*(Na/1.0d24)
  Keo = (K0eo)*(1.0d24/Na) !!!!!!!!!!!!!!!!!!!!!!!
  

!  cHplusbulk=xHplusbulk
!  cOHminbulk=xOHminbulk
!csalt=xsalt/(vsal*vs*6.02/10)
expmupos=xposbulk/xsolbulk**vsal
expmuneg=xnegbulk/xsolbulk**vsal
expmuHplus=xHplusbulk/xsolbulk ! vHplus=vsol
expmuOHmin=xOHminbulk/xsolbulk ! vOHminus=vsol

  call solve(flagcrash)
!call fe(cc, ccc)         ! calculates and saves free energy to disk
!  call salvar(flagcrash)
!  print*, 'Save OK',yes

enddo ! i

  call endall     ! clean up and terminate


end 
subroutine endall
stop
end subroutine
