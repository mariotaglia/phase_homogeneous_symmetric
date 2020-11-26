module system
!real*16 delta   ! delta is the discretization lenght in z direction
!real*16 sigmaA
real*16 vsal
real*16 Ma
real*16 rsal
real*16 vpol
real*16 phi_sal
real*16 phi_pol
integer yes
real*16 Mb
integer npasos
integer npasosgrid
real*16 phimin, phimax
real*16 sigmaB
real*16 csalt
real*16 nor
real*16 KaA
real*16 KaB
real*16 Knew
real*16 pKaA
real*16 pkaB
real*16 pKaNa
real*16 pKaCl
real*16 pkEo
real*16 Keo,csal,csalini,csalfin
integer ncsal
real*16 KaHplus
real*16 KaOHmin
real*16 pKaHplus,pKHbulk,cHplusbulk,pHbulk
real*16 pKaOHmin,pOHbulk,cOHminbulk
real*16 cNaplus,cClmin,pClbulk,pNabulk
real*16 xNaplus, xClmin,xmNaplus,xmClmin
real*16 KaNa,KaCl
real*16 expmupos, expmuneg, expmuHplus, expmuOHmin  
real*16 x2alpha,x3alpha,x2beta,x3beta

!real*16 pHbulk
!real*8 st
!integer VUELTA
endmodule

module results
real*16 x2alphafixed
!real*16 x3alphafixed
real*16 arrayalpha(2,400000)
real*16 arraybeta(2,400000)
real*16 finalalpha(2,400000)
real*16 finalbeta(2,400000)
real*16 cargaalpha(1,400000)
real*16 cargabeta(1,400000)
integer conteo
real*16 arrayalphaspi(2,400000)
real*16 arraybetaspi(2,400000)
real*16 array_fe_alpha(2,400000)
real*16 array_fe_beta(2,400000)
real*16 array_mu2_alpha(2,400000)
real*16 array_mu2_beta(2,400000)
real*16 array_mu3_alpha(2,400000)
real*16 array_mu3_beta(2,400000)
real*16 sumalphaNaCl(1,400000)
real*16 sumbetaNaCl(1,400000)
real*16 arraysolvalpha(1,400000)
real*16 arraysolvbeta(1,400000)

real*16 arrayNaplusalpha(1,400000)
real*16 arrayNaplusbeta(1,400000)

real*16 arrayClminalpha(1,400000)
real*16 arrayClminbeta(1,400000)

real*16 arrayHplusalpha(1,400000)
real*16 arrayHplusbeta(1,400000)

real*16 arrayOHminalpha(1,400000)
real*16 arrayOHminbeta(1,400000)

real*16 array_fchA_alpha(1,400000)
real*16 array_fucA_alpha(1,400000)
real*16 array_fIPA_alpha(1,400000)
real*16 array_fPPA_alpha(1,400000)

real*16 array_fchB_alpha(1,400000)
real*16 array_fucB_alpha(1,400000)
real*16 array_fIPB_alpha(1,400000)
real*16 array_fPPB_alpha(1,400000)

real*16 array_fchA_beta(1,400000)
real*16 array_fucA_beta(1,400000)
real*16 array_fIPA_beta(1,400000)
real*16 array_fPPA_beta(1,400000)

real*16 array_fchB_beta(1,400000)
real*16 array_fucB_beta(1,400000)
real*16 array_fIPB_beta(1,400000)
real*16 array_fPPB_beta(1,400000)


integer cont
endmodule

module const
real*16, parameter :: pi = 3.14159 ! pi 
real*16, parameter :: Na = 6.02d23 ! Avogadro's number
!real*16, parameter :: lb = 0.714   ! bjerrum lenght in water in nm
real*16, parameter :: vs = 0.03  ! bjerrum lenght in water in nm
!real*16, parameter :: vp =0.195 ! 0.11 ! 
real*16 xsolbulk
real*16 vp
!real*16 constq
!real*16 pKw
endmodule


module solver
!real*16, allocatable :: xflag(:)
!integer infile
!real*16, parameter :: error = 1.0d-6 ! maximum kinsol norm
real*16 norma
integer iter
integer linearsolver
endmodule


!module molecules
!real*16 zpos, zneg, zpolA, zpolB ! charges of cation, anions and polyelectrolyte segment
!real*16 vsalt, vpol   ! volume of salt and polyelectrolyte segments in units of vsol
!real*16 vsol             ! solvent volume 
!real*16 K0A, K0B ,K0ANa,K0BCl, K0Eo!K0
!endmodule

!module bulk
!real*16 xHplusbulk, xOHminbulk ! bulk volume fraction of H+ and OH-
!real*16 xposbulk, xnegbulk     ! bulk volume fraction of cation and anion
!real*16 xsolbulk               ! bulk volume fraction of solvent
!real*16 expmupos, expmuneg, expmuHplus, expmuOHmin  ! exp(-beta*mu)*(bulk volume fraction), where mu is the chemical potential
!endmodule

!module rand
!integer seed
!endmodule

module mkinsol
double precision, allocatable :: pp(:)
endmodule
