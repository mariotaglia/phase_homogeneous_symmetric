! $UWHPSC/codes/fortran/newton/functions.f90

module functions

contains


real(kind=8) function falpha(x)
    use system
    implicit none
    real(kind=8), intent(in) :: x

    falpha=1.0-x2alpha-x3alpha-vsal*(x**vsal)*(expmupos+expmuneg)-x*(1+expmuOHmin+expmuHplus)

end function falpha

real(kind=8) function fbeta(x)
    use system
    implicit none
    real(kind=8), intent(in) :: x

    fbeta=1.0-x2beta-x3beta-vsal*(x**vsal)*(expmupos+expmuneg)-x*(1+expmuOHmin+expmuHplus)

end function fbeta


real(kind=8) function fp(x)
    use system
    implicit none
    real(kind=8), intent(in) :: x
    
    fp=-vsal*(x**(vsal-1))*vsal*(expmupos+expmuneg)-(1+expmuOHmin+expmuHplus)
end function fp

end module functions
