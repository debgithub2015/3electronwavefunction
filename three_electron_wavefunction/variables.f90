!------------------------Created by DC 11/19/2012----------------------------------!
MODULE variables

!-----------------------------------------------------------------------------------!
! This module declares the global variable for the entire programme                 !
! BB(1:15):: Bernoulli numbers                                                      ! 
! CV(1:10)::MacLaurin coefficients                                                  !
! chireal -- real part of chi function                                              !
! chiimag -- imaginary part of the chi function                                     !
! clausen -- carries the numerical value of Clausen function                        !
! fdilog -- numerical value of dilogarithm under branch constraint                  !
! vclfun -- VCL function                                                            !
! vsifun -- Imaginary part of VS function                                           !
! vsfun -- Real part of VS function                                                 !
! vnfun -- VN function                                                              !
! vfunc -- v function                                                               !
!-----------------------------------------------------------------------------------!

USE kinds

IMPLICIT NONE 

REAL(dbl),ALLOCATABLE::BB,CV,aij,mu,g,q,gam
REAL(dbl)::chireal,chiimag,clausen,fdilog
REAL(dbl)::vclfun,vsifun,vsfun,vnfun,vfunc

END MODULE variables
