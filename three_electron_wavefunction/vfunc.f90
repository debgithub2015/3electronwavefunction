!-------------------------Created by DC 11/19/2012---------------------------------!
MODULE vfunction

IMPLICIT NONE
!----------------------------------------------------------------------------------!
! This module calculates the v(z) function                                         !
!----------------------------------------------------------------------------------!
CONTAINS

SUBROUTINE vmain(x)
USE kinds; USE constants
USE variables, ONLY:vfunc
USE function

REAL(dbl),INTENT(IN)::x
REAL(dbl)::fun1,fun2,fun3,fun4,fun5,vfunc1
REAL(dbl)::fout1,fout2,fout3,fout4
LOGICAL::debug

debug=.true.

IF (x < 0.0_dbl) THEN
    x=-x

IF (x == 0.0_dbl) THEN
  vfunc=0.0_dbl
ELSEIF (x == 1.0_dbl)THEN
  vfunc=-(Log(2.0_dbl)**2/4.0_dbl)-(pi**2/12.0_dbl)
ELSEIF (x > 0.0_dbl) THEN
  IF (0.0_dbl< x < 1.0_dbl) THEN
     fun1=(1.0_dbl-x)/2.0_dbl
     fun2=(1.0_dbl+x)/2.0_dbl
     CALL dilogarithm(fun1,fout1)
     CALL dilogarithm(fun2,fout2)
     vfunc=((fout1-fout2)/2.0_dbl)-((Log(fun1)**2-Log(fun2)**2)/4.0_dbl)
  ELSEIF (x > 1.0_dbl)THEN
     fun1=1.0_dbl/x
     fun2=(1.0_dbl-fun1)/2.0_dbl
     fun3=(1.0_dbl+fun1)/2.0_dbl
     CALL dilogarithm(fun2,fout1)
     CALL dilogarithm(fun3,fout2)
     vfunc1=((fout1-fout2)/2.0_dbl)-((Log(fun2)**2-Log(fun3)**2)/4.0_dbl)
     fun4=1.0_dbl-x
     fun5=-x
     CALL dilogarithm(fun4,fout3)
     CALL dilogarithm(fun5,fout4)
     vfunc=vfunc1+fout3+fout4+((Log(x)*Log(1.0_dbl+x))+(pi**2/12.0_dbl))
  END IF
ELSEIF (x < 0.0_dbl) THEN
   x= -x
   IF (0.0_dbl< x < 1.0_dbl) THEN
     fun1=(1.0_dbl-x)/2.0_dbl
     fun2=(1.0_dbl+x)/2.0_dbl
     CALL dilogarithm(fun1,fout1)
     CALL dilogarithm(fun2,fout2)
     vfunc=((fout1-fout2)/2.0_dbl)-((Log(fun1)**2-Log(fun2)**2)/4.0_dbl)
     vfunc=-vfunc
  ELSEIF (x > 1.0_dbl)THEN
     fun1=1.0_dbl/x
     fun2=(1.0_dbl-fun1)/2.0_dbl
     fun3=(1.0_dbl+fun1)/2.0_dbl
     CALL dilogarithm(fun2,fout1)
     CALL dilogarithm(fun3,fout2)
     vfunc1=((fout1-fout2)/2.0_dbl)-((Log(fun2)**2-Log(fun3)**2)/4.0_dbl)
     fun4=1.0_dbl-x
     fun5=-x
     CALL dilogarithm(fun4,fout3)
     CALL dilogarithm(fun5,fout4)
     vfunc=vfunc1+fout3+fout4+((Log(x)*Log(1.0_dbl+x))+(pi**2/12.0_dbl))
     vfunc=-vfunc
  END IF
END IF
   
WRITE(*,*)'The v function',vfunc

END SUBROUTINE vmain(x)

SUBROUTINE vn(x,y)
!------------------------------------------------------------------------------!
!  Calculates the VN function                                                 !
!------------------------------------------------------------------------------!
USE kinds
USE constants
USE variables, ONLY:CV,vnfun
USE function
USE bernoulli

REAL(dbl),INTENT(IN)::x,y
INTEGER(istd)::i
REAL(dbl)::func,fout

func=(1.0_dbl-ABS(x))/2.0_dbl

CALL dilogarithm(func,fout)

vnfun=SIGN(1.0_dbl,x)*(-((Log(ABS(y))**2)/4.0_dbl)-(pi**2/12.0_dbl)+  &
       fout+((Log((1.0_dbl+ABS(g))/2.0_dbl)**2)/2.0_dbl))

WRITE(*,*)'The VS function',vsfun

END SUBROUTINE vn(x)

SUBROUTINE vs(x)
!------------------------------------------------------------------------------!
!  Calculates the VS function                                                 !
!------------------------------------------------------------------------------!
USE kinds
USE variables, ONLY:CV,vsfun
USE bernoulli

REAL(dbl),INTENT(IN)::x
INTEGER(istd)::i
REAL(dbl)::sum

sum=0.0_dbl

DO i=2,10
  sum=sum+CV(i)*(x**(2*i-2))
END DO

vsifun=sum

WRITE(*,*)'The VS function',vsfun

END SUBROUTINE vs(x)

SUBROUTINE vsi(x)
!------------------------------------------------------------------------------!
!  Calculates the VSI function                                                 !
!------------------------------------------------------------------------------!
USE kinds
USE variables, ONLY:CV,vsifun
USE bernoulli

REAL(dbl),INTENT(IN)::x
INTEGER(istd)::i
REAL(dbl)::sum

sum=0.0_dbl

DO i=2,10
  sum=sum+CV(i)*(x**(2*i-2))*(-1.0_dbl)**(i-1)
END DO

vsifun=sum

WRITE(*,*)'The VSI function',vsifun

END SUBROUTINE vsi(x)

SUBROUTINE vcl(x)
!----------------------------------------------------------------------------------!
! Calculates VCL function
!----------------------------------------------------------------------------------!
USE kinds
USE variables, ONLY:BB,vclfun
USE bernoulli

REAL(dbl),INTENT(IN)::x
INTEGER(istd)::i
REAL(dbl)::w,sum

sum=0.0_dbl

CALL arctanh(x)

w=2.0_dbl*arctanh(x)

DO i=1,15
   sum=sum+((BB(i)*(2.0_dbl**(2*i)-1.0_dbl)*w**(2*i))/  &
            (2.0_dbl*i*(2.0_dbl*i+1.0_dbl))) 
END DO

vclfun=-(w*Log(2.0_dbl))+sum

WRITE(*,*)'The VCL function',vclfun

END SUBROUTINE vcl(x) 

REAL FUNCTION arctanh(x)

!--------Calculates the arctanh(x) function analytically--------------!

USE kinds

IMPLICIT NONE

REAL(dbl),INTENT(IN)::x

arctanh=(1.0_dbl/2.0_dbl)*Log((1.0_dbl+x)/(1.0_dbl-x))

WRITE(*,*)'The inverse-hyperbolic function', arctanh(x)

END FUNCTION arctanh(x)

END MODULE vfunction
