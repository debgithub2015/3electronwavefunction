!-----------------------Created by DC 11/19/2012-----------------------------------!
MODULE function
!----------------------------------------------------------------------------------!
!  This module calculates the function, like dilogarithms, Clausen function and chi!
!  function                                                                        !
!----------------------------------------------------------------------------------!
IMPLICIT NONE

CONTAINS

SUBROUTINE chi(x)
!-----------------------------Subroutine chi---------------------------------------!
! This subroutine calculates the chi function                                      !
!  Re[chi(x)]=Sum(y^(2*n)/(2*n+1)^2,n=1,10)                                        !
!  Im[chi(x)]=Sum(y^(2*n)*(-1)^n/(2*n+1)^2,n=1,10)                                 !
!----------------------------------------------------------------------------------!
USE kinds; USE variables,ONLY:chireal,chiimag

INTEGER(istd)::i
REAL(dbl),INTENT(IN)::x
REAL(dbl)::sumchiRe,sumchiIm,chiRe(1:10),chiIm(1:10)

sumchiRe=0.0_dbl
sumchiIm=0.0_dbl

DO i=1,10
  chiRe(i)=(x**(2*i))/((2.0_dbl*i+1.0_dbl)**2)
  sumchiRe=sumchiRe+chiRe(i)
  chiIm(i)=((x**(2*i))*((-1.0_dbl)**i))/((2.0_dbl*i+1.0_dbl)**2)
  sumchiIm=sumchiIm+chiIm(i)
END DO

chireal=sumchiRe
chiimag=sumchiIm

WRITE(*,*)"Sum of chi real and imaginary",chireal,chiimag

END SUBROUTINE chi(x)

SUBROUTINE Clausenfunc(x)
!--------------------------Subroutine Clausen-------------------------------------- !
!   This subroutine calculates the Clausen function, with added branch constraints  !
!-----------------------------------------------------------------------------------!
USE kinds
USE constants
USE variables, ONLY:BB,clausen ; USE Bernoulli

INTEGER(istd)::i
REAL(dbl),INTENT(IN)::x
REAL(dbl)::sumclausen

sumclausen=0.0_dbl

DO i=1,15
   IF (x.le.2.0944_dbl)THEN
       sumclausen=sumclausen+((-1.0_dbl**(i-1)*BB(i)*x**(2*i))/ &
                           (2.0_dbl*i*(2.0_dl*i+1.0_dbl)))
   ELSE IF (x.gt.4.18879)THEN
       sumclausen=sumclausen+((-1.0_dbl**(i-1)*BB(i)*(2.0_dbl*pi-x)**(2*i))/ &
                           (2.0_dbl*i*(2.0_dl*i+1.0_dbl)))      
   ELSE
       sumclausen=sumclausen+(-1.0_dbl**(i-1)*BB(i)*(2.0_dbl**(2*i)-1.0_dbl)*(x-pi)**(2*i))/  &
                           (2.0_dbl*i*(2.0_dbl*i+1.0_dbl))       
   END IF
END DO

IF (x.le.2.0944_dbl)THEN
     clausen=x*(1.0_dbl-Log(ABS(x)))+sumclausen
ELSE IF  (x.gt.4.18879)THEN
     clausen=-((2.0_dbl*pi-x)*(1.0_dbl-Log(ABS(x)))+sumclausen)
ELSE
     clausen=(x-pi)*Log(2.0_dbl)+sumclausen
END IF

WRITE(*,*)'Clausen function',clausen

END SUBROUTINE Clausenfunc(x)

SUBROUTINE dilogarithm(x,fout)
!------------------------Subroutine dilogarithm(x)---------------------------------!
!     This subroutine calculates dilogarithm according to branch condition         !
!----------------------------------------------------------------------------------!
USE kinds
USE constants
USE variables, ONLY:fdilog

REAL(dbl),INTENT(IN)::x
REAL(dbl),INTENT(OUT)::fout
REAL(dbl)::w,fx

IF (ABS(x).le. (1.0_dbl/2.0_dbl))THEN
     CALL dilogarithmseries(x,fdilog)
ELSEIF(x.le.-1.0_dbl)THEN
      w=1.0_dbl/(1.0_dbl-x)
      CALL dilogarithmseries(w,fx)
      fdilog=fx+(Log(w)*Log(w*x**2)/2.0_dbl)-(pi**2/6.0_dbl)
ELSEIF(x.lt.-(1.0_dbl/2.0_dbl))THEN
      w=x/(x-1.0_dbl)
      CALL dilogarithmseries(w,fx)
      fdilog=-fx-((Log(-w/x)**2)/2.0_dbl)
ELSEIF(x.eq.1.0_dbl)THEN
      fdilog=(pi**2)/6.0_dbl
ELSEIF(x.lt.1.0_dbl)THEN
      w=1.0_dbl-x
      CALL dilogarithmseries(w,fx)
      fdilog=-fx-(Log(x)*Log(w))+((pi**2)/6.0_dbl)
ELSEIF(x.le.2.0_dbl)THEN
      w=1.0_dbl-(1.0_dbl/x)
      CALL dilogarithmseries(w,fx)
      fdilog=fx-(Log(x)*Log(w))-(Log(x)**2)/2.0_dbl-(pi**2/6.0_dbl)
ELSE 
      CALL dilogarithmseries(1.0_dbl/x,fx)
      fdilog=fx-(Log(x)**2/2.0_dbl)+pi**(2.0_dbl/3) 
END IF

fout=fdilog

END SUBROUTINE dilogarithm(x,fout)

SUBROUTINE dilogarithmseries(x,fx)
!-------------------------Subroutine dilogarithmseries(x,n,fx)---------------------!
!    This subroutine calculates the dilogarithm series                             !
!----------------------------------------------------------------------------------!
USE kinds

INTEGER(istd)::i
REAL(dbl),INTENT(IN)::x
REAL(dbl),INTENT(OUT)::fx
REAL(dbl)::sum

sum=0.0_dbl

IF (ABS(x)>(1.0_dbl/2.0_dbl)) THEN
     WRITE(*,*) 'The dilogarithmic function is not defined'
ELSE
DO i=1,30
  sum=sum+(x**i)/(i**2)
END DO
END IF

fx=sum

WRITE(*,*)'dilogarithm series',fx

END SUBROUTINE dilogterm(x,fx)
