!--------------------Created by DC on 11/20/2012-----------------------------------!

MODULE build_twoelfunc

IMPLICIT NONE

CONTAINS

SUBROUTINE sigma

END SUBROUTINE sigma

SUBROUTINE remsingsigma

END SUBROUTINE remsingsigma

SUBROUTINE muandgamma
USE kinds
USE variables, ONLY: mu,g,q,gam

INTEGER(istd)::i,j,i1,i2,i3

ALLOCATE(mu(1:4,1:4),g(1:4,1:4),q(1:3),gam(1:3))

DO j=1,4
  i1=2
  i2=3
  i3=4
  IF (j==2) THEN
      i1=1
  ELSEIF(j==3)THEN
      i2=1
  ELSEIF(j==4)THEN
      i3=1
  END IF
  mu(j,j)=2.0_dbl*a(i1,i2)*a(i1,i3)*a(i2,i3)
  mu(j,i1)=a(i2,i3)*(a(i1,i2)**2+a(i1,i3)**2-a(j,i1)**2);
  mu(j,i2)=a(i1,i3)*(a(i2,i1)**2+a(i2,i3)**2-a(j,i2)**2);
  mu(j,i3)=a(i1,i2)*(a(i3,i1)**2+a(i3,i2)**2-a(j,i3)**2);
  g(j,j)=SUM('mu[j,i]','i'=0..3);


END SUBROUTINE muandgamma

SUBROUTINE inteledist
USE kinds
USE variables, ONLY:aij

INTEGER(istd)::i,j

ALLOCATE(aij(1:4,1:4))

a(1,2)=a1
a(1,3)=a2
a(1,4)=a3
a(2,3)=a12
a(3,4)=a23
a(2,4)=a31

DO i=2,4
 DO j=1,j-1
    a(j,i)=a(i,j)
 END DO
END DO  

WRITE(*,*)'Inter electronic distance',aij

END SUBROUTINE inteledist

END MODULE build_twoelfunc
