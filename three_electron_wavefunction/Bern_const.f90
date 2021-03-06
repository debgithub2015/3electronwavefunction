!-----------------------Created by DC 11/19/2012----------------------------------!
MODULE Bernoulli

USE kinds
USE variables,ONLY:BB,CV

IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assigning the Bernoulli Numbers and MacLaurin Series for v(z)                   !
! Depending on the complexity of the programme one can control the order of the polynomial !

REAL(dbl)::X

ALLOCATE(BB(1:15),CV(1:10)) 

BB(1)= 1.0_dbl/6.0_dbl;
BB(2)= -1.0_dbl/30.0_dbl;
BB(3)= 1.0_dbl/42.0_dbl;
BB(4)= -1.0_dbl/30.0_dbl;
BB(5)= 5.0_dbl/66.0_dbl;
BB(6)= -691.0_dbl/2730.0_dbl;
BB(7)= 7.0_dbl/6.0_dbl;
BB(8)= -3617.0_dbl/510.0_dbl;
BB(9)= 43867.0_dbl/798.0_dbl;
BB(10)= -174611.0_dbl/330.0_dbl;
BB(11)= 854513.0_dbl/138.0_dbl;
BB(12)= -236364091.0_dbl/2730.0_dbl;
BB(13)= 8553103.0_dbl/6.0_dbl;
BB(14)= -23749461029.0_dbl/870.0_dbl;
BB(15)= 8615841276005.0_dbl/14322.0_dbl;

X=Log(2.0_dbl)
CV(1)=-2.0_dbl*X;       !      CV[n] is coefficient of z^(2n-1)
CV(2)=-((2.0_dbl/3.0_dbl)*X) - (1.0_dbl/3.0_dbl);
CV(3)=-((2.0_dbl/5.0_dbl)*X) - (3.0_dbl/10.0_dbl);
CV(4)=-((2.0_dbl/7.0_dbl)*X) - (11.0_dbl/42.0_dbl);
CV(5)=-((2.0_dbl/9.0_dbl)*X) - (25.0_dbl/108.0_dbl);
CV(6)=-((2.0_dbl/11.0_dbl)*X) - (137.0_dbl/660.0_dbl);
CV(7)=-((2.0_dbl/13.0_dbl)*X) - (49.0_dbl/260.0_dbl);
CV(8)=-((2.0_dbl/15.0_dbl)*X) - (121.0_dbl/700.0_dbl);
CV(9)=-((2.0_dbl/17.0_dbl)*X) - (761.0_dbl/4760.0_dbl);
CV(10)=-((2.0_dbl/19.0_dbl)*X) - (7129.0_dbl/47880.0_dbl);
 
END MODULE Bernoulli
