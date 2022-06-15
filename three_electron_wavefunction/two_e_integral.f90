!------------------------Added by DC 11/19/2012------------------------------------!
MODULE twoeleint

USE kinds
USE variables

IMPLICIT NONE


SUBROUTINE fullint

k12=(-4.0_dbl*pi**2/b12)*log((b1-b12)/(b1+b12))*log((b2-b12)/(b2+b12))
CALL dilogarithm()
CALL dilogarithm()
k1=(4.0_dbl*pi**2/b12)*(f(-(b2-b12)/(b1+b12))-f(-(b2+b12)/(b1-b12)))
CALL dilogarithm()
CALL dilogarithm()
k2=(4.0_dbl*pi**2/b12)*(f(-(b1-b12)/(b2+b12))-f(-(b1+b12)/(b2-b12)))+ &
                       (1.0_dbl/2.0_dbl)*Log((b2+b12)/(b1+b12))**2-   &
                       (1.0_dbl/2.0_dbl)*Log((b2-b12)/(b1-b12)**2)
BH=k1+k2+k12

absb12=ABS(b12)
T=SIGN(1.0_dbl,B12)*(Log((b1+b12)*(b2+b12))**2-Log((b1+b12)*(b2-b12))**2  &
                   -Log((b1-b12)*(b2+b12))**2+Log((b1-b12)*(b2-b12))**2)

BFH=4.0_dbl*((pi**2)/absb12)*((pi**2/2.0_dbl)+2.0_dbl*V(b2/absb12)+        &
      2.0_dbl*V(b1/absb12)-2.0_dbl*V((b12**2+(b1*b2))/(absb12*(b1+b2)))-   &
                                                            (0.5_dbl)*T)

END SUBROUTINE fullint

SUBROUTINE remiddi

w1=(a2+a3)/a1
w2=(a1+a3)/a2
w3=(a1+a2)/a3

CALL dilogarithm
......6

RM=((32.0_dbl*pi**3)/(a1*a2*a3))*(Log(w1)*log(1+(1+1/w1)))

END SUBROUTINE remiddi
END MODULE twoeleint
