import numpy as np

## Function to compute a double factorial; see https://en.wikipedia.org/wiki/Double_factorial
def doublefactorial(n):
     if n <= 0:
        return 1
     else:
        return n * doublefactorial(n-2)

## Compute the Newtonian prefix m
def Newtonian_m(m1,m2,l,m):
    epsilon = (l+m)%2
    doubfact = doublefactorial(2*l+1)
    n = np.power(complex(0,m), l)*8*np.divide(np.pi,doubfact)
    if epsilon==0:
        n *= np.sqrt(np.divide((l+1)*(l+2),l*(l-1)))
    elif epsilon==1:
        num = (2*l + 1)*(l + 2)*(l*l - m*m)
        div = (2*l - 1)*(l + 1)*(l)*(l - 1)
        n *= -2.j*np.sqrt( np.divide( num , div ) )
    else:
        print("Epsilon must be 0 or 1")
        exit()
    return n

## Compute the Newtonian prefix c
def Newtonian_c(m1,m2,l,m):
    Mtot = m1 + m2
    m1hat = np.divide(m1,Mtot)
    m2hat = np.divide(m2,Mtot)
    epsilon = (l+m)%2
    if (m%2)==0:
        sign = 1
    else:
        sign = -1
    lpepm1 = l + epsilon - 1
    if (m1!=m2) or sign==1:
        c = np.power(m2hat,lpepm1) + sign*np.power(m1hat,lpepm1)
    else:
        if l==2 or l==3:
            c = -1.
        elif l==4 or l==5:
            c = -0.5
        else:
            c = 0.
    return c

## Compute the Associate Legendre Polynomial for input value 0
def AssociatedLegendre(l,m):
    if l==1:
        if m==1:
            return -1.
        else:
            print("You used a bad (l,m)")
    if l==2:
        if m==2:
            return 3.
        elif m==1:
            return 0.
        else:
            print("You used a bad (l,m)")
    if l==3:
        if m==3:
            return -15.
        elif m==2:
            return 0.
        elif m==1:
            return 1.5
        else:
            print("You used a bad (l,m)")
    if l==4:
        if m==4:
            return 105.
        elif m==3:
            return 0.
        elif m==2:
            return -7.5
        elif m==1:
            return 0.
        else:
            print("You used a bad (l,m)")
    if l==5:
        if m==5:
            return -945.
        elif m==4:
            return 0.
        elif m==3:
            return 52.5
        elif m==2:
            return 0.
        elif m==1:
            return -1.875
        else:
            print("You used a bad (l,m)")
    if l==6:
        if m==6:
            return 10395.
        elif m==5:
            return 0.
        elif m==4:
            return -472.5
        elif m==3:
            return 0.
        elif m==2:
            return 13.125
        elif m==1:
            return 0.
        else:
            print("You used a bad (l,m)")
    if l==7:
        if m==7:
            return -135135.
        elif m==6:
            return 0.
        elif m==5:
            return 5197.5
        elif m==4:
            return 0.
        elif m==3:
            return -118.125
        elif m==2:
            return 0.
        elif m==1:
            return 2.1875
        else:
            print("You used a bad (l,m)")
    if l==8:
        if m==8:
            return 2027025.
        elif m==7:
            return 0.
        elif m==6:
            return -67567.5
        elif m==5:
            return 0.
        elif m==4:
            return 1299.375
        elif m==3:
            return 0.
        elif m==2:
            return -19.6875
        elif m==1:
            return 0.
        else:
            print("You used a bad (l,m)")

## Compute the absolute value of Scalar Spherical Harmonics at the equatorial plane
def SphericalHarmonicAtPiOver2(l,m):
    absM = abs( m )
    