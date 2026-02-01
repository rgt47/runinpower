X0 = {{t,0}, {t,0}}
X1 = {{t,t}, {t,t}}
Y0 = {{C11},{C12}}
Y1 = {{C21},{C22}}

Z ={{t}, {t}}
G = {{b}}
R = {{2 a,-1 a}, {-1 a,2 a}}
Ri=Inverse[R]

X0RiX0 = Transpose[X0] . Ri . X0
X1RiX1 = Transpose[X1] . Ri . X1
XRiX = X0RiX0 + X1RiX1

X0RiZ = Transpose[X0] . Ri . Z
X1RiZ = Transpose[X1] . Ri . Z
XRiZ = X0RiZ + X1RiZ

X0RiY0 = Transpose[X0] . Ri . Y0
X1RiY1 = Transpose[X1] . Ri . Y1
XRiY = X0RiY0+ X1RiY1

ZRiY0 = Transpose[Z] . Ri . Y0
ZRiY1 = Transpose[Z] . Ri . Y1
ZRiY = ZRiY0+ ZRiY1

ZRiZ = Transpose[Z] . Ri . Z

H = Simplify[Inverse[(Inverse[G] + ZRiZ)]]

W0 =   X0RiZ . H .  Transpose[X0RiZ]
W1 =   X1RiZ . H .  Transpose[X1RiZ]
W2 = Simplify[W0 + W1]

XSigmaiXinv = Simplify[Inverse[XRiX - W2]]

T0 = X0RiY0 - X0RiZ . H . ZRiY0
T1 = X1RiY1 - X1RiZ . H . ZRiY1
T2 = Simplify[T1 + T0]

beta = XSigmaiXinv . T2]
EstGamma = Simplify[beta[[2,1]]]
VarGamma = Simplify[XSigmaiXinv[[2,2]]]

sigInv = Simplify[Ri - Ri . Z . H . Transpose[Z] . Ri]
ZZ = Simplify[Transpose[X0] . sigInv . Y0] + Simplify[Transpose[X1] . sigInv . Y1]
Simplify[ZZ]
Simplify[XSigmaiXinv . ZZ]


e = {{2,1}}
ee = Transpose[e]
Dimensions[ee]
ee . e
