X0 = {{t,0,0}, {0,t,0}}
X1 = {{t,0,0}, {0,t,t}}
Y0 = {{C11},{C12}}
Y1 = {{C21},{C22}}

Z ={{t}, {t}}
G = {{b}}
R = {{2 a,-1 a}, {-1 a,2 a}}
Ri=Inverse[R]

ZRiZ = Transpose[Z] . Ri . Z

H = Simplify[Inverse[(Inverse[G] + ZRiZ)]]
X0RiX0 = Transpose[X0] . Ri . X0
X1RiX1 = Transpose[X1] . Ri . X1
XRiX = X0RiX0 + X1RiX1

X0RiZHZRiX0 = Transpose[X0] . Ri . Z . H . Transpose[Z] . Ri . X0
X1RiZHZRiX1 = Transpose[X1] . Ri . Z . H . Transpose[Z] . Ri . X1
XRiZHZRiX =  (X0RiZHZRiX0 + X1RiZHZRiX1)
XRiZHZRiX = Simplify[XRiZHZRiX]
XSigmaiX = Simplify[XRiX - XRiZHZRiX]
VVV = Simplify[Inverse[XSigmaiX]]
VV2 = Simplify[VVV[[3,3]]]

X0RiY0 = Transpose[X0] . Ri . Y0
X1RiY1 = Transpose[X1] . Ri . Y1
XRiY = X0RiY0+ X1RiY1

ZRiY0 = Transpose[Z] . Ri . Y0
ZRiY1 = Transpose[Z] . Ri . Y1
ZRiY = ZRiY0+ ZRiY1
mat = Simplify[VVV . (XRiY - XRiZ . H . ZRiY) ]
Dimensions[mat]
Dimensions[XRiY]
mat[[3,1]]

% H = ReplaceAll[H,   a -> (d - c)/3]
% H = ReplaceAll[H,   b -> (d + (2 c))/(3 t t)]
H = Simplify[H]
XRiZHZRiX = ReplaceAll[XRiZHZRiX ,   a -> (d - c)/3]
XRiZHZRiX = ReplaceAll[XRiZHZRiX ,   b -> (d + (2 c))/(3 t t)]
XRiZ = X0RiZ + X1RiZ
X1RiZ = Transpose[X1] . Ri . Z
X0RiZ = Transpose[X0] . Ri . Z
XRiX= ReplaceAll[XRiX,   a -> (d - c)/3]
Sigmai = Simplify[Ri - (Ri . Z . H . Transpose[Z] . Ri)]

% Sigmai = ReplaceAll[Sigmai,   2 a + b t t -> c]
% Sigmai = ReplaceAll[Sigmai,    a - b t t -> d]
% Sigmai = ReplaceAll[Sigmai,    a - b t t -> d]
Sigmai = ReplaceAll[Sigmai,   2 a + b t t -> d]
% Sigmai = ReplaceAll[Sigmai,     3 a a + 6 a b t t -> d d - c c  ]
Sigmai = ReplaceAll[Sigmai,     3 a a + 6 a b t t -> 3 a (d - c)  ]
Sigmai = ReplaceAll[Sigmai,    a - b t t -> c]


X0SigmaiX0 = Transpose[X0] . Sigmai . X0
X1SigmaiX1 = Transpose[X1] . Sigmai . X1
XSigmaiX=  X0SigmaiX0 +  X1SigmaiX1
Simplify[XSigmaiX]
K =  (3 a (d - c))/(t t) 
KK = Simplify[K XSigmaiX]
KKK = Simplify[Inverse[KK]]
XSigmaiXi = Simplify[Inverse[XSigmaiX]]
XSigmaiXi[[2,2]]
XSigmaiXi[[3,3]]

VV2 = ReplaceAll[VV2 ,   a -> (d - c)/3]
VV2 = ReplaceAll[VV2 ,   b -> (d + (2 c))/(3 t t)]
Simplify[VV2]

ZRiY0 = Transpose[Z] . Ri . Y0
ZRiY1 = Transpose[Z] . Ri . Y1
ZRiY = ZRiY0+ ZRiY1


W0 =   X0RiZ . H .  Transpose[X0RiZ]
W0 = ReplaceAll[W0,   a -> (d - c)/3]
W0 = ReplaceAll[W0,   b -> (d + (2 c))/(3 t t)]
Simplify[W0]
W1 =   X1RiZ . H .  Transpose[X1RiZ]
W1 = ReplaceAll[W1,   a -> (d - c)/3]
W1 = ReplaceAll[W1,   b -> (d + (2 c))/(3 t t)]
Simplify[W1]
W2 = Simplify[W0 + W1]
XSigmaiX = Simplify[XRiX - W2]

XSigmaiXinv = FullSimplify[Inverse[XRiX - W2]]
zz = XSigmaiXinv[[2,2]]
T0 = X0RiY0 - X0RiZ . H . ZRiY0
T1 = X1RiY1 - X1RiZ . H . ZRiY1
T2 = Simplify[T1 + T0]

beta = XSigmaiXinv . T2
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
