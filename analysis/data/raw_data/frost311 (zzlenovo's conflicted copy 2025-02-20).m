R = {{2 a,-1 a}, {-1 a,2 a}}
Print[R]
Rinv=Inverse[R]
Rinv = Simplify[Rinv]
Print[Rinv]

X0 = {{t,0}, {t,0}}
X1 = {{t,t}, {t,t}}

Z ={{t}, {t}}
G = {{b}}
ZprimeRinvZ = Simplify[Transpose[Z] . Rinv . Z]

H = Inverse[(Inverse[G] + ZprimeRinvZ)]
H = Simplify[H]

XprimeRinvX0 = Simplify[Transpose[X0] . Rinv . X0]
XprimeRinvX1 = Simplify[Transpose[X1] . Rinv . X1]
XprimeRinvX = XprimeRinvX0 + XprimeRinvX1

XprimeRinvZ0 = Simplify[Transpose[X0] . Rinv . Z]
XprimeRinvZ1 = Simplify[Transpose[X1] . Rinv . Z]

W0 =   XprimeRinvZ0 . H .  Transpose[XprimeRinvZ0]

W1 =   XprimeRinvZ1 . H .  Transpose[XprimeRinvZ1]

W = W0 + W1

XprimeSigmainvXinv = Simplify[Inverse[XprimeRinvX - W]]

VarGamma = XprimeSigmainvXinv[[2,2]]
