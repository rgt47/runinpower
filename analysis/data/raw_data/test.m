R = {{2 a,-1 a}, {-1 a,2 a}}
Rinv=Inverse[R]
Rinv = Simplify[Rinv]
X0 = {{t,0}, {t,0}}
X1 = {{t,t}, {t,t}}
Z ={{t}, {t}}
G = {{b}}
W = Simplify[Inverse[R + Z . Inverse[G] . Transpose[Z]]]
Print[W]
ZprimeRinvZ = Simplify[Transpose[Z] . Rinv . Z]
H = Inverse[(Inverse[G] + ZprimeRinvZ)]

H = Simplify[H]
F =   Rinv - Rinv . Z . H . Transpose[Z]  . Rinv
F = Simplify[F]
Print[F]
(*
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
Print[VarGamma]
*)
