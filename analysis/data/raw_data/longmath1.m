Clear["Global`*"]


rr = {{a, 0},{0, a}};
tt = Array[t, 2];
v =rr + b Transpose[{tt}] . {tt};
vinv = Inverse[v];
vinv2 = Simplify[vinv,  1 == t[1]^2+t[2]^2 ] ;
vinv3 = ReplaceAll[vinv2, b -> r (a + b)]
vinv4 = ReplaceAll[vinv3, t[2]^2 -> 1 - t[1]^2]
vinv4[[2,2]] = ReplaceAll[vinv4[[2,2]], t[1]^2 -> 1 - t[2]^2]
vinv5 = Simplify[vinv4, b == r (a + b)]
vinv6 = ReplaceAll[vinv5, a^2 + a b -> a (a + b)]
vinv7 = ReplaceAll[vinv6, a + b -> b/r]
vinv8 = Simplify[vinv7]

xxt = Array[f, {2, 3}]
xxt[[All,1]] = 1
xxt[[All,{2,3}]] = tt
xxc = xxt
xxc[[All,3]] = 0
xpx = Simplify[Transpose[xxt] . xxt, 1 == t[1]^2+t[2]^2]
xpx = Simplify[xpx, 0 == t[1]+t[2]]
xt = Simplify[Transpose[xxt] . tt, 1 == t[1]^2+t[2]^2]
xt = Simplify[xt, 0 == t[1]+t[2]]; xt
xttx = Simplify[xt . Transpose[xt], 0 == t[1]+t[2]]
xpvinvxt = Simplify[Transpose[xxt] . vinv8 . xxt];
xpvinvxt2 = ReplaceAll[xpvinvxt, t[1]+t[2] -> 0];
xpvinvxt3 = ReplaceAll[xpvinvxt2, t[1]^2+t[2]^2 -> 1]
xpvinvxc = Simplify[Transpose[xxc] . vinv8 . xxc];
xpvinvxc2 = ReplaceAll[xpvinvxc, t[1]+t[2] -> 0];
xpvinvxc3 = ReplaceAll[xpvinvxc2, t[1]^2+t[2]^2 -> 1]
xpvinvx = (m xpvinvxt3)  + (m xpvinvxc3)
xpvinvxinv = Simplify[Inverse[xpvinvx]]
var = xpvinvxinv
yt = Array[y, 2];
yc = Array[z, 2];
xpvinvyt = Simplify[Transpose[xxt] . vinv8 . yt]
xpvinvyt[[1]] = ReplaceAll[xpvinvyt[[1]], t[1] + t[2] -> 0]
xpvinvyt = ReplaceAll[xpvinvyt, t[1]^2+t[2]^2 -> 1]

xpvinvyc = Simplify[Transpose[xxc] . vinv8 . yc]
xpvinvyc[[1]] = ReplaceAll[xpvinvyc[[1]], t[1] + t[2] -> 0]
xpvinvyc = ReplaceAll[xpvinvyc, t[1]^2+t[2]^2 -> 1]


xpvinvy = (m xpvinvyt)  + (m xpvinvyc)
hat = Simplify[xpvinvxinv . xpvinvy]



