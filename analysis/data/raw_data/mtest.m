d = {{v0, v01},{v01, v1}}
Print[d]
zpz={{n,0},{0,s}}
xpx={{n,0},{0,s}}
j = zpz+(ve Inverse[d])
jinv = Inverse[j]
xvx=(xpx-(xpx . jinv . xpx))/ve
xvx= Simplify[xvx]
xvxinv=Inverse[m xvx]
Print[d]
aa = Simplify[xvxinv]
Print[aa];
Quit[];
