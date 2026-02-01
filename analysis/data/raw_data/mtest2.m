ClearAll[]
r = {{a, 0},{0, a}}
Print[r]
u = {{w1, w2}}
v =r + b Transpose[u] . u
Print[v]
vinv = Inverse[v]
Print[vinv]
uu = {w1, w2}
u2 = Mean[uu] 
u3 = uu - u2 
u4 = Power[u3, 2]
u5 = Total[u4]
u6 = Power[u5, .5]
u7 = u3/u6
u8 = Total[Power[u7, 2]]

uu = {w1, w2, w3}
u2 = Mean[uu] 
u3 = uu - u2 
u4 = Power[u3, 2]
u5 = Total[u4]
u6 = Power[u5, .5]
u7 = u3/u6
Total[Power[u7, 2]]

aa = Array[a, 3]
a2 = Mean[aa] 
a3 = aa - a2 
a4 = Power[a3, 2]
a5 = Total[a4]
a6 = Power[a5, .5]
a7 = a3/a6
a8 = Total[Power[a7, 2]]
