restart
R = QQ[x,y,z]
curve = ideal( x^4-y^5, x^3-y^7 )
gens gb curve
dim curve
codim curve
degree curve
curve1 = saturate(curve,ideal(x))
curve2 = saturate(curve,curve1)
curve == radical curve
curve = curve1
curve == radical curve
degree curve


 surface = ideal( x^5 + y^5 + z^5 - 1)
theirunion = intersect(curve,surface)
 curve*surface == theirunion
ourpoints = curve + surface
dim ourpoints
degree ourpoints
degree radical ourpoints

 R’ = ZZ/101[x,y,z]
 ourpoints’ = substitute(ourpoints,R’)
 decompose ourpoints’
 oo / degree
rationalPoints ourpoints’ 
loadPackage "RationalPoints"
rationalPoints ourpoints’

restart
R = QQ[vars(0..17)];
M = coker genericMatrix(R,a,3,6)
isHomogeneous M
codim M
degree M
hilbertSeries M
reduceHilbert hilbertSeries M
C = resolution M
C.dd_3
betti C

 restart
 ringP3 = QQ[x_0..x_3]
ringP1 = QQ[s,t]
 cubicMap = map(ringP1,ringP3,{s^3, s^2*t, s*t^2, t^3})
 idealCubic = kernel cubicMap
idealCubic2 = monomialCurveIdeal(ringP3,{1,2,3})
 M = matrix{{x_0,x_1,x_2},{x_1,x_2,x_3}}
 idealCubic3 = minors(2, M)
idealCubic==idealCubic2
idealCubic==idealCubic3
 codim idealCubic
degree idealCubic 
dim idealCubic

restart
cubic = Proj(QQ[x_0..x_3]/ideal(x_0*x_2-x_1^2,x_1*x_3-x_2^2,x_0*x_3-x_1*x_2))
codim cubic
dim cubic
tangentSheaf(cubic)
cotangentSheaf(cubic)

restart
Quintic = Proj(QQ[x_0..x_4]/ideal(x_0^5+x_1^5+x_2^5+x_3^5+x_4^5-101*x_0*x_1*x_2*x_3*x_4))
codim Quintic
dim Quintic
singularLocus(Quintic)
HH^1(cotangentSheaf(Quintic))
HH^2(cotangentSheaf(Quintic))
HH^1(cotangentSheaf(2,Quintic))
HH^2(cotangentSheaf(2,Quintic))
hh^(1,1)(Quintic)
hh^(1,2)(Quintic)
hh^(2,1)(Quintic)
hh^(2,2)(Quintic)
euler(Quintic)

restart
R = QQ[x_0..x_3];
F = sum(4,i->x_i^4)
S = R/F;
matrix table(3,3,(p,q)-> hh^(p,q)(Proj(S)))

