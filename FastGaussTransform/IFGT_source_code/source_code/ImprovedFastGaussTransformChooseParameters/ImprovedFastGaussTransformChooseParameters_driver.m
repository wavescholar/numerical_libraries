d=2
h=0.3
epsilon=1e-6
Klimit=round(40*sqrt(d)/h)

[K,p_max,r]=ImprovedFastGaussTransformChooseParameters(d,h,epsilon,Klimit)

