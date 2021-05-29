
from matplotlib.pylab import *
from mpl_toolkits import mplot3d

xi  = linspace(-1, 1, 20)    #chi
eta = linspace(-1, 1, 20)    #n


XI, ETA = meshgrid(xi,eta)
#print(f"Xi = {XI}")
#print(f"ETA= {ETA}")


N0 = (XI-1)*(ETA-1)*XI*ETA/4
N1 = (XI+1)*(ETA-1)*XI*ETA/4
N2 = (XI+1)*(ETA+1)*XI*ETA/4
N3 = (XI-1)*(ETA+1)*XI*ETA/4
N4 = (1-XI**2)*(ETA-1)*ETA/2
N5 = (1-ETA**2)*(XI+1)*XI/2
N6 = (1-XI**2)*(ETA+1)*ETA/2
N7 = (1-ETA**2)*(XI-1)*XI/2
N8 = (1-ETA**2)*(1-XI**2)


fig = figure()
ax = axes(projection='3d')
ax.plot_surface(XI,ETA,N7)
show()

'''

from sympy import *
from sympy.matrices import Matrix


x0 = Symbol('x0')
y0 = Symbol('y0')
x1 = Symbol('x1')
y1 = Symbol('y1')
x2 = Symbol('x2')
y2 = Symbol('y2')
x3 = Symbol('x3')
y3 = Symbol('y3')
x4 = Symbol('x4')
y4 = Symbol('y4')
x5 = Symbol('x5')
y5 = Symbol('y5')
x6 = Symbol('x6')
y6 = Symbol('y6')
x7 = Symbol('x7')
y7 = Symbol('y7')
x8 = Symbol('x8')
y8 = Symbol('y8')

xi = Symbol('xi')
eta = Symbol('eta')

#Shape functions
N0 = (xi-1)*(eta-1)*xi*eta/4
N1 = (xi+1)*(eta-1)*xi*eta/4
N2 = (xi+1)*(eta+1)*xi*eta/4
N3 = (xi-1)*(eta+1)*xi*eta/4
N4 = (1-xi**2)*(eta-1)*eta/2
N5 = (1-eta**2)*(xi+1)*xi/2
N6 = (1-xi**2)*(eta+1)*eta/2
N7 = (1-eta**2)*(xi-1)*xi/2
N8 = (1-eta**2)*(1-xi**2)


dN0_dxi = N0.diff(xi)
dN0_deta = N0.diff(eta)

dN1_dxi = N1.diff(xi)
dN1_deta = N1.diff(eta)

dN2_dxi = N2.diff(xi)
dN2_deta = N2.diff(eta)

dN3_dxi = N3.diff(xi)
dN3_deta = N3.diff(eta)

dN4_dxi = N4.diff(xi)
dN4_deta = N4.diff(eta)

dN5_dxi = N5.diff(xi)
dN5_deta = N5.diff(eta)

dN6_dxi = N6.diff(xi)
dN6_deta = N6.diff(eta)

dN7_dxi = N7.diff(xi)
dN7_deta = N7.diff(eta)

dN8_dxi = N8.diff(xi)
dN8_deta = N8.diff(eta)


print(f"dN0_dxi= {dN0_dxi}")
print("")
print(f"dN0_deta= {dN0_deta}")
print("")

print(f"dN1_dxi= {dN1_dxi}")
print("")
print(f"dN1_deta= {dN1_deta}")
print("")

print(f"dN2_dxi= {dN2_dxi}")
print("")
print(f"dN2_deta= {dN2_deta}")
print("")

print(f"dN3_dxi= {dN3_dxi}")
print("")
print(f"dN3_deta= {dN3_deta}")
print("")

print(f"dN4_dxi= {dN4_dxi}")
print("")
print(f"dN4_deta= {dN4_deta}")
print("")

print(f"dN5_dxi= {dN5_dxi}")
print("")
print(f"dN5_deta= {dN5_deta}")
print("")

print(f"dN6_dxi= {dN6_dxi}")
print("")
print(f"dN6_deta= {dN6_deta}")
print("")

print(f"dN7_dxi= {dN7_dxi}")
print("")
print(f"dN7_deta= {dN7_deta}")
print("")

print(f"dN8_dxi= {dN8_dxi}")
print("")
print(f"dN8_deta= {dN8_deta}")
print("")


x = x0*N0 + x1*N1 + x2*N2 + x3*N3 + x4*N4 + x5*N5 + x6*N6 + x7*N7 + x8*N8 
y = y0*N0 + y1*N1 + y2*N2 + y3*N3 + y4*N4 + y5*N5 + y6*N6 + y7*N7 + y8*N8

print(f"x={x}")
print(f"y={y}")


dx_dxi = x.diff(xi)
dx_deta = x.diff(eta)

dy_dxi = y.diff(xi)
dy_deta = y.diff(eta)

print("")
print(f"dx_dxi= {dx_dxi}")
print("")
print(f"dx_deta= {dx_deta}")
print("")
print(f"dy_dxi= {dy_dxi}")
print("")
print(f"dy_deta= {dy_deta}")

'''