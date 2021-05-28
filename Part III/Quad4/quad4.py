from numpy import array, sqrt, zeros, ix_
from scipy.linalg import det, inv

def quad4(xy, properties):

    E = properties["E"]
    ν = properties["nu"]
    bx = properties["bx"]
    by = properties["by"]
    t = properties["t"]

    Eσ = E / (1-ν**2) * array(
        [
        [1 , ν , 0       ]       ,
        [ν , 1 , 0       ]       ,
        [0 , 0 , (1-ν)/2 ]
        ])

    x0 = xy[0,0]
    x1 = xy[1,0]
    x2 = xy[2,0]
    x3 = xy[3,0]

    y0 = xy[0,1]
    y1 = xy[1,1]
    y2 = xy[2,1]
    y3 = xy[3,1]

    ke = zeros((8,8))
    fe = zeros((8,1))

    #Primer punto de Gauss de la regla 2x2
    # xi = 1.0 / sqrt(3)
    # eta = -1.0 / sqrt(3)
    # wi = 1.0
    # wj = 1.0

    gauss_rule = [
        (-1.0 / sqrt(3), -1.0 / sqrt(3), 1.0, 1.0),
        ( 1.0 / sqrt(3), -1.0 / sqrt(3), 1.0, 1.0),
        ( 1.0 / sqrt(3),  1.0 / sqrt(3), 1.0, 1.0),
        (-1.0 / sqrt(3),  1.0 / sqrt(3), 1.0, 1.0),
    ]

    for xi, eta, wi, wj in gauss_rule:

        # print(f"xi = {xi} eta = {eta}")

        x = x0*(1 - eta)*(1 - xi)/4 + x1*(1 - eta)*(xi + 1)/4 + x2*(eta + 1)*(xi + 1)/4 + x3*(1 - xi)*(eta + 1)/4
        y = y0*(1 - eta)*(1 - xi)/4 + y1*(1 - eta)*(xi + 1)/4 + y2*(eta + 1)*(xi + 1)/4 + y3*(1 - xi)*(eta + 1)/4
        dx_dxi = -x0*(1 - eta)/4 + x1*(1 - eta)/4 + x2*(eta + 1)/4 - x3*(eta + 1)/4
        dx_deta = -x0*(1 - xi)/4 - x1*(xi + 1)/4 + x2*(xi + 1)/4 + x3*(1 - xi)/4
        dy_dxi = -y0*(1 - eta)/4 + y1*(1 - eta)/4 + y2*(eta + 1)/4 - y3*(eta + 1)/4
        dy_deta = -y0*(1 - xi)/4 - y1*(xi + 1)/4 + y2*(xi + 1)/4 + y3*(1 - xi)/4
        
        dN0_dxi = eta/4. - 1/4.
        dN0_deta = xi/4. - 1/4.
        dN1_dxi = 1/4. - eta/4.
        dN1_deta = -xi/4. - 1/4.
        dN2_dxi = eta/4. + 1/4.
        dN2_deta = xi/4. + 1/4.
        dN3_dxi = -eta/4. - 1/4.
        dN3_deta = 1/4. - xi/4.

        # print(f"x = {x} y = {y}")

        J = array([
        [dx_dxi, dx_deta],
        [dy_dxi, dy_deta]
        ]).T

        detJ = det(J)

        if detJ <= 0.:
            print(f"FATAL! detJ <= 0...")
            exit(-1)

        Jinv = inv(J)

        # print(f"J = {J}")
        # print(f"detJ = {detJ}")

        dN0_dxy = Jinv@array([ dN0_dxi, dN0_deta ])
        dN1_dxy = Jinv@array([ dN1_dxi, dN1_deta ])
        dN2_dxy = Jinv@array([ dN2_dxi, dN2_deta ])
        dN3_dxy = Jinv@array([ dN3_dxi, dN3_deta ])

        # ε = B ue
        B = zeros((3, 8))
        B[0,0] = dN0_dxy[0]
        B[1,1] = dN0_dxy[1]
        B[2,0] = dN0_dxy[1]
        B[2,1] = dN0_dxy[0]
        B[0,2] = dN1_dxy[0]
        B[1,3] = dN1_dxy[1]
        B[2,2] = dN1_dxy[1]
        B[2,3] = dN1_dxy[0]
        B[0,4] = dN2_dxy[0]
        B[1,5] = dN2_dxy[1]
        B[2,4] = dN2_dxy[1]
        B[2,5] = dN2_dxy[0]
        B[0,6] = dN3_dxy[0]
        B[1,7] = dN3_dxy[1]
        B[2,6] = dN3_dxy[1]
        B[2,7] = dN3_dxy[0]

        # print(f"B = {B}")

        x0 = xy[0,0]
        x1 = xy[1,0]
        x2 = xy[2,0]
        x3 = xy[3,0]

        y0 = xy[0,1]
        y1 = xy[1,1]
        y2 = xy[2,1]
        y3 = xy[3,1]

        l0 = sqrt((x1-x2)**2 + (y1+y2)**2)
        l1 = sqrt((x0-x2)**2 + (y0+y2)**2)
        l2 = sqrt((x1-x0)**2 + (y1+y0)**2)
        l3 = sqrt((x3-x0)**2 + (y3+y0)**2)
        sp = (l0+l1+l2+l3) /2
        Ae = sqrt((sp-l0)*(sp-l1)*(sp-l2)*(sp-l3)*sp)

        ke += t * wi * wj * B.T @ Eσ @ B * detJ
        fe += t * wi * wj * Ae * (t/3) * by * detJ
        

    return ke, fe


def quad4_post(xy, u_e, properties):

    E = properties["E"]
    ν = properties["nu"]
    bx = properties["bx"]
    by = properties["by"]
    t = properties["t"]

    #Podemos pasarle otros valores de xi y eta
    if "xi" in properties:
        xi = properties["xi"]
    else:
        xi = 0.0

    if "eta" in properties:
        eta = properties["eta"]
    else:
        eta = 0.0

    Eσ = E / (1-ν**2) * array(
        [
        [1 , ν , 0       ]       ,
        [ν , 1 , 0       ]       ,
        [0 , 0 , (1-ν)/2 ]
        ])

    x0 = xy[0,0]
    x1 = xy[1,0]
    x2 = xy[2,0]
    x3 = xy[3,0]

    y0 = xy[0,1]
    y1 = xy[1,1]
    y2 = xy[2,1]
    y3 = xy[3,1]

    x = x0*(1 - eta)*(1 - xi)/4 + x1*(1 - eta)*(xi + 1)/4 + x2*(eta + 1)*(xi + 1)/4 + x3*(1 - xi)*(eta + 1)/4
    y = y0*(1 - eta)*(1 - xi)/4 + y1*(1 - eta)*(xi + 1)/4 + y2*(eta + 1)*(xi + 1)/4 + y3*(1 - xi)*(eta + 1)/4
    dx_dxi = -x0*(1 - eta)/4 + x1*(1 - eta)/4 + x2*(eta + 1)/4 - x3*(eta + 1)/4
    dx_deta = -x0*(1 - xi)/4 - x1*(xi + 1)/4 + x2*(xi + 1)/4 + x3*(1 - xi)/4
    dy_dxi = -y0*(1 - eta)/4 + y1*(1 - eta)/4 + y2*(eta + 1)/4 - y3*(eta + 1)/4
    dy_deta = -y0*(1 - xi)/4 - y1*(xi + 1)/4 + y2*(xi + 1)/4 + y3*(1 - xi)/4
    
    dN0_dxi = eta/4. - 1/4.
    dN0_deta = xi/4. - 1/4.
    dN1_dxi = 1/4. - eta/4.
    dN1_deta = -xi/4. - 1/4.
    dN2_dxi = eta/4. + 1/4.
    dN2_deta = xi/4. + 1/4.
    dN3_dxi = -eta/4. - 1/4.
    dN3_deta = 1/4. - xi/4.

    # print(f"x = {x} y = {y}")

    J = array([
    [dx_dxi, dx_deta],
    [dy_dxi, dy_deta]
    ]).T

    detJ = det(J)

    if detJ <= 0.:
        print(f"FATAL! detJ <= 0...")
        exit(-1)

    Jinv = inv(J)

    # print(f"J = {J}")
    # print(f"detJ = {detJ}")

    dN0_dxy = Jinv@array([ dN0_dxi, dN0_deta ])
    dN1_dxy = Jinv@array([ dN1_dxi, dN1_deta ])
    dN2_dxy = Jinv@array([ dN2_dxi, dN2_deta ])
    dN3_dxy = Jinv@array([ dN3_dxi, dN3_deta ])

    # ε = B ue
    B = zeros((3, 8))
    B[0,0] = dN0_dxy[0]
    B[1,1] = dN0_dxy[1]
    B[2,0] = dN0_dxy[1]
    B[2,1] = dN0_dxy[0]
    B[0,2] = dN1_dxy[0]
    B[1,3] = dN1_dxy[1]
    B[2,2] = dN1_dxy[1]
    B[2,3] = dN1_dxy[0]
    B[0,4] = dN2_dxy[0]
    B[1,5] = dN2_dxy[1]
    B[2,4] = dN2_dxy[1]
    B[2,5] = dN2_dxy[0]
    B[0,6] = dN3_dxy[0]
    B[1,7] = dN3_dxy[1]
    B[2,6] = dN3_dxy[1]
    B[2,7] = dN3_dxy[0]


    ε = B @ u_e
    σ = Eσ @ ε

    return ε, σ


def quad4_nodal_stress_averaging(xy, u_e, properties, W):
    E = properties["E"]
    ν = properties["nu"]
    bx = properties["bx"]
    by = properties["by"]
    t = properties["t"]

    #Podemos pasarle otros valores de xi y eta
    if "xi" in properties:
        xi = properties["xi"]
    else:
        xi = 0.0

    if "eta" in properties:
        eta = properties["eta"]
    else:
        eta = 0.0

    Eσ = E / (1-ν**2) * array(
        [
        [1 , ν , 0       ]       ,
        [ν , 1 , 0       ]       ,
        [0 , 0 , (1-ν)/2 ]
        ])

    x0 = xy[0,0]
    x1 = xy[1,0]
    x2 = xy[2,0]
    x3 = xy[3,0]

    y0 = xy[0,1]
    y1 = xy[1,1]
    y2 = xy[2,1]
    y3 = xy[3,1]


    gauss_rule = [
        (-1.0 / sqrt(3), -1.0 / sqrt(3), -1.0, -1.0),
        ( 1.0 / sqrt(3), -1.0 / sqrt(3), +1.0, -1.0),
        ( 1.0 / sqrt(3),  1.0 / sqrt(3), +1.0, +1.0),
        (-1.0 / sqrt(3),  1.0 / sqrt(3), -1.0, +1.0),]

    for xi, eta, wi, wj in gauss_rule:

        chi = xi*(3**0.5)
        n = eta*(3**0.5)

        N0 = (1/4)*(1-chi)*(1-n)
        N1 = (1/4)*(1+chi)*(1-n)
        N2 = (1/4)*(1+chi)*(1+n)
        N3 = (1/4)*(1-chi)*(1+n)

        N = zeros((4, 4))
        N[0,0] = N2
        N[0,1] = N3
        N[0,2] = N0
        N[0,3] = N1
        N[1,0] = N1
        N[1,1] = N2
        N[1,2] = N3
        N[1,3] = N0
        N[2,0] = N0
        N[2,1] = N1
        N[2,2] = N2
        N[2,3] = N3
        N[3,0] = N3
        N[3,1] = N0
        N[3,2] = N1
        N[3,3] = N2

        w = N @ W

    return w


def quad4_line_load(xy,properties_load,ni,nj):
    t = properties_load["t"]
    tx = properties_load["tx"]
    ty = properties_load["ty"]

    x0 = xy[0,0]
    x1 = xy[1,0]
    y0 = xy[0,1]
    y1 = xy[1,1]

    xi = -1.0 / sqrt(3)
    eta =  -1.0 / sqrt(3)

    N0 = 1-xi/2
    N1 = 1+xi/2

    N = zeros((2, 4))
    N[0,0] = 1
    N[0,1] = 0
    N[0,2] = 1
    N[0,3] = 0
    N[1,0] = 0
    N[1,1] = 1
    N[1,2] = 0
    N[1,3] = 1

    # dN1 en dchi es 1 y para dN0 es 1      

    x = N0*x0 + N1*x1
    y = N0*y0 + N1*y1

    dN0_dxi = -1/2
    dN1_dxi = 1/2

    dx = dN0_dxi* x0 + dN1_dxi* x1
    dy = dN0_dxi* y0 + dN1_dxi* y1

    d2_cuadrado = dx**2 + dy**2

    Fe = t*(dx**2 + dy**2)* N.T * [tx,ty]

    fe = [Fe[0,0],Fe[1,0],Fe[2,0],Fe[3,0]]

    return fe


#xy = array([
# [-1,-1],
# [1,-1],
# [1,1],
# [-1,1],
#   ])

'''
xy = array([
[0,0],
[1,0],
[1,1],
[0,1],  ])

properties = {}
properties["E"] = 1.
properties["nu"] = 0.25
properties["bx"] = 0
properties["by"] = 1.
properties["t"] = 1.

ke, fe = quad4(xy, properties)

print(f"ke = {ke}")
print(f"fe = {fe}")

fixed_dofs = [0, 1, 2, 3]
free_dofs = [4, 5, 6, 7]

ke_ff = ke[ix_(free_dofs, free_dofs)]
fe_ff = array([0, -1, 0, -1])

#print(f"ke_ff = {ke_ff}")

from scipy.linalg import solve

u = zeros((8,1))
uf = solve(ke_ff, fe_ff)

print(f"uf = {uf}")
'''