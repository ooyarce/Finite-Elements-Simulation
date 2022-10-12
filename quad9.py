from numpy import array, sqrt, zeros, ix_
from scipy.linalg import det, inv
from sys import exit
import pprint

def quad9(xy, properties):

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
    x4 = xy[4,0]
    x5 = xy[5,0]
    x6 = xy[6,0]
    x7 = xy[7,0]
    x8 = xy[8,0]
    
    y0 = xy[0,1]
    y1 = xy[1,1]
    y2 = xy[2,1]
    y3 = xy[3,1]
    y4 = xy[4,1]
    y5 = xy[5,1]
    y6 = xy[6,1]
    y7 = xy[7,1]
    y8 = xy[8,1]

    ke = zeros((18,18))
    fe = zeros((18,1))

    #Primer punto de Gauss de la regla 2x2
    # xi = 1.0 / sqrt(3)
    # eta = -1.0 / sqrt(3)
    # wi = 1.0
    # wj = 1.0

    """
    gauss_rule = [
        (-1.0 / sqrt(3), -1.0 / sqrt(3), 1.0, 1.0),
        ( 1.0 / sqrt(3), -1.0 / sqrt(3), 1.0, 1.0),
        ( 1.0 / sqrt(3),  1.0 / sqrt(3), 1.0, 1.0),
        (-1.0 / sqrt(3),  1.0 / sqrt(3), 1.0, 1.0),
    ]
    """
    gauss_rule = [
        (-sqrt(0.6), -sqrt(0.6), 25/81 , 25/81 ),
        ( sqrt(0.6), -sqrt(0.6), 25/81 , 25/81 ),
        ( sqrt(0.6),  sqrt(0.6), 25/81 , 25/81 ),
        (-sqrt(0.6),  sqrt(0.6), 25/81 , 25/81 ),
        (         0,          0, 64/81 , 64/81 ),
        (         0, -sqrt(0.6), 40/81 , 40/81 ),
        ( sqrt(0.6),          0, 40/81 , 40/81 ),
        (         0,  sqrt(0.6), 40/81 , 40/81 ),
        (-sqrt(0.6),          0, 40/81 , 40/81 ),
    ]
    
    A = 0.0

    for xi, eta, wi, wj in gauss_rule:

        # print(f"xi = {xi} eta = {eta}")

        x = eta*x0*xi*(eta - 1)*(xi - 1)/4 + eta*x1*xi*(eta - 1)*(xi + 1)/4 + eta*x2*xi*(eta + 1)*(xi + 1)/4 + eta*x3*xi*(eta + 1)*(xi - 1)/4 + eta*x4*(1 - xi**2)*(eta - 1)/2 + eta*x6*(1 - xi**2)*(eta + 1)/2 + x5*xi*(1 - eta**2)*(xi + 1)/2 + x7*xi*(1 - eta**2)*(xi - 1)/2 + x8*(1 - eta**2)*(1 - xi**2)
        y = eta*xi*y0*(eta - 1)*(xi - 1)/4 + eta*xi*y1*(eta - 1)*(xi + 1)/4 + eta*xi*y2*(eta + 1)*(xi + 1)/4 + eta*xi*y3*(eta + 1)*(xi - 1)/4 + eta*y4*(1 - xi**2)*(eta - 1)/2 + eta*y6*(1 - xi**2)*(eta + 1)/2 + xi*y5*(1 - eta**2)*(xi + 1)/2 + xi*y7*(1 - eta**2)*(xi - 1)/2 + y8*(1 - eta**2)*(1 - xi**2)
        
        dx_dxi = eta*x0*xi*(eta - 1)/4 + eta*x0*(eta - 1)*(xi - 1)/4 + eta*x1*xi*(eta - 1)/4 + eta*x1*(eta - 1)*(xi + 1)/4 + eta*x2*xi*(eta + 1)/4 + eta*x2*(eta + 1)*(xi + 1)/4 + eta*x3*xi*(eta + 1)/4 + eta*x3*(eta + 1)*(xi - 1)/4 - eta*x4*xi*(eta - 1) - eta*x6*xi*(eta + 1) + x5*xi*(1 - eta**2)/2 + x5*(1 - eta**2)*(xi + 1)/2 + x7*xi*(1 - eta**2)/2 + x7*(1 - eta**2)*(xi - 1)/2 - 2*x8*xi*(1 - eta**2)
        dx_deta = eta*x0*xi*(xi - 1)/4 + eta*x1*xi*(xi + 1)/4 + eta*x2*xi*(xi + 1)/4 + eta*x3*xi*(xi - 1)/4 + eta*x4*(1 - xi**2)/2 - eta*x5*xi*(xi + 1) + eta*x6*(1 - xi**2)/2 - eta*x7*xi*(xi - 1) - 2*eta*x8*(1 - xi**2) + x0*xi*(eta - 1)*(xi - 1)/4 + x1*xi*(eta - 1)*(xi + 1)/4 + x2*xi*(eta + 1)*(xi + 1)/4 + x3*xi*(eta + 1)*(xi - 1)/4 + x4*(1 - xi**2)*(eta - 1)/2 + x6*(1 - xi**2)*(eta + 1)/2
        dy_dxi = eta*xi*y0*(eta - 1)/4 + eta*xi*y1*(eta - 1)/4 + eta*xi*y2*(eta + 1)/4 + eta*xi*y3*(eta + 1)/4 - eta*xi*y4*(eta - 1) - eta*xi*y6*(eta + 1) + eta*y0*(eta - 1)*(xi - 1)/4 + eta*y1*(eta - 1)*(xi + 1)/4 + eta*y2*(eta + 1)*(xi + 1)/4 + eta*y3*(eta + 1)*(xi - 1)/4 + xi*y5*(1 - eta**2)/2 + xi*y7*(1 - eta**2)/2 - 2*xi*y8*(1 - eta**2) + y5*(1 - eta**2)*(xi + 1)/2 + y7*(1 - eta**2)*(xi - 1)/2
        dy_deta = eta*xi*y0*(xi - 1)/4 + eta*xi*y1*(xi + 1)/4 + eta*xi*y2*(xi + 1)/4 + eta*xi*y3*(xi - 1)/4 - eta*xi*y5*(xi + 1) - eta*xi*y7*(xi - 1) + eta*y4*(1 - xi**2)/2 + eta*y6*(1 - xi**2)/2 - 2*eta*y8*(1 - xi**2) + xi*y0*(eta - 1)*(xi - 1)/4 + xi*y1*(eta - 1)*(xi + 1)/4 + xi*y2*(eta + 1)*(xi + 1)/4 + xi*y3*(eta + 1)*(xi - 1)/4 + y4*(1 - xi**2)*(eta - 1)/2 + y6*(1 - xi**2)*(eta + 1)/2
        
        
        dN0_dxi = eta*xi*(eta - 1)/4 + eta*(eta - 1)*(xi - 1)/4
        dN1_dxi = eta*xi*(eta - 1)/4 + eta*(eta - 1)*(xi + 1)/4
        dN2_dxi = eta*xi*(eta + 1)/4 + eta*(eta + 1)*(xi + 1)/4
        dN3_dxi = eta*xi*(eta + 1)/4 + eta*(eta + 1)*(xi - 1)/4
        dN4_dxi = -eta*xi*(eta - 1)
        dN5_dxi = xi*(1 - eta**2)/2 + (1 - eta**2)*(xi/2 + 1/2)
        dN6_dxi = -eta*xi*(eta + 1)
        dN7_dxi = xi*(1 - eta**2)/2 + (1/2 - eta**2/2)*(xi - 1)
        dN8_dxi = -2*xi*(1 - eta**2)
        
        dN0_deta = eta*xi*(xi - 1)/4 + xi*(eta - 1)*(xi - 1)/4
        dN1_deta = eta*xi*(xi + 1)/4 + xi*(eta - 1)*(xi + 1)/4
        dN2_deta = eta*xi*(xi + 1)/4 + xi*(eta + 1)*(xi + 1)/4
        dN3_deta = eta*xi*(xi - 1)/4 + xi*(eta + 1)*(xi - 1)/4
        dN4_deta = eta*(1 - xi**2)/2 + (1/2 - xi**2/2)*(eta - 1)
        dN5_deta = -eta*xi*(xi + 1)
        dN6_deta = eta*(1 - xi**2)/2 + (1 - xi**2)*(eta/2 + 1/2)
        dN7_deta = -eta*xi*(xi - 1)
        dN8_deta = -2*eta*(1 - xi**2)
        

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
        dN4_dxy = Jinv@array([ dN4_dxi, dN4_deta ])
        dN5_dxy = Jinv@array([ dN5_dxi, dN5_deta ])
        dN6_dxy = Jinv@array([ dN6_dxi, dN6_deta ])
        dN7_dxy = Jinv@array([ dN7_dxi, dN7_deta ])
        dN8_dxy = Jinv@array([ dN8_dxi, dN8_deta ])
        

        # ε = B ue
        B = zeros((3, 18))
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
        
        B[0,8] = dN4_dxy[0]
        B[1,9] = dN4_dxy[1]
        B[2,8] = dN4_dxy[1]
        B[2,9] = dN4_dxy[0]
        B[0,10] = dN5_dxy[0]
        B[1,11] = dN5_dxy[1]
        B[2,10] = dN5_dxy[1]
        B[2,11] = dN5_dxy[0]
        B[0,12] = dN6_dxy[0]
        B[1,13] = dN6_dxy[1]
        B[2,12] = dN6_dxy[1]
        B[2,13] = dN6_dxy[0]
        B[0,14] = dN7_dxy[0]
        B[1,15] = dN7_dxy[1]
        B[2,14] = dN7_dxy[1]
        B[2,15] = dN7_dxy[0]
        B[0,16] = dN8_dxy[0]
        B[1,17] = dN8_dxy[1]
        B[2,16] = dN8_dxy[1]
        B[2,17] = dN8_dxy[0]
        

        # print(f"B = {B}")


        ke += t * wi * wj * B.T @ Eσ @ B * detJ
        
        A += wi * wj * detJ

    return ke, fe, A


def quad9_post(xy, u_e, properties):

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
    x4 = xy[4,0]
    x5 = xy[5,0]
    x6 = xy[6,0]
    x7 = xy[7,0]
    x8 = xy[8,0]
    
    y0 = xy[0,1]
    y1 = xy[1,1]
    y2 = xy[2,1]
    y3 = xy[3,1]
    y4 = xy[4,1]
    y5 = xy[5,1]
    y6 = xy[6,1]
    y7 = xy[7,1]
    y8 = xy[8,1]
    
    x = eta*x0*xi*(eta - 1)*(xi - 1)/4 + eta*x1*xi*(eta - 1)*(xi + 1)/4 + eta*x2*xi*(eta + 1)*(xi + 1)/4 + eta*x3*xi*(eta + 1)*(xi - 1)/4 + eta*x4*(1 - xi**2)*(eta - 1)/2 + eta*x6*(1 - xi**2)*(eta + 1)/2 + x5*xi*(1 - eta**2)*(xi + 1)/2 + x7*xi*(1 - eta**2)*(xi - 1)/2 + x8*(1 - eta**2)*(1 - xi**2)
    y = eta*xi*y0*(eta - 1)*(xi - 1)/4 + eta*xi*y1*(eta - 1)*(xi + 1)/4 + eta*xi*y2*(eta + 1)*(xi + 1)/4 + eta*xi*y3*(eta + 1)*(xi - 1)/4 + eta*y4*(1 - xi**2)*(eta - 1)/2 + eta*y6*(1 - xi**2)*(eta + 1)/2 + xi*y5*(1 - eta**2)*(xi + 1)/2 + xi*y7*(1 - eta**2)*(xi - 1)/2 + y8*(1 - eta**2)*(1 - xi**2)
    
    dx_dxi = eta*x0*xi*(eta - 1)/4 + eta*x0*(eta - 1)*(xi - 1)/4 + eta*x1*xi*(eta - 1)/4 + eta*x1*(eta - 1)*(xi + 1)/4 + eta*x2*xi*(eta + 1)/4 + eta*x2*(eta + 1)*(xi + 1)/4 + eta*x3*xi*(eta + 1)/4 + eta*x3*(eta + 1)*(xi - 1)/4 - eta*x4*xi*(eta - 1) - eta*x6*xi*(eta + 1) + x5*xi*(1 - eta**2)/2 + x5*(1 - eta**2)*(xi + 1)/2 + x7*xi*(1 - eta**2)/2 + x7*(1 - eta**2)*(xi - 1)/2 - 2*x8*xi*(1 - eta**2)
    dx_deta = eta*x0*xi*(xi - 1)/4 + eta*x1*xi*(xi + 1)/4 + eta*x2*xi*(xi + 1)/4 + eta*x3*xi*(xi - 1)/4 + eta*x4*(1 - xi**2)/2 - eta*x5*xi*(xi + 1) + eta*x6*(1 - xi**2)/2 - eta*x7*xi*(xi - 1) - 2*eta*x8*(1 - xi**2) + x0*xi*(eta - 1)*(xi - 1)/4 + x1*xi*(eta - 1)*(xi + 1)/4 + x2*xi*(eta + 1)*(xi + 1)/4 + x3*xi*(eta + 1)*(xi - 1)/4 + x4*(1 - xi**2)*(eta - 1)/2 + x6*(1 - xi**2)*(eta + 1)/2
    dy_dxi = eta*xi*y0*(eta - 1)/4 + eta*xi*y1*(eta - 1)/4 + eta*xi*y2*(eta + 1)/4 + eta*xi*y3*(eta + 1)/4 - eta*xi*y4*(eta - 1) - eta*xi*y6*(eta + 1) + eta*y0*(eta - 1)*(xi - 1)/4 + eta*y1*(eta - 1)*(xi + 1)/4 + eta*y2*(eta + 1)*(xi + 1)/4 + eta*y3*(eta + 1)*(xi - 1)/4 + xi*y5*(1 - eta**2)/2 + xi*y7*(1 - eta**2)/2 - 2*xi*y8*(1 - eta**2) + y5*(1 - eta**2)*(xi + 1)/2 + y7*(1 - eta**2)*(xi - 1)/2
    dy_deta = eta*xi*y0*(xi - 1)/4 + eta*xi*y1*(xi + 1)/4 + eta*xi*y2*(xi + 1)/4 + eta*xi*y3*(xi - 1)/4 - eta*xi*y5*(xi + 1) - eta*xi*y7*(xi - 1) + eta*y4*(1 - xi**2)/2 + eta*y6*(1 - xi**2)/2 - 2*eta*y8*(1 - xi**2) + xi*y0*(eta - 1)*(xi - 1)/4 + xi*y1*(eta - 1)*(xi + 1)/4 + xi*y2*(eta + 1)*(xi + 1)/4 + xi*y3*(eta + 1)*(xi - 1)/4 + y4*(1 - xi**2)*(eta - 1)/2 + y6*(1 - xi**2)*(eta + 1)/2
    
    dN0_dxi = eta*xi*(eta - 1)/4 + eta*(eta - 1)*(xi - 1)/4
    dN1_dxi = eta*xi*(eta - 1)/4 + eta*(eta - 1)*(xi + 1)/4
    dN2_dxi = eta*xi*(eta + 1)/4 + eta*(eta + 1)*(xi + 1)/4
    dN3_dxi = eta*xi*(eta + 1)/4 + eta*(eta + 1)*(xi - 1)/4
    dN4_dxi = -eta*xi*(eta - 1)
    dN5_dxi = xi*(1 - eta**2)/2 + (1 - eta**2)*(xi/2 + 1/2)
    dN6_dxi = -eta*xi*(eta + 1)
    dN7_dxi = xi*(1 - eta**2)/2 + (1/2 - eta**2/2)*(xi - 1)
    dN8_dxi = -2*xi*(1 - eta**2)
    
    dN0_deta = eta*xi*(xi - 1)/4 + xi*(eta - 1)*(xi - 1)/4
    dN1_deta = eta*xi*(xi + 1)/4 + xi*(eta - 1)*(xi + 1)/4
    dN2_deta = eta*xi*(xi + 1)/4 + xi*(eta + 1)*(xi + 1)/4
    dN3_deta = eta*xi*(xi - 1)/4 + xi*(eta + 1)*(xi - 1)/4
    dN4_deta = eta*(1 - xi**2)/2 + (1/2 - xi**2/2)*(eta - 1)
    dN5_deta = -eta*xi*(xi + 1)
    dN6_deta = eta*(1 - xi**2)/2 + (1 - xi**2)*(eta/2 + 1/2)
    dN7_deta = -eta*xi*(xi - 1)
    dN8_deta = -2*eta*(1 - xi**2)

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
    dN4_dxy = Jinv@array([ dN4_dxi, dN4_deta ])
    dN5_dxy = Jinv@array([ dN5_dxi, dN5_deta ])
    dN6_dxy = Jinv@array([ dN6_dxi, dN6_deta ])
    dN7_dxy = Jinv@array([ dN7_dxi, dN7_deta ])
    dN8_dxy = Jinv@array([ dN8_dxi, dN8_deta ])

    # ε = B ue
    B = zeros((3, 18))
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
    
    B[0,8] = dN4_dxy[0]
    B[1,9] = dN4_dxy[1]
    B[2,8] = dN4_dxy[1]
    B[2,9] = dN4_dxy[0]
    B[0,10] = dN5_dxy[0]
    B[1,11] = dN5_dxy[1]
    B[2,10] = dN5_dxy[1]
    B[2,11] = dN5_dxy[0]
    B[0,12] = dN6_dxy[0]
    B[1,13] = dN6_dxy[1]
    B[2,12] = dN6_dxy[1]
    B[2,13] = dN6_dxy[0]
    B[0,14] = dN7_dxy[0]
    B[1,15] = dN7_dxy[1]
    B[2,14] = dN7_dxy[1]
    B[2,15] = dN7_dxy[0]
    B[0,16] = dN8_dxy[0]
    B[1,17] = dN8_dxy[1]
    B[2,16] = dN8_dxy[1]
    B[2,17] = dN8_dxy[0]


    ε = B @ u_e
    σ = Eσ @ ε

    return ε, σ

"""
xy = array([
[-1,-1],
[1,-1],
[1,1],
[-1,1],
    ])

# xy = array([
# [0,0],
# [1,0],
# [1,1],
# [0,1],
#     ])

properties = {}
properties["E"] = 1.
properties["nu"] = 0.25
properties["bx"] = 0
properties["by"] = 1.
properties["t"] = 1.

ke, fe = quad9(xy, properties)

print(f"ke = {ke}")

fixed_dofs = [0, 1, 2, 3]
free_dofs = [4, 5, 6, 7]

ke_ff = ke[ix_(free_dofs, free_dofs)]
fe_ff = array([0, -1, 0, -1])

print(f"ke_ff = {ke_ff}")

from scipy.linalg import solve

u = zeros((8,1))
uf = solve(ke_ff, fe_ff)

print(f"uf = {uf}")
"""