#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 20:23:22 2021

@author: Tomy
"""


#from numpy import  zeros , ix_ , array , int
import matplotlib.pylab as plt
from quad4 import quad4
densidad_hormigon = 2500 #kg/m3

import numpy as np

g = 9.8

fid = open("malla_con_elipse_5_OFICIAL2.msh","r")

LINE_ELEMENT = 1
TRIANGLE_ELEMENT = 2
QUAD_ELEMENT = 3

Empotrado = 1
BordeNatural = 2
Placa = 3
Extremos = 4

while True:
    line = fid.readline()
    
    if line.find("$Nodes") >= 0:
        break

Nnodes = int(fid.readline())

xy = np.zeros([Nnodes,2])

for i in range(Nnodes):
    line = fid.readline()
    sl = line.split()
    xy[i,0] = float(sl[1])
    xy[i,1] = float(sl[2])
    
print(xy)

while True:
    line = fid.readline()
    
    if line.find("$Elements") >= 0:
        break

Nelem = int(fid.readline())

conect = np.zeros((Nelem,4), dtype = np.int32)

fixed_nodes = []

Nquads = 0
Quadrangles = []

for i in range(Nelem):
    line = fid.readline()
    sl = line.split()
    element_number = np.int32(sl[0]) -1
    element_type = np.int32(sl[1])
    physical_grp = np.int32(sl[3])
    entity_number = np.int32(sl[4])


    if element_type == LINE_ELEMENT and physical_grp == Empotrado:
        n1 = np.int32(sl[5]) -1
        n2 = np.int32(sl[6]) -1
        fixed_nodes += [n1,n2]
        print (n1)
        
    if element_type == QUAD_ELEMENT and (physical_grp == Placa or physical_grp == Extremos):
        n0 = np.int32(sl[5]) -1
        n1 = np.int32(sl[6]) -1
        n2 = np.int32(sl[7]) -1
        n3 = np.int32(sl[8]) -1
        print (n0,n1,n2)
        

        conect[element_number, :] = [ n0 , n1 , n2 , n3] 
        
        Quadrangles.append(element_number)
        Nquads += 1

print (conect)

h0 = 50e-2
properties_0 = {}
properties_0["E"] = 20e9 #Pa
properties_0["nu"] = 0.25
properties_0["bx"] = 0
properties_0["t"] = 4e-3
properties_0["by"] = -densidad_hormigon * g
#properties_0["by"] = -densidad_hormigon * h0 * g *properties_0["t"]

Ndofs_per_node = 2
Ndofs = Ndofs_per_node * Nnodes

K = np.zeros((Ndofs,Ndofs))
f = np.zeros((Ndofs,1))

for e in Quadrangles:
    ni = int(conect[e,0])
    nj = int(conect[e,1])
    nk = int(conect[e,2])
    nl = int(conect[e,3])
    
    xy_e = xy[[ ni , nj , nk , nl ],:] #reescribiendo xy
    ke , fe = quad4( xy_e , properties_0)
    d = [ 2*ni , 2*ni+1 , 2*nj , 2*nj+1 , 2*nk , 2*nk+1 , 2*nl , 2*nl+1 ]    
    print (f"Elemento {e}: {d}")
    
    #DIRECT STIFFNESS METHOD
    for i in range(4*Ndofs_per_node):
        p = d[i]
        #print (f"p = {p}")
        for j in range(4*Ndofs_per_node):
            q = d[j]
            K[p,q] += ke[i,j]
        f[p] += fe[i]
"""
plt.figure()
plt.matshow(K)
plt.colorbar()
plt.show()
"""

fixed_nodes = np.unique(fixed_nodes)

c_DOFs = []

for n in fixed_nodes:
    c_DOFs += [ 2*n , 2*n+1 ]    
        
free_DOFs = np.arange(Ndofs)
free_DOFs = np.setdiff1d(free_DOFs, c_DOFs)

print (free_DOFs,c_DOFs)

nodes = [3,4,11,12,13]
i = 0
for n in nodes:
    if i <=1:
        f[2*n] = 1000.0/8
        i +=1
    else:
        f[2*n] = 1000.0/4
        i +=1
Kff = K[np.ix_(free_DOFs, free_DOFs)]
Kfc = K[np.ix_(free_DOFs, c_DOFs)]
Kcf = K[np.ix_(c_DOFs, free_DOFs)]
Kcc = K[np.ix_(c_DOFs, c_DOFs)]

ff = f[free_DOFs]
fc = f[c_DOFs]



#SOLVE
from scipy.linalg import solve
u = np.zeros((Ndofs,1))
u[free_DOFs] = solve( Kff , ff )

#GET REACTION FORCES
R = Kcf @ u[free_DOFs] + Kcc @ u[c_DOFs] - fc 

print (f"u={u}")
print (f"R={R}")

factor = 1e4

uv = u.reshape([-1,2])

plt.plot(xy[:,0] + factor*uv[:,0] , xy[:,1] + factor*uv[:,1],".")

xy2 = xy
for e in Quadrangles:
    ni = int(conect[e,0])
    nj = int(conect[e,1])
    nk = int(conect[e,2])
    nl = int(conect[e,3])

    xy_e = xy[[ ni , nj, nk, nl, ni],:] + factor*uv[[ ni , nj, nk, nl, ni],:] 
    
    plt.plot(xy_e[:,0] , xy_e[:,1], "k")


plt.axis("equal")
plt.show()        
        
#ESCRITURA DE ARCHIVOS PARA QUE EL MESH VEA LAS DEFORMACIONES

from gmesh_post import write_node_data,write_node_data_2,write_element_data
from quad4 import *
nodes = np.arange(1,Nnodes+1)
write_node_data("ux.msh",nodes,uv[:,0],"Despl. X")
write_node_data("uy.msh",nodes,uv[:,1],"Despl. Y")
write_node_data_2("desplazamientos.msh",nodes,uv[:,0],uv[:,1],"Despl")

#-----------------------------------------------------------------------------
#CALCULO DE TENSIONES

sigma_xx = np.zeros((Nquads+1))
sigma_xy = np.zeros((Nquads+1))
sigma_yy = np.zeros((Nquads+1))
i=0
Sigma_max = []
for e in Quadrangles:
    ni = int(conect[e,0])
    nj = int(conect[e,1])
    nk = int(conect[e,2])
    nl = int(conect[e,3])
    xy_e = xy2[[ ni , nj , nk , nl, ni ],:]
    uv_e = uv[[ni,nj,nk,nl],:]
    u_e = uv_e.reshape((-1))
    epsilon,sigma = quad4_post(xy_e,u_e, properties_0)
    sigma_xx[i] = sigma[0] 
    sigma_yy[i] = sigma[1]
    sigma_xy[i] = sigma[2]
    i+=1
    smax = max([abs(min(sigma_xx)),abs(max(sigma_xx))],[abs(min(sigma_xy)),abs(max(sigma_xy))],[abs(min(sigma_yy)),abs(max(sigma_yy))])
    Sigma_max.append(smax)
    
print (Sigma_max)
elementos = np.array(Quadrangles)+1
write_element_data("sigma_x.msh",elementos,sigma_xx,"Sigma_x")





        