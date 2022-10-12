#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 20:23:22 2021

@author: Tomy
"""
import matplotlib.pylab as plt
from gmesh_post import write_node_data,write_node_data_2,write_element_data, write_element_node_data
from quad9 import *
import numpy as np
def optimizar(espesores,iteration):
    print(f"----------------------------------------------------------------------------------")
    print(f"----------------------------------------------------------------------------------")
    print(f"------------------------INICIO ITERACION {iteration}------------------------------\n")
    infill = 0.98
    densidad_PLA = 1.240*infill #gr/cm3
    
    g = 9.8
    
    fid = open("malla_con_elipse_5_OFICIAL2.msh","r")
    
    LINE_ELEMENT = 8
    TRIANGLE_ELEMENT = 9
    QUAD_ELEMENT = 10
    
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
        
    
    
    while True:
        line = fid.readline()
        
        if line.find("$Elements") >= 0:
            break
    
    Nelem = int(fid.readline())
    
    conect = np.zeros((Nelem,9), dtype = np.int32)
    
    fixed_nodes = []
    load_nodes = []
    
    Nquads = 0
    Quadrangles = []
    tipo = []
    
    for i in range(Nelem):
        line = fid.readline()
        sl = line.split()
        element_number = np.int32(sl[0]) -1
        element_type = np.int32(sl[1])
        physical_grp = np.int32(sl[3])
        tipo.append(physical_grp)
        entity_number = np.int32(sl[4])
    
        if element_type == LINE_ELEMENT and physical_grp == Empotrado:
            n1 = np.int32(sl[5]) -1
            n2 = np.int32(sl[6]) -1
            n3 = np.int32(sl[7]) -1
            
            fixed_nodes += [n1,n2,n3]
            
        
        if (element_type == LINE_ELEMENT or element_type == QUAD_ELEMENT) and physical_grp == BordeNatural:
            n1 = np.int32(sl[5]) -1
            n2 = np.int32(sl[6]) -1
            n3 = np.int32(sl[7]) -1
            
            load_nodes += [n1,n2]
            load_nodes += [n2,n3]
    
        if element_type == QUAD_ELEMENT and (physical_grp == Placa or physical_grp == Extremos):
            n0 = np.int32(sl[5]) -1
            n1 = np.int32(sl[6]) -1
            n2 = np.int32(sl[7]) -1
            n3 = np.int32(sl[8]) -1
            n4 = np.int32(sl[9]) -1
            n5 = np.int32(sl[10]) -1
            n6 = np.int32(sl[11]) -1
            n7 = np.int32(sl[12]) -1
            n8 = np.int32(sl[13]) -1
            
    
            conect[element_number, :] = [ n0 , n1 , n2 , n3, n4, n5, n6, n7, n8] 
            
            Quadrangles.append(element_number)
            Nquads += 1
    
    
    load_nodes.sort()
    final_load_nodes = []
    for i in load_nodes:
        if i not in final_load_nodes:
            final_load_nodes.append(i)
    
    
    h0 = 50e-2
    alpha = 1
    properties_0 = {}
    properties_0["E"] = 20e9*alpha #Pa REVISAR
    properties_0["nu"] = 0.25
    properties_0["bx"] = 0
    properties_0["t"] = 4e-3
    properties_0["by"] = -densidad_PLA * g
    #properties_0["by"] = -densidad_hormigon * h0 * g *properties_0["t"]
    
    Ndofs_per_node = 2
    Ndofs = Ndofs_per_node * Nnodes
    
    K = np.zeros((Ndofs,Ndofs))
    f = np.zeros((Ndofs,1))
    
    peso = 0
    A0 = 1.
    indice = 0
    for e in Quadrangles:
        ni = int(conect[e,0])
        nj = int(conect[e,1])
        nk = int(conect[e,2])
        nl = int(conect[e,3])
        nm = int(conect[e,4])
        nn = int(conect[e,5])
        no = int(conect[e,6])
        npe = int(conect[e,7])
        nq = int(conect[e,8])
        
        xy_e = xy[[ ni , nj , nk , nl , nm , nn , no , npe , nq ],:] #reescribiendo xy
        
        if tipo[e] == Placa:
            espesor = espesores[indice]
        
        else:
            espesor = 5.2e-1
        
        properties_0["t"] = espesor
        
        ke , fe , Ae = quad9( xy_e , properties_0)
        
        peso += Ae * espesor * densidad_PLA
        d = [ 2*ni , 2*ni+1 , 2*nj , 2*nj+1 , 2*nk , 2*nk+1 , 2*nl , 2*nl+1 , 2*nm , 2*nm+1 , 2*nn , 2*nn+1 , 2*no , 2*no+1 , 2*npe , 2*npe+1 , 2*nq , 2*nq+1 ]    
        #print (f"Elemento {e}: {d}")
        
        #DIRECT STIFFNESS METHOD
        for i in range(9*Ndofs_per_node):
            p = d[i]
            #print (f"p = {p}")
            for j in range(9*Ndofs_per_node):
                q = d[j]
                K[p,q] += ke[i,j]
            f[p] += fe[i]
    
        indice += 1
    
    
    fixed_nodes = np.unique(fixed_nodes)
    
    c_DOFs = []
    
    for n in fixed_nodes:
        c_DOFs += [ 2*n , 2*n+1 ]    
            
    free_DOFs = np.arange(Ndofs)
    free_DOFs = np.setdiff1d(free_DOFs, c_DOFs)
    
    
    q = 1000.0/(len(final_load_nodes)-1)
    i = 0
    for n in final_load_nodes:
        if i <=1:
            f[2*n] = 5.2e-3*q/8
            i +=1
        else:
            f[2*n] = properties_0["t"] *q/4
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
    
    #print (f"u={u}")
    #print (f"R={R}")
    
    factor = 2.5e5
    
    uv = u.reshape([-1,2])
    
    plt.plot(xy[:,0] + factor*uv[:,0] , xy[:,1] + factor*uv[:,1],".")
    
    xy2 = xy
    for e in Quadrangles:
        ni = int(conect[e,0])
        nj = int(conect[e,1])
        nk = int(conect[e,2])
        nl = int(conect[e,3])
        nm = int(conect[e,4])
        nn = int(conect[e,5])
        no = int(conect[e,6])
        npe = int(conect[e,7])
        nq = int(conect[e,8])
    
        xy_e = xy[[ ni, nm, nj, nn, nk, no, nl, npe, ni],:] + factor*uv[[ ni, nm, nj, nn, nk, no, nl, npe, ni],:] 
        
        plt.plot(xy_e[:,0] , xy_e[:,1], "k")
    
    
    plt.axis("equal")
    plt.show()        
    
    
    #ESCRITURA DE ARCHIVOS PARA QUE EL MESH VEA LAS DEFORMACIONES
    
    
    nodes = np.arange(1,Nnodes+1)
    write_node_data("ux.msh",nodes,uv[:,0],"Despl. X")
    write_node_data("uy.msh",nodes,uv[:,1],"Despl. Y")
    write_node_data_2("desplazamientos.msh",nodes,uv[:,0],uv[:,1],"Despl")
    
    #-----------------------------------------------------------------------------
    #CALCULO DE TENSIONES
    
    sigma_xx = np.zeros((Nquads+1))
    sigma_xy = np.zeros((Nquads+1))
    sigma_yy = np.zeros((Nquads+1))
    
    #sigma_xx = np.zeros(((Nquads+1,9)))
    #sigma_xy = np.zeros(((Nquads+1,9)))
    #sigma_yy = np.zeros(((Nquads+1,9)))
    
    i=0
    peso_it = 0
    for e in Quadrangles:
        ni = int(conect[e,0])
        nj = int(conect[e,1])
        nk = int(conect[e,2])
        nl = int(conect[e,3])
        nm = int(conect[e,4])
        nn = int(conect[e,5])
        no = int(conect[e,6])
        npe = int(conect[e,7])
        nq = int(conect[e,8])
        
        properties_0["t"] = espesores[i]
     
        xy_e = xy2[[ ni , nj , nk , nl, nm , nn , no , npe , nq , ni ],:]
        
        uv_e = uv[[ni,nj,nk,nl,nm,nn,no,npe,nq],:]
        u_e = uv_e.reshape((-1))
        epsilon,sigma = quad9_post(xy_e,u_e, properties_0)
        
        ke , fe , Ae = quad9( xy_e , properties_0)
        
        peso_it += Ae * properties_0["t"] * densidad_PLA      
        
        xy_e = xy2[[ ni , nj , nk , nl, nm , nn , no , npe , nq , ni ],:]
        
        uv_e = uv[[ni,nj,nk,nl,nm,nn,no,npe,nq],:]
        u_e = uv_e.reshape((-1))
        epsilon,sigma_it = quad9_post(xy_e,u_e, properties_0)
        sigma_xx[i] = sigma[0]
        
        i+=1    
    smax = max(sigma_xx)
    smin = min(sigma_xx)
    i=0
    for e in Quadrangles:
        properties_0["t"] = espesores[i]
        a = properties_0["t"]
        #print (f"elemento : {e} espesor original : {a}")
        
        #MODIFICACION DE ESPESORES
        if tipo[e] == Placa:
            
            if smin <= sigma_xx[i] < 0:
                properties_0["t"] = properties_0["t"] - 0.01
                espesores[i] = properties_0["t"]
                
            elif 0 <= sigma_xx[i] < -smin:
                properties_0["t"] = properties_0["t"] - 0.02
                espesores[i] = properties_0["t"]
                
            elif -smin <= sigma_xx[i] < -2*smin:
                properties_0["t"] = properties_0["t"] - 0.03
                espesores[i] = properties_0["t"]
                
            elif -smin*2 <= sigma_xx[i] < -smin*4:
                properties_0["t"] = properties_0["t"] + 0.005
                espesores[i] = properties_0["t"]
                
            elif -smin*4 <= sigma_xx[i] < smax/2:
                properties_0["t"] = properties_0["t"] + 0.008
                espesores[i] = properties_0["t"]
                
            elif smax/3 <= sigma_xx[i] < 2*smax/3:
                properties_0["t"] = properties_0["t"] + 0.01
                espesores[i] = properties_0["t"]
                
            elif 2*smax/3 <= sigma_xx[i] <= smax:
                properties_0["t"] = properties_0["t"] + 0.027
                espesores[i] = properties_0["t"]
                
        else:
            properties_0["t"] = 0.52
            espesores[i] = 0.52
        #print (f"elemento : {e} espesor nuevo : {espesores[i]}\n")
        i+=1
    peso_it = peso_it*2*10
    print(f"-----> peso_it {iteration} = {peso_it}gr")
    elementos = np.array(Quadrangles)+1
    write_element_data("sigma_x_OFICIAL_"+str(iteration)+".msh",elementos,sigma_xx,"Sigma"+str(iteration))
    write_element_data("espesores_OFICIAL_"+str(iteration)+".msh",elementos, espesores,"Espesores"+str(iteration))
    print(f"----------------------------------------------------------------------------------")
    print(f"----------------------------------------------------------------------------------")
    print(f"------------------------FIN ITERACION {iteration}---------------------------------\n")
    return espesores

from hola import espesores
"""
espesores = espesores()[0]
espesores = optimizar(espesores,1)
espesores = optimizar(espesores,2)
espesores = optimizar(espesores,3)
espesores = optimizar(espesores,4)
espesores = optimizar(espesores,5)
espesores = optimizar(espesores,6)
espesores = optimizar(espesores,7)
espesores = optimizar(espesores,8)
espesores = optimizar(espesores,9)
espesores = optimizar(espesores,10)
espesores = optimizar(espesores,11)
espesores = optimizar(espesores,12)
espesores = optimizar(espesores,13)
espesores = optimizar(espesores,14)
"""







        