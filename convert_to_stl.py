import numpy as np

#Proposito: convertir una malla de cuadrilateros, con espesor constante por
# cada quad, en una malla STL para impresion 3D. 
#
# Argumentos:
#   name : string con el nombre del modelo
#   xy   : numpy-array con las coordenadas de los nodos
#   Quadrangles : numero (id) de cada elemento en gmsh, llamado 'e' abajo
#   conec : conectividad de cada elemento 'e'. 
#   conec[e,:] contiene las 4 esquinas del elemento e
#   nodos_borde : arreglo con los elementos linea que pertenecen al borde
#   nodos_borde[i,:] : contiene los nodos inicial y final de una linea
#   que pertenece a la frontera del modelo
#
# Esta funcion realiza el promedio nodal de los espesores antes de escribir
# el STL. 

from hola import espesores
from optimization import optimizar
name = "Placa_intento_3"

espesores,xy,Quadrangles,conect,nodos_borde,Nelem = espesores()

espesores = optimizar(espesores,1)
espesores = optimizar(espesores,2)
espesores = optimizar(espesores,3)
espesores = optimizar(espesores,4)
espesores = optimizar(espesores,5)
espesores = optimizar(espesores,6)
espesores = optimizar(espesores,7)
espesores = optimizar(espesores,8)



#Funcion para escribir 4 facetas (triangulos) dadas las esquinas de un
  #quad y su espesor el z. 
def write_facet_from_quad_coords(fid, x,y,z):
	fid.write(f"""facet normal 0. 0. 0.
  outer loop
    vertex {x[0]} {y[0]} {z[0]/2}
    vertex {x[1]} {y[1]} {z[1]/2}
    vertex {x[2]} {y[2]} {z[2]/2}
  endloop
endfacet
""")
	fid.write(f"""facet normal 0. 0. 0.
  outer loop
    vertex {x[0]} {y[0]} {z[0]/2}
    vertex {x[2]} {y[2]} {z[2]/2}
    vertex {x[3]} {y[3]} {z[3]/2}
  endloop
endfacet
""")
	fid.write(f"""facet normal 0. 0. 0.
  outer loop
    vertex {x[0]} {y[0]} {-z[0]/2}
    vertex {x[1]} {y[1]} {-z[1]/2}
    vertex {x[2]} {y[2]} {-z[2]/2}
  endloop
endfacet
""")
	fid.write(f"""facet normal 0. 0. 0.
  outer loop
    vertex {x[0]} {y[0]} {-z[0]/2}
    vertex {x[2]} {y[2]} {-z[2]/2}
    vertex {x[3]} {y[3]} {-z[3]/2}
  endloop
endfacet
""")

#Funcion para escribir 2 facetas (triangulos) dadas las esquinas de un
#elemento linea perteneciente a la frontera y su espesor el z
def write_facet_from_line_coords(fid, x,y,z):
	fid.write(f"""facet normal 0. 0. 0.
  outer loop
    vertex {x[0]} {y[0]} {z[0]/2}
    vertex {x[1]} {y[1]} {z[1]/2}
    vertex {x[0]} {y[0]} {-z[0]/2}
  endloop
endfacet
""")
	fid.write(f"""facet normal 0. 0. 0.
  outer loop
    vertex {x[0]} {y[0]} {-z[0]/2}
    vertex {x[1]} {y[1]} {-z[1]/2}
    vertex {x[1]} {y[1]} { z[1]/2}
  endloop
endfacet
""")

def convert_to_stl(name, espesores, xy, Quadrangles, conect, nodos_borde, Nelem):

  fid = open(f"{name}.stl","w")

  fid.write(f"solid {name}\n")



  #Promediar espesor en los nodos y determinar conectividad
  Nnodes = xy.shape[0]
  z = np.zeros(Nnodes)
  nconec = np.zeros(Nnodes,dtype=int)
  
  elementoslinea = len(conect)-len(Quadrangles)
  Nelem = elementoslinea + len(Quadrangles*4)
  
  conec = np.zeros((Nelem,4),dtype=int)
  
  index = Quadrangles[0]
  print(index)
  j = 0
  j2 = 0
  espesores2 = np.zeros((len(Quadrangles))*4)
  
  for e in Quadrangles:       
    if index <= Nelem-4:
        ni = int(conect[e,0])
        nj = int(conect[e,1])
        nk = int(conect[e,2])
        nl = int(conect[e,3])
        nm = int(conect[e,4])
        nn = int(conect[e,5])
        no = int(conect[e,6])
        npe = int(conect[e,7])
        nq = int(conect[e,8])
        
        
        conec1 = [ni,nm,nq,npe] #0487
        conec2 = [nm,nj,nn,nq]  #4158
        conec3 = [nq,nn,nk,no]  #8526
        conec4 = [npe,nq,no,nl] #7863
        
        conec[index] = conec1
        conec[index+1] = conec2
        conec[index+2] = conec3
        conec[index+3] = conec4
       
    if j2 <= len(espesores):
        espesor = espesores[j2] 
        espesores2[j] = espesor
        espesores2[j+1] = espesor
        espesores2[j+2] = espesor
        espesores2[j+3] = espesor
        #print(f"espesor = {espesores2[j]}")
        j +=4
        j2+=1
    index +=4
    
  
  print(f"len espesores = {len(espesores)}")
  espesores = espesores2
  dim_ini = len(Quadrangles)
  Quadrangles2 = np.zeros((dim_ini*4),dtype=int)
  j = 0  
  eu = Quadrangles[0]
  for e in Quadrangles:
      Quadrangles2[j] = eu
      Quadrangles2[j+1] = eu + 1
      Quadrangles2[j+2] = eu + 2
      Quadrangles2[j+3] = eu + 3
      j+=4
      eu+=4
      
  Quadrangles = Quadrangles2
  print(conec.shape)
  print(f"len espesores = {len(espesores)}")
  print(f"len Quadrangles = {len(Quadrangles)}")
  
  for i,e in enumerate(Quadrangles):
    nodes = conec[e,:]
    x = xy[nodes,0]
    y = xy[nodes,1]
    z[nodes] += espesores[i]
    #print(espesores2[i])
    nconec[nodes] += 1
  z = z/nconec
  #Ahora crear las facetas del STL a partir de los cuadrilateros
  for i,e in enumerate(Quadrangles):
      nodes = conec[e,:]
      x = xy[nodes,0]
      y = xy[nodes,1]
      zz = z[nodes]
      write_facet_from_quad_coords(fid, 10*x,10*y,10*zz)
 
  #Crear las facetas de los bordes
  for nodes in nodos_borde:
      print(nodes)
      x = xy[nodes,0]
      y = xy[nodes,1]
      zz = z[nodes]
      write_facet_from_line_coords(fid, 10*x,10*y,10*zz)
  
  #Cerrar el archivo
  fid.write(f"endsolid {name}\n")

  fid.close()
  
convert_to_stl(name,espesores,xy,Quadrangles,conect,nodos_borde,Nelem)
