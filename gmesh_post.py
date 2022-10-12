"""
$NodeData
number-of-string-tags
< "string-tag" >
…
number-of-real-tags
< real-tag >
…
number-of-integer-tags
< integer-tag >
…
node-number value …
…
$EndNodeData
"""
 
def write_node_data(file_name,nodes,data, nombre_datos=""): #dame un número de nodos y dame un valor de desplazamientos, se pueden escribir dezpl. en x en un archivo, en y en un archivo y juntos
    fid= open(file_name,"w")
    Ndata = len(nodes)
    #HEADER
    fid.write("$MeshFormat\n")
    fid.write("2.2 0 8\n")
    fid.write("$EndMeshFormat\n")
    fid.write("$NodeData\n")
    fid.write("1\n")
    fid.write(f"\"{nombre_datos}\"\n")
    fid.write("1\n")
    fid.write("0\n") #puede indicar el tiempo atmosférico                                                               (bromas XD)
    fid.write("3\n")
    fid.write("0\n") #puede indicar el tiempo (temporal)
    fid.write("1\n") #dimension de datos 1 porque grafico escalares
    fid.write(f"{Ndata}\n") #numero de datos
    
    for i in range(Ndata):
        fid.write(f"{nodes[i]} {data[i]}\n")
    
    fid.write("$EndNodeData\n")
    fid.close()
    return     
    

def write_node_data_2(file_name,nodes,data1,data2, nombre_datos=""): #dame un número de nodos y dame un valor de desplazamientos, se pueden escribir dezpl. en x en un archivo, en y en un archivo y juntos
    fid= open(file_name,"w")
    Ndata = len(nodes)
    #HEADER
    fid.write("$MeshFormat\n")
    fid.write("2.2 0 8\n")
    fid.write("$EndMeshFormat\n")
    fid.write("$NodeData\n")
    fid.write("1\n")
    fid.write(f"\"{nombre_datos}\"\n")
    fid.write("1\n")
    fid.write("0\n") #puede indicar el tiempo atmosférico                                                               (bromas XD)
    fid.write("3\n")
    fid.write("0\n") #puede indicar el tiempo (temporal)
    fid.write("3\n") #dimension de datos 1 porque grafico escalares
    fid.write(f"{Ndata}\n") #numero de datos
    
    for i in range(Ndata):
        fid.write(f"{nodes[i]} {data1[i]} {data2[i]} 0.0 \n")

    fid.write("$EndNodeData\n")
    fid.close()
    return     

def write_element_data(file_name,elements,data, nombre_datos=""): #dame un número de nodos y dame un valor de desplazamientos, se pueden escribir dezpl. en x en un archivo, en y en un archivo y juntos
    fid= open(file_name,"w")
    Ndata = len(elements)
    #HEADER
    fid.write("$MeshFormat\n")
    fid.write("2.2 0 8\n")
    fid.write("$EndMeshFormat\n")
    fid.write("$ElementData\n")
    fid.write("1\n")
    fid.write(f"\"{nombre_datos}\"\n")
    fid.write("1\n")
    fid.write("0\n") #puede indicar el tiempo atmosférico                                                               (bromas XD)
    fid.write("3\n")
    fid.write("0\n") #puede indicar el tiempo (temporal)
    fid.write("1\n") #dimension de datos 1 porque grafico escalares
    fid.write(f"{Ndata}\n") #numero de datos
    
    for i in range(Ndata):
        fid.write(f"{elements[i]} {data[i]} \n")

    fid.write("$EndElementData\n")
    fid.close()
    return   
  
def write_element_node_data(file_name,elements,data, nombre_datos=""): #dame un número de nodos y dame un valor de desplazamientos, se pueden escribir dezpl. en x en un archivo, en y en un archivo y juntos
    fid= open(file_name,"w")
    Ndata = len(elements)
    #HEADER
    fid.write("$MeshFormat\n")
    fid.write("2.2 0 8\n")
    fid.write("$EndMeshFormat\n")
    fid.write("$ElementNodeData\n")
    fid.write("1\n")
    fid.write(f"\"{nombre_datos}\"\n")
    fid.write("1\n")
    fid.write("0\n") #puede indicar el tiempo atmosférico                                                               (bromas XD)
    fid.write("3\n")
    fid.write("0\n") #puede indicar el tiempo (temporal)
    fid.write("1\n") #dimension de datos 1 porque grafico escalares
    fid.write(f"{Ndata}\n") #numero de datos
    print (f"pene: {data.shape}")
    for i in range(Ndata):
        
        #print (f"{data[i,0]} {data[i,1]} {data[i,2]} {data[i,3]} {data[i,4]} {data[i,5]} {data[i,6]} {data[i,7]} {data[i,8]}")
        #fid.write(f"{elements[i]} 9 {data[i][0]} {data[i][1]} {data[i][2]} {data[i][3]} {data[i][4]} {data[i][5]} {data[i][6]} {data[i][7]} {data[i][8]}\n")
        fid.write(f"{elements[i]} 9 {data[i,0]} {data[i,1]} {data[i,2]} {data[i,3]} {data[i,4]} {data[i,5]} {data[i,6]} {data[i,7]} {data[i,8]}\n")

    fid.write("$EndElementNodeData\n")
    fid.close()
    return     
