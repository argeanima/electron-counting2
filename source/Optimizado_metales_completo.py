# -*- coding: utf-8 -*-

"""
Created on Wed Jan 26 18:47:35 2022

@author: Uriel Alejandro
"""

from collections import defaultdict
import math
import networkx as nx
import MatrizBO_completo as MBOC
import numpy as np
import itertools
import copy

global __ATOM_LIST__
__ATOM_LIST__ = \
    ['h',  'he',
     'li', 'be', 'b',  'c',  'n',  'o',  'f',  'ne',
     'na', 'mg', 'al', 'si', 'p',  's',  'cl', 'ar',
     'k',  'ca', 'sc', 'ti', 'v ', 'cr', 'mn', 'fe', 'co', 'ni', 'cu',
     'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr',
     'rb', 'sr', 'y',  'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag',
     'cd', 'in', 'sn', 'sb', 'te', 'i',  'xe',
     'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy',
     'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w',  're', 'os', 'ir', 'pt',
     'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn',
     'fr', 'ra', 'ac', 'th', 'pa', 'u',  'np', 'pu']


global atomic_valence

atomic_valence = defaultdict(list)
atomic_valence[1] = [1]
atomic_valence[5] = [3]
atomic_valence[6] = [4]
atomic_valence[7] = [3,4]
atomic_valence[8] = [3,2]
atomic_valence[9] = [1]
atomic_valence[14] = [4]
atomic_valence[15] = [5,3] #[5,4,3]
atomic_valence[16] = [6,3,2] #[6,4,2]
atomic_valence[17] = [1]
atomic_valence[32] = [4]
atomic_valence[35] = [1]
atomic_valence[53] = [1]

global atomic_valence_electrons

atomic_valence_electrons = {}
atomic_valence_electrons[1] = 1
atomic_valence_electrons[5] = 3
atomic_valence_electrons[6] = 4
atomic_valence_electrons[7] = 5
atomic_valence_electrons[8] = 6
atomic_valence_electrons[9] = 7
atomic_valence_electrons[14] = 4
atomic_valence_electrons[15] = 5
atomic_valence_electrons[16] = 6
atomic_valence_electrons[17] = 7
atomic_valence_electrons[32] = 4
atomic_valence_electrons[35] = 7
atomic_valence_electrons[53] = 7
radios_cov = [] # Radios covalentes(Encontrados en Wikipedia)
simbolos_cov = [] # Símbolos de átomos asociados a los radios

#Variables para generar la molecula
matriz_adj = []
simbolos_molecula = []
molecula = nx.Graph()
vertices = []
coordx = []
coordy = []
coordz = []
simbolos_cov = []
radios_cov = []

def simbolo_atomo(atom): #convert integer atom to string atom

    global __ATOM_LIST__
    atom = __ATOM_LIST__[atom - 1]
    
    return atom


def entero_atomo(atom): # convert str atom to integer atom
    global __ATOM_LIST__
    atom = atom.lower()
    
    return __ATOM_LIST__.index(atom) + 1
        

def distancia_R3(n,m,x,y,z): #Función de distancia en R3
    
    d = math.sqrt((x[n]-x[m])**2+(y[n]-y[m])**2+(z[n]-z[m])**2)
    
    return d

def numero_de_enlaces(n,ver,m_ad,enl_c):

    conect = 0
        
    for j in range(len(ver)):
        conect = conect+m_ad[j][n]
 
    #for par in enl_c:
     #   if (par[1] == ver[n]):
      #      conect = +1
    
    return conect

def id_radio_cov(mol):#Función que asocia un radio covalente a un símbolo 
    
    m = 0
    
    for i in range(len(simbolos_cov)):
        
        if(simbolos_cov[i]==mol):
            m = i
            break
            
    return radios_cov[m]

def crea_ma(mol): #Función que crea enlaces y matriz de adjacencia

    matriz_np = nx.to_numpy_matrix(mol, nodelist=mol.nodes) #Crea matriz de adjacencia
    
    return matriz_np.tolist()

def crea_info_t(ver): #Para cada vertice de la subgrafica se asocia información
    
    simbolos_t = []
    
    for atomo in ver:
        simbolos_t.append(simbolos_molecula[atomo])

    return simbolos_t
    
def detecta_coord_completos(mol,ver,sim_mol):
    
    C_completos = []
    enl_c_par = []
    enl_c = []
    
    for atomo in mol.nodes:
        coordinaciones=mol.degree(atomo)
        
        if (coordinaciones==4 and sim_mol[atomo]=="C"):
            C_completos.append(atomo)
            enl_c_par +=list(molecula.edges(atomo))
    
    for par in enl_c_par:
        enl_c.append(par[1])
    
    return enl_c, C_completos


#Parte molecula
def parte_molecula(mol,Cs):
    
    mol_c = mol
    for atomo in Cs:
        mol_c.remove_node(atomo)

    return mol

def componentes_conexas(mol):
    
    S = [mol.subgraph(c).copy() for c in nx.connected_components(mol)]
    
    return S

def corrige_matrizBO(vec_ver,vec_bo,ma_ad):
    
    BO = ma_ad

    for k in range(len(vec_bo)):
        for i in range(len(vec_ver[k])):
            for j in range(i,len(vec_ver[k])):
                BO[vec_ver[k][i]][vec_ver[k][j]] = vec_bo[k][i][j]
                BO[vec_ver[k][j]][vec_ver[k][i]] = vec_bo[k][j][i]
    
    return BO

def reconstruye_lista(vert_in,lista,nodos_eliminados,aristas_eliminadas):

    lista_ord = []
    pares_lista = []
    lista_el = []
  
    
    for nodo in nodos_eliminados:
        lista_el.append(4)
        
    vert_in = vert_in + [nodos_eliminados]
    lista = lista + [lista_el]
    
    for i in range(len(vert_in)):
        for j in range(len(vert_in[i])):
            pares_lista.append([vert_in[i][j],lista[i][j]])
    pares_lista=sorted(pares_lista, key=lambda x: x[0])
    
    for i in range(len(pares_lista)):
        lista_ord.append(pares_lista[i][1])
    
    return lista_ord

def reconstruye_cargas(vertices_ligantes,cargas_factibles):
    
    cargas_finales = []
    vertices = (np.array(vertices_ligantes)).sum()
    
    for lista_cargas in cargas_factibles:
        
        pares = []
        pares_ordenados = []
        cargas_ordenadas =[]
        
        cargas = (np.array(lista_cargas)).sum()
        
        for i in range(len(vertices)):
            pares.append([vertices[i],cargas[i]])
            
        pares_ordenados = sorted(pares, key=lambda x: x[0])
        
        for i in range(len(pares)):
            cargas_ordenadas.append(pares_ordenados[i][1])

        cargas_finales.append(cargas_ordenadas)  
          
    return cargas_finales
    

def list_enl_c(l_ver,vert,en_c):
    
    lista_carbonos = []
    lista_vecinos = []
    l_ordenes = []
    
    for i in range(len(vert)):
        if vert[i] in en_c:
            lista_carbonos.append(i)
        for j in range(len(vecinos_metal)):
            if (vert[i] == vecinos_metal[j]):
                lista_vecinos.append(i)
                l_ordenes.append(ordenes_vecinos_metal[j])

    return lista_carbonos, lista_vecinos, l_ordenes
                
def detecta_enl_coord(cargas):
    
    pares = []
    
    for i in range(len(cargas)):
        for j in range(i+1,len(cargas)):
            if (cargas[i]==1 and cargas[j] == -1):
                pares.append([i,j])
            elif (cargas[i]==-1 and cargas[j] == 1):
                pares.append([i,j])
    
    return pares

def enumera_grafica(mol):
    
    rmol = mol.copy()
    n = len(rmol.nodes)
    
    mapping = dict(zip(rmol, range(0, n)))
    emol = nx.relabel_nodes(rmol, mapping)
    
    return emol

def corrige_carga(carbonos,cargas_lista):

    cargas_completas = []
    for cargas_f in cargas_lista:
        cargas_carb = []
        for cargas in cargas_f:
            cargas_carb += cargas
        for carbono in carbonos:
            cargas_carb.insert(carbono,0)
        cargas_completas.append(cargas_carb)
        
    return cargas_completas

def corrige_BO(soluciones_bos,vector_vertices):
    
    BO = matriz_adj.copy()

    for k in range(len(soluciones_bos)):
        for i in range(len(vector_vertices[k])):
            for j in range(len(vector_vertices[k])):
                BO[vector_vertices[k][i]][vector_vertices[k][j]] = soluciones_bos[k][i][j]
 
    return BO

def corrige_matricesBOs (soluciones_bos,q_lista,vertices_f,carbonos_completos):    

    BOS = []
    soluciones_BO = []
    combinaciones_lista = []
    cargas = []
    
    for i in range(len(soluciones_bos)):
        combinaciones_lista.append([ n for n in range(len(soluciones_bos[i]))])   
    combinaciones_soluciones = list(itertools.product(*combinaciones_lista))  
    
    for combinacion in combinaciones_soluciones:
        
        posible_carga = []
        posibles_BOS = []
        
        for i in range(len(list(combinacion))):
            posible_carga.append(q_lista[i][combinacion[i]])
            posibles_BOS.append(soluciones_bos[i][combinacion[i]])
    
    cargas.append(posible_carga)
    BOS.append(posibles_BOS)
    
    for i in range(len(BOS)):
        BO = corrige_BO(BOS[i],vertices_f)
        BO_copy = copy.deepcopy(BO)
        soluciones_BO.append(BO_copy)
    
    sol_carga = corrige_carga(carbonos_completos,cargas)
        
    return soluciones_BO, sol_carga
                
def optimiza():
    
    carga = 0 
    BOS_lista_f = [] #Conjunto de matrices cuadradas de posibles soluciones
    q_lista_f = []
    vertices_f = []

    
    enl_car,carbonos_completos = detecta_coord_completos(molecula,vertices,simbolos_molecula) #Se detectan los atomos de Carbono con 4 coordinaciones y los atomos coordinados a ellos
    molecula_copia = molecula.copy() #Se copia la molecula para no perder la infomación original
    molecula_copia = parte_molecula(molecula_copia,carbonos_completos) #Parte_molecula
    comp = componentes_conexas(molecula_copia) 
    
    for fragmento in comp:

    
        vert_f = list(fragmento.nodes)#
        fragmento_en = enumera_grafica(fragmento)#
        l_ver = list(fragmento_en.nodes)#
        sim_mol_t = crea_info_t(vert_f)#
        l_enl_c, l_vecinos, lista_ordenes = list_enl_c(l_ver,vert_f,enl_car)
        ma_ad_t = crea_ma(fragmento_en) #Crea enlaces y matriz de adjacencia
        posibles_soluciones ,posibles_cargas = MBOC.crea_BO(l_ver,sim_mol_t, ma_ad_t,carga,l_enl_c,lista_ordenes,l_vecinos) 
        if len(posibles_soluciones) == 0:   
            return [],[]
        
        BOS_lista_f.append(posibles_soluciones)
        q_lista_f.append(posibles_cargas)
        vertices_f.append(vert_f)
    soluciones_BO,soluciones_cargas = corrige_matricesBOs(BOS_lista_f,q_lista_f,vertices_f,carbonos_completos)  

    return soluciones_BO, soluciones_cargas

def crea_molecula(molecula_in,vertices_in,simbolos_molecula_in,matriz_adj_in,carga,ordenes_in,vecinos_in):
    
    global simbolos_molecula #Declaramos globales los datos recibidos
    global molecula 
    global vertices
    global matriz_adj
    global ordenes_vecinos_metal
    global vecinos_metal
    
    simbolos_molecula = simbolos_molecula_in.copy() #Declaramos globales los datos recibidos
    molecula = molecula_in.copy()
    vertices = vertices_in.copy()
    matriz_adj = matriz_adj_in.copy()
    ordenes_vecinos_metal = list(ordenes_in).copy()
    vecinos_metal = vecinos_in.copy()
    
    
    if len(molecula.nodes) <= 100:    #Si el total de atomos es mejor o igual a n se aplica el algoritmo sin optimizar

        posibles_soluciones, posibles_cargas = MBOC.crea_BO(vertices,simbolos_molecula,matriz_adj,carga,[],ordenes_vecinos_metal,vecinos_metal)            

        return posibles_soluciones, posibles_cargas
    
        
    
    else: #En otro caso se aplica el algoritmo de optimización
        posibles_soluciones, posibles_cargas = optimiza()

        return posibles_soluciones, posibles_cargas

        

