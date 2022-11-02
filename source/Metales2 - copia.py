# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 16:44:30 2022

@author: Uriel Alejandro
"""
# Hola
from timeit import default_timer
from collections import defaultdict
import math
import networkx as nx
import numpy as np
import Optimizado_metales_completo as OM
import itertools
import copy
import os
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

metales_transicion = [21,22,25,26,41,42,44,45,46,73,74,75,76,77,78]
global atomic_valence

atomic_valence = defaultdict(list)
atomic_valence[1] = [1]
atomic_valence[5] = [3,4]
atomic_valence[6] = [4]
atomic_valence[7] = [4,3]
atomic_valence[8] = [3,2]
atomic_valence[9] = [1]
atomic_valence[14] = [4]
atomic_valence[15] = [4,3] #[5,3]
atomic_valence[16] = [6,3,2] #[6,4,2]
atomic_valence[17] = [1]
atomic_valence[21] = [15,14,13,12,11,10]
atomic_valence[22] = [14,13,12,11,10]
atomic_valence[25] = [7,8,9,10,11]
atomic_valence[26] = [5,6,7,8,9,10]
atomic_valence[32] = [4]
atomic_valence[35] = [1]
atomic_valence[41] = [13,12,11,10,9]
atomic_valence[42] = [12,11,10,9]
atomic_valence[44] = [5,6,7,8,9,10]
atomic_valence[45] = [5,6,7,8,9]
atomic_valence[46] = [4,5,6,7,8]
atomic_valence[53] = [1]
atomic_valence[73] = [13,12,11,10,9]
atomic_valence[74] = [12,11,10,9]
atomic_valence[75] = [11,10,9,8]
atomic_valence[76] = [10,9,8,7,6]
atomic_valence[77] = [9,8,7,6,5]
atomic_valence[78] = [8,7,6,5,4]

global numeros_oxidacion

numeros_oxidacion = defaultdict(list)

numeros_oxidacion[21] = [3]
numeros_oxidacion[22] = [4]
numeros_oxidacion[23] = [2,3,4,5]
numeros_oxidacion[24] = [2,3,6]
numeros_oxidacion[25] = [2,3,4,6,7]
numeros_oxidacion[26] = [2,3]
numeros_oxidacion[27] = [2,3]
numeros_oxidacion[28] = [2]
numeros_oxidacion[29] = [1,2]
numeros_oxidacion[30] = [2]
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
atomic_valence_electrons[21] = 3
atomic_valence_electrons[22] = 4
atomic_valence_electrons[25] = 7
atomic_valence_electrons[26] = 8
atomic_valence_electrons[32] = 4
atomic_valence_electrons[35] = 7
atomic_valence_electrons[41] = 5
atomic_valence_electrons[42] = 6
atomic_valence_electrons[44] = 8
atomic_valence_electrons[45] = 9
atomic_valence_electrons[46] = 10
atomic_valence_electrons[53] = 7
atomic_valence_electrons[73] = 5
atomic_valence_electrons[74] = 6
atomic_valence_electrons[75] = 7
atomic_valence_electrons[76] = 8
atomic_valence_electrons[77] = 9
atomic_valence_electrons[78] = 10
metales_transicion = ["Sc","Ti","Mn","Fe","Nb","Mb","Ru","Rh","Pd","Ta","W","Re"
                      "Os","Ir","Pt","Mo"]
radios_cov = [] # Radios covalentes(Encontrados en Wikipedia)
simbolos_cov = [] # Símbolos de átomos asociados a los radios

#Variables para generar la molecula
matriz_adj = []
simbolos_moleculas = []
moleculas = nx.Graph()
vertices_completos = []
coordx = []
coordy = []
coordz = []
# Imprimir
atomos_mol2 = []
enlaces_mol2 = []

with open('Radios_Covalentes.txt') as rad: #Lee documento de radios covalentes
    lineas = rad.readlines()
rad.close()

for i in range (len(lineas)):
    k=lineas[i].split()
    simbolos_cov.append(k[0]) 
    radios_cov.append(float(k[1]))

def simbolo_atomo(atom): #convert integer atom to string atom

    global __ATOM_LIST__
    atom = __ATOM_LIST__[atom - 1]
    
    return atom


def entero_atomo(atom): # convert str atom to integer atom
    global __ATOM_LIST__
    atom = atom.lower()
    
    return __ATOM_LIST__.index(atom) + 1
    
def lee_xyz(nombre):#función que lee documento .xyz 
    
    #nodo_sh=0   
    #print(os.path.splitext(file_name)[0])
    with open(nombre) as xyz:  #Abre archivo xyz
        lineas = xyz.readlines()
        n = int(lineas[0])
    
    for i in range (2,2+n): #Lee las coordenadas y atomos
        k = lineas[i].split()
        simbolos_moleculas.append(k[0]) #Simbolo de cada átomo
        moleculas.add_node(len(coordx))
        vertices_completos.append(len(coordx)) #Número de etiqueta del vertice
        coordx.append(float(k[1]))
        coordy.append(float(k[2]))
        coordz.append(float(k[3]))
    xyz.close()
    

def distancia_R3(n,m): #Función de distancia en R3
    
    d = math.sqrt((coordx[n]-coordx[m])**2+(coordy[n]-coordy[m])**2+(coordz[n]-coordz[m])**2)
    
    return d

def numero_de_enlaces(n,ver,m_ad):

    conect = 0
        
    for j in range(len(ver)):
        conect = conect+m_ad[j][n]
 
    return conect

def id_radio_cov(mol):#Función que asocia un radio covalente a un símbolo 
    
    m = 0
    
    for i in range(len(simbolos_cov)):
        
        if(simbolos_cov[i]==mol):
            m = i
            break
            
    return radios_cov[m]


def crea_ma(mol): #Función que crea enlaces y matriz de adjacencia
    
    for i in range(len(vertices_completos)): #Creación de enlaces 
        r1 = id_radio_cov(simbolos_moleculas[vertices_completos[i]])
        
        for j in range(i+1,len(vertices_completos)):
            r2 = id_radio_cov(simbolos_moleculas[vertices_completos[j]])
            dis = distancia_R3(i,j)                
            if(dis<=(r1+r2)*1.2): #Tolerancia de 0.3, valor arb.
                mol.add_edge(i, j)
                
    matriz_np = nx.to_numpy_matrix(mol, nodelist=mol.nodes)#Crea matriz de adjacencia
    
    return matriz_np.tolist()

def crea_sub_ma(mol): #Función que crea enlaces y matriz de adjacencia

    matriz_np = nx.to_numpy_matrix(mol, nodelist=mol.nodes) #Crea matriz de adjacencia
    
    return matriz_np.tolist()

def crea_simbolos_individuales(vertices): #Para cada vertice de la subgrafica se el símbolo atómico
    
    simbolos_ind = []
    
    for atomo in vertices:
        simbolos_ind.append(simbolos_moleculas[atomo])

    return simbolos_ind

def crea_enlaces_print(BO_final):
    
    
    enlaces=[]
    for i in vertices_completos:
        for j in range(i+1,len(vertices_completos)):
            if BO_final[i][j]!=0:
                k=[i,j,int(BO_final[i][j])]
                enlaces.append(k)
    
    return enlaces

def detecta_ciclos(mol):
    
    ciclos=nx.cycle_basis(mol, 0)
    
    return ciclos  

def colinealidad_buena(ciclo):

    x1 = coordx[ciclo[0]]
    x2 = coordx[ciclo[1]]
    x3 = coordx[ciclo[2]]
    y1 = coordy[ciclo[0]]
    y2 = coordy[ciclo[1]]
    y3 = coordy[ciclo[2]]
    z1 = coordz[ciclo[0]]
    z2 = coordz[ciclo[1]]
    z3 = coordz[ciclo[2]]

    for i in range(3,len(ciclo)):
        a1 = coordx[ciclo[i]]-x1
        a2 = coordy[ciclo[i]]-y1
        a3 = coordz[ciclo[i]]-z1
        b1 = x2-x1
        b2 = y2-y1
        b3 = z2-z1
        c1 = x3-x1
        c2 = y3-y1
        c3 = z3-z1
        M = np.array([[a1, a2, a3], [b1, b2, b3], [c1, c2 ,c3]])
        det = np.linalg.det(M)
        if(0.5<abs(det)):
            
            return False
            
    return True
    
def pares_buenos(ciclo):

    #e_pi=0

    #for atomo in ciclo:
     #   entero = entero_atomo(simbolos_moleculas[atomo])
      #  coordinaciones = moleculas.degree(atomo)
      #  eletrones_valencia = atomic_valence_electrons[entero]
      #  valencia = valencias_finales[atomo]
      #  pi = valencia-coordinaciones-eletrones_valencia
      #  e_pi+=pi
    #if(e_pi%4==2):
        #return True
   # return False
   return False
def ciclos_mol2(atomos_p,enlaces_p,ciclos_ar):
    
    for ciclo in ciclos_ar:
        for atomo in ciclo:
            sintaxis = atomos_p[atomo]+".ar"
            atomos_p[atomo]=sintaxis
            for enlace in range(len(enlaces_p)):
                if(enlaces_p[enlace][0]==atomo or enlaces_p[enlace][1]==atomo):
                    if(enlaces_p[enlace][0]==atomo):
                        for atomo2 in ciclo:
                            if (enlaces_p[enlace][1]==atomo2):
                                enlaces_p[enlace][2] = 'ar'
                                break
                    if(enlaces_p[enlace][1]==atomo):
                        for atomo2 in ciclo:
                            if (enlaces_p[enlace][2]==atomo2):
                                enlaces_p[enlace][2] = 'ar'

    return atomos_p,enlaces_p

def identifica_ciclos_ar(atomos_p,enlaces_p,ciclos):
    
    ciclos_ar = []
    
    for ciclo in ciclos:
        verifica_colinealidad = colinealidad_buena(ciclo)

        verficia_pares = pares_buenos(ciclo)
        
        if verifica_colinealidad and verficia_pares:
            ciclos_ar.append(ciclo)
            
    atomos_p,enlaces_p = ciclos_mol2(atomos_p,enlaces_p, ciclos_ar)
    
    return atomos_p,enlaces_p

def valencias(BO_final):
    
    val = []
    
    for i in range(len(BO_final)):
        suma = 0
        for j in range(len(BO_final)):
            suma = + BO_final[i][j]
        val.append(suma)
        
    return val
        
def crea_archivo_mol(BO_final,nombre,ciclos,cargas,electrones_libres): # Función que crea del documento .mol2  

    atomos_p = simbolos_moleculas.copy()
    enlaces_p = crea_enlaces_print(BO_final)
    atomos_p,enlaces_p = identifica_ciclos_ar(atomos_p,enlaces_p,ciclos)
    nombre_xyz = nombre+'.mol2'
    f = open(nombre_xyz, "w") 
    f.write("@<TRIPOS>MOLECULE\n")
    f.write(" Molden generated mol2\n")
    n1=len(simbolos_moleculas)
    n2=len(enlaces_p)
    f.write("    %d     %d      1\n"%(n1,n2))
    f.write(" SMALL\n")
    f.write(" CHARGES\n")
    f.write(" ****\n")
    f.write(" ****\n")
    f.write("@<TRIPOS>ATOM\n")
    for i in range (len(simbolos_moleculas)): #Escribe información de los átomos de las moléculas
        f.write("{:>6}  {:>1}    {: .4f}    {: .4f}    {: .4f}  {:<3}   {:}    {:}\n".format(int(cargas[i]),simbolos_moleculas[i],coordx[i],coordy[i],coordz[i],atomos_p[i],"2 RES1",cargas[i]))
        #f.write("{:>6}  {:>1}    {: .4f}    {: .4f}    {: .4f}  {:<3}   {:}    {:}\n".format(i+1,simbolos_moleculas[i],coordx[i],coordy[i],coordz[i],atomos_p[i],"2 RES1",electrones_libres[i]/2))
    f.write("@<TRIPOS>BOND\n")
    for i in range (len(enlaces_p)): #Escribe información de los enlaces de las moléculas
        f.write("{:>6}{:>6}{:>6}{:>7}\n".format(i+1,enlaces_p[i][0]+1,enlaces_p[i][1]+1,enlaces_p[i][2])) 
    f.write("@<TRIPOS>SUBSTRUCTURE\n")
    f.write("      1 RES1       1\n")
    f.close()

def componentes_conexas(mol):
    
    S = [mol.subgraph(c).copy() for c in nx.connected_components(mol)]
    
    return S


def corrige_BO(vec_ver,vector_bos,ma_ad):
    
    BO = ma_ad.copy()
    
    for k in range(len(vector_bos)):
        for i in range(len(vec_ver[k])):
            for j in range(len(vec_ver[k])):
                BO[vec_ver[k][i]][vec_ver[k][j]] = vector_bos[k][i][j]
    
    return BO        

def corrige_matricesBOs(vector_vertices,soluciones_bos,ma_ad):
    
    soluciones_BO = []

    for soluciones in soluciones_bos:
        
        solucion_BO = corrige_BO(vector_vertices,soluciones,ma_ad)
        solucion_BO_copy = copy.deepcopy(solucion_BO)
        soluciones_BO.append(solucion_BO_copy)
        
    return soluciones_BO
  
def corrige_BO_metal(vecinos,BO,lista_ordenes,metal_molecula):
    

            
    for i in range(len(lista_ordenes)):
        BO[vecinos[i]][metal_molecula] = lista_ordenes[i]
        BO[metal_molecula][vecinos[i]] = lista_ordenes[i]   

    return BO

def corrige_BO_ligante(vector_vertices,soluciones_bos,ma_ad,vecinos_metales,ordenes,metales_molecula):
    
    BO = ma_ad.copy()

    for k in range(len(soluciones_bos)):
        for i in range(len(vector_vertices[k])):
            for j in range(len(vector_vertices[k])):
                BO[vector_vertices[k][i]][vector_vertices[k][j]] = soluciones_bos[k][i][j]

    for i in range(len(metales_molecula)):
        BO_final = corrige_BO_metal(vecinos_metales[i],BO.copy(),ordenes[i],metales_molecula[i])
    
    return BO_final

def corrige_matricesBOs_ligante (vector_vertices,soluciones_bos,ma_ad,vecinos_metal,ordenes,metales_molecula):    

    soluciones_BO = []

    for i in range(len(soluciones_bos)):
        
        solucion_BO = corrige_BO_ligante(vector_vertices,soluciones_bos[i],ma_ad,vecinos_metal,ordenes[i] ,metales_molecula)
        solucion_BO_copy = copy.deepcopy(solucion_BO)
        soluciones_BO.append(solucion_BO_copy)

        
    return soluciones_BO

def corrige_matrizBOs_moleculas(vec_ver,vector_bos,ma_ad):
    
    BOS = []
     
    for k in range(len(vector_bos[0])):
        BO = ma_ad
        for i in range(len(vec_ver[k])):
            for j in range(i,len(vec_ver[k])):
                BO[vec_ver[k][i]][vec_ver[k][j]] = 0#vector_bos[k][i][j]
                BO[vec_ver[k][j]][vec_ver[k][i]] = 0 #vector_bos[k][i][j]
        BOS.append(BO)
    
    return BOS

def ajustar_ciclos(ciclos_enumerados,vertices_originales,metales_molecula):
    
    ciclos_originales = []
    
    for ciclo in ciclos_enumerados:
        ciclo_original = []
        for atomo in ciclo:
            ciclo_original.append(vertices_originales[atomo])
        for metal in metales_molecula:
            if metal not in ciclo_original: 
                ciclos_originales.append(ciclo_original)
    
    return ciclos_originales

def enumera_grafica(mol):
    
    rmol = mol.copy()
    n = len(rmol.nodes)
    
    mapping = dict(zip(rmol, range(0, n)))
    emol = nx.relabel_nodes(rmol, mapping)
    
    return emol


def detecta_metales(lista_vertices,ma_ad,sim_mol): #Detecta los metales en un conjunto de vertices
    
    metales = []
    for atomo in lista_vertices:
        for metal in metales_transicion:
            if(simbolos_moleculas[atomo] == metal):
                metales.append(atomo)           

    return metales

def suprime_metales(molecula,metales): #Elimina metales del la molecula organometálica
    
    vecinos_metales = []
    enlaces_metales = []
    
    molecula_sn_metal = molecula.copy()
    #Guarda los vecinos de cada metal detectado
    for metal in metales:  
        enlaces_metales.append(moleculas.edges(metal)) #Guarda los vertices del metal
    for enlaces_metal in enlaces_metales:
        vecinos_metal = []
        for enlace in enlaces_metal:
            vecinos_metal.append(enlace[1])
        #Guarda conjuntos de vecinos, uno por metal.
        vecinos_metales.append(vecinos_metal) 
    #Elimina los metales de la molecula organometálcia
    for metal in metales:    
        molecula_sn_metal.remove_node(metal)
    
    ligantes = componentes_conexas(molecula_sn_metal)
    
    return vecinos_metales, ligantes

def ordenes_enlaces(lista_vecinos,vecinos_metal):
    
    ordenes = [1,2,3] #Posibles ordenes de los enlaces con el metal
    ordenes_lista_metales = []
    combinaciones_lista = []
    combinaciones_ordenes_finales = []
    combinaciones_ordenes_finales_pmetal = []
   # ordenes_lista = []  
    #combinaciones_ordenes_test = []
    for lista in lista_vecinos:
        ordenes_metal = []
        for i in range(len(lista)):
            ordenes_metal.append(ordenes)
        combinaciones_ordenes = list(itertools.product(*ordenes_metal))
        ordenes_lista_metales.append(combinaciones_ordenes)  

    for i in range(len(ordenes_lista_metales)):
        combinaciones_lista.append([ n for n in range(len(ordenes_lista_metales[i]))])
    
    combinaciones_soluciones = list(itertools.product(*combinaciones_lista))
    
    for combinacion in combinaciones_soluciones:
        
        ordenes_pmetal = []
        ordenes_pmetal_ind = []
        
        for i in range(len(list(ordenes_lista_metales))):
            ordenes_pmetal += list(ordenes_lista_metales[i][combinacion[i]])
            ordenes_pmetal_ind.append(list(ordenes_lista_metales[i][combinacion[i]]))
        combinaciones_ordenes_finales.append(ordenes_pmetal)
        combinaciones_ordenes_finales_pmetal.append(ordenes_pmetal_ind)

    #for i in range(len(vecinos_metal)):
     #   ordenes_lista.append(ordenes)
    #combinaciones_ordenes_test = list(itertools.product(*ordenes_lista))
    return combinaciones_ordenes_finales, combinaciones_ordenes_finales_pmetal,
    #return combinaciones_ordenes_test , combinaciones_ordenes_finales_pmetal
        
def vecinos_ligante(vertices_ligante,vecinos_metales): 

    # Identifica los vecinos que se encuentran en cada ligante
    vec_enumerados = [] #Vecinos reenumerados
    vec_metales = [] # Vecinos originales por metal
    
    for vecinos in vecinos_metales:
        vec_metal = []
        for i in range(len(vertices_ligante)):
            if vertices_ligante[i] in vecinos:
                vec_enumerados.append(i)
                vec_metal.append(vertices_ligante[i])
            
        vec_metales.append(vec_metal)
    
    return vec_metales, vec_enumerados

def corrige_cargas(vertices_ligantes,cargas_factibles):
    
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

def neutralidad_ligante(cargas, ordenes):
    
    for i in range(len(cargas)):
        
        carga_ligante = abs(np.array(cargas[i]).sum())-np.array(ordenes[i]).sum()
        
        if (carga_ligante != 0):
            
            return False
    
    return True

def vefifica_no_oxidacion (cargas_posibles,metales_ligantes,vecinos_ligante,vertices_ligantes,ordenes_posibles_pmetal):  

    for i in range (len(metales_ligantes)):

        carga = 0
        no_enlaces = 0
        no_coordinaciones = 0
        #prueba = []
        for j in range(len(ordenes_posibles_pmetal)):
            no_coordinaciones += len(ordenes_posibles_pmetal[j][i])

        for j in range(len(vertices_ligantes)):
            if (len(list(set(vecinos_ligante[j][i]) & set(vertices_ligantes[j])))!=0):
                
                carga += abs(np.array(cargas_posibles[j]).sum())
                no_enlaces += np.array(ordenes_posibles_pmetal[j]).sum()
        no_oxidacion = no_enlaces-carga
        
        numero_atomico = entero_atomo(simbolos_moleculas[metales_ligantes[i]])
        if (no_oxidacion in numeros_oxidacion[numero_atomico]):
            
            continue 
        else:
            return False, no_oxidacion
    return True, no_oxidacion
def carga_molecula_metal(cargas_posibles,ordenes_posibles,ordenes_posibles_pmetal,metales_ligantes,vecinos_ligante,vertices_ligantes):

    numero_e = [12,14,16,18]
    cargas_lista = []

    no_oxidacion_ok, no_oxidacion = vefifica_no_oxidacion (cargas_posibles,metales_ligantes,vecinos_ligante,vertices_ligantes,ordenes_posibles_pmetal)
    #if no_oxidacion_ok == False:
        #return False, []
    for i in range(len(metales_ligantes)):
        cargas_lista.append(numero_e)
    combinaciones_cargas = list(itertools.product(*cargas_lista))
    cargas_ligantes =np.array(cargas_posibles).sum()
    enlaces_totales = np.array(ordenes_posibles).sum()
    enlaces = np.array(enlaces_totales).sum()
    for combinacion in combinaciones_cargas:
        carga_metales = enlaces-np.array(combinacion).sum()
        for metal in metales_ligantes:
            numero_atomico_metal = entero_atomo(simbolos_moleculas[metal])
            carga_metales += atomic_valence_electrons[numero_atomico_metal]
        carga_total = np.array(cargas_ligantes).sum()+ carga_metales
        if carga_total == carga_molecula:
            #print("Número de Oxidación:", no_oxidacion)
            return True, combinacion
    return False, []

def carga_metal_ligante(cargas_metal,vertices_molecula,cargas_ligante):
    

    for i in range(len(cargas_ligante)):
        for j in range(len(metales)):
            posicion_metal = vertices_molecula.index(metales[j])   
            cargas_ligante[i].insert(posicion_metal, cargas_metal[i][j])


    return cargas_ligante  
    
                
def soluciones_ligantes(vertices_ligantes,simbolos_ligantes,vecinos_ligantes,BOS_lista_ligantes,q_lista_ligantes,ordenes_vecinos_ligantes,ordenes_vecinos_ligantes_pmetal,metales_ligantes):   
    
    combinaciones_lista = []
    BO_soluciones = []
    carga_soluciones = []
    ordenes_soluciones = []
    electrones_metal = []
    
    for i in range(len(vertices_ligantes)):
        combinaciones_lista.append([ n for n in range(len(ordenes_vecinos_ligantes[i]))])
    
    combinaciones_soluciones = list(itertools.product(*combinaciones_lista))
    
    for combinacion in combinaciones_soluciones:
        
        posible_carga = []
        posible_orden = []
        posible_orden_pmetal = []
        posibles_BOS = []
        
        for i in range(len(list(combinacion))):
            posible_carga.append(q_lista_ligantes[i][combinacion[i]])
            posible_orden.append(ordenes_vecinos_ligantes[i][combinacion[i]])
            posible_orden_pmetal.append(ordenes_vecinos_ligantes_pmetal[i][combinacion[i]])
            posibles_BOS.append(BOS_lista_ligantes[i][combinacion[i]])
        ligantes_neutros, e_metales = carga_molecula_metal(posible_carga,posible_orden,posible_orden_pmetal,metales_ligantes,vecinos_ligantes,vertices_ligantes)

        if (ligantes_neutros == True):
            
            BO_soluciones.append(posibles_BOS)
            carga_soluciones.append(posible_carga)
            ordenes_soluciones.append(posible_orden)
            electrones_metal.append(e_metales)

    return BO_soluciones,carga_soluciones, ordenes_soluciones, electrones_metal

        
def crea_molecula_coordinada(molecula, vertices_molecula,carga):
        
    ciclos_lista = [] 
    vertices_ligantes = []
    vertices_ligantes_enumerados = []
    simbolos_ligantes = []
    vecinos_ligantes = []
    BOS_lista_ligantes = [] #Conjunto de matrices cuadradas de posibles soluciones
    q_lista_ligantes = []
    ordenes_vecinos_ligantes = []
    ordenes_vecinos_ligantes_pmetal = []
    
    matriz_adj_mol = crea_sub_ma(molecula) #Crea matriz de adj. para la molecula
    simbolos_molecula = crea_simbolos_individuales(vertices_molecula)
    metales_molecula = detecta_metales(molecula,matriz_adj_mol,simbolos_molecula) #Metales de la molecula
    vecinos_metales,ligantes = suprime_metales(molecula,metales_molecula) #Detecta los vecinos de cada metal

    for ligante in ligantes:        
        
        BOS_lista = []
        q_lista = []
        ordenes_vecinos = []
        ordenes_vecinos_pmetal = []
        
        vertices_ligante = list(ligante.nodes).copy()
        simbolos_ligante = crea_simbolos_individuales (vertices_ligante)
        vecinos,vecinos_enumerados = vecinos_ligante(vertices_ligante,vecinos_metales.copy())
        posible_ord_enlaces, posible_ord_enlaces_pmetal = ordenes_enlaces(vecinos.copy(),vecinos_enumerados.copy())
        ligante_enumerado = enumera_grafica(ligante)
        vertices_ligante_enumerado = list(ligante_enumerado.nodes).copy()
        ma_adj_ligante = crea_sub_ma(ligante_enumerado)
        ciclos_ligante = detecta_ciclos(ligante_enumerado)
        ciclos_ajustados = ajustar_ciclos (ciclos_ligante,vertices_ligante,metales_molecula)
        
        vertices_ligantes_enumerados.append(vertices_ligante_enumerado)
        vertices_ligantes.append(vertices_ligante)
        simbolos_ligantes.append(simbolos_ligante)
        vecinos_ligantes.append(vecinos)
        ciclos_lista += ciclos_ajustados
        
        for i in range(len(posible_ord_enlaces)):
            
            posibles_soluciones_ligante ,posibles_cargas_ligante = OM.crea_molecula(ligante_enumerado,vertices_ligante_enumerado,simbolos_ligante,ma_adj_ligante,0,posible_ord_enlaces[i],vecinos_enumerados) 

            if len(posibles_soluciones_ligante) == 0:   
                continue
            for j in range(len(posibles_soluciones_ligante)):
                BOS_lista.append(posibles_soluciones_ligante[j])
                q_lista.append(posibles_cargas_ligante[j])
                ordenes_vecinos.append(posible_ord_enlaces[i])
                ordenes_vecinos_pmetal.append(posible_ord_enlaces_pmetal[i])

        BOS_lista_ligantes.append(BOS_lista)
        q_lista_ligantes.append(q_lista)
        ordenes_vecinos_ligantes.append(ordenes_vecinos)
        ordenes_vecinos_ligantes_pmetal.append(ordenes_vecinos_pmetal)
    soluciones_factibles,cargas_factibles, ordenes_factibles, cargas_metal = soluciones_ligantes(vertices_ligantes,simbolos_ligantes,vecinos_ligantes,BOS_lista_ligantes,q_lista_ligantes,ordenes_vecinos_ligantes.copy(),ordenes_vecinos_ligantes_pmetal,metales_molecula)
    soluciones_BO = corrige_matricesBOs_ligante(vertices_ligantes,soluciones_factibles,matriz_adj_mol,vecinos_metales,ordenes_factibles,metales_molecula)
    cargas_lingate = corrige_cargas(vertices_ligantes,cargas_factibles)
    
    return soluciones_BO, cargas_lingate, ciclos_lista, cargas_metal

def crea_MO_moleculas(moleculas_ind):
    
    BOS_lista = []
    vert_lista = []
    q_lista = []
    ciclos_lista = [] 
       
    for molecula in moleculas_ind:
               
        vertices_molecula = list(molecula.nodes).copy()
     
        if (len(list(set(metales) & set(vertices_molecula)))!=0):   
            
            posibles_soluciones_molecula,posibles_cargas_ligante, ciclos_ajustados,cargas_metal = crea_molecula_coordinada(molecula,vertices_molecula,0)  
            cargas_ligante = carga_metal_ligante(cargas_metal,vertices_molecula,posibles_cargas_ligante)
            BOS_lista.append(posibles_soluciones_molecula)    
            vert_lista.append(vertices_molecula)
            q_lista.append(cargas_ligante)
            ciclos_lista += ciclos_ajustados
            
        
        else:
            simbolos_molecula = crea_simbolos_individuales (vertices_molecula)
            molecula_enumerada = enumera_grafica(molecula) #Reasigna etiqueta a cada molecula individual
            vertices_molecula_enumerados = list(molecula_enumerada.nodes).copy()
            ma_adj_submolecula = crea_sub_ma(molecula_enumerada)
            ciclos_molecula = detecta_ciclos(molecula_enumerada)
            ciclos_ajustados = ajustar_ciclos (ciclos_molecula,vertices_molecula,[]) #Se reasigna la etiqueta original a los atómos del ciclo         
            posibles_soluciones_molecula, posibles_cargas_molecula = OM.crea_molecula(molecula_enumerada,vertices_molecula_enumerados,simbolos_molecula,ma_adj_submolecula,0,[],[])  
            BOS_lista.append(posibles_soluciones_molecula)
            vert_lista.append(vertices_molecula)
            q_lista.append(posibles_cargas_molecula)
            ciclos_lista = ciclos_lista + ciclos_ajustados
    
    return vert_lista,BOS_lista,ciclos_lista,q_lista

def pares_libres(cargas,BO):
    
    pares_libres = []
    enlaces = list(np.array(BO).sum(axis=1))
    
    for atomo in vertices_completos:           
        numero_atomico = entero_atomo(simbolos_moleculas[atomo]) 
        e_libres = atomic_valence_electrons[numero_atomico] - enlaces[atomo] -cargas[atomo]
        pares_libres.append(e_libres)

    return pares_libres
            
            
            
            
dir_arc = input("Directorio del archivo archivo:") #Ingesa directorio
nombre_arc =(os.path.splitext(os.path.basename(dir_arc))[0])
carga_molecula = int(input("Carga molecular:"))
inicio = default_timer()
#Pide datos
carga = 0 #la carga se supuso, pero puede ser generaizado
lee_xyz(dir_arc) #Lee archivo xyz
matriz_adj = crea_ma(moleculas) #Crea matriz de adjacencia principal
print("Átomos:",len(simbolos_moleculas))
metales = detecta_metales(moleculas,matriz_adj,simbolos_moleculas)
moleculas_individuales = componentes_conexas(moleculas)

if (len(metales) == 0):
    print("Molecula sin metal de coordinación") 
else:
    print("Metales en las moleculas:", [simbolos_moleculas[metal] for metal in metales])
    
vert_lista,BOS_lista,ciclos,q_lista = crea_MO_moleculas(moleculas_individuales)

for i in range(len(BOS_lista[0])):
    nombre_archivo = nombre_arc +"_" +str(i)
    electrones_libres = pares_libres(q_lista[0][i],BOS_lista[0][i])
    crea_archivo_mol(BOS_lista[0][i],nombre_archivo,ciclos,q_lista[0][i],electrones_libres)     

        
final = default_timer()
print("Tiempo de ejecución:",final-inicio, "s")




#Por hacer
# 1. Detectar más de una molecula
# 2. Fragmentar la molecula quitando el metal y corregir enlaces
# 3. Corregir MA_AD De varias moleculas ligadas a Fe