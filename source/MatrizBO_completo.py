# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 16:10:23 2022

@author: Uriel Alejandro
"""

#Importación librelias

import copy
import networkx as nx
import numpy as np 
from collections import defaultdict
import itertools

#Variables

#Variable de bases de datos

#Lista de valencias y atomos de valencia 

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

metales_transicion = [21,22,24,25,26,27,41,42,44,45,46,73,74,75,76,77,78]
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
atomic_valence[24] = [8,9,10,11,12]
atomic_valence[25] = [7,8,9,10,11]
atomic_valence[26] = [5,6,7,8,9,10]
atomic_valence[27] = [5,6,7,8,9]
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
atomic_valence_electrons[24] = 6
atomic_valence_electrons[25] = 7
atomic_valence_electrons[26] = 8
atomic_valence_electrons[27] = 9
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
def simbolo_atomo(atom): #convert integer atom to string atom

    global __ATOM_LIST__
    atom = __ATOM_LIST__[atom - 1]
    
    return atom


def entero_atomo(atom): # convert str atom to integer atom
    
    global __ATOM_LIST__
    atom = atom.lower()
    
    return __ATOM_LIST__.index(atom) + 1



def numero_de_enlaces(n):

    conect = 0
        
    for i in vertices:
        conect += matriz_adj[i][n]
    
    for i in vecinos_carbono:
        if i == n:
            conect += 1
    for i in range(len(vecinos_metal)):
       if vecinos_metal[i] == n:
           conect += ordenes_vecinos_metal[i]
       
    return conect


def valencias_buenas(BO,valencias):#

    numero_enlaces = np.array(BO).sum(axis=1)        
    for i in vecinos_carbono:
        numero_enlaces[i] += 1
    for i in range(len(vecinos_metal)):
        numero_enlaces[vecinos_metal[i]] += ordenes_vecinos_metal[i]
    for i in range(len(valencias)): 
        if numero_enlaces[i] > valencias[i]:
            return False

    return True

def crea_enlaces(UA, AC):
    
    enlaces = []
    
    for k, i in enumerate(UA):
        for j in UA[k + 1:]:
          
            if AC[i][j]==1 or AC[i][j]==2:
                enlaces.append(tuple(sorted([i, j])))
                
    return enlaces

def obten_pares_UA(UA,m_ad):
    
    enlaces = crea_enlaces(UA,m_ad)

    if len(enlaces) == 0:
        
        return[()]
    
    G = nx.Graph()
    G.add_edges_from(enlaces)
    UA_pares = [list(nx.max_weight_matching(G))]
    
    return UA_pares
    
def calcula_carga_atomica(atomo,no_enlaces):
    
    if atomo == 1:
        carga = 1 - no_enlaces
    elif atomo == 5:
        carga = 3 - no_enlaces
    elif atomo == 15 and no_enlaces == 5:
        carga = 0
    elif atomo == 16 and no_enlaces == 6:
        carga = 0
    elif atomo in metales_transicion:
        carga = atomic_valence_electrons[atomo] - 18 + no_enlaces
    else:
        carga = atomic_valence_electrons[atomo] - 8 + no_enlaces

    return carga

def carga_buena(BO,du,ua):
    # Carga total
    Q = 0
    # Carga fragmentos
    no_enlaces = list(np.array(BO).sum(axis=1))
    for i in vecinos_carbono:
        no_enlaces[i] += 1
    for i in range(len(vecinos_metal)):
        no_enlaces[vecinos_metal[i]] += ordenes_vecinos_metal[i]       
    
    for i in vertices:
        simbolo_atomo = simbolos_molecula[i]
        atomo = entero_atomo(simbolo_atomo)
        carga = calcula_carga_atomica(atomo,no_enlaces[i])
        Q += carga
        if atomo==6:
            numero_enlaces_unicos_C = list(np.array(BO)[i, :]).count(1)
            if numero_enlaces_unicos_C == 2 and no_enlaces[i] == 2:
                Q += 1
                carga = 2
            if numero_enlaces_unicos_C == 3 and Q + 1 < carga: #Cambiar por CARGA
                Q += 2
                carga = 1                
    return (0 == Q) #Cambiar por CARGA

def carga_final(BO,du,ua):

    q_lista = []
    no_enlaces = list(np.array(BO).sum(axis=1))

    for i in vecinos_carbono:
        no_enlaces[i] += 1
    for i in range(len(vecinos_metal)):
        no_enlaces[vecinos_metal[i]] += ordenes_vecinos_metal[i]   
    
    for i in vertices:
        simbolo_atomo = simbolos_molecula[i]
        atomo = entero_atomo(simbolo_atomo)
        carga = calcula_carga_atomica(atomo,no_enlaces[i])         
        q_lista.append(carga)

    return q_lista


def calcula_ui(fila,ver,vale,sim_mol,m_ad,enl_c):
    
    enlaces = numero_de_enlaces(fila)
    ui = vale[fila]-enlaces
    
    return ui

def crea_matriz_vs():
    
    posibles_valencias = []
    enlaces = []
    for i in vertices:
        numero_atomico = entero_atomo(simbolos_molecula[i])
        enlaces_atomo = numero_de_enlaces(i)
        pos_val_atom=[]
        enlaces.append(enlaces_atomo)
        for j in atomic_valence[numero_atomico]:
            
            if(enlaces_atomo<=j):
                pos_val_atom.append(j)
                
        
        if(pos_val_atom==[]):
            
            return [],[],False
        
        else:
            posibles_valencias.append(pos_val_atom)
    
    matriz_vs = list(itertools.product(*posibles_valencias)).copy()
   
    return matriz_vs, enlaces,True      

def BO_buena(BO,du,ua,valencia):
    
    
    verifica_valencias = valencias_buenas(np.array(BO),valencia)#
    if verifica_valencias:
        return False
    
    numero_enlaces_BO = (np.array(BO)).sum()
    for i in vecinos_carbono:
        numero_enlaces_BO[i] += 1
    for i in range(len(vecinos_metal)):
        numero_enlaces_BO[vecinos_metal[i]] += ordenes_vecinos_metal[i]
    verifica_suma = numero_enlaces_BO - (np.array(matriz_adj)).sum() == sum(du)#
    verifica_carga = carga_buena(np.array(BO),du,ua)#
    if verifica_suma and verifica_carga:#
        return True
    
    return False

def calcular_DU_UA(valencias,valencias_maximas):
    
    ua = []
    du = []
    
    for i, (valencia_maxima, valencia) in enumerate(zip(valencias_maximas, valencias)):
        if valencia_maxima - valencia > 0:

            ua.append(i)
            du.append(valencia_maxima - valencia)
    

    return ua, du

def construye_BO(ua,du,valencia,par):
    
    BO = np.array(matriz_adj).copy()
    DU_guardar = []

    while DU_guardar != du :
        
        for i, j in par:
            BO[i,j] += 1
            BO[j,i] += 1
        
        valencias_BO = list(BO.sum(axis=1))
        DU_guardar = copy.copy(du)
        for i in vecinos_carbono:
            valencias_BO[i] += 1
        for i in range(len(vecinos_metal)):
            valencias_BO[vecinos_metal[i]] += ordenes_vecinos_metal[i]
        
        ua, du = calcular_DU_UA(valencias_BO,valencia)
        par = obten_pares_UA(ua,matriz_adj)[0]
        
    return BO

def crea_q_lista(valencias,sim_mol):
    
    q_lista = []
    
    for i in range(len(valencias)):
        sim_atomo = sim_mol[i]
        atomo = entero_atomo(sim_atomo)
        q = calcula_carga_atomica(atomo,valencias[i])
        q_lista.append(q)
    
    return q_lista
                                                 
def valencias_suficientes(BO_posible,valencias):
        
    enlaces = list(np.array(BO_posible).sum(axis=1))
    for i in vecinos_carbono:
        enlaces[i] += 1
    for i in range(len(vecinos_metal)):
        enlaces[vecinos_metal[i]] += ordenes_vecinos_metal[i]
    
    for i in range(len(enlaces)):
        numero_atomico = entero_atomo(simbolos_molecula[i])
        if (enlaces[i]<min(atomic_valence[numero_atomico])):
            return False
    return True        

def neutralidad_ligante(cargas):
    
    carga_ligante = abs(np.array(cargas).sum())-np.array(ordenes_vecinos_metal).sum()
    if (carga_ligante != 0):
            return False
    
    return True

def repeticion_solucion (lista_soluciones, nueva_solucion):
    
    for solucion in lista_soluciones:
        if (np.array_equal(solucion,nueva_solucion) == True):
            return True
    return False
def crea_BO(vertices_in,simbolos_molecula_in,matriz_adj_in,carga_in,vecinos_carbono_in,ordenes_in,vecinos_metal_in):
    
    global simbolos_molecula #Declaramos globales los datos recibidos
    global vertices
    global matriz_adj
    global vecinos_carbono
    global ordenes_vecinos_metal
    global vecinos_metal
    
    soluciones = []
    cargas = []
    simbolos_molecula = simbolos_molecula_in.copy() #Declaramos globales los datos recibidos
    vertices = vertices_in.copy()
    matriz_adj = matriz_adj_in.copy()
    vecinos_carbono = vecinos_carbono_in.copy()
    ordenes_vecinos_metal = ordenes_in.copy()
    vecinos_metal = vecinos_metal_in.copy()
    
    
    matriz_pos_val,valencias_mat_adj, existen_posibles_soluciones= crea_matriz_vs()#
    if existen_posibles_soluciones == False:
        return [],[]
    mejor_BO = matriz_adj.copy()#

    for valencias in matriz_pos_val:#
        
        
        UA,DU = calcular_DU_UA(valencias_mat_adj,valencias)#
        pares_UA= obten_pares_UA(UA,matriz_adj)

        for pares in pares_UA:
            
            BO_posible = construye_BO(UA,DU,valencias,pares)
            #ligante_neutro = neutralidad_ligante(cargas_finales)
            verifica_valencias = valencias_buenas(BO_posible,valencias)
            valencias_suf_grandes = valencias_suficientes(BO_posible,valencias)
            
            if np.array(mejor_BO).sum()<=np.array(BO_posible).sum() and verifica_valencias and valencias_suf_grandes:
                mejor_BO = BO_posible
                if repeticion_solucion(soluciones,mejor_BO) == False:
                    soluciones.append(mejor_BO.copy())
                    cargas.append(carga_final(mejor_BO.copy(), DU, UA))
             
    return soluciones,cargas

