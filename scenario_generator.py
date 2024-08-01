# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 16:51:58 2024

@author: Trabajador
"""
import networkx as nx
import random as rd
import auxiliar as aux

# =============================================================================
def iterador(node, lista, tree):
    """
    iterador auxiliar function, do not remember what is does
    ----
    :param node: int
    :param lista: list
    :param tree: graph
    ----
    :return lista: lista
    """ 
    if len(lista) != len(tree.nodes):
        for node in tree.adj[node].keys():
            if node not in set(lista):
                lista.append(node)
                iterador(node, lista, tree)
    return lista
# =============================================================================


# =============================================================================
def generateBasicTree(nodesInTree):
    """
    generateBasicTree auxiliar function, creates a random tree of nodesInTree number of nodes.
    The nodes are properly enumerated and the tree is directed.
    ----
    :param nodesInTree: int
    ----
    :return newTree: graph
    """ 
    # Generate random Tree
    tree = nx.random_tree(nodesInTree)
    
    # Enumerate nodes tree properly
    lista = [0]
    node = 0
    lista = iterador(node, lista, tree)
    dictRelacion = dict(zip(lista, list(tree.nodes)))
   
    # Generate new graph directed
    newTree = nx.DiGraph()
    for i in tree.nodes:
        newTree.add_node(dictRelacion[i]) 
    for i in tree.edges:
        if dictRelacion[i[0]] < dictRelacion[i[1]]:
            newTree.add_edge(dictRelacion[i[0]], dictRelacion[i[1]])
        else:
            newTree.add_edge(dictRelacion[i[1]], dictRelacion[i[0]])
    
    return newTree
# =============================================================================


# =============================================================================
def generateTree(n, b, s):
    """
    generateTree creates a random tree of n nodes, b spaces for node and the spaces
    are the (2*s)^2 u^2 dimensions. The nodes are properly enumerated and the tree
    is directed.
    ----
    :param n: int
    :param b: int
    :param s: int
    ----
    :return T: graph
    """ 
    # Generate tree as directed graph
    T = generateBasicTree(n)

    # Generate attributes for each node
    xMax = 2000
    yMax = 2000
    dictAttributes = {}
    for i in T.nodes:
        # positions = []
        spaces = {}
        centers = {}
        epsilons = {}
        for j in range(b):
            # Position for the node
            position = (rd.randint(0, xMax - s), rd.randint(3, yMax - s))
            # positions.append(position)
            # Box space for the node
            space = ((position[0], position[1]),
                     (position[0] + 2*s, position[1] + 2*s))
            spaces[j] = space
            center = (position[0] + s, position[1] + s)
            centers[j] = center
            epsilon = (s, s)
            epsilons[j] = epsilon
        # Attribute
        dictAttribute = {"center": centers,
                              "epsilon": epsilons,
                              "space": spaces}
        dictAttributes[i] = dictAttribute
    
    # Set attributes for each node
    nx.set_node_attributes(T, dictAttributes)
    
    return T
# =============================================================================


# =============================================================================
def modifyGraph(G):
    """
    modifyGraph modify the graph so it can be solved the continous steiner tree
    model
    ----
    :param G: graph
    ----
    :return G: graph
    """ 
    listNodes = list(G.nodes)
    newArcs = []
    for i in G.nodes:
        sucesores = [i for i in G.successors(i)]
        if len(sucesores) > 1:
            # Número de nodos de steiner por crear
            nSteiner = len(sucesores) + 1 - 2
            # Lista de nuevos nodos steiner
            listNSteiner = range(len(listNodes), len(listNodes) + nSteiner)
            for k in listNSteiner:
                listNodes.append(k)
                newArcs.append((i, k))
                for l in listNSteiner:
                    if k != l:
                        newArcs.append((k, l))
                for l in sucesores:
                    newArcs.append((k, l))
            for j in sucesores:
                G.remove_edge(i, j)
    
    for i in listNodes:
        if i not in set(G.nodes):
            G.add_node(i)
    
    for i in newArcs:
        G.add_edge(i[0], i[1])
        
    return G
# =============================================================================

  
# =============================================================================
def createDiscreteSquareScenario(minX, maxX, minY, maxY, T, dx, dy):
    # Generate Grid
    # -------------------------------------------------------------------------
    # xs = []
    # ys = []
    # for i in T.nodes:
    #     for j in T.nodes[i]["space"].values():
    #         for k in j:
    #             xs.append(k[0])
    #             ys.append(k[1])
    
    # minX = min(xs)
    # maxX = max(xs)
    # minY = min(ys)
    # maxY = max(ys)
    

    # while (maxX - minX) % dx != 0:
    #     maxX = maxX + 1
        
    # while (maxY - minY) % dy != 0:
    #     maxY = maxY + 1

    # Nodes
    nodes = aux.nodesWithCoordinates((minX, maxX + dx, dx), (minY, maxY + dy, dy))
    
    # Edges
    edges = aux.generateEdges(len(range(minX, maxX + dx, dx)), len(range(minY, maxY + dy, dy)))
    
    edges = list(set(edges))
    
    # Create Grid
    grid = nx.DiGraph()
    
    for key, coordinates in nodes.items():
        grid.add_node(key, coor = coordinates, typ = 'grid')
        
    for edge in edges:
        grid.add_edge(edge[0],
                      edge[1],
                      weight = rd.randrange(100),
                      typ = 'grid')
        grid.add_edge(edge[1],
                      edge[0],
                      weight = rd.randrange(100),
                      typ = 'grid')
    # -------------------------------------------------------------------------

    # Add tree nodes to the grid
    # -------------------------------------------------------------------------
    count = grid.number_of_nodes()
    dictAyudaNodes = {}
    for i in T.nodes:
        number = grid.number_of_nodes() + i
        spaces = T.nodes[i]['space']
        
        nodes = []
        for space in spaces.values():
            for i in grid.nodes:
                coor = grid.nodes[i]["coor"]
                if (coor[0] >= space[0][0]
                    and coor[1] >= space[0][1]
                    and coor[0] <= space[1][0]
                    and coor[1] <= space[1][1]
                    ):
                    nodes.append(i)
        dictAyudaNodes[number] = {"space": spaces,
                                  "nodes": nodes,
                                  "typ": 'tree'}
        
    for number, value in dictAyudaNodes.items():
        grid.add_node(number, space = value["space"],
                      nodes = value["nodes"],
                      typ = value["typ"])
    # -------------------------------------------------------------------------
    
    # Add tree edges to the grid
    # -------------------------------------------------------------------------
    # dictAyudaEdges = {}
    for edge in T.edges():
        new_edge = (edge[0]+ count, edge[1]+ count)
        aa = T.nodes[edge[0]]["center"][0]
        partida = (aa[0], aa[1] - T.nodes[edge[0]]["epsilon"][0][1])
        bb = T.nodes[edge[1]]["center"][0]
        final =  (bb[0], bb[1] + T.nodes[edge[1]]["epsilon"][0][1])
        grid.add_edge(new_edge[0],
                      new_edge[1],
                      link = (partida, final),
                      typ = 'tree')
    # -------------------------------------------------------------------------

    return grid
# =============================================================================


# ========================================================================= 
def dameMalladoTriangularCompleto(minX, maxX, minY, maxY, d, T):
    G = aux.malladoTriangular(minX, maxX, minY, maxY, d)
            
    G = aux.formarMalladoFinal(G)
    
    G = aux.anadirArbol(G, T)

    return G
# ========================================================================= 


# =============================================================================
def dameVariosMalladosCuadrados(T):

    listDictCoor = []
    listDictMallado = []
    for i in T.nodes:
        nodosPadreEHijos, listaPadreHijos = aux.giveMeFatherAndChildren(T, i)
        if nodosPadreEHijos != []:
            puntos = []
            for j in nodosPadreEHijos:
                for key, center in T.nodes[j]["center"].items():
                    punto0 = (center[0] - T.nodes[j]["epsilon"][key][0],
                              center[1] - T.nodes[j]["epsilon"][key][1])
                    punto1 = (center[0] + T.nodes[j]["epsilon"][key][0],
                              center[1] - T.nodes[j]["epsilon"][key][1])
                    punto2 = (center[0] - T.nodes[j]["epsilon"][key][0],
                              center[1] + T.nodes[j]["epsilon"][key][1])
                    punto3 = (center[0] + T.nodes[j]["epsilon"][key][0],
                              center[1] + T.nodes[j]["epsilon"][key][1])
                    puntos.append(punto0)
                    puntos.append(punto1)
                    puntos.append(punto2)
                    puntos.append(punto3)
                    
            hananGridEdges, hananGridPoints = aux.gridHanan(puntos)
            
            dictMallado = {"edges": hananGridEdges, "nodes": hananGridPoints,
                           "listaPadreHijos": listaPadreHijos }
            listDictMallado.append(dictMallado)
            
            dictCoor = aux.damePuntosEnEspacios(T, nodosPadreEHijos, hananGridPoints)
            listDictCoor.append(dictCoor)
    
    dictCoor = aux.joinDictionaries(listDictCoor)
        
    grids = aux.formarGrafo(T, listDictMallado, listDictCoor)
        
    return grids
# =============================================================================


# =============================================================================
def dameVariosMalladosTriangulares(T):

    listDictCoor = []
    listDictMallado = []
    for i in T.nodes:
        nodosPadreEHijos, listaPadreHijos = aux.giveMeFatherAndChildren(T, i)
        if nodosPadreEHijos != []:
            puntos = []
            for j in nodosPadreEHijos:
                for key, center in T.nodes[j]["center"].items():
                    punto0 = (center[0] - T.nodes[j]["epsilon"][key][0],
                              center[1] - T.nodes[j]["epsilon"][key][1])
                    punto1 = (center[0] + T.nodes[j]["epsilon"][key][0],
                              center[1] - T.nodes[j]["epsilon"][key][1])
                    punto2 = (center[0] - T.nodes[j]["epsilon"][key][0],
                              center[1] + T.nodes[j]["epsilon"][key][1])
                    punto3 = (center[0] + T.nodes[j]["epsilon"][key][0],
                              center[1] + T.nodes[j]["epsilon"][key][1])
                    # punto4 = center
                    puntos.append(punto0)
                    puntos.append(punto1)
                    puntos.append(punto2)
                    puntos.append(punto3)
                    # puntos.append(punto4)
                    
            cajaLimite = ((-99999, -99999), (99999, 99999))
            hananGridPoints, hananGridEdges = aux.malladoTriangularAlternativo(puntos, cajaLimite)
            dictMallado = {"edges": hananGridEdges, "nodes": hananGridPoints,
                           "listaPadreHijos": listaPadreHijos }
            listDictMallado.append(dictMallado)
            
            dictCoor = aux.damePuntosEnEspacios(T, nodosPadreEHijos, hananGridPoints)
            listDictCoor.append(dictCoor)
    
    dictCoor = aux.joinDictionaries(listDictCoor)
        
    grids = aux.formarGrafo(T, listDictMallado, listDictCoor)
        
    return grids
# =============================================================================









