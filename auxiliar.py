# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 13:45:06 2023

@author: Trabajador
"""
import math
import networkx as nx
import random as rd
from collections import Counter



# =============================================================================
def nodesBetweenSourceAndDestination(point1, point2):
    """
    nodesBetweenSourceAndDestination auxiliar function, do not remember what it
    does but is important to plot the L_1 results
    ----
    :param point1: dupla or list 
    :param point2: dupla or list
    """ 
    if (point1[0] != point2[0]
        and point1[1] != point2[1]):
        aux0 = [point2[0], point1[1]]
        aux1 = [point2[0], point2[1]]
        chain = [point1, aux0, aux1, point2]
    elif (point1[0] != point2[0] 
          and point1[1] == point2[1]):
        aux0 = [point2[0], point1[1]]
        chain = [point1, aux0, point2]
    elif (point1[0] == point2[0]
          and point1[1] != point2[1]):
        aux0 = [point1[0], point2[1]]
        chain = [point1, aux0, point2]
    else:
        chain = [point1, point2]
    return chain
# =============================================================================


# =============================================================================
def saveData(dictData, name):
    """
    saveData save result
    ----
    :param dictData: dict 
    :param name: string 
    """ 
    f = open(name + ".txt", "a")
    for key, value in dictData.items():
        f.write(str(key) + ":" + str(value) + "\n")
    f.close()
# =============================================================================
        

# =============================================================================
# Create a dictionary that associates a numeric label to each x, y, z 
# coordinate in order.
def nodesWithCoordinates(minMaxDistX, minMaxDistY):
    """
    nodesWithCoordinates create a dictionary that associates a numeric label 
    to each x, y, z coordinate in order.
    ----
    :param minMaxDistX: int
    :param minMaxDistY: int
    ----
    :return nodes: list
    """ 
    minX = minMaxDistX[0]
    maxX = minMaxDistX[1]
    dx = minMaxDistX[2]
    
    minY = minMaxDistY[0]
    maxY = minMaxDistY[1]
    dy = minMaxDistY[2]
    
    nodes = {}
    counter = 0
    for x in range(minX, maxX, dx):
        for y in range(minY, maxY, dy):
            coordinate = (x, y)
            nodes[counter] = coordinate
            counter = counter + 1 
    return nodes
# =============================================================================


# =============================================================================
# Generates the edges of the graph with elbows.
def generateEdges(nx, ny):
    """
    generateEdges creates a dictionary that associates a numeric label 
    to each x, y, z coordinate in order.
    ----
    :param minMaxDistX: int
    :param minMaxDistY: int
    ----
    :return lista_aristas: list
    """ 
    # Lists to generate the edges.
    # -------------------------------------------------------------------------
    vertical_list = []
    for j in range(nx):
        ape2 = j*ny
        vertical_list.append(ape2)

    horizontal_lists = []
    for i in vertical_list:
        horizontal_list = []
        for j in range(ny):
            ape3 = j + i
            horizontal_list.append(ape3)
        horizontal_lists.append(horizontal_list)
    # -------------------------------------------------------------------------
     
    # Edges        
    lista_aristas = []
          
    # vertical edges
    # -------------------------------------------------------------------------
    for lista in horizontal_lists:
        for a in lista[:-1]:
            ariy = (a, a + 1)
            lista_aristas.append(ariy)
    # -------------------------------------------------------------------------

    # horizontal edges
    # -------------------------------------------------------------------------
    for lista in horizontal_lists[:-1]:
        for a in lista:
            arix = (a, ny + a)
            lista_aristas.append(arix)
    # -------------------------------------------------------------------------
    
    return lista_aristas
# =============================================================================


# =========================================================================
def malladoTriangular(minX, maxX, minY, maxY, d):
    """
    malladoTriangular create a triangular grid
    ----
    :param minX: int or float
    :param maxX: int or float
    :param minY: int or float
    :param maxY: int or float
    :param d: int or float
    ----
    :return H: graph
    """ 
    distancia_horizontal = d
    distancia_vertical = d * math.sqrt(3) / 2

    # n : int
    # The number of columns in the lattice.
    # print(dimensionX, distancia_horizontal, int(dimensionX // distancia_horizontal))
    n = int( (maxX - minX) // distancia_horizontal)

    # m : int
    # The number of rows in the lattice.
    # print(dimensionY, distancia_vertical, int(dimensionY // distancia_vertical))
    m = int( (maxY - minY) // distancia_vertical)
    
    H = nx.empty_graph(0)
    if n == 0 or m == 0:
        return H
    
    rows = range(m + 1)
    cols = range(n + 1)
    
    # Make grid
    H.add_edges_from(((i, j), (i + 1, j))
                      for j in rows 
                      for i in cols[:n])
    H.add_edges_from(((i, j), (i, j + 1)) 
                      for j in rows[:m] 
                      for i in cols)
    # # add diagonals
    H.add_edges_from(((i, j), (i + 1, j + 1)) 
                      for j in rows[1:m:2] 
                      for i in cols[:n])
    H.add_edges_from(((i + 1, j), (i, j + 1)) 
                      for j in rows[:m:2] 
                      for i in cols[:n])
        
    # Add position node attributes
    ii = (i for i in cols for j in rows)
    jj = (j for i in cols for j in rows)
    xx = [0.5 *  d * (j % 2) + d*i 
          for i in cols 
          for j in rows]
    h = d * math.sqrt(3) / 2
    yy = (h * j for i in cols for j in rows)
        
    pos = {(i, j): (x, y) 
           for i, j, x, y in zip(ii, jj, xx, yy) 
           if (i, j) in H}
    nx.set_node_attributes(H, pos, "pos")
    
    list_to_remove = []
    for i, d in H.nodes(data = True):
        if (d['pos'][0] < minX or
            d['pos'][0] > maxX or
            d['pos'][1] < minY or
            d['pos'][1] > maxY):
            list_to_remove.append(i)
            
    H.remove_nodes_from(list_to_remove)
        
    return H
# =========================================================================



# =========================================================================
def formarMalladoFinal(G):
    """
    malladoTriangular create a grid
    ----
    :param G: graph
    ----
    :return H: graph
    """ 
    listaNodos = list(G.nodes())
    nombreNodos = range(len(listaNodos))
    
    dictNombreNodos = dict(zip(listaNodos, nombreNodos))
        
    H = nx.DiGraph()
    
    for key, value in dictNombreNodos.items():
        # pass
        H.add_node(value, coor = G.nodes[key]['pos'], typ = "grid")
        
    for a in G.edges:
        H.add_edge(dictNombreNodos[a[0]],
                   dictNombreNodos[a[1]],
                   weight = rd.randrange(1, 100),
                   typ = "grid")
        H.add_edge(dictNombreNodos[a[1]],
                   dictNombreNodos[a[0]],
                   weight = rd.randrange(1, 100),
                   typ = "grid")
    
    return H
# =========================================================================


  
# ========================================================================= 
def anadirArbol(grid, T):
    """
    anadirArbol add T to the grid
    ----
    :param grid: graph
    ----
    :return T: graph
    """ 
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
# ========================================================================= 


# =============================================================================
def giveMeFatherAndChildren(graph, node):
    """
    giveMeFatherAndChildren auxiliar function, do not remember what is does
    ----
    :param graph: graph
    ----
    :return : list
    """ 
    if [i for i in list(graph.successors(node))] == []:
        return [], []
    else:
        return [node] + [i for i in list(graph.successors(node))], [(node, i) for i in list(graph.successors(node))]
# =============================================================================


# =============================================================================
def gridHanan(points):
    """
    giveMeFatherAndChildren auxiliar function, do not remember what is does
    ----
    :param points: list
    ----
    :return hananGridEdges: list
    :return hananGridPoints: list
    """ 
    # Edges in Hannan Grid
    # -------------------------------------------------------------------------
    xs = []
    ys = []
    
    for i in points:
        xs.append(i[0])
        ys.append(i[1])
        
    xs = sorted(list(dict.fromkeys(xs)))
    ys = sorted(list(dict.fromkeys(ys)))
    
    edgesXs = zip(xs[:-1], xs[1:])
    edgesYs = zip(ys[:-1], ys[1:])
    
    # Por alguna razón que no consigo comprender hay que poner el for de las 
    # aristas antes que los demás, de otra forma los valores de x o z se
    # vuelven 0
    hananGridEdges = []
    for egdeX in edgesXs:
        for y in ys:
            edge = ((egdeX[0], y), (egdeX[1], y))
            # ([egdeX[0], y],[egdeX[1], y])
            hananGridEdges.append(edge)
    
    for edgeY in edgesYs:
       for x in xs:
           edge = ((x, edgeY[0]), (x, edgeY[1]))
           # ([x, edgeY[0]],[x, edgeY[1]])
           hananGridEdges.append(edge)
                
    # Remove duplicates
    hananGridEdges = [entry for idx, entry in enumerate(hananGridEdges)
                      if entry not in hananGridEdges[:idx]]
    # -------------------------------------------------------------------------
    
    # Points in Hannan Grid
    # -------------------------------------------------------------------------
    hananGridPoints = []
    for i in hananGridEdges:
        hananGridPoints.append(i[0])
        hananGridPoints.append(i[1])
    
    # Remove duplicates
    hananGridPoints = [entry for idx, entry in enumerate(hananGridPoints)
                      if entry not in hananGridPoints[:idx]]
    # -------------------------------------------------------------------------
    
    return hananGridEdges, hananGridPoints
# =============================================================================


# =============================================================================
def damePuntosEnEspacios(G, nodosPadreEHijos, hananGridPoints):
    """
    damePuntosEnEspacios auxiliar function, do not remember what is does
    ----
    :param G: graph
    :param nodosPadreEHijos: list
    :param hananGridPoints: list
    ----
    :return dictCoor: dict
    """ 
    dictCoor = {}
    for i in nodosPadreEHijos:
        lista = []
        dictSpacios = G.nodes[i]["space"]
        for point in hananGridPoints:
            for space in dictSpacios.values():
                if (point[0] >= space[0][0]
                    and point[1] >= space[0][1]
                    and point[0] <= space[1][0]
                    and point[1] <= space[1][1]
                    ):
                    lista.append(point)
                    break
        dictCoor[i] = lista
    return dictCoor
# =============================================================================


# =============================================================================
def joinDictionaries(listDictCoor):
    """
    joinDictionaries auxiliar function, join dictionaries
    ----
    :param listDictCoor: list
    ----
    :return dictSol: dict
    """ 
    keys = []
    for dictMini in listDictCoor:
        for i in dictMini.keys():
            keys.append(i)
    keys = list(set(keys))
    keys.sort()
    
    dictSol = {}
    for i in keys:
        listaDeListas = []
        for dictMini in listDictCoor:
            for key, value in dictMini.items():
                if key == i:
                    listaDeListas.append(value)
        newLista = list(set.intersection(*map(set, listaDeListas)))
        dictSol[i] = newLista
    
    return dictSol
# =============================================================================


# =============================================================================
def formarGrafo(G, listDictMallado, listDictCoor):
    """
    formarGrafo auxiliar function, add data to the graph I think
    ----
    :param G: graph
    :param listDictCoor: list
    :param listDictCoor: list
    ----
    :return newG: graph
    """ 
    listaCoor = []
    for i in listDictMallado:
        for j in i["nodes"]:
            listaCoor.append(j)
    listaCoor = list(set(listaCoor))
    listaCoor.sort()
    
    # Enumerate nodes
    dictCoorNum = dict(zip(listaCoor, range(len(listaCoor))))
    dictNumCoor = dict(zip(range(len(listaCoor)), listaCoor))
    
    newG = nx.DiGraph()
    
    for i in listDictMallado:
        for j in i["edges"]:
            newG.add_edge(dictCoorNum[j[0]], dictCoorNum[j[1]], 
                          coor = j, 
                          padreHijo = i['listaPadreHijos'])
            newG.add_edge(dictCoorNum[j[1]], dictCoorNum[j[0]], 
                          coor = j, 
                          padreHijo = i['listaPadreHijos'])
    
    for i in listDictCoor:
        for key, value in i.items():
            for j in value:
                newG.add_node(dictCoorNum[j], coor = j, node = key, typ = "grid")
    
    for i in newG.nodes:
        if "coor" not in newG.nodes[i]:
            nx.set_node_attributes(newG, {i:{'coor': dictNumCoor[i], 'node':"none"}})

    return newG
# =============================================================================
 

# =============================================================================
def malladoTriangularAlternativo(puntos, cajaLimite):
    """
    malladoTriangularAlternativo no sé de qué va lo de alternativo
    ----
    :param puntos: list
    :param cajaLimite: list
    ----
    :return todosLosPuntos: list
    :return todasLasAristas: list
    """     
    rectas = dameTodaLasRectas(puntos, cajaLimite)    
    todasLasAristas = []
    todosLosPuntos = []
    for idx, i in enumerate(rectas):
        restoDeRectas = rectas.copy()
        del restoDeRectas[idx]
        puntosCorte = puntosCorteUnaLineaConElResto(i, restoDeRectas, cajaLimite)
        for j in puntosCorte:
            todosLosPuntos.append(j)
        aristas = dameAristas(puntosCorte)
        for j in aristas:
            todasLasAristas.append(j)
    return todosLosPuntos, todasLasAristas
# =============================================================================


# =============================================================================
def dameTodaLasRectas(puntos, cajaLimite):
    """
    dameTodaLasRectas function auxiliar, no recuerdo bien qué hace
    ----
    :param puntos: list
    :param cajaLimite: list
    ----
    :return rectas: list
    """    
    # =============================================================================
    def dameRectas(punto, cajaLimite):
        """
        dameRectas function auxiliar, no recuerdo bien qué hace
        ----
        :param punto: tupla
        :param cajaLimite: list
        ----
        :return rectas: list
        """    
        dimGrande = max(cajaLimite[1][0], cajaLimite[1][1])

        rectaHorizontal = ((punto[0] - 2*dimGrande, punto[1]),
                           (punto[0] + 2*dimGrande, punto[1]))
        
        rectaDiagonalUno = ((punto[0] - dimGrande , punto[1] - math.sqrt(3)*dimGrande), 
                            (punto[0] + dimGrande, punto[1] + math.sqrt(3)*dimGrande))
        
        rectaDiagonalDos = ((punto[0] - dimGrande , punto[1] + math.sqrt(3)*dimGrande), 
                            (punto[0] + dimGrande, punto[1] - math.sqrt(3)*dimGrande))
        
        rectas = [rectaHorizontal, rectaDiagonalUno, rectaDiagonalDos]
        
        return rectas
    # =============================================================================
    rectas = []
    for punto in puntos:
        algunasRectas = dameRectas(punto, cajaLimite)
        for recta in algunasRectas:
            rectas.append(recta)
    return rectas
# =============================================================================


# =============================================================================
def puntosCorteUnaLineaConElResto(linea, lineas, cajaLimite):
    """
    puntosCorteUnaLineaConElResto gives the cut point between one line and the
    others
    ----
    :param linea: list
    :param lineas: list
    :param cajaLimite: list
    ----
    :return puntosCorte: list
    """  
    # -------------------------------
    def det(a, b):
        """
        Dont know
        """
        return a[0] * b[1] - a[1] * b[0]
    # -------------------------------
    def preguntaLineasCruzar(line1, line2):
        """
        preguntaLineasCruzar tell if two lines cut each or not
        ----
        :param line1: list
        :param line2: list
        ----
        :return bool:
        """  
        xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
        ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])
        div = det(xdiff, ydiff)
        if div == 0:
            return False
        else:
            return True
    # -------------------------------
    def lineIntersection(line1, line2):
        """
        lineIntersection tell where two lines cut
        ----
        :param line1: list
        :param line2: list
        ----
        :return x: int
        :return y: int
        """  
        xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
        ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])
    
        div = det(xdiff, ydiff)
        if div == 0:
           raise Exception('lines do not intersect')
    
        d = (det(*line1), det(*line2))
        x = det(d, xdiff) / div
        y = det(d, ydiff) / div
        return x, y
    # -------------------------------
    puntosCorte = []
    for i in lineas:
        if preguntaLineasCruzar(linea, i):
            punto = lineIntersection(linea, i)
            if (punto[0] >= cajaLimite[0][0]
                and punto[1] >= cajaLimite[0][1]
                and punto[0] <= cajaLimite[1][0]
                and punto[1] <= cajaLimite[1][1]
                ):
                puntosCorte.append(punto)
                    
    puntosCorte = eliminaPuntosRepetidos(puntosCorte)
             
    return puntosCorte
# =============================================================================


# =============================================================================
def dameAristas(puntos):
    """
    dameAristas get edges from points
    ----
    :param puntos: list
    ----
    :return aristas: list
    """  
    puntos.sort()
    aristas = list(zip(puntos[:-1], puntos[1:]))
    return aristas
# =============================================================================


# =============================================================================
def eliminaPuntosRepetidos(puntos):
    """
    eliminaPuntosRepetidos deletes repeted points
    ----
    :param puntos: list
    ----
    :return puntosFinales: list
    """  
    duplicates = []

    for idx in range(len(puntos)):
        if idx not in set(duplicates):
            punto1 = puntos[idx]
            nuevaLista = puntos[:idx] + puntos[idx+1:]
            for indi, i in enumerate(nuevaLista):
                if indi not in set(duplicates):
                    xx = math.isclose(punto1[0], i[0], rel_tol = 0.001)
                    yy = math.isclose(punto1[0], i[0], rel_tol = 0.001)
                    if xx == True and yy == True:
                        duplicates.append(indi)
     
    puntosFinales = [i for idx, i in enumerate(puntos) 
                     if idx not in set(duplicates)]
    
    return puntosFinales
# =============================================================================


# =============================================================================
def giveMeParameters(G):
    """
    giveMeParameters get parameters from the graph
    ----
    :param G: graph
    ----
    :return dictParameters: dict
    """  
    dictParameters = {}
    
    # --------------------
    arcsInTree = []
    arcsInGrid = []

    for i in G.edges:
        if G.edges[i]["typ"] == "tree":
            arcsInTree.append(i)
        elif G.edges[i]["typ"] == "grid":
            arcsInGrid.append(i)
    # --------------------

    # --------------------
    listValveNode = []
    for i in G.nodes:
        if G.nodes[i]["typ"] == "tree":
            for j in G.nodes[i]["nodes"]:
                listValveNode.append((i, j))
    # --------------------

    # --------------------            
    listValves = [i for i in G.nodes if G.nodes[i]["typ"] == "tree"]
    raiz = min(listValves)
    
    nodesInTree = []
    for i in arcsInTree:
        nodesInTree.append(i[0])
        nodesInTree.append(i[1])
    dictRepe = Counter(nodesInTree)
    
    dictValves = {}
    for key, value in dictRepe.items():
        if key == raiz:
            dictValves[key] = "root"
        elif value >= 2:
            dictValves[key] = "valve"
        else:
            dictValves[key] = "leaf"
    # --------------------            

    # combinacionesMismaFuente = []
    # for a in arcsInTree:
    #     fuentes = [a[0]]
    #     combinacion = [a]
    #     for b in arcsInTree:
    #         if b[0] not in set(fuentes):
    #             combinacion.append(b)
    #             fuentes.append(b[0])
    #     combinacionesMismaFuente.append(combinacion)    
          
    dictParameters["arcsInTree"] = arcsInTree
    dictParameters["arcsInGrid"] = arcsInGrid
    dictParameters["listValveNode"] = listValveNode
    dictParameters["dictValves"] = dictValves
    # dictParameters["combinacionesMismaFuente"] = combinacionesMismaFuente
    
    return dictParameters
# =============================================================================











