from codigo import *
import networkx as nx
import cvxpy as cp
import numpy as np
import random

class MaxCutIch(AbstractMaxCut):

    def __init__(self, graph, solver='scs'):
        super().__init__(graph)
        solver = solver.upper()
        if solver not in cp.installed_solvers():
            raise KeyError("Solver '%s' is not installed." % solver)
        self.solver = getattr(cp, solver)


    def solve(self,grafo, k, tol):
        self.graph = grafo
        parts = self._solve_aux(self.graph)
        aux_grafo = self.cambiar_grafo(parts, self.graph)
        print("Hecho")
        #while len(parts) > k:
        #    aux_grafo = self.cambiar_grafo(parts, self.graph)
            #parts = self._solve_aux(aux_grafo, k, tol, len(parts))
        return parts

    def cambiar_grafo(self, part, grafo):
        adj = nx.adjacency_matrix(grafo)
        n_nodes = len(part)
        nueva_adj = np.zeros((n_nodes, n_nodes))
        edges = grafo.edges()
        for u, v in edges:
            ind1 = self._search_part(part, u)
            ind2 = self._search_part(part, v)
            nueva_adj[ind1][ind2] += edges[u, v]['weight']

        n_aristas = len(np.nonzero(nueva_adj[0]))
        with open('./temp', 'w') as f:
            f.write('{} {}'.format(int(len(part)), int(n_aristas)))
            f.write('\n')
            for i in range(len(part)):
                for j in range(len(part)):
                    if nueva_adj[i][j] == 0:
                        continue
                    f.write('{} {} {}'.format(i+1, j+1, int(nueva_adj[i][j])))
                    f.write('\n')
            f.close()
        return cargar_grafos('./temp')

    def _solve_aux(self,grafo, k=2, tol=1.7,r=0):
        #r = 0
        n_nodes = len(grafo)
        m = n_nodes
        tol = tol
        adjacent = nx.adjacency_matrix(grafo).toarray()
        matrix = cp.Variable((n_nodes, n_nodes), PSD=True)
        tam=cp.multiply(adjacent, 1 - matrix).shape
        suma=0
        for i in range(tam[0]):
            for j in range(tam[1]):
                if i<j:
                    suma+=cp.multiply(adjacent, 1 - matrix)[i][j] #se multipican solamente las entradas donde i es menor que j
        cut=(k-1)/k*suma
        constraints=[]
        constraints+=[cp.diag(matrix)==1]
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                if i is not j:
                    constraints+=[matrix[i][j] >=-1/(k-1)]
                elif i is j:
                    constraints+=[matrix[i][j]==1]
        #cut = (k-1)/k * cp.sum(cp.multiply(adjacent, 1 - matrix))
        problem = cp.Problem(cp.Maximize(cut), constraints)
        problem.solve(getattr(cp, self.solver))
        print("Problema resuelto ...")
        T = []
        for i in range(n_nodes):
            for j in range(i, n_nodes):
                if i == j:
                    continue
                for h in range(j, n_nodes):
                    if h == i or h ==j:
                        continue
                    val = matrix.value[i][j]
                    val += matrix.value[i][h]
                    val += matrix.value[j][h]
                    T.append((val, i, j, h))
        T = sorted(T, reverse=True) ##checar reverse
        parts = [] #vector de sets
        for valor in T:
            if valor[0] > tol:
                print(parts)
                ind1 = self._search_part(parts, valor[1])
                ind2 = self._search_part(parts, valor[2])
                ind3 = self._search_part(parts, valor[3])

                print(ind1, ind2, ind3)

                if ind1 == -1 and ind2 == -1 and ind3 == -1:
                    parts.append(set([valor[1], valor[2], valor[3]]))
                    continue

                if ind1 == -1:
                    parts.append(set([valor[1]]))
                    ind1 = len(parts)-1

                if ind2 == -1:
                    parts.append(set([valor[2]]))
                    ind2 = len(parts)-1
                    print('@@', parts)

                if ind3 == -1:
                    parts.append(set([valor[3]]))
                    ind3 = len(parts)-1
                #paso 7
                parts = self._set_handle(parts, ind1, ind2, ind3)

        return parts



    def _set_handle(self, parts, ind1, ind2, ind3):
        if ind1 != ind2:
            parts[min(ind1, ind2)] |= parts[max(ind1, ind2)]
            if ind3 == ind1 or ind3 == ind2:
                parts.pop(max(ind1, ind2))
                return parts
            parts[min(ind1, ind2)] |= parts[ind3]
            parts[min(ind1, ind2, ind3)] = parts[min(ind1, ind2)]
            aux = sorted([ind1, ind2, ind3], reverse=True) #checar el orden de mayor a menor
            print(aux)
            parts.pop(aux[2])
            parts.pop(aux[1])
            return parts
        if ind1 == ind2:
            if ind1 == ind3:
                return parts
            parts[min(ind1, ind3)] |= parts[max(ind1, ind3)]
            parts.pop(max(ind1, ind3))
            return parts



    def _search_part(self, parts, val):
        for i in range(len(parts)):
            if val in parts[i]:
                return i
        return -1








############################################################
file = "./G1"


grafo=cargar_grafos(file)
print(type(grafo))

import matplotlib.pyplot as plt
nx.draw(grafo, with_labels=True)  # networkx draw()
plt.draw()
plt.savefig("G.png")

# In[53]:


sdp = MaxCutIch(grafo)
print("Solver")
print(sdp.solve(grafo, 3, 1.7))
