
# ## Proyecto final de optimización 1
# ### Max K cut par k=2 usando semi definite relaxation

# # Parte 1: Importar el conjunto de grafos que se van a probar

import networkx as nx
import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt

def cargar_grafos(path):
    graph = nx.Graph()
    with open(path) as file:
        n_nodes = int(file.readline().split(' ', 1)[0])
        graph.add_nodes_from(range(n_nodes))
        for row in file:
            start, end, weight = [int(e) for e in row.strip('\n').split()]
            graph.add_edge(start - 1, end - 1, weight=weight)

    return graph


# ## Parte 2: Creación de la clase abstracta para la solución de máximo corte

#creación de la clase para resolver el problema de la SDP
#Una clase abstracta importada de la librería abc
import cvxpy as cp
from abc import ABCMeta, abstractmethod


# In[50]:


class AbstractMaxCut(metaclass=ABCMeta):
    def __init__(self,graph):
        self.graph=graph
        self._solution=None

    def get_results(self, item='cut', verbose=False): #Regresa los lazy evaluated max-cut alcanzados, regresa el corte o su valor en la matriz inicial resolviendo el programa SDP
        if self._results is None:
            self.solve(verbose)
        if item not in self._results:

            valid = ', '.join(["'%s'" % key for key in self._results.keys()])
            raise KeyError(
                "No se encuentra la opción que se solicitó" % valid
            )
        return self._results.get(item)
    @abstractmethod
    def solve(self, verbose=True):
        #Resolver el problema BM formulado del max cut usando RTR
        return NotImplemented


def get_partition(vectors):
	#Usa el redondeo Goemans-Williamson
    random = np.random.normal(size=vectors.shape[1])
    # Calcula las probabilidades de partición y redondea el resultado
    random/=la.norm(random,2)
    return np.sign(np.dot(vectors, random))

def get_3_partition(vectors):
    random_1=np.random.normal(size=vectors.shape[1])
    random_2=np.random.normal(size=vectors.shape[1])
    random_3=np.random.normal(size=vectors.shape[1])
    random_1/=la.norm(random_1)
    random_2/=la.norm(random_2)
    random_3/=la.norm(random_3)
    cut=[]
    for i in range(vectors.shape[0]):
        cantidades=[la.norm(np.dot(vectors[i].T,random_1)),la.norm(np.dot(vectors[i].T,random_2)),la.norm(np.dot(vectors[i].T,random_3))]
        cut.append(cantidades.index(max(cantidades))+1)
    return cut


def get_cut_value(graph, partition):
    # Regresa el costo de la partición del grafo
    in_cut = sum(attr['weight'] for u, v, attr in graph.edges(data=True) if partition[u] != partition[v])
    total = .5 * nx.adjacency_matrix(graph).sum()
    return in_cut / total

class MaxCutSDP(AbstractMaxCut):
    "solución del problema de los k cortes, para el caso k=3"
    def __init__(self, graph, solver='scs'):
        super().__init__(graph)
        solver = solver.upper()
        if solver not in cp.installed_solvers():
            raise KeyError("Solver '%s' no instalado." % solver)
        self.solver = getattr(cp, solver)

    def solve(self, verbose=True):
        matrix = self._solve_sdp()
        matrix = nearest_psd(matrix)
        #print("resuelto")
        # Tenemos el corte definido por la matriz
        vectors = np.linalg.cholesky(matrix)
        cut = get_3_partition(vectors)
        print(cut)
        # Tenemos el valor del corte
        value = get_cut_value(self.graph, cut)
        self._results = {'matrix': matrix, 'cut': cut, 'value': value}
        # Optionally be verbose about the results.
        if verbose:
            print(
                "Problema SDP-relaxed max-cut resuelto.\n"
                "Peso total de la solución %f." % value
            )
        #    print(self._results['cut'])
        #    print(self._results['matrix'])

    def _solve_ich_heuristic(self):
        r=0 #Numero actual de particiones
        m=self.graph.number_of_nodes()
        tol=1.7
        matrix=self._solve_sdp()
        matrix=nearest_psd(matrix)
        vectors=la.cholesky(matrix)
        cut=get_3_partitio
        return m

    def lower_bound(self,verbose=True):
        k=self._solve_ich_heuristic()
    def _solve_sdp(self,graph=self.graph):
        """Resuelve el problema SDP del maximo corte.
        regresa la matriz que maximiza <C, 1 - X>
        """
        # Propiedades del grafo a cortar
        n_nodes = len(graph)
        adjacent = nx.adjacency_matrix(graph).toarray() #Convertir la matriz de adyacencia a los
        matrix = cp.Variable((n_nodes, n_nodes), PSD=True)
        tam=cp.multiply(adjacent, 1 - matrix).shape
        suma=0
        for i in range(tam[0]):
            for j in range(tam[1]):
                if i<j:
                    suma+=cp.multiply(adjacent, 1 - matrix)[i][j] #se multipican solamente las entradas donde i es menor que j
        cut=2/3*suma
        constraints=[]
        constraints+=[cp.diag(matrix)==1]
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                if i is not j:
                    constraints+=[matrix[i][j] >=-1/2]
                elif i is j:
                    constraints+=[matrix[i][j]==1]
        problem = cp.Problem(cp.Maximize(cut), constraints)
        # resolver el programa
        problem.solve(getattr(cp, self.solver))
        #problem.solve(solver=cp.MOSEK)
        return matrix.value


def nearest_psd(matrix):
	#la matriz regresa la matriz positiva definida más cercana a la matriz original, se checa que sea semipositiva definida si no se crea una matriz semi positiva definida
    if is_psd(matrix):
        return matrix
    spacing = np.spacing(np.linalg.norm(matrix))
    identity = np.identity(len(matrix))
    k = 1
    while not is_psd(matrix):
        min_eig = np.min(np.real(np.linalg.eigvals(matrix)))
        matrix += identity * (- min_eig * (k ** 2) + spacing)
        k += 1
    return matrix


def is_psd(matrix):
    #Checamos si una matriz es semi definida positiva con la factorización de Cholesky
    try:
        _ = np.linalg.cholesky(matrix)
        return True
    except np.linalg.LinAlgError:
        return False


import progressbar
from time import sleep
bar = progressbar.ProgressBar(maxval=20, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])

# ## Parte de la prueba del algoritmo de Max Cut
file="./G1"
grafo=cargar_grafos(file)
print("Grafo cargado del archivo : {} ".format(file))
sdp = MaxCutSDP(grafo)
print("Resolviendo el problema de corte:")
bar.start()
for i in range(20):
    bar.update(i+1)
    sleep(0.1)
bar.finish()
sdp.lower_bound(grafo)
"""cut = sdp.get_results('cut')

A = [i for i in range(len(cut)) if cut[i] == 1]
B = [i for i in range(len(cut)) if cut[i] == 2]
C = [i for i in range(len(cut)) if cut[i] == 3]
print(len(A), len(B),len(C))


colors = []
edges = grafo.edges()
for u, v in edges:
    if (u in A) and (v in A):
        colors.append('blue')
        continue
    if (u in B) and (v in B):
        colors.append('green')
        continue
    if (u in C) and (v in C):
        colors.append('yellow')
        continue
    #colors.append('red')
    grafo.remove_edge(u, v)

print(colors)

plt.figure()
nx.draw(grafo, with_labels=True, edge_color=colors)  # networkx draw()
plt.draw()
plt.savefig("G_3_cut.png")"""
