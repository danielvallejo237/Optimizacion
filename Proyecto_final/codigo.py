
# ## Proyecto final de optimización 1
# ### Max K cut par k=2 usando semi definite relaxation

# # Parte 1: Importar el conjunto de grafos que se van a probar

import networkx as nx
import numpy as np

def cargar_grafos(path):
    graph = nx.Graph()
    with open(path) as file:
        n_nodes = int(file.readline().split(' ', 1)[0])
        graph.add_nodes_from(range(n_nodes))
        for row in file:
            start, end, weight = [int(e) for e in row.strip('\n').split()]
            graph.add_edge(start - 1, end - 1, weight=weight)
    '''graph.add_nodes_from([1, 2, 3, 4])
    graph.add_weighted_edges_from([(1, 2, 1), (2, 3, 4), (1, 4, 3), (3, 4, 2)])'''
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
                "In valid 'item' keyword: should be one of {{%s}}." % valid
            )
        return self._results.get(item)

    @abstractmethod
    def solve(self, verbose=True):
        #Resolver el problema BM formulado del max cut usando RTR
        return NotImplemented


def get_partition(vectors):
	#Usa el redondeo Goemans-Williamson
    random = np.random.normal(size=vectors.shape[1])
    random /= np.linalg.norm(random, 2)
    # Calcula las probabilidades de partición y redondea el resultado
    return np.sign(np.dot(vectors, random))


def get_cut_value(graph, partition):
    # Regresa el costo de la partición del grafo
    in_cut = sum(attr['weight'] for u, v, attr in graph.edges(data=True) if partition[u] != partition[v])
    total = .5 * nx.adjacency_matrix(graph).sum()
    return in_cut / total

class MaxCutSDP(AbstractMaxCut):
    """Solver Semi-Definite Programming para la solución del Max-Cut Problem, el resultado final es un conjunto de {+1,-1} que indica la pertenencia al grafo """

    def __init__(self, graph, solver='scs'):
        super().__init__(graph)
        solver = solver.upper()
        if solver not in cp.installed_solvers():
            raise KeyError("Solver '%s' no instalado." % solver)
        self.solver = getattr(cp, solver)

    def solve(self, verbose=True):
        matrix = self._solve_sdp()
        matrix = nearest_psd(matrix)
        # Tenemos el corte definido por la matriz
        vectors = np.linalg.cholesky(matrix)
        cut = get_partition(vectors)
        # Tenemos el valor del crte
        value = get_cut_value(self.graph, cut)
        self._results = {'matrix': matrix, 'cut': cut, 'value': value}
        # Optionally be verbose about the results.
        if verbose:
            print(
                "Problema SDP-relaxed max-cut resuelto.\n"
                "Cortes solución de %f Comparten un peso total." % value
            )
            print(self._results['cut'])
            print(self._results['matrix'])

    def _solve_sdp(self):
        """Resuelve el problema SDP del maximo corte.
        regresa la matriz que maximiza <C, 1 - X>
        """
        # Gather properties of the graph to cut.
        n_nodes = len(self.graph)
        adjacent = nx.adjacency_matrix(self.graph).toarray() #Convertir la matriz de adyacencia a los
        matrix = cp.Variable((n_nodes, n_nodes), PSD=True)
        cut = .25 * cp.sum(cp.multiply(adjacent, 1 - matrix))
        problem = cp.Problem(cp.Maximize(cut), [cp.diag(matrix) == 1])
        # resolver el programa
        problem.solve(getattr(cp, self.solver))
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


# ## Parte de la prueba del algoritmo de Max Cut

import random


random.seed(0)
n = 100
p = 0.8
file = "G1"
MAX = 10000
vec = []
edges = []
set = set()
for i in range(n):
    for j in range(n):
        if i == j:
            continue
        if (j, i) in edges:
            print(0)
            continue
        if random.randint(1, MAX)/MAX < p:
            continue
        edges.append((i, j))
        vec.append((i+1, j+1, random.randint(1, 100)))
        set.add(i)
        set.add(j)
with open(file, "w") as f:
    f.write('{} {}'.format(len(set), len(vec)))
    f.write('\n')
    for i in vec:
        f.write('{} {} {}'.format(*i))
        f.write('\n')
    f.close()


grafo=cargar_grafos(file)
print(type(grafo))

import matplotlib.pyplot as plt
nx.draw(grafo, with_labels=True)  # networkx draw()
plt.draw()
plt.savefig("G.png")

# In[53]:


sdp = MaxCutSDP(grafo)
print("Solver")
sdp.solve(grafo)


cut = sdp.get_results('cut')
print(cut)

A = [i for i in range(len(cut)) if cut[i] == 1]
B = [i for i in range(len(cut)) if cut[i] == -1]
print(A, B)


colors = []
edges = grafo.edges()
for u, v in edges:
    if (u in A) and (v in A):
        colors.append('blue')
        continue
    if (u in B) and (v in B):
        colors.append('green')
        continue
    #colors.append('red')
    grafo.remove_edge(u, v)

print(colors)

plt.figure()
nx.draw(grafo, with_labels=True, edge_color=colors)  # networkx draw()
plt.draw()
plt.savefig("G_cut.png")


for i in range(1, n):
    for j in range(1, n):
        pass
