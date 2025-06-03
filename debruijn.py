# Parte 1 - Configuração do Ambiente
# Parte 2 - Preparação de Dados
# Aqui iremos importar as bibliotecas necessárias para manipular os dados
import copy
import string
import random
from itertools import product
from collections import defaultdict
import numpy as np

# Função para gerar aleatoriamente uma sequência de DNA
def generate_dna(length):
    return ''.join(random.choices('ACGT', k=length))

# Função para gerar reads a partir do DNA, com o número e comprimento de reads específicos
def generate_reads(genome, num_reads, read_length):
    reads = []
    for _ in range(num_reads):
        start = random.randint(0, len(genome) - read_length)
        reads.append(genome[start:start + read_length])
    return reads

# Parâmetros de simulação
genome_size = 50000  # tamanho do genoma
read_length = 100  # comprimento de cada read
num_reads = 1000  # número de reads

# Gerando o DNA e os reads
genome = generate_dna(genome_size)
reads = generate_reads(genome, num_reads, read_length)

# Salvando a sequência original em um .txt para comparar dpeois
with open("sequencia_original.txt", "w") as f:
  f.write(genome)

# Parte 3 - Implementação do Grafo de De Bruijn
# Função para construir o grafo de Bruijin
def build_de_bruijn_graph(reads, k):
    edges = defaultdict(list)
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            left = kmer[:-1]
            right = kmer[1:]
            if left not in edges:
                edges[left] = []
            edges[left].append(right)
    return edges

# Construir grafo com k-mer de tamanho 20
de_bruijn_graph = build_de_bruijn_graph(reads, 20)
print("Visualização dos 5 primeiros nodos gerados: ", list(de_bruijn_graph.keys())[:5])

# Parte 4 - Montagem de Sequências
# Funções para encontrar um caminho euleriano e montar a sequência
def edges(graph):
    for node in graph:
        for target in graph[node]:
            yield (node, target)

def follow_tour(tour, graph):
    edges_ = list(edges(graph))
    for start, end in zip(tour, tour[1:]):
        try:
            edges_.remove((start, end))
        except:
            return False

    if edges_:
        return False
    else:
        return True

def check_tour(start, graph):
    our_tour = tour(start, graph)
    valid_tour = follow_tour(our_tour, graph)
    return valid_tour, "".join(s[0:-1] for s in our_tour) + our_tour[-1][1:]

def tour(start_node, graph):
    graph = copy.deepcopy(graph)
    return _tour(start_node, graph)

def _tour(start_node, graph, end=None):
    tour = [start_node]
    finish_on = end if end is not None else start_node
    while True:
        options = graph[tour[-1]]
        if not options:
            break

        tour.append(options.pop())
        if tour[-1] == finish_on:
            break

    offset = 0
    for n,step in enumerate(tour[:]):
        options = graph[step]
        if options:
            t = _tour(options.pop(), graph, step)
            n += offset
            tour = tour[:n+1] + t + tour[n+1:]
            offset += len(t)

    return tour

# Seleciona um nodo para iniciar o caminho euleriano
start_node = list(de_bruijn_graph.keys())[0]
print("Nodo escolhido: ", start_node)
valid, genome_assembly = check_tour(start_node, de_bruijn_graph)

# Mostra o tamanho das sequências para fins de comparação
print("\nTamanho do genoma montado: ", len(genome_assembly))
print("Tamanho do genoma original: ", len(genome))

# Salvando a sequência montada em um .txt para comparar dpeois
with open("sequencia_montada.txt", "w") as f:
  f.write(genome_assembly)
