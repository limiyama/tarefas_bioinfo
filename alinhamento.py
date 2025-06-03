import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

#Definição de Parâmetros e Sequências
seq1 = "ACTGGGTCAAC"
seq2 = "ATTGGCCAC"
match = 2
mismatch = -1
gap = -2

#Inicialização da Matriz de Pontuação
def initialize_matrices(seq1, seq2, gap):
  n = len(seq1) + 1
  m = len(seq2) + 1
  score_matrix = np.zeros((n, m))
  for i in range(n):
    score_matrix[i][0] = i * gap
  for j in range(m):
    score_matrix[0][j] = j * gap
  return score_matrix

score_matrix1 = initialize_matrices(seq1, seq2, gap)

print(score_matrix1)

# Preenchimento da Matriz de Pontuação
def fill_score_matrix(seq1, seq2, score_matrix, match, mismatch, gap):
  n, m = score_matrix.shape
  for i in range(1, n):
    for j in range(1, m):
      if seq1[i-1] == seq2[j-1]:
        score_diag = score_matrix[i-1][j-1] + match
      else:
        score_diag = score_matrix[i-1][j-1] + mismatch
      score_up = score_matrix[i-1][j] + gap
      score_left = score_matrix[i][j-1] + gap
      score_matrix[i][j] = max(score_diag, score_up, score_left)
  return score_matrix

score_matrix = fill_score_matrix(seq1, seq2, score_matrix1, match, mismatch, gap)

print(score_matrix)

def traceback(score_matrix, seq1, seq2, gap):
  alignment1, alignment2 = "", ""
  i, j = len(seq1), len(seq2)
  while i > 0 and j > 0:
    score_current = score_matrix[i][j]
    score_diag = score_matrix[i-1][j-1]
    score_up = score_matrix[i-1][j]
    score_left = score_matrix[i][j-1]

    if seq1[i-1] == seq2[j-1]:
      score_diag += match
    else:
      score_diag += mismatch

    if score_current == score_diag:
      alignment1 += seq1[i-1]
      alignment2 += seq2[j-1]
      i -= 1
      j -= 1
    elif score_current == score_up + gap:
      alignment1 += seq1[i-1]
      alignment2 += "-"
      i -= 1
    elif score_current == score_left + gap:
      alignment1 += "-"
      alignment2 += seq2[j-1]
      j -= 1

  while i > 0:
    alignment1 += seq1[i-1]
    alignment2 += "-"
    i -= 1

  while j > 0:
    alignment1 += "-"
    alignment2 += seq2[j-1]
    j -= 1

  return alignment1[::-1], alignment2[::-1]

alignment1, alignment2 = traceback(score_matrix, seq1, seq2, gap)
print("Alinhamento 1:", alignment1)
print("Alinhamento 2:", alignment2)


# Utilizando o BioPython para alinhamento global com os mesmos parâmetros
bio_alignments = pairwise2.align.globalms(seq1, seq2, match, mismatch, gap, gap)

# Exibir o melhor alinhamento obtido pelo BioPython
print("Melhor alinhamento obtido com BioPython:")
print(format_alignment(*bio_alignments[0]))
