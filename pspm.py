# Parte 2 - Carregamento de Sequência de DNA:
from Bio import Entrez
from Bio import SeqIO
import numpy as np

Entrez.email = "ligia.mizuyama@usp.br"
# Utilizando o SARS-CoV-2
handle = Entrez.efetch(db="nucleotide", id="NC_045512.2", rettype="fasta", retmode="text")
seq_record = SeqIO.read(handle, "fasta")
handle.close()
dna_sequence = str(seq_record.seq)

# Parte 3 - Definição da Matriz de Probabilidade Específica da Posição (PSPM):
pspm = {
    'A': [1.0, 0.0, 0.33, 1.0, 1.0, 0.0],
    'C': [0.0, 1.0, 0.0, 0.0, 0.0, 0.67],
    'G': [0.0, 0.0, 0.67, 0.0, 0.0, 0.33],
    'T': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
}
# Parte 4 - Implementação da Busca de Motivos Usando PSPM:
def calculate_pspm_score(subseq, pspm):
  score = 1.0
  for i, nucleotide in enumerate(subseq):
    if nucleotide in pspm:
      score *= pspm[nucleotide][i]
    else:
      score *= 0.0  # Penaliza nucleotídeos não reconhecidos
  return score

def find_motifs_with_pspm(sequence, pspm, motif_length):
  candidates = []
  for i in range(len(sequence) - motif_length + 1):
    subseq = sequence[i:i+motif_length]
    score = calculate_pspm_score(subseq, pspm)
    candidates.append((i, subseq, score))
  return candidates

motif_length = 6
candidates = find_motifs_with_pspm(dna_sequence, pspm, motif_length)
candidates = sorted(candidates, key=lambda x: x[2], reverse=True)

# Exibir os melhores candidatos
for i in range(10):
  print(f"Posição: {candidates[i][0]}, Motivo: {candidates[i][1]}, Pontuação: {candidates[i][2]}")
