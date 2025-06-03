# 4. Limpeza da Estrutura Proteica e 5. Visualização da Estrutura Proteica com Py3Dmol
from Bio.PDB import Select, PDBIO
import py3Dmol

pdbl = PDBList()
parser = PDBParser(QUIET=True)

# link do canal iônico de potássio: https://www.rcsb.org/structure/3UKM
pdbl.retrieve_pdb_file("3UKM", pdir='.', file_format='pdb')

structure = parser.get_structure("3UKM", "pdb{}.ent".format("3UKM".lower()))

# 3. Exploração da Estrutura Proteica:
for model in structure:
      for chain in model:
       print("Cadeia:", chain.id)
       for residue in chain:
          if residue.id[0] != " ":
           print("Composto:", residue.id)
           atoms = list(residue.get_atoms())
           print("Átomos :", atoms)

# 4. Limpeza da Estrutura Proteica e 5. Visualização da Estrutura Proteica com Py3Dmol
class NonWaterSelect(Select):
  def accept_residue(self, residue):
      return residue.get_resname() != "HOH"

io = PDBIO()
io.set_structure(structure)
io.save("clean_structure_potassio.pdb", NonWaterSelect())

with open("clean_structure_potassio.pdb", "r") as pdb_file:
  pdb_content = pdb_file.read()

view = py3Dmol.view(width=800, height=400)
view.addModel(pdb_content, "pdb")
# usando view.setStyle para visualizar melhor as cadeias
view.setStyle({'chain':'A'},{'cartoon': {'color':'red'}})
view.setStyle({'chain':'B'},{'cartoon': {'color':'green'}})
view.setStyle({'chain':'C'},{'cartoon': {'color':'purple'}})
view.setStyle({'chain':'D'},{'cartoon': {'color':'blue'}})
view.zoomTo()
view.show()
