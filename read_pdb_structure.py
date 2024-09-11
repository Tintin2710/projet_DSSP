from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

class ReadPdbStructure:
    """Represent and extract the main chain atoms (N, CA, C, O) in the PDB structure."""

    def __init__(self, pdb_filename):
        self.pdb_filename = pdb_filename
        self.structure = self._read_pdb_structure()
        self.main_chain_atoms = self._get_main_chain_atoms()

    def _read_pdb_structure(self):
        """Lire le fichier PDB et renvoyer la structure"""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", self.pdb_filename)
        return structure
    
    def _get_main_chain_atoms(self):
        """Extract the main chain atomic coordinates of each residue (N, CA, C, O)"""
        main_chain_atoms = []
        previous_res_id = None
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue) and residue.has_id('N') and residue.has_id('CA') and residue.has_id('C') and residue.has_id('O'):
                        res_id = residue.get_id()[1]
                        chain_id = chain.get_id()
                        if previous_res_id is not None and abs(res_id - previous_res_id) > 1:
                            print(f"Non-continuous residue IDs: previous_res_id={previous_res_id}, current_res_id={res_id}, chain={chain_id}")    
                        previous_res_id = res_id

                        main_chain_atoms.append({
                            'residue': residue.get_resname(),
                            'res_id': res_id,
                            'chain': chain_id,
                            'N': residue['N'].get_coord(),
                            'CA': residue['CA'].get_coord(),
                            'C': residue['C'].get_coord(),
                            'O': residue['O'].get_coord()
                        })
        return main_chain_atoms
    
    def get_main_chain_atoms(self):
        """Obtenir des informations atomiques sur la chaîne principale"""
        return self.main_chain_atoms
        