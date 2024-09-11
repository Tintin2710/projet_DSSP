from Bio.PDB import PDBParser, DSSP

class Compare:
    """Classes for comparing custom secondary structures with DSSP allocation results"""

    def __init__(self, pdb_filename, predicted_structure):
        self.pdb_filename = pdb_filename
        self.predicted_structure = predicted_structure
        self.true_structure = self._dssp_structure()

    def _dssp_structure(self):
        """Extract the secondary structure of the DSSP allocation from the PDB file."""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", self.pdb_filename)
        
        try:
            dssp = DSSP(structure[0], self.pdb_filename)
        except Exception as e:
            print(f"DSSP parsing error: {e}")
            return []

        dssp_structure = []
        for key in dssp.keys():
            dssp_code = dssp[key][2]  # Secondary structure code for DSSP
            if dssp_code in ('H', 'G', 'I'):  # Helix
                dssp_structure.append('H')
            elif dssp_code in ('E', 'B'):  # β-sheet
                dssp_structure.append('E')
            else:
                dssp_structure.append('C')  # Coli
        return dssp_structure

    def compare(self):
        """Compare the predicted secondary structure with the DSSP results and calculate the accuracy and sensitivity."""
        if len(self.predicted_structure) != len(self.true_structure):
            print("Predicted and true structure lengths do not match and cannot be compared.")
            return None

        print(f"Projected structure: {self.predicted_structure}")
        print(f"True structure: {self.true_structure}")

        # Calculation accuracy
        tp = sum([1 for i in range(len(self.predicted_structure)) if self.predicted_structure[i] == self.true_structure[i]])
        accuracy = tp / len(self.predicted_structure)

        # Calculating Helix sensitivity
        total_helices_in_dssp = sum([1 for i in range(len(self.true_structure)) if self.true_structure[i] == 'H'])
        helix_sensitivity = (sum([1 for i in range(len(self.predicted_structure)) 
                                  if self.predicted_structure[i] == 'H' and self.true_structure[i] == 'H']) 
                             / total_helices_in_dssp) if total_helices_in_dssp > 0 else 0
        
        # Calculation of β-sheet sensitivity
        total_beta_in_dssp = sum([1 for i in range(len(self.true_structure)) if self.true_structure[i] == 'E'])
        beta_sensitivity = (sum([1 for i in range(len(self.predicted_structure)) 
                                 if self.predicted_structure[i] == 'E' and self.true_structure[i] == 'E']) 
                            / total_beta_in_dssp) if total_beta_in_dssp > 0 else 0

        print(f"Helix sensitivity: {helix_sensitivity}")
        print(f"β sheet sensitivity: {beta_sensitivity}")

        return accuracy, helix_sensitivity, beta_sensitivity
