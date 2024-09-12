from Bio.PDB import PDBParser, DSSP

class Compare:
    """
    Class for comparing custom secondary structure predictions with DSSP-derived secondary structures.

    Attributes:
    pdb_filename (str): The path to the PDB file to be analyzed.
    predicted_structure (list): The predicted secondary structure ('H', 'E', 'C') for each residue.
    true_structure (list): The DSSP-derived secondary structure for each residue.
    """

    def __init__(self, pdb_filename, predicted_structure):
        """
        Initialize the Compare class with the PDB file and the predicted structure.

        Parameters:
        pdb_filename (str): Path to the PDB file.
        predicted_structure (list): List of predicted secondary structures ('H', 'E', 'C').
        """
        self.pdb_filename = pdb_filename
        self.predicted_structure = predicted_structure
        self.true_structure = self._dssp_structure()  # Get true structure from DSSP

    def _dssp_structure(self):
        """
        Extract the DSSP secondary structure assignment from the PDB file.

        This method runs DSSP on the structure and translates DSSP's secondary structure
        codes into simplified forms: 'H' for Helix, 'E' for β-Sheet, and 'C' for Coil.

        Returns:
        list: A list of secondary structure labels ('H', 'E', 'C').
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", self.pdb_filename)
        
        try:
            dssp = DSSP(structure[0], self.pdb_filename)  # DSSP analysis on the first model
        except Exception as e:
            print(f"DSSP parsing error: {e}")
            return []

        dssp_structure = []
        for key in dssp.keys():
            dssp_code = dssp[key][2]  # Extract the DSSP secondary structure code
            if dssp_code in ('H', 'G', 'I'):  # Helix types in DSSP
                dssp_structure.append('H')
            elif dssp_code in ('E', 'B'):  # Beta-sheet types in DSSP
                dssp_structure.append('E')
            else:
                dssp_structure.append('C')  # Coil (everything else)
        return dssp_structure

    def compare(self):
        """
        Compare the predicted structure to the DSSP-derived structure and compute accuracy and sensitivity.

        This method checks whether the predicted and DSSP-derived secondary structures match, then
        calculates overall accuracy, as well as helix and β-sheet sensitivity.

        Returns:
        tuple: accuracy (float), helix sensitivity (float), and β-sheet sensitivity (float).
        """
        if len(self.predicted_structure) != len(self.true_structure):
            print("Predicted and true structure lengths do not match and cannot be compared.")
            return None

        print(f"Predicted structure: {self.predicted_structure}")
        print(f"True structure: {self.true_structure}")

        # Calculate overall accuracy (True Positive matches)
        tp = sum([1 for i in range(len(self.predicted_structure)) if self.predicted_structure[i] == self.true_structure[i]])
        accuracy = tp / len(self.predicted_structure)

        # Calculate Helix Sensitivity
        # Step 1: Count the total number of helices ('H') in the true structure (Total Actual Positives for helices)
        total_helices_in_dssp = sum([1 for i in range(len(self.true_structure)) if self.true_structure[i] == 'H'])

        # Step 2: Calculate True Positives (TP) for helices, i.e., when both the predicted and true structure are 'H'
        tp_helix = sum([1 for i in range(len(self.predicted_structure)) 
                        if self.predicted_structure[i] == 'H' and self.true_structure[i] == 'H'])

        # Step 3: Calculate False Negatives (FN) for helices, i.e., when the true structure is 'H' but predicted structure is not 'H'
        fn_helix = sum([1 for i in range(len(self.predicted_structure)) 
                        if self.predicted_structure[i] != 'H' and self.true_structure[i] == 'H'])

        # Step 4: Calculate helix sensitivity (recall) = TP / (TP + FN), but handle the case where no helices are present in the true structure
        helix_sensitivity = (tp_helix / (tp_helix + fn_helix)) if (tp_helix + fn_helix) > 0 else 0


        # Calculate β-Sheet Sensitivity
        # Step 1: Count the total number of β-sheets ('E') in the true structure (Total Actual Positives for β-sheets)
        total_beta_in_dssp = sum([1 for i in range(len(self.true_structure)) if self.true_structure[i] == 'E'])

        # Step 2: Calculate True Positives (TP) for β-sheets, i.e., when both the predicted and true structure are 'E'
        tp_beta = sum([1 for i in range(len(self.predicted_structure)) 
                    if self.predicted_structure[i] == 'E' and self.true_structure[i] == 'E'])

        # Step 3: Calculate False Negatives (FN) for β-sheets, i.e., when the true structure is 'E' but predicted structure is not 'E'
        fn_beta = sum([1 for i in range(len(self.predicted_structure)) 
                    if self.predicted_structure[i] != 'E' and self.true_structure[i] == 'E'])

        # Step 4: Calculate β-sheet sensitivity (recall) = TP / (TP + FN), but handle the case where no β-sheets are present in the true structure
        beta_sensitivity = (tp_beta / (tp_beta + fn_beta)) if (tp_beta + fn_beta) > 0 else 0


        return accuracy, helix_sensitivity, beta_sensitivity

