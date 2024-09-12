from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

class ReadPdbStructure:
    """
    Represent and extract the main chain atoms (N, CA, C, O) in the PDB structure.

    Attributes:
    pdb_filename (str): The filename of the PDB structure to be read.
    structure: The parsed structure from the PDB file.
    main_chain_atoms (list): A list containing atomic coordinates of the main chain atoms for each residue.
    """

    def __init__(self, pdb_filename):
        """
        Initialize the ReadPdbStructure class.

        Parameters:
        pdb_filename (str): The path to the PDB file to read.
        """
        self.pdb_filename = pdb_filename
        self.structure = self._read_pdb_structure()  # Parse the PDB structure
        self.main_chain_atoms = self._get_main_chain_atoms()  # Extract main chain atoms

    def _read_pdb_structure(self):
        """
        Read the PDB file and return the structure.

        Returns:
        structure: Parsed structure from the PDB file using PDBParser.
        """
        parser = PDBParser(QUIET=True)  # Suppress warnings during parsing
        structure = parser.get_structure("protein", self.pdb_filename)  # Parse the PDB file
        return structure
    
    def _get_main_chain_atoms(self):
        """
        Extract the main chain atomic coordinates of each residue (N, CA, C, O).

        This method iterates through the structure and collects the main chain atoms
        of residues that are amino acids (is_aa check). It also ensures that residues
        are continuous in sequence by checking the residue IDs and logs discontinuities.

        Returns:
        list: A list of dictionaries, each containing the residue name, ID, chain ID,
              and coordinates of the N, CA, C, and O atoms, along with chain direction.
        """
        main_chain_atoms = []  # List to store atomic coordinates of main chain
        previous_res_id = None  # Track previous residue ID to check sequence continuity
        direction = "forward"  # Direction of chain traversal

        # Iterate through each model in the structure
        for model in self.structure:
            for chain in model:
                # Iterate through residues in the chain
                for residue in chain:
                    # Check if the residue is an amino acid and has the main chain atoms
                    if is_aa(residue) and residue.has_id('N') and residue.has_id('CA') and residue.has_id('C') and residue.has_id('O'):
                        res_id = residue.get_id()[1]  # Get residue ID
                        chain_id = chain.get_id()  # Get chain ID

                        # Check for discontinuities in residue sequence
                        if previous_res_id is not None and abs(res_id - previous_res_id) > 1:
                            print(f"Non-continuous residue IDs: previous_res_id={previous_res_id}, current_res_id={res_id}, chain={chain_id}")    
                        previous_res_id = res_id  # Update previous residue ID

                        # Append main chain atom coordinates and information
                        main_chain_atoms.append({
                            'residue': residue.get_resname(),  # Residue name
                            'res_id': res_id,  # Residue ID
                            'chain': chain_id,  # Chain ID
                            'N': residue['N'].get_coord(),  # N atom coordinates
                            'CA': residue['CA'].get_coord(),  # CA atom coordinates
                            'C': residue['C'].get_coord(),  # C atom coordinates
                            'O': residue['O'].get_coord(),  # O atom coordinates
                            'direction': direction  # Direction of the chain
                        })
        return main_chain_atoms
    
    def get_main_chain_atoms(self):
        """
        Get the extracted main chain atom information.

        Returns:
        list: A list of dictionaries containing the main chain atom information.
        """
        return self.main_chain_atoms
