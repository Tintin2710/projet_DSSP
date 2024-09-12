import math
import numpy as np

class SecondaryStructureAssigner:
    """
    Secondary Structure Assigner for detecting α-Helices, β-Folds, and β-Bridges 
    based on hydrogen bonds and geometric features following the Kabsch & Sander method.
    """

    def __init__(self, main_chain_atoms, hbonds):
        """
        Initialize with main chain atom data and hydrogen bonds.

        Parameters:
        main_chain_atoms (list): List of dictionaries containing atomic coordinates of main chain atoms.
        hbonds (list): List of hydrogen bonds with their respective energies and atom information.
        """
        self.hbonds = hbonds
        self.main_chain_atoms = main_chain_atoms

    def calculate_dihedral_angle(self, atom1, atom2, atom3, atom4):
        """
        Calculate the dihedral angle between four atoms using their coordinates.

        Parameters:
        atom1, atom2, atom3, atom4 (numpy array): Coordinates of the four atoms.

        Returns:
        float: The dihedral angle in degrees.
        """
        # Vector calculations for the dihedral angle
        v1 = atom2 - atom1
        v2 = atom3 - atom2
        v3 = atom4 - atom3

        # Normal vectors to the planes
        n1 = np.cross(v1, v2)
        n2 = np.cross(v2, v3)

        # Calculate the angle between the normals
        x = np.dot(n1, n2)
        y = np.dot(np.cross(n1, n2), v2 / np.linalg.norm(v2))

        return math.atan2(y, x) * 180 / math.pi

    def assign_secondary_structure(self):
        """
        Assign secondary structure (α-Helix, β-Fold, or β-Bridge) based on hydrogen bonds and dihedral angles.

        This method evaluates each hydrogen bond and checks the dihedral angles (φ, ψ) to classify 
        residues as being part of an α-helix or β-sheet. β-sheets are further classified as parallel 
        or antiparallel based on bond direction.

        Returns:
        list: A list containing the secondary structure assignment for each residue, 
              where 'H' indicates α-Helix, 'E' indicates β-Sheet, and 'C' indicates Coil.
        """
        helix = set()
        beta_parallel = set()
        beta_antiparallel = set()

        # Map residue IDs to their indices in the main chain list
        res_id_to_index = {res['res_id']: idx for idx, res in enumerate(self.main_chain_atoms)}

        for bond in self.hbonds:
            donor_residue_id = bond['donor']['res_id']
            acceptor_residue_id = bond['acceptor']['res_id']
            diff = abs(donor_residue_id - acceptor_residue_id)

            # Calculate dihedral angles for donor and acceptor residues
            if donor_residue_id in res_id_to_index and acceptor_residue_id in res_id_to_index:
                donor_idx = res_id_to_index[donor_residue_id]
                acceptor_idx = res_id_to_index[acceptor_residue_id]

                # Ensure valid indices for dihedral angle calculations
                if 1 <= donor_idx < len(self.main_chain_atoms) - 2 and 1 <= acceptor_idx < len(self.main_chain_atoms) - 2:
                    donor_phi = self.calculate_dihedral_angle(
                        self.main_chain_atoms[donor_idx - 1]['C'],  
                        self.main_chain_atoms[donor_idx]['N'],      
                        self.main_chain_atoms[donor_idx]['CA'],     
                        self.main_chain_atoms[donor_idx + 1]['C']   
                    )
                    donor_psi = self.calculate_dihedral_angle(
                        self.main_chain_atoms[donor_idx]['N'],      
                        self.main_chain_atoms[donor_idx]['CA'],     
                        self.main_chain_atoms[donor_idx + 1]['C'],  
                        self.main_chain_atoms[donor_idx + 1]['N']   
                    )
                    acceptor_phi = self.calculate_dihedral_angle(
                        self.main_chain_atoms[acceptor_idx - 1]['C'],  
                        self.main_chain_atoms[acceptor_idx]['N'],      
                        self.main_chain_atoms[acceptor_idx]['CA'],     
                        self.main_chain_atoms[acceptor_idx + 1]['C']   
                    )
                    acceptor_psi = self.calculate_dihedral_angle(
                        self.main_chain_atoms[acceptor_idx]['N'],      
                        self.main_chain_atoms[acceptor_idx]['CA'],     
                        self.main_chain_atoms[acceptor_idx + 1]['C'],  
                        self.main_chain_atoms[acceptor_idx + 1]['N']   
                    )

                    # Helix detection: Check φ and ψ angles and distance
                    if diff in [3,4,5] and (-90 <= donor_phi <= -40 and -70 <= donor_psi <= -20) and (-90 <= acceptor_phi <= -40 and -70 <= acceptor_psi <= -20):
                        helix.add(donor_residue_id)
                        helix.add(acceptor_residue_id)

                    # β-sheet detection based on φ and ψ angles
                    elif (-160 <= donor_phi <= -110 and 110 <= donor_psi <= 160) or (-140 <= donor_phi <= -90 and 90 <= donor_psi <= 140):
                        donor_index = res_id_to_index[donor_residue_id]
                        acceptor_index = res_id_to_index[acceptor_residue_id]
                        if donor_index < acceptor_index:
                            direction = self.main_chain_atoms[donor_index]['direction']
                            if direction == self.main_chain_atoms[acceptor_index]['direction']:
                                beta_parallel.add(donor_residue_id)
                                beta_parallel.add(acceptor_residue_id)
                            else:
                                beta_antiparallel.add(donor_residue_id)
                                beta_antiparallel.add(acceptor_residue_id)

        # Initialize all residues to 'C' (coil)
        secondary_structure = ['C'] * len(self.main_chain_atoms)

        # Assign α-Helices ('H')
        for res_id in helix:
            if res_id in res_id_to_index:
                index = res_id_to_index[res_id]
                secondary_structure[index] = 'H'  

        # Assign β-sheets ('E')
        for res_id in beta_parallel.union(beta_antiparallel):
            if res_id in res_id_to_index:
                index = res_id_to_index[res_id]
                if secondary_structure[index] != 'H':
                    secondary_structure[index] = 'E'

        # Debugging information: print the number of atoms in the main chain and the allocation result
        print("Length of main_chain_atoms:", len(self.main_chain_atoms))
        print("Helix indices:", helix)
        print("Parallel Beta indices:", beta_parallel)
        print("Antiparallel Beta indices:", beta_antiparallel)

        return secondary_structure




