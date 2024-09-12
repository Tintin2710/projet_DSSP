import numpy as np

class HydrogenBondDetector:
    """
    Identify hydrogen bonds based on energy between donor and acceptor residues.
    
    Attributes:
    main_chain_atoms (list): A list of main chain atom information, including atomic coordinates.
    """

    def __init__(self, main_chain_atoms):
        """
        Initialize with main chain atom data.

        Parameters:
        main_chain_atoms (list): List containing atomic coordinates of main chain atoms.
        """
        self.main_chain_atoms = main_chain_atoms
    
    def calculate_hydrogen_bonds(self, energy_threshold=-0.5):
        """
        Calculate hydrogen bonds between residues based on energy threshold.

        Parameters:
        energy_threshold (float): Energy cutoff for valid hydrogen bonds, default is -0.5 kcal/mol.

        Returns:
        list: List of hydrogen bonds where energy is below the threshold.
        """
        hbonds = []
        num_residues = len(self.main_chain_atoms)

        for i in range(num_residues):  # Donor residue loop
            donor_residue = self.main_chain_atoms[i]
            donor_chain_id = donor_residue.get('chain_id', None)  # Get donor's chain ID
            
            for j in range(i + 1, num_residues):  # Acceptor residue loop
                acceptor_residue = self.main_chain_atoms[j]
                acceptor_chain_id = acceptor_residue.get('chain_id', None)  # Get acceptor's chain ID
                
                if donor_chain_id != acceptor_chain_id:
                    continue  # Skip bonds between different chains

                # Calculate hydrogen bond energy
                energy, r_ON = self._calculate_hydrogen_bond_energy(donor_residue, acceptor_residue)
                
                if energy < energy_threshold:  # Check if bond meets energy criteria
                    hbonds.append({
                        'donor': donor_residue,
                        'acceptor': acceptor_residue,
                        'energy': energy,
                        'r_ON': r_ON,
                    })

        return hbonds

    def _calculate_hydrogen_bond_energy(self, donor_residue, acceptor_residue):
        """
        Calculate hydrogen bond energy using a simplified equation.

        Parameters:
        donor_residue (dict): Donor residue atomic data.
        acceptor_residue (dict): Acceptor residue atomic data.

        Returns:
        tuple: Hydrogen bond energy and O-N distance (r_ON).
        """
        if 'O' not in acceptor_residue or 'N' not in donor_residue or 'C' not in donor_residue:
            return float('inf'), None  # Return high energy if required atoms are missing
        
        # Get atom coordinates
        O = acceptor_residue['O']
        N = donor_residue['N']
        C = donor_residue['C']
        
        # Calculate distances
        r_ON = np.linalg.norm(O - N)
        r_CN = np.linalg.norm(C - N)
        
        # Check if O-N distance is within reasonable range
        if r_ON > 4.5 or r_ON < 2.0:
            return float('inf'), r_ON  # Return high energy if distance is unrealistic

        # Simplified energy calculation
        energy = 0.084 * (1 / r_ON - 1 / r_CN) * 332
        return energy, r_ON




