import numpy as np

class HydrogenBondDetector:
    """Identifier les liaisons hydrogène en fonction de l'énergie"""

    def __init__(self, main_chain_atoms):
        self.main_chain_atoms = main_chain_atoms
    
    def calculate_hydrogen_bonds(self, energy_threshold= -0.5):
        hbonds = []
        num_residues = len(self.main_chain_atoms)
        
        for i in range(num_residues):  # Iterate over donor residues
            donor_residue = self.main_chain_atoms[i]
            donor_chain_id = donor_residue.get('chain_id', None)  # 获取 donor 的链号
            for j in range(i + 1, num_residues):  # Iterate over acceptor residues
                acceptor_residue = self.main_chain_atoms[j]
                acceptor_chain_id = acceptor_residue.get('chain_id', None)  # 获取 acceptor 的链号
                
                if donor_chain_id != acceptor_chain_id:
                    #print(f"跨链氢键跳过: donor_chain_id={donor_chain_id}, acceptor_chain_id={acceptor_chain_id}")
                    continue

                # Calculate hydrogen bond energy
                energy, r_ON = self._calculate_hydrogen_bond_energy(donor_residue, acceptor_residue)
                
                # Only consider valid hydrogen bonds (energy < threshold)
                if energy < energy_threshold:
                    hbonds.append({
                        'donor': donor_residue,
                        'acceptor': acceptor_residue,
                        'energy': energy,
                        'r_ON': r_ON,
                    })

        return hbonds

    def _calculate_hydrogen_bond_energy(self, donor_residue, acceptor_residue):
        """Calculer l'énergie d'une liaison hydrogène à partir de l'équation énergétique"""
        
        # Ensure that necessary atoms are present in both residues
        if 'O' not in acceptor_residue or 'N' not in donor_residue or 'C' not in donor_residue:
            return float('inf'), None
        
        # Get atomic coordinates
        O = acceptor_residue['O']  # Oxygen atom in acceptor residue
        N = donor_residue['N']     # Nitrogen atom in donor residue
        C = donor_residue['C']     # Carbon atom in donor residue
        
        # Calculate distances between atoms
        r_ON = np.linalg.norm(O - N)  
        r_CN = np.linalg.norm(C - N)  
        
        # Simplified energy calculation (considering O-N and C-N distances)
        if r_ON > 4.5 or r_ON < 2.0:
            return float('inf'), r_ON

        energy = 0.084 * (1 / r_ON - 1 / r_CN) * 332
        return energy, r_ON



