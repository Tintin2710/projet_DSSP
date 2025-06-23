import math
import subprocess

# Constants for electrostatic interactions
INTERACTION_ATOMS = ["N", "C", "O", "H", "H1", "CA"]
Q1, Q2, F, EMIN = 0.42, 0.20, 332, -0.5

class HydrogenBondFinder:
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.reduced_pdb_file = self.add_hydrogens()
        self.dico_res = self.parse_pdb(self.reduced_pdb_file)

    def add_hydrogens(self):
        """Adds hydrogens to the PDB file using the reduce tool and returns the modified file name."""
        reduced_pdb_file = "reduced_" + self.pdb_file
        try:
            with open(reduced_pdb_file, 'w') as output_file:
                subprocess.run(['reduce', '-HIS', self.pdb_file], stdout=output_file, check=True)
            print(f"Hydrogens added to: {reduced_pdb_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error adding hydrogens: {e}")
            return None
        return reduced_pdb_file

    def parse_pdb(self, pdb_file):
        """Parse the PDB file, storing the coordinates of each residue by chain."""
        dico_res = {}
        try:
            with open(pdb_file, "r") as file:
                for line in file:
                    if line.startswith("ATOM") and line[12:16].strip() in INTERACTION_ATOMS:
                        chain_id = line[21].strip()
                        res_num = int(line[22:26].strip())
                        atom_name = line[12:16].strip()
                        atom = {
                            'res_num': res_num,
                            'res_name': line[17:20].strip(),
                            'x': float(line[30:38]),
                            'y': float(line[38:46]),
                            'z': float(line[46:54])
                        }
                        if chain_id not in dico_res:
                            dico_res[chain_id] = {}
                        if res_num not in dico_res[chain_id]:
                            dico_res[chain_id][res_num] = {}
                        dico_res[chain_id][res_num][atom_name] = atom
        except Exception as e:
            print(f"Error parsing PDB file: {e}")
        return dico_res

    def calculate_distance(self, resA, resB, atomA, atomB):
        """Calculate the distance between two atoms, returning infinity if any atom is missing."""
        if resA not in self.dico_res or resB not in self.dico_res:
            return float('inf')
        if atomA not in self.dico_res[resA] or atomB not in self.dico_res[resB]:
            return float('inf')

        a, b = self.dico_res[resA][atomA], self.dico_res[resB][atomB]
        distance = math.sqrt((b['x'] - a['x']) ** 2 + (b['y'] - a['y']) ** 2 + (b['z'] - a['z']) ** 2)
        return distance

    def is_hydrogen_bond(self, resA, resB):
        """Determines if there is a hydrogen bond between residues A and B."""
        rON = self.calculate_distance(resA, resB, 'O', 'N')
        rOH = self.calculate_distance(resA, resB, 'O', 'H') or self.calculate_distance(resA, resB, 'O', 'H1')
        rCN = self.calculate_distance(resA, resB, 'C', 'N')
        rCH = self.calculate_distance(resA, resB, 'C', 'H') or self.calculate_distance(resA, resB, 'C', 'H1')

        if float('inf') in {rON, rOH, rCN, rCH}:
            return False  # Return false if any required atom is missing for bond calculation

        # Calculate electrostatic energy
        E_elec = Q1 * Q2 * F * (1 / rON + 1 / rCH - 1 / rOH - 1 / rCN)
        return E_elec < EMIN

    def find_all_hydrogen_bonds(self):
        """Iterates through residues to identify all potential hydrogen bonds and returns them."""
        bonds = {}
        for chain_id, residues in self.dico_res.items():
            for resA in residues:
                for resB in residues:
                    if resA < resB:  # Avoid redundant checks
                        if self.is_hydrogen_bond(resA, resB):
                            bonds.setdefault(chain_id, []).append((resA, resB))
        return bonds
