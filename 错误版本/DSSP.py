import sys, getopt, math
import numpy as np
from itertools import groupby
from operator import itemgetter
from Bio.PDB import PDBParser, DSSP
import subprocess


class HydrogenBondFinder:
    def __init__(self, pdbfile):
        self.pdbfile = pdbfile
        self.interac_atoms = ["N", "C", "O", "H", "H1", "CA"]
        self.q1, self.q2, self.f, self.Emin = 0.42, 0.20, 332, -0.5
        self.dico_res = self.parse_pdb(self.add_hydrogens_with_reduce(pdbfile))

    def add_hydrogens_with_reduce(self, pdb_file):
        """Add hydrogens to PDB file using reduce tool."""
        reduced_pdb_file = "reduced_" + pdb_file
        with open(reduced_pdb_file, 'w') as output_file:
            subprocess.run(['reduce', '-HIS', pdb_file], stdout=output_file)
        print(f"Hydrogens added to: {reduced_pdb_file}")
        return reduced_pdb_file

    def parse_pdb(self, pdbfile):
        """Parse PDB file and store coordinates of each residue by chain."""
        dico_res = {}
        with open(pdbfile, "r") as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM") and line[12:16].strip() in self.interac_atoms:
                    chain_id = line[21].strip()
                    if chain_id not in dico_res:
                        dico_res[chain_id] = {}
                    res_num = int(line[22:26])
                    atom_name = line[12:16].strip()
                    if res_num not in dico_res[chain_id]:
                        dico_res[chain_id][res_num] = {}
                    dico_res[chain_id][res_num][atom_name] = self.extract_atom_coord(line)
        return dico_res

    @staticmethod
    def extract_atom_coord(line):
        """Extract atom coordinates, residue number, and name from a PDB line."""
        atom = {
            'res_num': int(line[22:26]),
            'res_name': line[17:20].strip(),
            'x': float(line[30:38]),
            'y': float(line[38:46]),
            'z': float(line[46:54])
        }
        return atom

    def dist2atoms(self, numresA, numresB, atomA, atomB):
        """Calculate distance between two atoms; return inf if atoms do not exist."""
        if numresA not in self.dico_res or numresB not in self.dico_res:
            return float('inf')
        if atomA not in self.dico_res[numresA] or atomB not in self.dico_res[numresB]:
            return float('inf')
        
        xA, yA, zA = self.dico_res[numresA][atomA]['x'], self.dico_res[numresA][atomA]['y'], self.dico_res[numresA][atomA]['z']
        xB, yB, zB = self.dico_res[numresB][atomB]['x'], self.dico_res[numresB][atomB]['y'], self.dico_res[numresB][atomB]['z']
        return math.sqrt((xB - xA) ** 2 + (yB - yA) ** 2 + (zB - zA) ** 2)

    def is_Hbond(self, numresA, numresB):
        """Determine if there is a hydrogen bond."""
        rON = self.dist2atoms(numresA, numresB, 'O', 'N')
        rOH = self.dist2atoms(numresA, numresB, 'O', 'H') or self.dist2atoms(numresA, numresB, 'O', 'H1')
        rCN = self.dist2atoms(numresA, numresB, 'C', 'N')
        rCH = self.dist2atoms(numresA, numresB, 'C', 'H') or self.dist2atoms(numresA, numresB, 'C', 'H1')
        if rON == float('inf') or rOH == float('inf') or rCN == float('inf') or rCH == float('inf'):
            return False
        E_elec = self.q1 * self.q2 * self.f * (1 / rON + 1 / rCH - 1 / rOH - 1 / rCN)
        return E_elec < self.Emin


class SecondaryStructureAssigner:
    def __init__(self, hydrogen_bond_finder):
        self.dico_res = hydrogen_bond_finder.dico_res
        self.dico_res_struct = self.dict_resid_structure()

    def dict_resid_structure(self):
        """Create a structure dictionary for all residues."""
        res_struct = {}
        for chain_id, residues in self.dico_res.items():
            res_struct[chain_id] = {}
            for res_num, atoms in residues.items():
                struct = {
                    'res_num': atoms['N']['res_num'] if 'N' in atoms else res_num,
                    'res_name': atoms['N']['res_name'] if 'N' in atoms else "UNK",
                    'type': "",
                    'BP1': 0,
                    'BP2': 0,
                    '3T': "",
                    '4T': "",
                    '5T': ""
                }
                res_struct[chain_id][res_num] = struct
        return res_struct

    def assign_secondary_structure(self):
        """Assign secondary structure types to residues."""
        # Placeholder: Add logic for assigning secondary structures
        pass

    def write_predicted_structure(self, filename="predicted_structure.txt"):
        """Write predicted secondary structures to a file by chain."""
        predicted_structure = self.extract_predicted_structure()
        with open(filename, 'w') as file:
            file.write("Predicted Secondary Structure by Chain:\n")
            for chain_id, structure in predicted_structure.items():
                file.write(f"Chain {chain_id}:\n")
                file.write(" ".join(structure) + "\n\n")
        print(f"Predicted secondary structures saved to: {filename}")

    def extract_predicted_structure(self):
        """Extract predicted structure labels."""
        predicted_structure = {}
        for chain_id, residues in self.dico_res_struct.items():
            chain_structure = []
            for res_id, res in residues.items():
                if 'type' in res:
                    if res['type'] == 'H':
                        chain_structure.append('H')
                    elif res['type'] == 'E':
                        chain_structure.append('E')
                    else:
                        chain_structure.append('C')
                else:
                    print(f"Warning: Residue {res_id} in chain {chain_id} is not formatted correctly.")
            predicted_structure[chain_id] = chain_structure
        return predicted_structure


class StructureComparator:
    def __init__(self, pdbfile, predicted_structure):
        self.pdbfile = pdbfile
        self.predicted_structure = predicted_structure
        self.true_structure = self.get_true_structure_dssp()

    def get_true_structure_dssp(self):
        """Retrieve true structure from DSSP for each chain independently."""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", self.pdbfile)
        true_structure_per_chain = {}

        try:
            dssp = DSSP(structure[0], self.pdbfile)  # Using the first model in the structure
        except Exception as e:
            print(f"DSSP parsing error: {e}")
            return {}

        for key in dssp.keys():
            chain_id = key[0]
            if chain_id not in true_structure_per_chain:
                true_structure_per_chain[chain_id] = []

            dssp_code = dssp[key][2]
            if dssp_code in ('H', 'G', 'I'):  # Helix
                true_structure_per_chain[chain_id].append('H')
            elif dssp_code in ('E', 'B'):  # Beta sheet
                true_structure_per_chain[chain_id].append('E')
            else:
                true_structure_per_chain[chain_id].append('C')
        
        return true_structure_per_chain

    def compare_structures(self):
        """Compare predicted and true structures and calculate metrics."""
        overall_predicted = []
        overall_true = []

        for chain_id, chain_pred in self.predicted_structure.items():
            if chain_id not in self.true_structure:
                print(f"Missing chain {chain_id} in true structure data.")
                continue
            true_chain = self.true_structure[chain_id]
            if len(chain_pred) != len(true_chain):
                print(f"Length mismatch in chain {chain_id}")
                continue

            overall_predicted.extend(chain_pred)
            overall_true.extend(true_chain)

            tp = sum(1 for i in range(len(chain_pred)) if chain_pred[i] == true_chain[i])
            accuracy = tp / len(chain_pred)

            tp_helix = sum(1 for i in range(len(chain_pred)) if chain_pred[i] == 'H' and true_chain[i] == 'H')
           
            fn_helix = sum(1 for i in range(len(chain_pred)) if chain_pred[i] != 'H' and true_chain[i] == 'H')
            helix_sensitivity = tp_helix / (tp_helix + fn_helix) if (tp_helix + fn_helix) > 0 else 0

            tp_beta = sum(1 for i in range(len(chain_pred)) if chain_pred[i] == 'E' and true_chain[i] == 'E')
            fn_beta = sum(1 for i in range(len(chain_pred)) if chain_pred[i] != 'E' and true_chain[i] == 'E')
            beta_sensitivity = tp_beta / (tp_beta + fn_beta) if (tp_beta + fn_beta) > 0 else 0

            print(f"Chain {chain_id} - Accuracy: {accuracy}, Helix Sensitivity: {helix_sensitivity}, Beta-Sheet Sensitivity: {beta_sensitivity}")

        # Calculate overall metrics
        total_tp = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] == overall_true[i])
        total_accuracy = total_tp / len(overall_predicted)

        total_tp_helix = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] == 'H' and overall_true[i] == 'H')
        total_fn_helix = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] != 'H' and overall_true[i] == 'H')
        total_helix_sensitivity = total_tp_helix / (total_tp_helix + total_fn_helix) if (total_tp_helix + total_fn_helix) > 0 else 0

        total_tp_beta = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] == 'E' and overall_true[i] == 'E')
        total_fn_beta = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] != 'E' and overall_true[i] == 'E')
        total_beta_sensitivity = total_tp_beta / (total_tp_beta + total_fn_beta) if (total_tp_beta + total_fn_beta) > 0 else 0

        print(f"\nOverall Protein - Accuracy: {total_accuracy}, Helix Sensitivity: {total_helix_sensitivity}, Beta-Sheet Sensitivity: {total_beta_sensitivity}")


# ********************* MAIN ********************* #

def main(argv):
    pdbfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('dssp_like.py -i <pdbfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('dssp_like.py -i <pdbfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            pdbfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    print('Input file is', pdbfile)
    print('Output file is', outputfile)

    # Initialize hydrogen bond finder
    hb_finder = HydrogenBondFinder(pdbfile)
    
    # Assign secondary structures
    structure_assigner = SecondaryStructureAssigner(hb_finder)
    structure_assigner.assign_secondary_structure()
    structure_assigner.write_predicted_structure(outputfile)

    # Compare structures
    comparator = StructureComparator(pdbfile, structure_assigner.extract_predicted_structure())
    comparator.compare_structures()


if __name__ == "__main__":
   main(sys.argv[1:])
