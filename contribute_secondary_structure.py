import sys
import getopt
import math
import numpy as np
import subprocess
from itertools import groupby
from operator import itemgetter
from Bio.PDB import PDBParser, DSSP



class PDBProcessor:
    """Handles PDB file parsing and hydrogen addition."""

    interac_atoms = ["N", "C", "O", "H", "H1", "CA"]

    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.reduced_pdb_file = self.add_hydrogens_with_reduce()

    def add_hydrogens_with_reduce(self):
        """Adds hydrogens to the PDB file using the Reduce tool and returns the modified file name."""
        reduced_pdb_file = "reduced_" + self.pdb_file
        with open(reduced_pdb_file, 'w') as output_file:
            subprocess.run(['reduce', '-HIS', self.pdb_file], stdout=output_file)
        print(f"Hydrogens added to: {reduced_pdb_file}")
        return reduced_pdb_file

    def parse_pdb(self):
        """Parses the PDB file, storing atomic coordinates per residue and chain."""
        dico_res = {}
        with open(self.reduced_pdb_file, "r") as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM") and line[12:16].strip() in self.interac_atoms:
                    chain_id = line[21].strip()
                    res_num = int(line[22:26])
                    atom_name = line[12:16].strip()
                    if chain_id not in dico_res:
                        dico_res[chain_id] = {}
                    if res_num not in dico_res[chain_id]:
                        dico_res[chain_id][res_num] = {}
                    atom = {
                        'res_num': res_num,
                        'res_name': line[17:20].strip(),
                        'x': float(line[30:38]),
                        'y': float(line[38:46]),
                        'z': float(line[46:54])
                    }
                    dico_res[chain_id][res_num][atom_name] = atom
        return dico_res


import math

class HydrogenBondCalculator:
    """Calculates hydrogen bonds between residues."""

    q1, q2, f, Emin = 0.42, 0.20, 332, -0.5

    def __init__(self):
        self.hydrogen_bond_count = 0 

    @staticmethod
    def distance(atom_a, atom_b):
        """Calculates the distance between two atoms, returns a large value if not found."""
        if not atom_a or not atom_b:
            return float('inf')
        return math.sqrt((atom_b['x'] - atom_a['x']) ** 2 +
                         (atom_b['y'] - atom_a['y']) ** 2 +
                         (atom_b['z'] - atom_a['z']) ** 2)

    def is_hydrogen_bond(self, atom_a, atom_b):
        """
        Determines if a hydrogen bond exists based on distance and energy calculation.
        This version uses robust logic to select the correct hydrogen atom.
        """

        # First, try to get the 'H' atom.
        h_atom_coords = atom_b.get('H')
        #    If the 'H' atom does not exist (get returns None), try to get the 'H1' atom as a backup.
        if h_atom_coords is None:
            h_atom_coords = atom_b.get('H1')

        # If neither 'H' nor 'H1' exists, h_atom_coords will be None.
        # The distance method will correctly return float('inf').
        r_oh = self.distance(atom_a.get('O'), h_atom_coords)
        r_ch = self.distance(atom_a.get('C'), h_atom_coords)

        # Other distance calculations remain unchanged
        r_on = self.distance(atom_a.get('O'), atom_b.get('N'))
        r_cn = self.distance(atom_a.get('C'), atom_b.get('N'))

        # Check if any distance calculation failed (i.e., atom does not exist)
        if float('inf') in {r_on, r_oh, r_cn, r_ch}:
            return False

        # Energy calculation formula remains unchanged
        E_elec = self.q1 * self.q2 * self.f * (1 / r_on + 1 / r_ch - 1 / r_oh - 1 / r_cn)
        if E_elec < self.Emin:
            self.hydrogen_bond_count += 1
            return True
            
        return False


class SecondaryStructureAssigner:
    def __init__(self, dico_res, bond_calculator):
        self.dico_res = dico_res
        self.bond_calculator = bond_calculator
        self.dico_res_struct = self.dict_resid_structure()

    def dict_resid_structure(self):
        """Creates a dictionary with structure labels for each residue, ensuring structural information for each."""
        res_struct = {}
        for chain_id, residues in self.dico_res.items():
            res_struct[chain_id] = {}
            for res_num, atoms in residues.items():
                struct = {
                    'res_num': res_num,
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

    def list_nturn(self, nturn):
        """Identifies n-turns in each chain using hydrogen bond patterns."""
        chain_nturns = {}
        for chain_id, residues in self.dico_res.items():
            list_turn, list_partner = [], []
            res_keys = sorted(residues.keys())
            for res in res_keys:
                if res + nturn in residues and 'H' in residues[res + nturn]:
                    if self.bond_calculator.is_hydrogen_bond(residues[res], residues[res + nturn]):
                        list_turn.extend([res, res + nturn])
                        list_partner.append((res, res + nturn))
            chain_nturns[chain_id] = (list(set(list_turn)), list_partner)
        return chain_nturns

    @staticmethod
    def successives_nb(liste):
        """Finds successive numbers in a list and groups them."""
        succ_nb = []
        for _, group in groupby(enumerate(sorted(liste)), lambda x: x[0] - x[1]):
            group = list(map(itemgetter(1), group))
            if len(group) > 1:
                succ_nb.extend(group)
        return succ_nb

    def helix_alpha(self, list_4turn):
        """Identifies residues in alpha helices for each chain."""
        return {chain_id: self.successives_nb(turns[0]) for chain_id, turns in list_4turn.items()}

    def helix_310(self, list_3turn):
        """Identifies residues in 3,10 helices for each chain."""
        return {chain_id: self.successives_nb(turns[0]) for chain_id, turns in list_3turn.items()}

    def helix_pi(self, list_5turn):
        """Identifies residues in pi helices for each chain."""
        return {chain_id: self.successives_nb(turns[0]) for chain_id, turns in list_5turn.items()}

    def partner_bridge(self, list_ij):
        """Finds bridge partners for given residue pairs."""
        partner = []
        list_i, list_j = list_ij[0], list_ij[1]
        for i in range(len(list_i)):
            partner.append((list_i[i], list_j[i]))
        for j in range(len(list_j)):
            partner.append((list_j[j], list_i[j]))
        return partner

    def parallel_bridge(self):
        """Finds residues in parallel beta-sheet bridges for each chain."""
        bridge_per_chain = {}
        for chain_id, residues in self.dico_res.items():
            list_i, list_j = [], []
            res_keys = sorted(residues.keys())
            for i in res_keys:
                for j in res_keys:
                    if i < j and j >= i + 3:
                        # Check if i-1, i+1, and j are present in residues before accessing
                        if (i - 1 in residues and i + 1 in residues and j in residues) and \
                            self.bond_calculator.is_hydrogen_bond(residues[i - 1], residues[j]) and \
                            self.bond_calculator.is_hydrogen_bond(residues[j], residues[i + 1]):
                            list_i.append(i)
                            list_j.append(j)
            bridge_per_chain[chain_id] = [list_i, list_j]
        return bridge_per_chain


    def antiparallel_bridge(self):
        """Finds residues in antiparallel beta-sheet bridges for each chain."""
        bridge_per_chain = {}
        for chain_id, residues in self.dico_res.items():
            list_i, list_j = [], []
            res_keys = sorted(residues.keys())
            for i in res_keys:
                for j in res_keys:
                    if i < j and j >= i + 3:
                        if self.bond_calculator.is_hydrogen_bond(residues[i], residues[j]) and \
                           self.bond_calculator.is_hydrogen_bond(residues[j], residues[i]):
                            list_i.append(i)
                            list_j.append(j)
            bridge_per_chain[chain_id] = [list_i, list_j]
        return bridge_per_chain

    def vect2atoms(self, numresA, numresB, atomA, atomB):
        """Calculates the vector between atoms of two residues."""
        try:
            xA, yA, zA = self.dico_res[numresA][atomA]['x'], self.dico_res[numresA][atomA]['y'], self.dico_res[numresA][atomA]['z']
            xB, yB, zB = self.dico_res[numresB][atomB]['x'], self.dico_res[numresB][atomB]['y'], self.dico_res[numresB][atomB]['z']
            return (xB - xA, yB - yA, zB - zA)
        except KeyError as e:
            return (0.0, 0.0, 0.0)
    
    @staticmethod
    def unit_vector(vector):
        """Returns the normalized form of the input vector."""
        return vector / np.linalg.norm(vector)
    
    def angle_between(self, v1, v2):
        """Returns the angle in degrees between vectors v1 and v2."""
        unit_v1 = self.unit_vector(v1)  # updated call, now as a static method
        unit_v2 = self.unit_vector(v2)
        return math.degrees(np.arccos(np.clip(np.dot(unit_v1, unit_v2), -1.0, 1.0)))




    def bend(self):
        """Finds residues involved in bends for each chain based on vector angles."""
        list_bend_per_chain = {}
        for chain_id, residues in self.dico_res.items():
            list_bend = []
            res_keys = sorted(residues.keys())
            for i in range(2, len(res_keys) - 2):
                v1 = self.vect2atoms(res_keys[i - 2], res_keys[i], 'CA', 'CA')
                v2 = self.vect2atoms(res_keys[i], res_keys[i + 2], 'CA', 'CA')
                if self.angle_between(v1, v2) > 70:
                    list_bend.append(res_keys[i])
            list_bend_per_chain[chain_id] = list_bend
        return list_bend_per_chain

    def res_bend ( self, dico_res ,  dico_res_struct , list_bend) :
        """Refills dict with structural informations for all residues

        Parameters
        ----------
        dico_res : dict
            dict with all residues
        dico_res_struct : dict
            dict with all residues structures informations
        list_bend : list
            residues involved in bend

        Returns
        -------
        dict
            dico_res_struct

        """
        for res in dico_res:
            if res in list_bend  :
                dico_res_struct[res]['type'] = 'S'
        return dico_res_struct
    
    @staticmethod
    def res_turnT(dico_res, dico_res_struct, list_4turn, list_3turn, list_5turn):
        """逐链处理，将结构信息填入所有残基字典中"""
        for chain_id in dico_res_struct.keys():
            if chain_id not in list_4turn or chain_id not in list_3turn or chain_id not in list_5turn:
                continue  
            
            chain_struct = dico_res_struct[chain_id]
            
            for res in chain_struct.keys():
                # 处理4-turn, 3-turn, 和5-turn结构
                if res in set(list_4turn[chain_id][0]) or res in set(list_5turn[chain_id][0]) or res in set(list_3turn[chain_id][0]):
                    chain_struct[res]['type'] = 'T'
                
                # 处理每个结构类型的起始和终止残基
                if res in [x[0] for x in list_3turn[chain_id][1]]:
                    chain_struct[res]["3T"] = '>'
                    if res + 1 in chain_struct:
                        chain_struct[res + 1]["3T"] = '3'
                    if res + 2 in chain_struct:
                        chain_struct[res + 2]["3T"] = '3'
                if res in [x[1] for x in list_3turn[chain_id][1]]:
                    chain_struct[res]['3T'] = '<'
                    
                if res in [x[0] for x in list_4turn[chain_id][1]]:
                    chain_struct[res]["4T"] = '>'
                    if res + 1 in chain_struct:
                        chain_struct[res + 1]["4T"] = '4'
                    if res + 2 in chain_struct:
                        chain_struct[res + 2]["4T"] = '4'
                    if res + 3 in chain_struct:
                        chain_struct[res + 3]["4T"] = '4'
                if res in [x[1] for x in list_4turn[chain_id][1]]:
                    chain_struct[res]['4T'] = '<'
                    
                if res in [x[0] for x in list_5turn[chain_id][1]]:
                    chain_struct[res]["5T"] = '>'
                    if res + 1 in chain_struct:
                        chain_struct[res + 1]["5T"] = '5'
                    if res + 2 in chain_struct:
                        chain_struct[res + 2]["5T"] = '5'
                    if res + 3 in chain_struct:
                        chain_struct[res + 3]["5T"] = '5'
                    if res + 4 in chain_struct:
                        chain_struct[res + 4]["5T"] = '5'
                if res in [x[1] for x in list_5turn[chain_id][1]]:
                    chain_struct[res]['5T'] = '<'
                    
        return dico_res_struct


    def res_betasheetE(self, dico_res, dico_res_struct, par, antipar):
        """Add beta-sheet information to all residues in each chain.
        
        Parameters
        ----------
        dico_res : dict
            Includes all residues organized by chain.
        dico_res_struct : dict
            Includes structural information for all residues.
        par : dict
            Contains parallel bridge information for each chain.
        antipar : dict
            Contains antiparallel bridge information for each chain.

        Returns
        -------
        dict
            Updated `dico_res_struct` dictionary.
        """
        for chain_id in dico_res_struct.keys():
            if chain_id not in par or chain_id not in antipar:
                continue  

            chain_struct = dico_res_struct[chain_id]
            
            # Treat parallel bridges
            for res in chain_struct.keys():
                if res in par[chain_id][0] or res in par[chain_id][1]:
                    chain_struct[res]['type'] = 'E'
                    partner = self.partner_bridge(par[chain_id])
                    list_partner = [item for item in partner if item[0] == res]
                    chain_struct[res]['BP1'] = list_partner[0][1]
                    if len(list_partner) == 2:
                        chain_struct[res]['BP2'] = list_partner[1][1]
                        
            # Treat antiparallel bridges
            for res in chain_struct.keys():
                if res in antipar[chain_id][0] or res in antipar[chain_id][1]:
                    chain_struct[res]['type'] = 'E'
                    partner = self.partner_bridge(antipar[chain_id])
                    list_partner = [item for item in partner if item[0] == res]
                    chain_struct[res]['BP1'] = list_partner[0][1]
                    if len(list_partner) == 2:
                        chain_struct[res]['BP2'] = list_partner[1][1]

        return dico_res_struct

    @staticmethod
    def res_helixpi(dico_res, dico_res_struct, h_pi):
        """Add pi helix information to all residues in each chain.
        
        Parameters
        ----------
        dico_res : dict
            Includes all residues organized by chain.
        dico_res_struct : dict
            Includes structural information for all residues.
        h_pi : dict
            Each chain's residues in pi helices.

        Returns
        -------
        dict
            Updated `dico_res_struct` dictionary.
        """
        for chain_id, chain_residues in dico_res.items():
            if chain_id not in h_pi:
                continue
            for res in chain_residues.keys():
                if res in h_pi[chain_id]:
                    dico_res_struct[chain_id][res]['type'] = 'I'
        return dico_res_struct

    @staticmethod
    def res_helixH(dico_res, dico_res_struct, h_alpha):
        """逐链填充所有残基的结构信息，处理 alpha 螺旋信息。
        
        Parameters
        ----------
        dico_res : dict
            包含所有残基的字典，按链组织。
        dico_res_struct : dict
            包含所有残基结构信息的字典。
        h_alpha : dict
            每条链中 alpha 螺旋的残基列表。

        Returns
        -------
        dict
            更新后的 `dico_res_struct` 字典。
        """
        for chain_id, chain_residues in dico_res.items():
            if chain_id not in h_alpha:
                continue
            for res in chain_residues.keys():
                if res in h_alpha[chain_id]:
                    dico_res_struct[chain_id][res]['type'] = 'H'
        return dico_res_struct

    @staticmethod
    def res_helix310(dico_res, dico_res_struct, h_310):
        """Add 3,10 helix information to all residues in each chain.
        
        Parameters
        ----------
        dico_res : dict
            Includes all residues organized by chain.
        dico_res_struct : dict
            Includes structural information for all residues.
        h_310 : dict
            Each chain's residues in 3,10 helices.

        Returns
        -------
        dict
            Updated `dico_res_struct` dictionary.
        """
        for chain_id, chain_residues in dico_res.items():
            if chain_id not in h_310:
                continue
            for res in chain_residues.keys():
                if res in h_310[chain_id]:
                    dico_res_struct[chain_id][res]['type'] = 'G'
        return dico_res_struct

    def assign_structure(self):
        """Assigns secondary structure types based on hydrogen bonds and geometric properties."""
        list_4turn = self.list_nturn(4)
        list_3turn = self.list_nturn(3)
        list_5turn = self.list_nturn(5)

        # Assign helices based on turn data
        h_310 = self.helix_310(list_3turn)
        h_alpha = self.helix_alpha(list_4turn)
        h_pi = self.helix_pi(list_5turn)

        # Identify beta-sheets and bends
        par_bridge = self.parallel_bridge()
        anti_par_bridge = self.antiparallel_bridge()
        list_bend = self.bend()

        # Update structure dictionary with identified structures
        self.dico_res_struct = self.res_turnT(self.dico_res, self.dico_res_struct, list_4turn, list_3turn, list_5turn)
        self.dico_res_struct = self.res_bend(self.dico_res, self.dico_res_struct, list_bend)
        self.dico_res_struct = self.res_betasheetE(self.dico_res, self.dico_res_struct, par_bridge, anti_par_bridge)
        self.dico_res_struct = self.res_helix310(self.dico_res, self.dico_res_struct, h_310)
        self.dico_res_struct = self.res_helixH(self.dico_res, self.dico_res_struct, h_alpha)
        self.dico_res_struct = self.res_helixpi(self.dico_res, self.dico_res_struct, h_pi)

        return self.dico_res_struct



class CompareStructures:
    
    @staticmethod
    def extract_predicted_structure(dico_res_struct):
        """
        Extracts predicted secondary structure results from `dico_res_struct` for each chain.

        Parameters
        ----------
        dico_res_struct : dict
            Dictionary containing structural information for each chain

        Returns
        -------
        dict
            Dictionary with simplified secondary structure labels ('H', 'E', 'C') for each chain
        """
        predicted_structure = {}
        for chain_id, residues in dico_res_struct.items():
            chain_structure = []
            for res_id, res in residues.items():
                if isinstance(res, dict) and 'type' in res:
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

    @staticmethod
    def write_results(predicted_structure, total_accuracy, helix_sensitivity, beta_sensitivity, outputfile):
        """
        Write predicted secondary structures and final accuracy/sensitivity results to a specified file.

        Parameters:
        predicted_structure (dict): Dictionary with chain IDs as keys and predicted secondary structures as values.
        total_accuracy (float): The overall accuracy of the secondary structure prediction.
        helix_sensitivity (float): The sensitivity of helix detection.
        beta_sensitivity (float): The sensitivity of beta-sheet detection.
        outputfile (str): The full path to the output file where results will be saved.
        """
        try:
            with open(outputfile, 'w') as file:
                # Write the predicted structures by chain
                file.write("Predicted Secondary Structure by Chain:\n")
                for chain_id, structure in predicted_structure.items():
                    file.write(f"Chain {chain_id}:\n")
                    file.write(" ".join(structure) + "\n\n")
                
                # Write the overall comparison results
                file.write("Overall Comparison Results\n")
                file.write("-" * 30 + "\n")
                file.write(f"Overall Accuracy: {total_accuracy:.3f}\n")
                file.write(f"Helix Sensitivity: {helix_sensitivity:.3f}\n")
                file.write(f"Beta-Sheet Sensitivity: {beta_sensitivity:.3f}\n")
                
            print(f"Results successfully written to {outputfile}")
        except IOError as e:
            print(f"An error occurred while writing to {outputfile}: {e}")

    @staticmethod
    def get_true_structure_dssp(pdb_filename):
        """Retrieve true structure from DSSP for each chain independently."""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_filename)
        true_structure_per_chain = {}

        try:
            dssp = DSSP(structure[0], pdb_filename)
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

    @staticmethod
    def compare_structures(predicted_structure, true_structure):
        """Compare predicted and true structures per chain and calculate accuracy and sensitivity."""
        overall_predicted = []
        overall_true = []
        
        for chain_id in predicted_structure.keys():
            if chain_id not in true_structure:
                print(f"Chain {chain_id} missing in true structure data.")
                continue
            
            pred = predicted_structure[chain_id]
            true = true_structure[chain_id]

            if len(pred) != len(true):
                print(f"Length mismatch in chain {chain_id}: predicted={len(pred)}, true={len(true)}")
                continue

            overall_predicted.extend(pred)
            overall_true.extend(true)

            tp = sum(1 for i in range(len(pred)) if pred[i] == true[i])
            accuracy = tp / len(pred)

            tp_helix = sum(1 for i in range(len(pred)) if pred[i] == 'H' and true[i] == 'H')
            fn_helix = sum(1 for i in range(len(pred)) if pred[i] != 'H' and true[i] == 'H')
            helix_sensitivity = tp_helix / (tp_helix + fn_helix) if (tp_helix + fn_helix) > 0 else 0

            tp_beta = sum(1 for i in range(len(pred)) if pred[i] == 'E' and true[i] == 'E')
            fn_beta = sum(1 for i in range(len(pred)) if pred[i] != 'E' and true[i] == 'E')
            beta_sensitivity = tp_beta / (tp_beta + fn_beta) if (tp_beta + fn_beta) > 0 else 0

            print(f"Chain {chain_id} - Accuracy: {accuracy:.3f}, Helix Sensitivity: {helix_sensitivity:.3f}, Beta-Sheet Sensitivity: {beta_sensitivity:.3f}")

        total_tp = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] == overall_true[i])
        total_accuracy = total_tp / len(overall_predicted)

        total_tp_helix = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] == 'H' and overall_true[i] == 'H')
        total_fn_helix = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] != 'H' and overall_true[i] == 'H')
        total_helix_sensitivity = total_tp_helix / (total_tp_helix + total_fn_helix) if (total_tp_helix + total_fn_helix) > 0 else 0

        total_tp_beta = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] == 'E' and overall_true[i] == 'E')
        total_fn_beta = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] != 'E' and overall_true[i] == 'E')
        total_beta_sensitivity = total_tp_beta / (total_tp_beta + total_fn_beta) if (total_tp_beta + total_fn_beta) > 0 else 0

        print(f"\nOverall Protein - Accuracy: {total_accuracy:.3f}, Helix Sensitivity: {total_helix_sensitivity:.3f}, Beta-Sheet Sensitivity: {total_beta_sensitivity:.3f}")

        # Return the overall metrics for recording
        return total_accuracy, total_helix_sensitivity, total_beta_sensitivity




# ********************* MAIN ********************* #

import sys
import getopt

def main(argv):
    pdbfile = ''
    outputfile = ''
    
    # Parse command-line arguments
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('Usage: dssp_like.py -i <pdbfile> -o <outputfile>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print('Usage: dssp_like.py -i <pdbfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            pdbfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    print(f'Input file is {pdbfile}')
    print(f'Output file is {outputfile}')

    # Instantiate PDBProcessor and process the PDB file
    pdb_processor = PDBProcessor(pdbfile)
    dico_res = pdb_processor.parse_pdb()

    # Instantiate HydrogenBondCalculator and SecondaryStructureAssigner
    bond_calculator = HydrogenBondCalculator()
    structure_assigner = SecondaryStructureAssigner(dico_res, bond_calculator)

    # Assign secondary structures based on the parsed PDB data
    dico_res_struct = structure_assigner.assign_structure()

    # Extract and save the predicted structure
    predicted_structure = CompareStructures.extract_predicted_structure(dico_res_struct)

    # Retrieve true structure from DSSP and compare with predictions
    true_structure = CompareStructures.get_true_structure_dssp(pdb_processor.reduced_pdb_file)
    print(f"Total hydrogen bonds detected: {bond_calculator.hydrogen_bond_count}")
    print("Predicted Structure:", predicted_structure)
    print("True Structure:", true_structure)

    # Perform comparison
    CompareStructures.compare_structures(predicted_structure, true_structure)
    total_accuracy, helix_sensitivity, beta_sensitivity = CompareStructures.compare_structures(predicted_structure, true_structure) 
    CompareStructures.write_results(predicted_structure, total_accuracy, helix_sensitivity, beta_sensitivity, outputfile)  
    


if __name__ == "__main__":
    main(sys.argv[1:])
