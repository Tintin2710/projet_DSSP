import math
import numpy as np
from itertools import groupby
from operator import itemgetter

class SecondaryStructureAssigner:
    def __init__(self, dico_res, hydrogen_bond_finder):
        self.dico_res = dico_res
        self.hydrogen_bond_finder = hydrogen_bond_finder
        self.dico_res_struct = self.dict_resid_structure(dico_res)

    def dict_resid_structure(self, dico_res):
        """Creates a dictionary for each residue to hold structural information."""
        res_struct = {}
        for chain_id, residues in dico_res.items():
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
        """Identifies n-turns in chains based on hydrogen bonds."""
        chain_nturns = {}
        for chain_id, residues in self.dico_res.items():
            list_turn = []
            list_partner = []
            res_keys = sorted(residues.keys())
            for res in res_keys:
                if res + nturn in residues:
                    if self.hydrogen_bond_finder.is_hydrogen_bond(res, res + nturn):
                        list_turn.append(res)
                        list_turn.append(res + nturn)
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
        """Identifies alpha helices."""
        return {chain_id: self.successives_nb(turns[0]) for chain_id, turns in list_4turn.items()}

    def helix_310(self, list_3turn):
        """Identifies 3,10 helices."""
        return {chain_id: self.successives_nb(turns[0]) for chain_id, turns in list_3turn.items()}

    def helix_pi(self, list_5turn):
        """Identifies pi helices."""
        return {chain_id: self.successives_nb(turns[0]) for chain_id, turns in list_5turn.items()}

    def bridge(self, is_parallel=True):
        """Finds residues in parallel or antiparallel bridges."""
        bridge_per_chain = {}
        for chain_id, residues in self.dico_res.items():
            list_i, list_j = [], []
            res_keys = sorted(residues.keys())
            for i in res_keys:
                for j in res_keys:
                    if i < j and j >= i + 3:
                        conditions = [
                            self.hydrogen_bond_finder.is_hydrogen_bond(i, j) and self.hydrogen_bond_finder.is_hydrogen_bond(j, i)
                        ]
                        if is_parallel:
                            conditions.extend([
                                self.hydrogen_bond_finder.is_hydrogen_bond(i - 1, j) and self.hydrogen_bond_finder.is_hydrogen_bond(j, i + 1)
                            ])
                        else:
                            conditions.extend([
                                self.hydrogen_bond_finder.is_hydrogen_bond(i - 1, j + 1) and self.hydrogen_bond_finder.is_hydrogen_bond(j - 1, i + 1)
                            ])
                        if any(conditions):
                            list_i.append(i)
                            list_j.append(j)
            bridge_per_chain[chain_id] = [list_i, list_j]
        return bridge_per_chain

    @staticmethod
    def partner_bridge(list_ij):
        """Returns bridge partners for each residue pair in a list."""
        return [(i, j) for i, j in zip(list_ij[0], list_ij[1])] + [(j, i) for i, j in zip(list_ij[0], list_ij[1])]

    @staticmethod
    def vect2atoms(dico_res, numresA, numresB, atomA, atomB):
        """Calculates vector AB from coordinates of atom A and atom B."""
        xA = dico_res[numresA][atomA]['x']
        yA = dico_res[numresA][atomA]['y']
        zA = dico_res[numresA][atomA]['z']
        xB = dico_res[numresB][atomB]['x']
        yB = dico_res[numresB][atomB]['y']
        zB = dico_res[numresB][atomB]['z']
        vectAB = (xB - xA, yB - yA, zB - zA)
        return vectAB

    @staticmethod
    def unit_vector(vector):
        """Returns the unit vector from the coordinates of a vector."""
        return vector / np.linalg.norm(vector)

    @staticmethod
    def angle_between(v1, v2):
        """Calculates the angle in degrees between vectors v1 and v2."""
        v1_u = SecondaryStructureAssigner.unit_vector(v1)
        v2_u = SecondaryStructureAssigner.unit_vector(v2)
        return math.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))

    @staticmethod
    def bend(dico_res):
        """Identifies bends based on the angle between CA atoms."""
        list_bend_per_chain = {}
        for chain_id, residues in dico_res.items():
            list_bend = []
            res_keys = sorted(residues.keys())
            for i in range(2, len(res_keys) - 2):
                res = res_keys[i]
                v1 = SecondaryStructureAssigner.vect2atoms(residues, res_keys[i - 2], res, 'CA', 'CA')
                v2 = SecondaryStructureAssigner.vect2atoms(residues, res, res_keys[i + 2], 'CA', 'CA')
                if SecondaryStructureAssigner.angle_between(v1, v2) > 70:
                    list_bend.append(res)
            list_bend_per_chain[chain_id] = list_bend
        return list_bend_per_chain

    def assign_structure(self):
        """Assigns secondary structure types to residues based on hydrogen bonds and geometry."""
        list_3turn = self.list_nturn(3)
        list_4turn = self.list_nturn(4)
        list_5turn = self.list_nturn(5)
        
        h_alpha = self.helix_alpha(list_4turn)
        h_310 = self.helix_310(list_3turn)
        h_pi = self.helix_pi(list_5turn)

        for chain_id in self.dico_res_struct.keys():
            for res in h_alpha.get(chain_id, []):
                self.dico_res_struct[chain_id][res]['type'] = 'H'
            for res in h_310.get(chain_id, []):
                self.dico_res_struct[chain_id][res]['type'] = 'G'
            for res in h_pi.get(chain_id, []):
                self.dico_res_struct[chain_id][res]['type'] = 'I'
        
        par_bridge = self.bridge(is_parallel=True)
        anti_par_bridge = self.bridge(is_parallel=False)

        for chain_id in self.dico_res_struct.keys():
            for res in par_bridge.get(chain_id, [])[0]:
                self.dico_res_struct[chain_id][res]['type'] = 'E'
            for res in anti_par_bridge.get(chain_id, [])[0]:
                self.dico_res_struct[chain_id][res]['type'] = 'E'
        
        return self.dico_res_struct
