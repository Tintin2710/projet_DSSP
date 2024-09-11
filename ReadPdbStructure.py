import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.DSSP import DSSP

class ReadPdbStructure:
    """Represent and extract the main chain atoms (N, CA, C, O) in the PDB structure."""

    def __init__(self, pdb_filename):
        self.pdb_filename = pdb_filename
        self.structure = self._read_pdb_structure()
        self.main_chain_atoms = self._get_main_chain_atoms()

    def _read_pdb_structure(self):
        """Lire le fichier PDB et renvoyer la structure"""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", self.pdb_filename)
        return structure
    
    def _get_main_chain_atoms(self):
        """Extract the main chain atomic coordinates of each residue (N, CA, C, O)"""
        main_chain_atoms = []
        previous_res_id = None
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue) and residue.has_id('N') and residue.has_id('CA') and residue.has_id('C') and residue.has_id('O'):
                        res_id = residue.get_id()[1]
                        chain_id = chain.get_id()
                        if previous_res_id is not None and abs(res_id - previous_res_id) > 1:
                            print(f"Non-continuous residue IDs: previous_res_id={previous_res_id}, current_res_id={res_id}, chain={chain_id}")    
                        previous_res_id = res_id

                        main_chain_atoms.append({
                            'residue': residue.get_resname(),
                            'res_id': res_id,
                            'chain': chain_id,
                            'N': residue['N'].get_coord(),
                            'CA': residue['CA'].get_coord(),
                            'C': residue['C'].get_coord(),
                            'O': residue['O'].get_coord()
                        })
        return main_chain_atoms
    
    def get_main_chain_atoms(self):
        """Obtenir des informations atomiques sur la chaîne principale"""
        return self.main_chain_atoms  
    

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
            donor_residue_id = donor_residue['res_id']
            for j in range(i + 1, num_residues):  # Iterate over acceptor residues
                acceptor_residue = self.main_chain_atoms[j]
                acceptor_chain_id = acceptor_residue.get('chain_id', None)  # 获取 acceptor 的链号
                acceptor_residue_id = acceptor_residue['res_id']
                
                if donor_chain_id != acceptor_chain_id:
                    #print(f"跨链氢键跳过: donor_chain_id={donor_chain_id}, acceptor_chain_id={acceptor_chain_id}")
                    continue

                diff = abs(donor_residue_id - acceptor_residue_id)

                if diff > 220:  
                    #print(f"远距离氢键跳过: donor_residue_id={donor_residue_id}, acceptor_residue_id={acceptor_residue_id}, diff={diff}")
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
                        'diff': diff
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


class SecondaryStructureAssigner:
    """Secondary Structure Assigner, detecting α-Helix, β-Folding, and β-Bridges based on hydrogen bonds"""

    def __init__(self, main_chain_atoms, hbonds):
        self.hbonds = hbonds
        self.main_chain_atoms = main_chain_atoms

    def assign_secondary_structure(self):
        helix = set()
        beta = set()
    
        res_id_to_index = {res['res_id']: idx for idx, res in enumerate(self.main_chain_atoms)}

        for bond in self.hbonds:

            donor_residue_id = bond['donor']['res_id']
            acceptor_residue_id = bond['acceptor']['res_id']
            diff = abs(donor_residue_id - acceptor_residue_id)

            print(f"Processing bond: donor_residue_id={donor_residue_id}, acceptor_residue_id={acceptor_residue_id}, diff={abs(donor_residue_id - acceptor_residue_id)}, energy = {bond['energy']}")

            # 优化后的螺旋判定条件，增加对氢键能量和角度的筛选
            if diff in [3,4,5]:
                helix.add(donor_residue_id)
                helix.add(acceptor_residue_id)

            # β 折叠判定条件
            elif 5 < diff <= 200 and bond['r_ON'] < 4.1 and bond['energy'] < -0.4:
                beta.add(donor_residue_id)
                beta.add(acceptor_residue_id)

        # 初始化所有残基为 'C'（环）
        secondary_structure = ['C'] * len(self.main_chain_atoms)

        # 分配 α-螺旋 ('H') 和 β-折叠 ('E')，通过res_id找到对应的索引
        for res_id in helix:
            if res_id in res_id_to_index:
                index = res_id_to_index[res_id]
                secondary_structure[index] = 'H'  # 先标记为螺旋

        for res_id in beta:
            if res_id in res_id_to_index:
                index = res_id_to_index[res_id]
                # 只有当残基没有被标记为螺旋时才标记为折叠
                if secondary_structure[index] != 'H':
                    secondary_structure[index] = 'E'

        # 调试信息：打印主链原子数量和分配结果
        print("Length of main_chain_atoms:", len(self.main_chain_atoms))
        print("Helix indices:", helix)
        print("Beta indices:", beta)

        return secondary_structure


from Bio.PDB import PDBParser, DSSP

class Compare:
    """比较自定义二级结构与 DSSP 分配结果的类"""

    def __init__(self, pdb_filename, predicted_structure):
        self.pdb_filename = pdb_filename
        self.predicted_structure = predicted_structure
        self.true_structure = self._dssp_structure()

    def _dssp_structure(self):
        """从 PDB 文件中提取 DSSP 分配的二级结构。"""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", self.pdb_filename)
        
        try:
            dssp = DSSP(structure[0], self.pdb_filename)
        except Exception as e:
            print(f"DSSP 解析错误: {e}")
            return []

        dssp_structure = []
        for key in dssp.keys():
            dssp_code = dssp[key][2]  # DSSP 的二级结构代码
            if dssp_code in ('H', 'G', 'I'):  # 螺旋
                dssp_structure.append('H')
            elif dssp_code in ('E', 'B'):  # β 折叠
                dssp_structure.append('E')
            else:
                dssp_structure.append('C')  # 环状结构
        return dssp_structure

    def compare(self):
        """比较预测的二级结构和 DSSP 结果，计算准确率和灵敏度。"""
        if len(self.predicted_structure) != len(self.true_structure):
            print("预测结构和真实结构长度不匹配，无法比较。")
            return None

        print(f"预测结构: {self.predicted_structure}")
        print(f"真实结构: {self.true_structure}")

        # 计算准确率
        tp = sum([1 for i in range(len(self.predicted_structure)) if self.predicted_structure[i] == self.true_structure[i]])
        accuracy = tp / len(self.predicted_structure)

        # 计算螺旋灵敏度
        total_helices_in_dssp = sum([1 for i in range(len(self.true_structure)) if self.true_structure[i] == 'H'])
        helix_sensitivity = (sum([1 for i in range(len(self.predicted_structure)) 
                                  if self.predicted_structure[i] == 'H' and self.true_structure[i] == 'H']) 
                             / total_helices_in_dssp) if total_helices_in_dssp > 0 else 0
        
        # 计算β折叠灵敏度
        total_beta_in_dssp = sum([1 for i in range(len(self.true_structure)) if self.true_structure[i] == 'E'])
        beta_sensitivity = (sum([1 for i in range(len(self.predicted_structure)) 
                                 if self.predicted_structure[i] == 'E' and self.true_structure[i] == 'E']) 
                            / total_beta_in_dssp) if total_beta_in_dssp > 0 else 0

        print(f"螺旋灵敏度: {helix_sensitivity}")
        print(f"β 折叠灵敏度: {beta_sensitivity}")

        return accuracy, helix_sensitivity, beta_sensitivity

# Main program to run the secondary structure assignment and evaluation
def main(pdb_filename):
    """Main function to run the secondary structure prediction and comparison with DSSP.
    
    Args:
        pdb_filename (str): The path to the PDB file.
    """
    # 1. Parse the PDB structure and extract main chain atoms
    pdb_structure = ReadPdbStructure(pdb_filename)
    main_chain_atoms = pdb_structure.get_main_chain_atoms()

    # Debug: print the number of atoms extracted
    print(f"Number of main chain atoms extracted: {len(main_chain_atoms)}")

    # 2. Detect hydrogen bonds
    hb_detector = HydrogenBondDetector(main_chain_atoms)
    hbonds = hb_detector.calculate_hydrogen_bonds()

    # Debug: print the number of hydrogen bonds detected
    print(f"Number of hydrogen bonds detected: {len(hbonds)}")

    # 3. Assign secondary structures based on hydrogen bonds
    # Correct parameter order: main_chain_atoms first, then hbonds
    ss_assigner = SecondaryStructureAssigner(main_chain_atoms, hbonds)
    predicted_structure = ss_assigner.assign_secondary_structure()

    # 4. Compare with DSSP results and calculate accuracy and sensitivity
    evaluator = Compare(pdb_filename, predicted_structure)
    accuracy, helix_sensitivity, beta_sensitivity = evaluator.compare()

    # Print the results
    print(f"Accuracy: {accuracy:.2f}")
    print(f"Helix Sensitivity: {helix_sensitivity:.2f}")
    print(f"Beta Sensitivity: {beta_sensitivity:.2f}")


if __name__ == "__main__":
    pdb_filename = 'pdbfile'  # Replace with your actual PDB file path
    main(pdb_filename)