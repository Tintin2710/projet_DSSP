import sys, getopt, math
import numpy as np
from itertools import groupby
from operator import itemgetter
import subprocess
from Bio.PDB import PDBParser, DSSP


#*********************GLOBAL VARIABLES*********************#

#List of targeted atoms ( electrostatic interaction )
interac_atoms = ["N", "C", "O", "H", "H1", "CA"]
q1, q2, f, Emin = 0.42, 0.20, 332, -0.5




#*********************FONCTIONS*********************#

def add_hydrogens_with_reduce(pdb_file):
    """Add hydrogens to PDB file using reduce tool."""
    reduced_pdb_file = "reduced_" + pdb_file
    with open(reduced_pdb_file, 'w') as output_file:
        subprocess.run(['reduce', '-HIS', pdb_file], stdout=output_file)
    print(f"Hydrogens added to: {reduced_pdb_file}")
    return reduced_pdb_file


def extract_atom_coord(line, dico_atom):
    """从 PDB 行中提取原子坐标、残基编号和名称"""
    atom = dico_atom
    atom['res_num'] = int(line[22:26])
    atom['res_name'] = line[17:20].strip()
    atom['x'] = float(line[30:38])
    atom['y'] = float(line[38:46])
    atom['z'] = float(line[46:54])
    return atom


def extract_first_resid(pdbfile):
    """提取第一个残基的编号"""
    with open(pdbfile, "r") as pdb_file:
        for line in pdb_file:
            atom_name = line[12:16].strip()
            if line.startswith("ATOM") and atom_name in interac_atoms:
                res_first = int(line[22:26])
                break
    return res_first


def parse_pdb(pdbfile):
    """Parse PDB file and store coordinates of each residue by chain."""
    with open(pdbfile, "r") as pdb_file:
        dico_res = {}
        for line in pdb_file:
            if line.startswith("ATOM") and line[12:16].strip() in interac_atoms:
                chain_id = line[21].strip()
                if chain_id not in dico_res:
                    dico_res[chain_id] = {}
                res_num = int(line[22:26])
                atom_name = line[12:16].strip()
                if res_num not in dico_res[chain_id]:
                    dico_res[chain_id][res_num] = {}
                atom = {}
                dico_res[chain_id][res_num][atom_name] = extract_atom_coord(line, atom)
    return dico_res

def dist2atoms(dico_res, numresA, numresB, atomA, atomB):
    """计算两个原子之间的距离，若不存在返回较大值"""
    if numresA not in dico_res or numresB not in dico_res:
        return float('inf')  # 返回无穷大，表示没有键
    if atomA not in dico_res[numresA] or atomB not in dico_res[numresB]:
        return float('inf')
    
    xA = dico_res[numresA][atomA]['x']
    yA = dico_res[numresA][atomA]['y']
    zA = dico_res[numresA][atomA]['z']
    xB = dico_res[numresB][atomB]['x']
    yB = dico_res[numresB][atomB]['y']
    zB = dico_res[numresB][atomB]['z']
    return math.sqrt((xB - xA) ** 2 + (yB - yA) ** 2 + (zB - zA) ** 2)


def is_Hbond(dico_res, numresA, numresB):
    """判断是否存在氢键，若不存在相关原子则忽略"""
    rON = dist2atoms(dico_res, numresA, numresB, 'O', 'N')
    rOH = dist2atoms(dico_res, numresA, numresB, 'O', 'H') or dist2atoms(dico_res, numresA, numresB, 'O', 'H1')
    rCN = dist2atoms(dico_res, numresA, numresB, 'C', 'N')
    rCH = dist2atoms(dico_res, numresA, numresB, 'C', 'H') or dist2atoms(dico_res, numresA, numresB, 'C', 'H1')
    if rON == float('inf') or rOH == float('inf') or rCN == float('inf') or rCH == float('inf'):
        return False  # 任何一个距离计算失败则认为没有氢键
    E_elec = q1 * q2 * f * (1 / rON + 1 / rCH - 1 / rOH - 1 / rCN)
    return E_elec < Emin


def dict_resid_structure(dico_res):
    """创建包含所有残基的结构标签的字典，确保每个残基有完整的结构信息."""
    res_struct = {}
    for chain_id, residues in dico_res.items():
        res_struct[chain_id] = {}
        for res_num, atoms in residues.items():
            struct = {
                'res_num': res_num,
                'res_name': atoms['N']['res_name'] if 'N' in atoms else "UNK",  # 如果'N'原子缺失，使用默认名称
                'type': "",
                'BP1': 0,
                'BP2': 0,
                '3T': "",
                '4T': "",
                '5T': ""
            }
            res_struct[chain_id][res_num] = struct
    return res_struct




def list_nturn(dico_res, nturn):
    """
    Find n-turns in each chain and return them separately for each chain.

    Parameters
    ----------
    dico_res : dict
        dict with all residues organized by chain
    nturn : int
        nturn = (3, 4, or 5)

    Returns
    -------
    dict
        A dictionary where each chain has a list of n-turns
    """
    chain_nturns = {}
    for chain_id, residues in dico_res.items():
        list_turn = []
        list_partner = []
        res_keys = sorted(residues.keys())
        for res in res_keys:
            if res + nturn in residues:
                if 'H' in residues[res + nturn]:
                    if is_Hbond(residues, res, res + nturn):
                        list_turn.append(res)
                        list_turn.append(res + nturn)
                        list_partner.append((res, res + nturn))
        chain_nturns[chain_id] = (list(set(list_turn)), list_partner)
    return chain_nturns



def successives_nb(liste):
    """ Finds succssives numbers in a given list

    Parameters
    ----------
    liste : list
        list of int

    Returns
    -------
    list
        list on int

    """
    succ_nb = []
    for k,g in groupby(enumerate(sorted(liste)),lambda x:x[0]-x[1]):
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        if len(group) > 1 :
            succ_nb = succ_nb + group
    return(succ_nb)


def helix_alpha(list_4turn):
    """Find residues in alpha helix for each chain."""
    helix_per_chain = {}
    for chain_id, turns in list_4turn.items():
        list_ha = successives_nb(turns[0])
        helix_per_chain[chain_id] = list_ha
    return helix_per_chain


def helix_310(list_3turn):
    """Find residues in 3,10 helix for each chain."""
    helix_310_per_chain = {}
    for chain_id, turns in list_3turn.items():
        list_h310 = successives_nb(turns[0])
        helix_310_per_chain[chain_id] = list_h310
    return helix_310_per_chain


def helix_pi(list_5turn):
    """Find residues in pi helix for each chain."""
    helix_pi_per_chain = {}
    for chain_id, turns in list_5turn.items():
        list_hpi = successives_nb(turns[0])
        helix_pi_per_chain[chain_id] = list_hpi
    return helix_pi_per_chain

def parallel_bridge(dico_res):
    """Find residues in parallel bridge for each chain."""
    parallel_per_chain = {}
    for chain_id, residues in dico_res.items():
        list_i, list_j = [], []
        res_keys = sorted(residues.keys())
        for i in res_keys:
            for j in res_keys:
                if i < j and j >= i + 3:
                    if i - 1 in residues and i + 1 in residues and j in residues:
                        cas1 = is_Hbond(residues, i - 1, j) and is_Hbond(residues, j, i + 1)
                    else:
                        cas1 = False
                    if i in residues and j - 1 in residues and j + 1 in residues:
                        cas2 = is_Hbond(residues, j - 1, i) and is_Hbond(residues, i, j + 1)
                    else:
                        cas2 = False
                    if cas1 or cas2:
                        list_i.append(i)
                        list_j.append(j)
        parallel_per_chain[chain_id] = [list_i, list_j]
    return parallel_per_chain


def antiparallel_bridge(dico_res):
    """Find residues in antiparallel bridge for each chain."""
    antiparallel_per_chain = {}
    for chain_id, residues in dico_res.items():
        list_i, list_j = [], []
        res_keys = sorted(residues.keys())
        for i in res_keys:
            for j in res_keys:
                if i < j and j >= i + 3:
                    if i in residues and j in residues:
                        cas1 = is_Hbond(residues, i, j) and is_Hbond(residues, j, i)
                    else:
                        cas1 = False
                    if i - 1 in residues and i + 1 in residues and j - 1 in residues and j + 1 in residues:
                        cas2 = is_Hbond(residues, i - 1, j + 1) and is_Hbond(residues, j - 1, i + 1)
                    else:
                        cas2 = False
                    if cas1 or cas2:
                        list_i.append(i)
                        list_j.append(j)
        antiparallel_per_chain[chain_id] = [list_i, list_j]
    return antiparallel_per_chain



def partner_bridge( list_ij ) :
    """Finds bridge partner(s) i,j

    Parameters
    ----------
    list_ij : list of list
        list of all residues i , list of all i bridge parteners j

    Returns
    -------
    list of tuples
        list of partners bridge (i,j)

    """
    partner = []
    list_i = list_ij[0]
    list_j = list_ij[1]
    for i in range(len(list_i)) :
        partner.append( (list_i[i] , list_j[i]))
    for j in range(len(list_j)) :
        partner.append(( list_j[j] , list_i[j]))
    return(partner)


def vect2atoms(dico_res, numresA, numresB, atomA, atomB) :
    """ Vector AB from coordinates of atom A and atom B

    Parameters
    ----------
    dico_res : dict
        dict with all residues
    numresA : int
        residue A number
    numresB : int
        residue B number
    atomA : str
        atom A name like 'C' , 'N' ...
    atomB : str
        atom B name like 'C' , 'N' ....

    Returns
    -------
    tuple
        vectAB( xB-xA , yB-yA ,zB-zA)

    """
    xA = dico_res[numresA][atomA]['x']
    yA = dico_res[numresA][atomA]['y']
    zA = dico_res[numresA][atomA]['z']
    xB = dico_res[numresB][atomB]['x']
    yB = dico_res[numresB][atomB]['y']
    zB = dico_res[numresB][atomB]['z']
    vectAB = ( xB-xA , yB-yA ,zB-zA)
    return(vectAB)



def unit_vector(vector):
    """Returns the norm of vector from the coordinates of vector


    Parameters
    ----------
    vector : tuple

    Returns
    -------
    numpy.ndarray
        norm of vector

    """
    return (vector / np.linalg.norm(vector))

def angle_between(v1, v2):
    """Returns the angle in degrees between vectors v1 and v2


    Parameters
    ----------
    v1 : tuple
    v2 : tuple

    Returns
    -------
    float
        angle in degrees

    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return (math.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))) )


def bend(dico_res):
    """Finds residues involved in bends for each chain."""
    list_bend_per_chain = {}
    for chain_id, residues in dico_res.items():
        list_bend = []
        res_keys = sorted(residues.keys())  # 获取该链的残基编号列表
        for i in range(2, len(res_keys) - 2):
            res = res_keys[i]
            v1 = vect2atoms(residues, res_keys[i - 2], res, 'CA', 'CA')
            v2 = vect2atoms(residues, res, res_keys[i + 2], 'CA', 'CA')
            if angle_between(v1, v2) > 70:
                list_bend.append(res)
        list_bend_per_chain[chain_id] = list_bend
    return list_bend_per_chain


def dict_resid_structure(dico_res):
    """创建包含所有残基的结构标签的字典，每条链单独处理。"""
    res_struct_per_chain = {}
    for chain_id, residues in dico_res.items():
        chain_struct = {}
        for res_num, atoms in residues.items():
            struct = {}
            struct['res_num'] = atoms['N']['res_num'] if 'N' in atoms else res_num
            struct['res_name'] = atoms['N']['res_name'] if 'N' in atoms else ""
            struct['type'] = ""
            struct['BP1'] = 0
            struct['BP2'] = 0
            struct['3T'] = ""
            struct['4T'] = ""
            struct['5T'] = ""
            chain_struct[res_num] = struct
        res_struct_per_chain[chain_id] = chain_struct
    return res_struct_per_chain



def res_bend ( dico_res ,  dico_res_struct , list_bend) :
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
    return( dico_res_struct)

def res_turnT(dico_res, dico_res_struct, list_4turn, list_3turn, list_5turn):
    """逐链处理，将结构信息填入所有残基字典中"""
    for chain_id in dico_res_struct.keys():
        if chain_id not in list_4turn or chain_id not in list_3turn or chain_id not in list_5turn:
            continue  # 跳过没有数据的链
        
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


def res_betasheetE(dico_res, dico_res_struct, par, antipar):
    """逐链填充所有残基的结构信息，处理平行和反平行的β片层信息。
    
    Parameters
    ----------
    dico_res : dict
        包含所有残基的字典。
    dico_res_struct : dict
        包含所有残基结构信息的字典。
    par : dict
        包含每条链的平行桥接信息。
    antipar : dict
        包含每条链的反平行桥接信息。

    Returns
    -------
    dict
        更新后的 `dico_res_struct` 字典。
    """
    for chain_id in dico_res_struct.keys():
        if chain_id not in par or chain_id not in antipar:
            continue  # 跳过没有桥接信息的链

        chain_struct = dico_res_struct[chain_id]
        
        # 处理平行桥接
        for res in chain_struct.keys():
            if res in par[chain_id][0] or res in par[chain_id][1]:
                chain_struct[res]['type'] = 'E'
                partner = partner_bridge(par[chain_id])
                list_partner = [item for item in partner if item[0] == res]
                chain_struct[res]['BP1'] = list_partner[0][1]
                if len(list_partner) == 2:
                    chain_struct[res]['BP2'] = list_partner[1][1]
                    
        # 处理反平行桥接
        for res in chain_struct.keys():
            if res in antipar[chain_id][0] or res in antipar[chain_id][1]:
                chain_struct[res]['type'] = 'E'
                partner = partner_bridge(antipar[chain_id])
                list_partner = [item for item in partner if item[0] == res]
                chain_struct[res]['BP1'] = list_partner[0][1]
                if len(list_partner) == 2:
                    chain_struct[res]['BP2'] = list_partner[1][1]

    return dico_res_struct


def res_helixpi(dico_res, dico_res_struct, h_pi):
    """逐链填充所有残基的结构信息，处理 pi 螺旋信息。
    
    Parameters
    ----------
    dico_res : dict
        包含所有残基的字典，按链组织。
    dico_res_struct : dict
        包含所有残基结构信息的字典。
    h_pi : dict
        每条链中 pi 螺旋的残基列表。

    Returns
    -------
    dict
        更新后的 `dico_res_struct` 字典。
    """
    for chain_id, chain_residues in dico_res.items():
        if chain_id not in h_pi:
            continue
        for res in chain_residues.keys():
            if res in h_pi[chain_id]:
                dico_res_struct[chain_id][res]['type'] = 'I'
    return dico_res_struct


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


def res_helix310(dico_res, dico_res_struct, h_310):
    """逐链填充所有残基的结构信息，处理 3,10 螺旋信息。
    
    Parameters
    ----------
    dico_res : dict
        包含所有残基的字典，按链组织。
    dico_res_struct : dict
        包含所有残基结构信息的字典。
    h_310 : dict
        每条链中 3,10 螺旋的残基列表。

    Returns
    -------
    dict
        更新后的 `dico_res_struct` 字典。
    """
    for chain_id, chain_residues in dico_res.items():
        if chain_id not in h_310:
            continue
        for res in chain_residues.keys():
            if res in h_310[chain_id]:
                dico_res_struct[chain_id][res]['type'] = 'G'
    return dico_res_struct



def affichage( dico_res_struct , outputfile ) :
    """Write structures assignation for each residue

    Parameters
    ----------
    dico_res_struct : dict
        dico_res_struct


    """
    with open(outputfile, 'w') as fillout:
        fillout.write("{:>3s}{:>4s}{:>12s}{:>4s}{:>4s}\n".format( "RES", "AA" , "STRUCTURE" , "BP1", "BP2"))
        for res in dico_res_struct :
            r = dico_res_struct[res]
            fillout.write("{:>3d}{:>4s}{:>3s}{:>3s}{:>3s}{:>3s}{:>4d}{:>4d}\n".format(res, r['res_name'] , r['type'], r['3T'] , r['4T'] , r['5T'] , r['BP1'] , r['BP2']))

def extract_predicted_structure(dico_res_struct):
    """
    从 dico_res_struct 中提取二级结构预测结果，支持每条链。

    Parameters
    ----------
    dico_res_struct : dict
        包含每条链的结构信息的字典

    Returns
    -------
    dict
        包含每条链简化二级结构标签的字典（'H', 'E', 'C'）
    """
    predicted_structure = {}
    for chain_id, residues in dico_res_struct.items():
        chain_structure = []
        for res_id, res in residues.items():
            if isinstance(res, dict) and 'type' in res:  # 确保 res 是字典且包含 'type'
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




def write_predicted_structure(predicted_structure, filename="predicted_structure.txt"):
    """
    Write predicted secondary structures to a file by chain.
    """
    with open(filename, 'w') as file:
        file.write("Predicted Secondary Structure by Chain:\n")
        for chain_id, structure in predicted_structure.items():
            file.write(f"Chain {chain_id}:\n")
            file.write(" ".join(structure) + "\n\n")
    print(f"Predicted secondary structures saved to: {filename}")


def get_true_structure_dssp(pdb_filename):
    """Retrieve true structure from DSSP for each chain independently."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_filename)
    true_structure_per_chain = {}

    try:
        dssp = DSSP(structure[0], pdb_filename)  # Using the first model in the structure
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


def compare_structures(predicted_structure, true_structure):
    """Compare predicted and true structures per chain and calculate accuracy and sensitivity."""

    # 用于累积所有链的预测和真实结构
    overall_predicted = []
    overall_true = []
    
    # 每条链的准确率和敏感度计算
    for chain_id in predicted_structure.keys():
        if chain_id not in true_structure:
            print(f"Chain {chain_id} missing in true structure data.")
            continue
        
        pred = predicted_structure[chain_id]
        true = true_structure[chain_id]

        if len(pred) != len(true):
            print(f"Length mismatch in chain {chain_id}: predicted={len(pred)}, true={len(true)}")
            continue

        # 累积链结果
        overall_predicted.extend(pred)
        overall_true.extend(true)

        # 计算每条链的准确性
        tp = sum(1 for i in range(len(pred)) if pred[i] == true[i])
        accuracy = tp / len(pred)

        # 计算螺旋灵敏度
        tp_helix = sum(1 for i in range(len(pred)) if pred[i] == 'H' and true[i] == 'H')
        fn_helix = sum(1 for i in range(len(pred)) if pred[i] != 'H' and true[i] == 'H')
        helix_sensitivity = tp_helix / (tp_helix + fn_helix) if (tp_helix + fn_helix) > 0 else 0

        # 计算β片层灵敏度
        tp_beta = sum(1 for i in range(len(pred)) if pred[i] == 'E' and true[i] == 'E')
        fn_beta = sum(1 for i in range(len(pred)) if pred[i] != 'E' and true[i] == 'E')
        beta_sensitivity = tp_beta / (tp_beta + fn_beta) if (tp_beta + fn_beta) > 0 else 0

        print(f"Chain {chain_id} - Accuracy: {accuracy}, Helix Sensitivity: {helix_sensitivity}, Beta-Sheet Sensitivity: {beta_sensitivity}")

    # 计算整个蛋白的总体准确率和灵敏度
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

    # 调用 add_hydrogens_with_reduce 函数
    reduced_pdb_file = add_hydrogens_with_reduce(pdbfile)

    # 使用 reduced PDB 文件解析和分析
    dico_res = parse_pdb(reduced_pdb_file)
    list_4turn = list_nturn(dico_res, 4)
    list_3turn = list_nturn(dico_res, 3)
    list_5turn = list_nturn(dico_res, 5)
    h_310 = helix_310(list_3turn)
    h_alpha = helix_alpha(list_4turn)
    h_pi = helix_pi(list_5turn)
    par = parallel_bridge(dico_res)
    antipar = antiparallel_bridge(dico_res)
    list_bend = bend(dico_res)
    dico_res_struct = dict_resid_structure(dico_res)
    dico_res_struct = res_turnT(dico_res, dico_res_struct, list_4turn, list_3turn, list_5turn)
    dico_res_struct = res_bend(dico_res, dico_res_struct, list_bend)
    dico_res_struct = res_betasheetE(dico_res, dico_res_struct, par, antipar)
    dico_res_struct = res_helix310(dico_res, dico_res_struct, h_310)
    dico_res_struct = res_helixH(dico_res, dico_res_struct, h_alpha)
    dico_res_struct = res_helixpi(dico_res, dico_res_struct, h_pi)

    # 提取自定义预测结果
    predicted_structure = extract_predicted_structure(dico_res_struct)

    # 保存预测的二级结构到文本文件
    write_predicted_structure(predicted_structure, "predicted_structure.txt")

    # 提取 DSSP 生成的真实结构
    true_structure = get_true_structure_dssp(reduced_pdb_file)
    print("Predicted Structure ", predicted_structure)
    print("True Structure ", true_structure)

    # 调用比较函数
    compare_structures(predicted_structure, true_structure)


if __name__ == "__main__":
   main(sys.argv[1:])