from read_pdb_structure import ReadPdbStructure
from hydrogen_bond_detector import HydrogenBondDetector
from secondary_structure_assigner import SecondaryStructureAssigner
from comparator import Compare

def main(pdb_filename):
    # Main program logic remains the same
    pdb_structure = ReadPdbStructure(pdb_filename)
    main_chain_atoms = pdb_structure.get_main_chain_atoms()
    print(f"Number of main chain atoms extracted: {len(main_chain_atoms)}")
    
    hb_detector = HydrogenBondDetector(main_chain_atoms)
    hbonds = hb_detector.calculate_hydrogen_bonds()
    print(f"Number of hydrogen bonds detected: {len(hbonds)}")
    
    ss_assigner = SecondaryStructureAssigner(main_chain_atoms, hbonds)
    predicted_structure = ss_assigner.assign_secondary_structure()
    
    evaluator = Compare(pdb_filename, predicted_structure)
    accuracy, helix_sensitivity, beta_sensitivity = evaluator.compare()
    print(f"Accuracy: {accuracy:.2f}")
    print(f"Helix Sensitivity: {helix_sensitivity:.2f}")
    print(f"Beta Sensitivity: {beta_sensitivity:.2f}")

if __name__ == "__main__":
    pdb_filename = 'pdbfile'  # Replace with your actual PDB file path
    main(pdb_filename)
