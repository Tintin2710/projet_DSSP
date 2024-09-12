from read_pdb_structure import ReadPdbStructure
from hydrogen_bond_detector import HydrogenBondDetector
from secondary_structure_assigner import SecondaryStructureAssigner
from comparator import Compare

# Main program to run the secondary structure assignment and evaluation
def main(pdb_filename):
    """
    Main function to run the secondary structure prediction and comparison with DSSP.

    Parameters:
    pdb_filename (str): The path to the PDB file to be processed.
    """
    # 1. Parse the PDB structure and extract main chain atoms
    pdb_structure = ReadPdbStructure(pdb_filename)
    main_chain_atoms = pdb_structure.get_main_chain_atoms()

    # Debug: print the number of main chain atoms extracted
    print(f"Number of main chain atoms extracted: {len(main_chain_atoms)}")

    # 2. Detect hydrogen bonds
    hb_detector = HydrogenBondDetector(main_chain_atoms)
    hbonds = hb_detector.calculate_hydrogen_bonds()

    # Debug: print the number of hydrogen bonds detected
    print(f"Number of hydrogen bonds detected: {len(hbonds)}")

    # 3. Assign secondary structures based on hydrogen bonds
    ss_assigner = SecondaryStructureAssigner(main_chain_atoms, hbonds)
    predicted_structure = ss_assigner.assign_secondary_structure()

    # 4. Compare with DSSP results and calculate accuracy and sensitivity
    evaluator = Compare(pdb_filename, predicted_structure)
    accuracy, helix_sensitivity, beta_sensitivity = evaluator.compare()

    # Print the evaluation results
    print(f"Accuracy: {accuracy:.2f}")
    print(f"Helix Sensitivity: {helix_sensitivity:.2f}")
    print(f"Beta-Sheet Sensitivity: {beta_sensitivity:.2f}")

if __name__ == "__main__":
    pdb_filename = pdb_file  # Replace with your actual PDB file path
    main(pdb_filename)
