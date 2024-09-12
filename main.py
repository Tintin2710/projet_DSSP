from read_pdb_structure import ReadPdbStructure
from hydrogen_bond_detector import HydrogenBondDetector
from secondary_structure_assigner import SecondaryStructureAssigner
from comparator import Compare

# Main program to run the secondary structure assignment and evaluation
def main(pdb_filename, output_filename="output.txt"):
    """
    Main function to run the secondary structure prediction and comparison with DSSP,
    and export the results to a file.

    Parameters:
    pdb_filename (str): The path to the PDB file to be processed.
    output_filename (str): The path to the output file where results will be saved.
    """
    with open(output_filename, 'w') as f:
        # 1. Parse the PDB structure and extract main chain atoms
        pdb_structure = ReadPdbStructure(pdb_filename)
        main_chain_atoms = pdb_structure.get_main_chain_atoms()

        # Write the number of main chain atoms to the file
        f.write(f"Number of main chain atoms extracted: {len(main_chain_atoms)}\n")

        # 2. Detect hydrogen bonds
        hb_detector = HydrogenBondDetector(main_chain_atoms)
        hbonds = hb_detector.calculate_hydrogen_bonds()

        # Write the number of hydrogen bonds detected to the file
        f.write(f"Number of hydrogen bonds detected: {len(hbonds)}\n")

        # 3. Assign secondary structures based on hydrogen bonds
        ss_assigner = SecondaryStructureAssigner(main_chain_atoms, hbonds)
        predicted_structure = ss_assigner.assign_secondary_structure()

        # 4. Compare with DSSP results and calculate accuracy and sensitivity
        evaluator = Compare(pdb_filename, predicted_structure)
        accuracy, helix_sensitivity, beta_sensitivity = evaluator.compare()

        # Write the evaluation results to the file
        f.write(f"Accuracy: {accuracy:.2f}\n")
        f.write(f"Helix Sensitivity: {helix_sensitivity:.2f}\n")
        f.write(f"Beta-Sheet Sensitivity: {beta_sensitivity:.2f}\n")

if __name__ == "__main__":
    pdb_filename = "path/to/file"  # Replace with your actual PDB file path
    output_filename = "file.txt"  # Define the output file name
    main(pdb_filename, output_filename)
