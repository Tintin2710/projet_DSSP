import sys
import getopt
from hydrogen_bond_detector import HydrogenBondFinder
from secondary_structure_assigner import SecondaryStructureAssigner
from comparator import StructureComparator

def main(argv):
    pdbfile = ''
    outputfile = ''

    # Parse command line arguments
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('Usage: main.py -i <pdbfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Usage: main.py -i <pdbfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            pdbfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    # Display input and output file names
    print('Input file is:', pdbfile)
    print('Output file is:', outputfile)

    # Step 1: Initialize hydrogen bond finder and get residues with hydrogen bond information
    hb_finder = HydrogenBondFinder(pdbfile)
    dico_res = hb_finder.dico_res  # Dictionary of residues parsed from the PDB file with hydrogen bonds

    # Step 2: Initialize the secondary structure assigner with hydrogen bond data and assign structures
    structure_assigner = SecondaryStructureAssigner(dico_res, hb_finder)
    assigned_structure = structure_assigner.assign_structure()

    # Step 3: Write the predicted structure to output file
    with open(outputfile, 'w') as file:
        for chain_id, residues in assigned_structure.items():
            file.write(f"Chain {chain_id}:\n")
            for res_num, res_info in residues.items():
                file.write(f"{res_num}: {res_info['type']}\n")

    # Step 4: Compare predicted structures with DSSP-generated true structures
    comparator = StructureComparator(pdbfile, structure_assigner.dico_res_struct)
    comparator.compare_structures()

if __name__ == "__main__":
    main(sys.argv[1:])

