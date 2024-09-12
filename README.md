# projet

## Description
This project is a Python-based tool for predicting and evaluating the secondary structure of proteins using hydrogen bond detection and energy calculations. The tool compares the predicted secondary structure with the DSSP-assigned secondary structure to calculate accuracy, helix sensitivity, and beta-sheet sensitivity.

## Features
- **Main chain Atom Extraction.** : The file `read_pdb_structure.py` contains the ReadPdbStructure class which can extract the main chain atoms (N, CA, C, O) from PDB file.
- **Hydrogen Bond Detection** : The file `hydrogen_bond_detector.py` contains the HydrogenBondDetector class identifies potential hydrogen bonds between main chain atoms based on energy thresholds and distances between donors and accpetors.
- **Secondary Structure Assignment** :The file `secondary_structure_assigner.py` contains the SecondaryStructureAssigner class, predicts secondary structures, helix and beta sheets, based on hydrogen bond patterns and geometric features following the Kabsch & Sander method.
- **DSSP Comparison** : The file `comparator.py` for compareing the predicted secondary structures with those provided by DSSP and calculates accuracy, helix sensitivity and beta-sheet sensitivity.
- **Main code**: The main script `main.py` that imports and utilizes all the classes mentioned above to run the complete process.

## Requirement
To run these scripts, the file `environment.yml` can create identical development environments and ensure that all dependencies are aligned with project requirements.

## Usage
This script processes a PDB file to extract main chain atoms, detect hydrogen bonds, assign secondary structures, and compare the results with DSSP. Follow the steps below to use the script:

**Install Dependencies**: Make sure the necessary Python packages are installed. Use either conda with the environment.yml file
for `conda`:
`conda env create -f environment.yml`

**Prepare the PDB File**: Ensure you have the PDB file you want to analyze. Place it in the appropriate directory.

**Run the Script**: To run the script, provide the path to the PDB file as an argument. You can modify the pdb_filename variable directly in `main.py` or pass it dynamically.

Here’s an example:
```python
if __name__ == "__main__":
    pdb_filename = "path/to/your/pdbfile.pdb"  # Replace with your actual PDB file path
    output_filename = "file.txt"
    main(pdb_filename, output_filename)
```

**Output**
The script will print the following information:

-Number of main chain atoms extracted.

-Number of hydrogen bonds detected.

-Predicted structre.

-Ture structure by DSSP.

-Accuracy

-Helix sensitivity

-Beta-sheet sensitivity

