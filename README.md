# DSSP 

## Description
This is a Python program to implement the DSSP **(Dictionary of Secondary Structures of Proteins)** method, implementinger its ability to use atomic coordinate data to identify hydrogen bonds and assign secondary structure elements to protein residues. Finally the predicted secondary structure is compared to the secondary structure assigned by DSSP to derive accuracy, helix sensitivity and β-sheet sensitivity.

## Features


## Requirement
- To run these scripts, the file `DSSPenv.yml` can create identical development environments and ensure that all dependencies are aligned with project requirements.

```
conda env create -f DSSPenv.yml
conda activate DSSPenv
```

- If you want to run DSSP (Dictionary of Secondary Structure of Proteins) locally, you need to download and install it on your computer!
for `conda`:
`conda install -c salilab dssp`

## Usage
This script processes a PDB file to extract main chain atoms, detect hydrogen bonds, assign secondary structures, and compare the results with DSSP. Follow the steps below to use the script:

**Install Dependencies**: Make sure the necessary Python packages are installed. Use either conda with the DSSPenv.yml file
for `conda`:
`conda env create -f DSSPenv.yml`

**Prepare the PDB File**: Ensure you have the PDB file you want to analyze. Place it in the appropriate directory.

**Run the Script**: To run the script, provide the path to the PDB file as an argument. You can modify the pdb_filename dynamically.

Here’s an example:
```
python contribute_secondary_structure.py -i ./data/4QR3.pdb -o ./result/4QR3_result.txt

```

**Output**:<br>
The script will print the following information in a file:

- Predicted Secondary Structure by Chain: <br>
- Overall Comparison Results: <br>
- Accuracy<br>
- Helix sensitivity<br>
- Beta-sheet sensitivity<br>


