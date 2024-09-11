# projet

## Description
This project is a Python-based tool for predicting and evaluating the secondary structure of proteins using hydrogen bond detection and energy calculations. The tool compares the predicted secondary structure with the DSSP-assigned secondary structure to calculate accuracy, helix sensitivity, and beta-sheet sensitivity.

## Features
- **Main chain Atom Extraction.** : The file `read_pdb_structure.py` extracts the main chain atoms (N, CA, C, O) from PDB file.
- **Hydrogen Bond Detection** : The file `hydrogen_bond_detector.py` identifies potential hydrogen bonds between main chain atoms based on energy thresholds and distances between donors and accpetors.
- **Secondary Structure Assignment** :Predicts secondary structures, helix and beta sheets, based on hydrogen bond patterns.
- **DSSP Comparison** : Compares the predicted secondary structures with those provided by DSSP and calculates accuracy, helix sensitivity and beta-sheet sensitivity.

## Requirement
To run these scripts, the file `environment.yml` can create identical development environments and ensure that all dependencies are aligned with project requirements.
