from Bio.PDB import PDBParser, DSSP

class StructureComparator:
    def __init__(self, pdbfile, predicted_structure):
        self.pdbfile = pdbfile
        self.predicted_structure = predicted_structure
        self.true_structure = self.get_true_structure_dssp()

    def get_true_structure_dssp(self):
        """Retrieve true structure from DSSP for each chain independently."""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", self.pdbfile)
        true_structure_per_chain = {}

        try:
            dssp = DSSP(structure[0], self.pdbfile)
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

    def compare_structures(self):
        """Compare predicted and true structures per chain and calculate accuracy and sensitivity."""

        # Used to accumulate overall structure data
        overall_predicted = []
        overall_true = []
        
        # Calculate accuracy and sensitivity per chain
        for chain_id, pred in self.predicted_structure.items():
            true = self.true_structure.get(chain_id, [])

            if not isinstance(pred, list) or not isinstance(true, list):
                print(f"Data for chain {chain_id} is not in the correct format (expected lists).")
                continue

            if len(pred) != len(true):
                print(f"Length mismatch in chain {chain_id}: predicted={len(pred)}, true={len(true)}")
                continue

            overall_predicted.extend(pred)
            overall_true.extend(true)

            tp = sum(1 for i in range(len(pred)) if pred[i] == true[i])
            accuracy = tp / len(pred) if len(pred) > 0 else 0

            tp_helix = sum(1 for i in range(len(pred)) if pred[i] == 'H' and true[i] == 'H')
            fn_helix = sum(1 for i in range(len(pred)) if pred[i] != 'H' and true[i] == 'H')
            helix_sensitivity = tp_helix / (tp_helix + fn_helix) if (tp_helix + fn_helix) > 0 else 0

            tp_beta = sum(1 for i in range(len(pred)) if pred[i] == 'E' and true[i] == 'E')
            fn_beta = sum(1 for i in range(len(pred)) if pred[i] != 'E' and true[i] == 'E')
            beta_sensitivity = tp_beta / (tp_beta + fn_beta) if (tp_beta + fn_beta) > 0 else 0

            print(f"Chain {chain_id} - Accuracy: {accuracy:.3f}, Helix Sensitivity: {helix_sensitivity:.3f}, Beta-Sheet Sensitivity: {beta_sensitivity:.3f}")

        # Overall metrics
        total_tp = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] == overall_true[i])
        total_accuracy = total_tp / len(overall_predicted) if len(overall_predicted) > 0 else 0

        total_tp_helix = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] == 'H' and overall_true[i] == 'H')
        total_fn_helix = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] != 'H' and overall_true[i] == 'H')
        total_helix_sensitivity = total_tp_helix / (total_tp_helix + total_fn_helix) if (total_tp_helix + total_fn_helix) > 0 else 0

        total_tp_beta = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] == 'E' and overall_true[i] == 'E')
        total_fn_beta = sum(1 for i in range(len(overall_predicted)) if overall_predicted[i] != 'E' and overall_true[i] == 'E')
        total_beta_sensitivity = total_tp_beta / (total_tp_beta + total_fn_beta) if (total_tp_beta + total_fn_beta) > 0 else 0

        print(f"\nOverall Protein - Accuracy: {total_accuracy:.3f}, Helix Sensitivity: {total_helix_sensitivity:.3f}, Beta-Sheet Sensitivity: {total_beta_sensitivity:.3f}")
