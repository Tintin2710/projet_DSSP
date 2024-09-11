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

            # Helix determination condition
            if diff in [3,4,5]:
                helix.add(donor_residue_id)
                helix.add(acceptor_residue_id)

            # β sheet determination condition
            elif 5 < diff <= 200 and bond['energy'] < -0.4:
                beta.add(donor_residue_id)
                beta.add(acceptor_residue_id)

        # Initialise all residues to ‘C’ (coli)
        secondary_structure = ['C'] * len(self.main_chain_atoms)

        # Assign helix (‘H’) and β-fold (‘E’), find the corresponding indexes by res_id
        for res_id in helix:
            if res_id in res_id_to_index:
                index = res_id_to_index[res_id]
                secondary_structure[index] = 'H'  

        for res_id in beta:
            if res_id in res_id_to_index:
                index = res_id_to_index[res_id]
                if secondary_structure[index] != 'H':
                    secondary_structure[index] = 'E'

        # Debugging information: print the number of atoms in the main chain and the allocation result
        print("Length of main_chain_atoms:", len(self.main_chain_atoms))
        print("Helix indices:", helix)
        print("Beta indices:", beta)

        return secondary_structure



