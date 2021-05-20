def find_b_turn(structure_name, file_name, chain, error):
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.vectors import calc_angle
    from Bio.PDB.DSSP import DSSP
    import pandas as pd
    import numpy as np

    parser = PDBParser()

    structure = parser.get_structure(structure_name, file_name)

    vectors = []
    phi = []
    psi = []
    exp_phi = np.array([60, -80])
    exp_psi = np.array([-120, 0])
    b_turns = []

    for atom in structure[0][chain].get_atoms():
        vectors.append(atom.get_vector())

    for i in range(len(vectors) - 2):
        if i % 2 == 0:
            phi.append(calc_angle(vectors[i], vectors[i + 1], vectors[i + 2]))
        else:
            psi.append(calc_angle(vectors[i], vectors[i + 1], vectors[i + 2]))

    df = pd.DataFrame(list(zip(phi, psi)), columns=['phi', 'psi'])

    for df_slice in df.rolling(window=2):
        if np.allclose(df_slice['phi'].values, exp_phi, atol=error) and np.allclose(df_slice['psi'].values, exp_psi, atol=error):
            b_turns.append(df_slice)

    return b_turns
