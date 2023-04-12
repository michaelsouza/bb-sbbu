import os
import pandas as pd
from tqdm import tqdm
from prody import parsePDB
from prody.proteins.localpdb import fetchPDB, pathPDBFolder


def extract_backbones(backbones_folder: str, pdb_id: str):
    # Parse the PDB file
    structure = parsePDB(pdb_id)
    
    # Get the list of chains
    chains = structure.getHierView().iterChains()
    
    # Extract the backbone atoms for each chain and save them in JSON format
    for chain in chains:
        chain_id = chain.getChid()
        backbone_file = os.path.join(backbones_folder, f'{pdb_id}_chain_{chain_id}.csv')

        # skip if file already exists
        if os.path.exists(backbone_file):
            print(f'File already exists for {backbone_file}')
            continue

        # Extract only the N-Ca-C atoms from the backbone
        backbone = chain.select('backbone and name N CA C')
        
        # skip if backbone is empty
        if backbone is None:
            print(f'Backbone is empty for {backbone_file}')
            continue
        
        # Prepare data for JSON
        backbone_data = []
        for atom in backbone:
            atom_data = {
                'chain_idx': atom.getChindex(),                
                'chain_name': atom.getChid(),
                'residue_idx': int(atom.getResnum()), # parse int64 to int
                'residue_name': atom.getResname(),
                'atom_idx': atom.getIndex(),
                'atom_name': atom.getName(),
                # split coordinates into x, y, z
                'x': atom.getCoords()[0],
                'y': atom.getCoords()[1],
                'z': atom.getCoords()[2],
            }
            backbone_data.append(atom_data)
        # save the backbone data in csv format
        df = pd.DataFrame(backbone_data)
        df.to_csv(backbone_file, index=False)


def create_nmr(nmr_folder, backbone_file: str, dmax: int):
    print('Processing file: ', backbone_file)
    
    # read the backbone file
    atoms = pd.read_csv(backbone_file)
    print('   number of atoms: ', len(atoms))

    # sort by residue number and N-Ca-C order
    atoms['atom_ord'] = atoms['atom_name'].map({'N': 0, 'CA': 1, 'C': 2})
    atoms = atoms.sort_values(['residue_idx', 'atom_ord'])
    atoms = atoms.drop(columns=['atom_ord'])

    # reset the index
    atoms.reset_index(drop=True, inplace=True)

    E = [] # list of edges
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            # compute the distance between the two atoms
            atom_i = atoms.loc[i] # get the i-th atom
            atom_j = atoms.loc[j] # get the j-th atom
            # compute the distance
            dist = ((atom_i['x'] - atom_j['x'])**2 + 
                    (atom_i['y'] - atom_j['y'])**2 + 
                    (atom_i['z'] - atom_j['z'])**2)**0.5
            # add the edge if the distance is less than dmax or the difference j - i is less than 3
            if (dist <= dmax) or (j - i <= 3):
                # add the edge to the list
                edge = {
                    'i': i,
                    'j': j,
                    'dist': dist,
                    'name_i': atom_i['atom_name'],
                    'name_j': atom_j['atom_name'],
                    'residue_i': atom_i['residue_name'],
                    'residue_j': atom_j['residue_name'],
                }
                E.append(edge)

    # create the dataframe with the edges
    nmr_file = os.path.basename(backbone_file).split('.')[0]
    nmr_file = os.path.join(nmr_folder, f'{nmr_file}_dmax_{dmax}.nmr')
    with open(nmr_file, 'w') as fd:
        fmt = '%4d %4d %.16e %.16e %4s %4s %5s %5s\n'
        for edge in E:
            edge = (
                edge["i"] + 1,
                edge["j"] + 1,
                edge["dist"],
                edge["dist"],
                edge["name_i"],
                edge["name_j"],
                edge["residue_i"],
                edge["residue_j"]
            )
            edge = fmt % edge
            fd.write(edge)


if __name__ == '__main__':
    # list of PDB IDs to download
    PDB = [
        '1n6t', '1fw5', '1adx', '1bdo', '1all', '6s61', '1fhl', '4wua',
        '6czf', '5ijn', '6rn2', '1cza', '6bco', '1epw', '5np0', '5nug',
        '4rh7', '3vkh'
    ]
    
    # create the pdb folder if it does not exist
    pdb_folder = os.path.join('data','pdb')
    if not os.path.exists(pdb_folder):
        os.makedirs(pdb_folder)    

    # create the backbones folder if it does not exist
    backbones_folder = os.path.join('data', 'backbones')
    if not os.path.exists(backbones_folder):
        os.makedirs(backbones_folder)

    # create the instances folder if it does not exist
    nmr_folder = os.path.join('data', 'nmr')
    if not os.path.exists(nmr_folder):
        os.makedirs(nmr_folder)

    # set the path to save the downloaded files
    pathPDBFolder(pdb_folder, divided=False)

    # download the PDB files
    print('Downloading PDB files')
    for pdb_id in tqdm(PDB):
        fetchPDB(pdb_id)  # download and save the PDB file
    
    # extract the backbones
    print('Extracting backbones')
    for pdb_id in tqdm(PDB):
        extract_backbones(backbones_folder, pdb_id)

    # create the NMR instances
    print('Creating NMR instances')
    DMAX = [5, 6] # list of dmax values
    for dmax in DMAX:
        print(f'Creating NMR instances with dmax = {dmax}')
        # list of backbone files
        backbone_files = [os.path.join(backbones_folder, fn) for fn in os.listdir(backbones_folder)]
        # keep only the backbone files with the A chain
        backbone_files = [fn for fn in backbone_files if fn.endswith('_A.csv')]
        # sort the backbone files by size
        backbone_files = sorted(backbone_files, key=os.path.getsize)
        # create the NMR instances
        for backbone_file in tqdm(backbone_files):
            create_nmr(nmr_folder, backbone_file, dmax)
