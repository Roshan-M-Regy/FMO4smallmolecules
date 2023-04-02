from rdkit.Chem import BRICS
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign


def get_BRICS_pairs(molecule):
    bonds = list(BRICS.FindBRICSBonds(molecule))
    atom_pairs = []
    for i,bond in enumerate(bonds):
        atom_pairs.append(bond[0])
    return atom_pairs

def gen_unique_conformers(mol,nconf):
    params = AllChem.ETKDGv3()
    params.numThreads = -1  # use all threads
    params.useSmallRingTorsions = True  # Use small ring torsions for speed
    params.pruneRmsThresh = 0.5  # prune conformers based on RMSD, 1 Angstrom
    params.useExpTorsionAnglePrefs = True  # Use experimental torsion preferences
    params.useBasicKnowledge = True  # Use basic knowledge from database
    params.enforceChirality = True  # Enforce chiral constraints
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=nconf, params=params)
    mp = AllChem.MMFFGetMoleculeProperties(mol,mmffVariant='MMFF94s')
    AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, mmffVariant='MMFF94s')
    res = []
    for cid in cids:
        print(cid)
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
        e = ff.CalcEnergy()
        res.append((cid,e))
    sorted_res = sorted(res, key=lambda x:x[1])
    rdMolAlign.AlignMolConformers(mol)
    return mol,sorted_res