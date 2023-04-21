from rdkit.Chem import BRICS
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import Draw
from rdkit import Chem
from copy import deepcopy


    
def draw_fragments_to_grid(fragments,name):
    frags_2d = []
    for frag in fragments:
        copy = deepcopy(frag)
        AllChem.Compute2DCoords(copy)
        frags_2d.append(copy)
    with open(f"{name}_fragments.svg",'w') as svgfile:
        svgfile.write(Draw._MolsToGridSVG(frags_2d))

def get_BRICS_fragments(molecule):
    return Chem.GetMolFrags(BRICS.BreakBRICSBonds(molecule),asMols=True),list(BRICS.FindBRICSBonds(molecule))

def get_BRICS_pairs(molecule):
    bonds = list(BRICS.FindBRICSBonds(molecule))
    atom_pairs = []
    for i,bond in enumerate(bonds):
        atom_pairs.append(bond[0])
    return atom_pairs

def display_sp3_carbons(molecule,name):
    for i,atom in enumerate(molecule.GetAtoms()):
        molecule.GetAtomWithIdx(atom.GetIdx()).SetProp('atomNote',f'{atom.GetIdx()+1}')
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 4 and atom.GetHybridization() == Chem.HybridizationType.SP3:
            molecule.GetAtomWithIdx(i).SetProp('atomNote', f'{atom.GetIdx()+1},sp3')
    
    copymol = deepcopy(molecule)
    AllChem.Compute2DCoords(copymol)
    Draw.MolToFile(copymol,f"{name}_sp3_shown_molecule.png",size=(600,600),dpi=300)
    return molecule


def fragment_selected_bonds(molecule,bda,baa):
    bondlist = []
    for i,bond in enumerate(molecule.GetBonds()):
        if (bond.GetBeginAtomIdx()+1 in bda and bond.GetEndAtomIdx()+1 in baa) or (bond.GetBeginAtomIdx()+1 in baa and bond.GetEndAtomIdx()+1 in bda):
            bondlist.append(bond.GetIdx())
    fragments = Chem.GetMolFrags(Chem.FragmentOnBonds(molecule,bondlist),asMols=True)
    return fragments
    

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
