from rdkit_utilties import *
from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np
from rdkit.Chem import BRICS

def get_basis_text(basis):
    basis_text = " $BASIS "
    basis_dict = {
        "6-31G": "GBASIS=N31 NGAUSS=6",
        "6-31G*": "GBASIS=N31 NGAUSS=6 NDFUNC=1",
        "6-31G**": "GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=1",
        "cc-pVDZ": "GBASIS=CCD"
    }
    try:
        basis_text += basis_dict[basis]
        basis_text += " $END\n"
        return basis_text
    except KeyError:
        print(f"Unknown basis {basis}, please use one of the following: ")
        print(list(basis_dict.keys()))
        exit()
        
def get_contrl_text(basis, theory):
    contrl_text = " $CONTRL "
    
    contrl_dict = {
        ("rhf","6-31G"):"SCFTYP=RHF RUNTYP=ENERGY",
        ("rhf","6-31G*"): "SCFTYP=RHF RUNTYP=ENERGY",
        ("rhf","6-31G**"): "SCFTYP=RHF RUNTYP=ENERGY",
        ("rhf","cc-pVDZ"): "SCFTYP=RHF RUNTYP=ENERGY ISPHER=1",
        ("mp2","6-31G"):"SCFTYP=RHF RUNTYP=ENERGY MPLEVL=2",
        ("mp2","6-31G*"):"SCFTYP=RHF RUNTYP=ENERGY MPLEVL=2",
        ("mp2","6-31G**"):"SCFTYP=RHF RUNTYP=ENERGY MPLEVL=2",
        ("mp2","cc-pVDZ"):"SCFTYP=RHF RUNTYP=ENERGY MPLEVL=2 ISPHER=1"
    }
    try:
        contrl_text += contrl_dict[(theory,basis)]
        contrl_text += " $END\n"
        return contrl_text
    except KeyError:
        print(f"Unknown theory/basis combination {theory}/{basis}! Please use one of the following:")
        print(list(contrl_dict.keys()))
        exit()

def get_data_text(molecule,name):
    data_text = " $DATA\n"
    data_text += f"{name}\n"
    # FMO doesn't support symmetry 
    data_text += "C1\n" 
    # find unique atoms and their atom numbers in a molecule
    atoms = molecule.GetAtoms()
    atomdict = {}
    for i in atoms:
        atomdict[i.GetSymbol()] = i.GetAtomicNum()
    for i in atomdict.keys():
        data_text += f"{i} {atomdict[i]}\n"
    data_text += " $END\n"
    return data_text

def get_fmoxyz_text(molecule,confID):
    conf = molecule.GetConformer(confID)
    atoms = molecule.GetAtoms()
    fmoxyz_text = " $FMOXYZ\n"
    for i,atom in enumerate(atoms):
        fmoxyz_text += f"{atom.GetSymbol()} {atom.GetAtomicNum():>10.1f} {conf.GetAtomPosition(i).x:>10.5f} {conf.GetAtomPosition(i).y:>10.5f} {conf.GetAtomPosition(i).z:>10.5f}\n"
    fmoxyz_text += " $END\n"
    return fmoxyz_text

def get_fmo_text(molecule,basis,fragmentation="BRICS"):
    for atom in molecule.GetAtoms():
    # For each atom, set the property "atomNote" to a index+1 of the atom
        atom.SetProp("atomNote", str(atom.GetIdx()))
    if fragmentation=='BRICS':
        atom_pairs = get_BRICS_pairs(molecule)

    fragmented = BRICS.BreakBRICSBonds(molecule)
    fragments = Chem.GetMolFrags(fragmented,asMols=True)
    nfrag = len(fragments)
    atomlist = [1 for i in molecule.GetAtoms()]
    for i,frag in enumerate(fragments):
        for atom in frag.GetAtoms():
            if atom.GetSymbol()!='*':
                atomlist[int(atom.GetProp("atomNote"))] = i+1
    atomlist = np.array(atomlist,dtype='str')

    fmo_text = " $FMO\n"
    fmo_text += f"NFRAG={nfrag}\n"
    fmo_text += f"INDAT(1)={','.join(atomlist)}\n"
    fmo_text += " $END\n"
    # Add initial and final bond terms 
    fmo_text += " $FMOBND\n"
    bonds = list(BRICS.FindBRICSBonds(molecule))
    bea = []
    bed = []
    for i,bondgrp in enumerate(bonds):
        bndatms = list(bondgrp[0])
        i = bndatms[0]+1
        j = bndatms[1]+1
        if i in bed:
            bed.append(j)
            bea.append(i)
        else:
            bed.append(i)
            bea.append(j)
    for i in range(len(bed)):
        fmo_text += f"{-bed[i]:>3} {bea[i]:>3} {basis}\n"
    fmo_text += " $END\n"

    return fmo_text 
     
def get_fmohyb_text(basis):
    fmohyb_text = " $FMOHYB\n"
    fmohyb_dict = {
        "mini":"MINI 5 5\n \
  1 0  -0.109772    0.515046    0.000000    0.000000    0.864512\n \
  0 1  -0.109775    0.515062    0.815062    0.000000   -0.288166\n \
  0 1  -0.109775    0.515062   -0.407531   -0.705864   -0.288166\n \
  0 1  -0.109775    0.515062   -0.407531    0.705864   -0.288166\n \
  0 1   0.996474    0.015610    0.000000    0.000000    0.000000\n",
        "hf-3c":"HF-3C 5 5\n \
  1 0  -0.066232  0.310761  0.000000  0.000000  0.521599\n \
  0 1  -0.066232  0.310761 -0.245884 -0.425884 -0.173866\n \
  0 1  -0.066232  0.310761 -0.245884  0.425884 -0.173866\n \
  0 1  -0.066232  0.310761  0.491769  0.000000 -0.173866\n \
  0 1   1.001501  0.015689  0.000000  0.000000  0.000000\n",
        "6-31G":"6-31G 9 5\n \
 1 0  -0.067724    0.300281    0.000000    0.000000    0.606750\n \
       0.306535    0.000000    0.000000    0.309793\n \
 0 1  -0.067730    0.300310    0.572037    0.000000   -0.20223\n \
       0.306552    0.292061    0.000000   -0.103255\n \
 0 1  -0.067730    0.300310   -0.286019   -0.495398   -0.202234\n \
       0.306552   -0.146031   -0.252933   -0.103255\n \
 0 1  -0.067730    0.300310   -0.286019    0.495398   -0.202234\n \
       0.306552   -0.146031    0.252933   -0.103255\n \
 0 1   1.011954   -0.016447    0.000000    0.000000    0.000000\n \
      -0.059374    0.000000    0.000000   -0.000001\n",
        "6-31G*":"6-31G* 15 5\n \
 1 0  -0.065034    0.288264    0.000000    0.000000    0.604413\n \
       0.290129    0.000000    0.000000    0.319045   -0.017106\n \
      -0.017106    0.057935    0.000000    0.000000    0.000000\n \
 0 1  -0.065041    0.288294    0.569833    0.000000   -0.201457\n \
       0.290147    0.300784    0.000000   -0.106342    0.049599\n \
      -0.017106   -0.008771    0.000000   -0.027223    0.000000\n \
 0 1  -0.065040    0.288293   -0.284917   -0.493490   -0.201456\n \
       0.290146   -0.150393   -0.260487   -0.106341   -0.000428\n \
       0.032923   -0.008771    0.033353    0.013612    0.023577\n \
 0 1  -0.065040    0.288293   -0.284917    0.493490   -0.201456\n \
       0.290146   -0.150393    0.260487   -0.106341   -0.000428\n \
       0.032923   -0.008771   -0.033353    0.013612   -0.023577\n \
 0 1   1.010938   -0.011976    0.000000    0.000000    0.000000\n \
      -0.054085    0.000000    0.000000   -0.000001   -0.003175\n \
      -0.003175   -0.003175    0.000000    0.000000    0.000000\n",
      "6-31G**":"6-31G** 15 5\n \
 1 0  -0.068254    0.305270    0.000003    0.000000    0.619132\n \
       0.287030    0.000002    0.000000    0.307201   -0.022701\n \
      -0.022701    0.042170    0.000000    0.000000    0.000000\n \
 0 1  -0.068257    0.305303    0.583705    0.000000   -0.206360\n \
       0.287057    0.289613    0.000000   -0.102393    0.034962\n \
      -0.022700   -0.015495    0.000000   -0.023534    0.000000\n \
 0 1  -0.068257    0.305306   -0.291851   -0.505502   -0.206358\n \
       0.287061   -0.144805   -0.250811   -0.102392   -0.008284\n \
       0.020546   -0.015495    0.028830    0.011767    0.020381\n \
 0 1  -0.068257    0.305306   -0.291851    0.505502   -0.206358\n \
       0.287061   -0.144805    0.250811   -0.102392   -0.008284\n \
       0.020546   -0.015495   -0.028830    0.011767   -0.020381\n \
 0 1   1.010732   -0.013164    0.000000    0.000000    0.000001\n \
      -0.052063    0.000000    0.000000    0.000000   -0.001621\n \
      -0.001621   -0.001620    0.000000    0.000000    0.000000\n",
      "cc-pVDZ":"cc-pVDZ 15\n \
 1 0  0.073001  0.327576  0.211428  0.000000  0.000000  0.659685  0.000000\n \
      0.000000  0.308077 -0.021065 -0.021065  0.042131  0.000000  0.000000\n \
      0.000000\n \
 0 1  0.073006  0.327600  0.211444  0.621943  0.000000 -0.219885  0.290451\n \
      0.000000 -0.102688  0.035108 -0.021065 -0.014043  0.000000 -0.022932\n \
      0.000000\n \
 0 1  0.073004  0.327592  0.211439 -0.310978 -0.538623 -0.219882 -0.145229\n \
     -0.251540 -0.102686 -0.007022  0.021065 -0.014043  0.028086  0.011466\n \
      0.019861\n \
 0 1  0.073004  0.327592  0.211439 -0.310978  0.538623 -0.219882 -0.145229\n \
      0.251540 -0.102686 -0.007022  0.021065 -0.014043 -0.028086  0.011466\n \
     -0.019861\n \
 0 1  0.996761 -0.035495 -0.034218  0.000000  0.000000  0.000000  0.000000\n \
      0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000\n \
      0.000000\n"
    }
    fmohyb_text += fmohyb_dict[basis]
    fmohyb_text += " $END\n"
    return fmohyb_text

def write_gamess_input_file(smiles,nconf,scftyp,basis,name):
    # write gamess FMO input files for each conformer 
    # only the FMOXYZ section would change for each conformation 
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    for atom in mol.GetAtoms():
        atom.SetProp("atomNote", str(atom.GetIdx()+1))
    Draw.MolToFile(mol,f"{name}_2Dimage.png",size=(300,300),dpi=300)
    fragmented=BRICS.BreakBRICSBonds(mol)
    pieces=Chem.GetMolFrags(fragmented,asMols=True)
    with open(f'{name}_fragments.svg','w') as svgfile:
        svgfile.write(Draw._MolsToGridSVG(pieces))
    basis_text = get_basis_text(basis)
    contrl_text = get_contrl_text(basis, scftyp)
    data_text = get_data_text(mol,name)
    fmohyb_text = get_fmohyb_text(basis)
    fmo_text = get_fmo_text(mol,basis,"BRICS")
    mol,sorted_res = gen_unique_conformers(mol,nconf)
    energyfile = open('conf_energies.txt','w')
    writer = Chem.SDWriter(f'{name}_conformers.sdf')
    for confID,e in sorted_res:
        mol.SetProp('ID',f"{name}_conformer_{confID}")
        writer.write(mol,confId=confID)
        energyfile.write('%s,%s\n'%(confID,e))
        fmoxyz_text = get_fmoxyz_text(mol,confID)
        filename = f"{name}_{confID}.inp"
        with open(filename,'w') as outfile:
            outfile.write(basis_text)
            outfile.write(contrl_text)
            outfile.write(data_text)
            outfile.write(fmoxyz_text)
            outfile.write(fmo_text)
            outfile.write(fmohyb_text)
    energyfile.close()