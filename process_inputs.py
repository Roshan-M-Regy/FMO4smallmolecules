import argparse 
from rdkit import Chem
def parse_inputs():
    parser = argparse.ArgumentParser(
            prog="FMO-GAMESS input file maker for small molecules",
            description="Generates an FMO (Fragment Molecular Orbital)  input file for GAMESS from a smiles string.\n Will print .inp files, .sdf files (conformers), .png for atom indexing, .svg for fragments and energy file with conformer energies."
    )
    parser.add_argument(
            '-s',
            '--smiles',
            required=False,
            type=str,
            help="Smiles string for small molecule"
    )
    parser.add_argument(
            '-sdf',
            '--sdffile',
            required=False,
            type=str,
            help="SDF file containing single conformer."
    )
    parser.add_argument(
            '-n',
            '--nconf',
            required=False,
            type=int,
            default=0,
            help="Number of unique conformations to generate. Actual number might be different than this depending on uniqueness."
    )
    parser.add_argument(
            '-st',
            '--scftyp',
            required=False,
            type=str,
            default='rhf',
            help="Type of SCF method"
    )
    parser.add_argument(
            '-b',
            '--basis',
            required=False,
            type=str,
            default='6-31G**',
            help="Basis set to use"
    )
    parser.add_argument(
            '-o',
            '--outprefix',
            required=False,
            type=str,
            default='AgentX',
            help="Generates outfile_confID.inp"
    )
    parser.add_argument(
            '-bda',
            '--bonddetached',
            required=False,
            nargs="+",
            type=int,
            default=[],
            help="Bond detached atoms for manual fragmentation"
    )
    parser.add_argument(
            "-baa",
            "--bondattached",
            required=False,
            nargs="+",
            type=int,
            default=[],
            help="Bond attached atoms for manual fragmentation"
    )
    parser.add_argument(
            "-fragstyle",
            "--fragmentation_style",
            required=False,
            type=str,
            default="BRICS",
            help="Fragmentation style for the molecule, BRICS or MANUAL.\n MANUAL requires -baa and -bda inputs to be provided."
    )
    parser.add_argument(
            "-showmol",
            "--show_molecule",
            required=False,
            type=bool,
            default=False,
            help="Save the 2D image of the molecule in a png file showing SP3 carbons and atom IDs.")

    args = parser.parse_args()

    if args.smiles==None:
        if args.sdffile==None:
            print ("Either give a SMILES string or an SDF file!")
            exit()
        else:
            molecule = Chem.AddHs(Chem.SDMolSupplier(args.sdffile)[0], addCoords=True)
    else:
        if args.nconf==0:
            args.nconf = 1
        molecule = Chem.AddHs(Chem.MolFromSmiles(args.smiles))

    if molecule is None:
        print("Can't make molecule from given inputs!")
        exit()
    return molecule,args
