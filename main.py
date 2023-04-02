import argparse 
parser = argparse.ArgumentParser(
        prog="FMO-GAMESS input file maker for small molecules",
        description="Generates an FMO (Fragment Molecular Orbital)  input file for GAMESS from a smiles string.\n Will print .inp files, .sdf files (conformers), .png for atom indexing, .svg for fragments and energy file with conformer energies."
)
parser.add_argument(
        '-s',
        '--smiles',
        required=True,
        type=str,
        help="Smiles string for small molecule"
)
parser.add_argument(
        '-n',
        '--nconf',
        required=False,
        type=int,
        default=1,
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
args = parser.parse_args()

from gamess_input import write_gamess_input_file

write_gamess_input_file(args.smiles,args.nconf,args.scftyp,args.basis,args.outprefix)