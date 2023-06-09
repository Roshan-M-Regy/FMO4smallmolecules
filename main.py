from process_inputs import parse_inputs
molecule,args = parse_inputs()

if args.show_molecule:
    from rdkit_utilities import display_sp3_carbons
    display_sp3_carbons(molecule,args.outprefix)
else:
    from gamess_input import write_gamess_input_file
    write_gamess_input_file(molecule,args.nconf,args.scftyp,args.basis,args.outprefix,args.fragmentation_style,args.bonddetached,args.bondattached)
