# FMO4smallmolecules
### FMO input file maker for small molecules
This set of scripts can be used to create input files for running Fragment Molecular Orbital calculations using GAMESS on small molecules. There are some closed sources packages like Facio which work very well on proteins but I don't find it very user friendly for small molecule calculations and some open source packages which aren't very well maintained hence why I created this particular package. It is still far from user friendly but it works very well for my use case where I can just throw a SMILES string or an sdf file from DiffDock or Vina and get it fragmented either using BRICS rules or by specifying manual fragmentation points. 


Usage when you only have smiles:
```
python main.py --smiles "C[C@H](c1nc2cnc3c(c2n1C4CCC(CC4)CC#N)cc[nH]3)O" --nconf 10 --basis "6-31G**" --scftyp rhf --outprefix 1J6
```
Usage when you only have an SDF file and want to try manual fragmentation you need to provide a list of bond attached and bond detached atoms:
```
python main.py -sdf conformer.sdf -bda 11 -baa 10 -fragstyle manual -o test_compound
```
Usage when you only want to draw the molecule with sp3 carbons and atoms highlighted
```
python main.py -sdf conformer.sdf -showmol True
```

if the -showmol True argument is passed, regardless of any other arguments being passed, the code will draw the molecule and stop.

The output consists of:
1. GAMESS .inp files for each conformer, 
2. .png file showing the atom numbering, 
3. .svg file showing the different fragments generated using BRICS for now, 
4. .sdf file with conformers generated,
5. conf_energies.txt with the energies of each generated conformer using the MMFF94s forcefield.


This code is under heavy development and should not be used for production calculations. Currently, following additions are being worked upon,
1. Adding hybrid orbitals for bond detached atoms other than SP3 hybridized Carbons
2. Better representation of bond detached and attached atoms on fragment images 
3. Adding SMIRNOFF forcefield support for conformer generation
4. Testing the code on macrocyclic ligands where I believe the FMO technique might be more useful. 
