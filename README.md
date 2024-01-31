# hydration

Amber-oriented scripts for adjusting water count on-the-fly in a molecular dynamics simulation

## Motivation

Getting the number of waters right in a constant-volume amber run is tricky. The pressure is
a traditional indicator, but cannot always be trusted if you are using restraints. Trouble is,
especially in the beginning of the run you probably want to have restraints to keep the molecule from 
straying too far from it starting point before everything settles down. Unfortunately, the 
number of waters is one of those things you need to decide at the very beginning, and then you 
are married to it.<br>
 The purpose of this project is to allow modification of the number of waters
in a simulation without having to start it all over again from scratch. Instead, you can edit
the waters on-the-fly.

## Description

 These scripts implement
a method for taking the constant N out of a constant NVT simulation. That is, you can safely add,
 remove or teleport waters in a simulation that is alredy flying. No need to run LEAP again, or to 
re-heat the system from absolute zero. Instead, you can keep your velocities and sprinkle in a few
more waters into any vacuum bubbles that have appeared. You may also strip out a few 
waters at a time, or even teleport a water from a place where there is
too much electron density to a place that needs more density, using a standard crystallographic
Fo-Fc difference map as a guide. <br>
The trick is to avoid any
"transporter accidents" where two molecules try to occupy the same space, and exploding the run.
The purpose of these scripts is to avoid such clashes and keep the simulation running smoothly.
Also, if you are using restraints on some of the waters (like I do), you will want to be sure that
these special waters dont get stripped out.

## Getting Started

### Dependencies

* These are linux c-shell scripts, so you will need a working `/bin/tcsh`
* `gemmi` to measure vacuum bubbles, which is distributed with the CCP4 suite
* the Phenix Suite to probe electron density map values from your mtz file
* `AMBER` and its tools installed, particularly `cpptraj` and `parmed`
* an `*.rst7` AMBER restart file for your system, and the accompanying `*.parm7` topology file
* an `orignames.pdb` file that contains all the original non-xyz information from the PDB file you initially gave to LEAP.
* an optional `current_restraints.pdb` file, containing the names and xyz positions of your restraint reference points
* You will also need to have prepared a `padded.parm7` topology file, which you can produce by adding ~50k dummy waters at the end of the PDB file you provide to LEAP. You only need to do this once. 

### Installing

* git clone this repo
* copy all the files into somewhere in your shell `$PATH`, and make them executable:
```
    chmod u+x *.com *.awk
```
Yes, I know the extension says `*.com`, but these are not Windows executables. The use of `.com` to denote shell scripts pre-dates Windows.

### Executing program

#### file format conversion 
The first program you may want to run is:
```
rst2pdb_runme.com amber.rst7 orignames=orignames.pdb parmfile=xtal.prmtop
```
This will serve as a good test to see if everything is working. Success is indicated by the output of a new pdb formatted file that will be called `amber.pdb`, which will have the same names as the atoms in `orignames.pdb`, but with the
xyz coordinates contained in the amber.rst7 file. The options shown for `orignames` and `parmfile` are the internally-stored defaults, so are not required if you already have those files with those names. I recommend using these names, but if you don't like those
names for some reason, I recommed editing the defaults at the top of the `rst2pdb_runme.com` script.<br>
Another internal default is `paddedparm=padded.parm7`, which is the topology file you will want to generate yourself
using LEAP, but with a few thousand extra dummy waters tacked on to the end. No need to use AddToBox, all the dummy waters can be at XYZ=0,0,0. LEAP won't care, and you are going to strip them out anyway. Remember, `*.parm7` files do not contain coordinates. If for some reason the `xtal.prmtop` is missing or otherwise incompatible with `amber.rst7`, the `rst2pdb_runme.com` script will look for `padded.parm7` and strip out waters until it can be used to read `amber.rst7`. The new, smaller, topology file will be written out
as `resized.parm7` as well as a shortened version of `orignames.pdb`. However, it is okay, and indeed recommended, for `orignames.pdb` to be longer than your active system. In fact, make it as big as `padded.parm7`.<br>

You can change the output file prefix using `outpreifx=whatever`, or by specifying a pdb file name on the command line. In general, for all these scripts, any variable you see at the beginning of the script can be changed on the command line using the `variable=value` syntax.

#### Abhor a Vacuum
The next thing you probably want to do is check to see if your current simulation has vacuum bubbles:
```
bubble_check_runme.com amber.pdb Vwater=96 minvoid=5
```
This will use `gemmi` to look for spaces between atoms and measure their volumes. In bulk water, there is about 96 A<sup>3</sup> of volume for each water molecule, which is the default for `Vwater`, so you don't have to specify it to take the default. The `minvoid` variable is the number of water molecules that a vacuum bubble can hold before it is considered a "void". How big are vacuum bubbles in normal bulk water at STP? I have no idea. But 5 strikes me as a bit much, so it is the default. If you have a vacuum bubble that can hold thousands of waters, then you might do well by filling it in.  This script will print out a reocmmended fill-in command, such as:<br>

#### Hydrate
Here is how you fill in those vacuum bubbles:
```
hydrarte_runme.com amber.pdb nadd=1000 outrst=wetter.rst7 outtop=xtal.prmtop
```
The output will be the result of running the AmberTools AddToBox program.  The last molecule in `amber.pdb` is extracted and taken as the "water" to add, thus preserving any extra points. In the end, you will have a new restart file, including velocities, and topology file to continue your simulation, but now with a few more waters. The filename values in the above command are the defaults. The key to adding in new waters to an AMBER parm file is the unmentioned `padded.parm7` file. The script will complain and exit if you don't have `padded.parm7`. It works by not adding extra waters to your previous `*.parm7` file, but rather by stripping an appropriate number of waters out of `padded.parm7` so that it is compatible with the new number of waters.<br> 
All new waters are given zero velocity, but will quickly heat up when the simulation restarts. I have not seen any need for an explicit re-heating protocol in my simulations so far.<br>
Also note: you can provide an `*.rst7` file on the command line instead of a PDB file and `hydrate_runme.com` will call `rst2pdb_runme.com` internally to generate the `*.pdb` file needed by AddToBox.<br>
If you provide a `restraints=current_restraints.pdb` file on the command line, this script will also update the reference points file called `ref.crd` so that it is also the same size as the new `wetter.rst7`, using the XYZ coordinates in `current_restraints.pdb`.

#### De-hydrate
Removing waters from an AMBER simulation is much easier than adding them. The existing and popular `cpptraj` and `parmed` programs can do this with the `strip` command. You can run those without my scripts.  However, a good question about rejecting waters is: which ones?
```
dehydrarte_amber_runme.com toowet.rst7 maxreject=10 outprefix=drier \
  mtzfile=phenixout.mtz mtzlabel=FOFCWT \
  notthese=current_restraints.pdb
```
An excellent way to prioritize your waters is with a density map, such as from a phenix or refmac refinement of a recent snapshot, or, even better, the average density of the simulation trajectory to date, subtracted from the original 2mFo-Fc map obtained from the original refmac or phenix refinement of the starting structure. With this map, you have an experimentaly-drived guide to which waters need to be there the least.  This script will probe the density map contained in the provided mtz file with all the water oxygen atoms from the `toowet.rst7` file.  The ones in the most-negative density will go first, with `maxreject` being the number of water molecules to discard.<br>
If you have restraints, or some other reason not to reject particular waters, then provide them, by name, in a PDB file using `notthese=`. The names should correspond to those in the underlying `orignames.pdb` file, where the ordinal position of each atom in `orignames.pdb` corresponds to the ordinal list of coordinates in `toowet.rst7`. The caveat here, however, is that not only will any atom named in `notthese=` be protected from rejection, so will every atom that comes before them. This is because a `strip` operation effectively re-names every atom that comes after it. If one of those low-on-the-list atoms is involved in a restraint, then you will get an explosion. This script avoids such disasters by first determining which waters are "disposable". But, if the last water in the `*.rst7` file is restrained, nothing can be done. Not until you run the `remap_waters_runme.com` script below.  In general, it may be a good idea to always remap before dehydration.


#### Teleportation
All this re-sizing of parm files is not neccesary if you want to add and subract the same number of waters. This is equivalent to teleportation. A water will vanish from one part of the run and appear somewhere else. This is a good way to leap over barriers in the simulation. It can often be the case that a buried water was not included at the start of the run, and never gets populated by random diffusion. The opposite can also be true: a water may have gotten jammed inside of the protein, unable to escape. An excellent way to detect these kinds of problems is to use a difference map. That is, subtract the average density from the simulation so far (Fcalc) from the observed density (Fobs) to obtain an Fobs-Fcalc map. Positive features in this map are usually waters that need to be populated, and negative features are waters that don't belong. The one big no-no, however, would be to teleport a water on top of another one. That would be a disaster. So, you want to check that the way is clear before engaging the transporter device. This script automates all that.
```
water_teleport_runme.com amber.rst7 destinations.pdb maxmoves=10 \
  outfile=teleported.rst7 \
  mtzfile=phenixout.mtz mtzlabel=FOFCWT \
  notthese=current_restraints.pdb
```
The value of `maxmoves` is the maximum number of water molecules to teleport in this run.<br>
The `destinations.pdb` file can be obtained from the original, refined and deposited structure, or perhaps a peak-pick of some density map. The Fobs-Fcalc map is one choice, but the Fobs map itself is another. Don't forget to symmetry-expand the list to cover your simulation supercell first. Think of `destinations.pdb` as a list of suggested destinations for teleporting waters. The script will internally score these sites by map density and screen them for clashes before moving anything.<br>
The `notthese` variable indicates a PDB file containing the names of existing atoms that should not be teleported. Those involved in restraints are a good example of things that are a bad idea to teleport across the unit cell. The resulting restraint energy will be gigantic, and crash the simulation. Best to leave those alone. So, this script does that.<br>
As a final check, this script will perform an energy minimization and display the total energy of the system before and after the teleport. Unless there is a bug, this step should be unexciting. You can disable these relatively time-consuming stages with `minimize=0` and `energycheck=0`<br>
Where do the teleported waters come from?  Like the dehydration script above, they are chosen based on the difference density. This actually creates a subtle problem, which is that if you restrain the waters after teleporting them, like I do, then you create restraints at the end of the list, and that prevents successful dehydration runs.
To fix all this, one must re-organize the water list before each dehydration run using the next script.


#### Reorganize
Which waters to reject? Turns out the pragmatic answer in AMBER is: the ones at the end of the file. Effectively, no matter which waters you pick to reject, everything after them gets re-named. You can try doing things like keeping track of waters by name instead of by slot number in the `*.rst7` file, but, trust me, there madness lies. It is better to just periodically re-sort the whole list.
```
remap_waters_runme.com amber.rst7 \
  mtzfile=phenixout.mtz mtzlabel=FOFCWT \
  restraints=current_restraints.pdb
```
This script works much like the `dehydrate_amber_runme.com` script above, except that it does not reject any waters. Rather, it prepares for a dehydration run by re-organizing the waters first. It probes the provided difference map to sort all the waters, and puts the most disposable at the bottom of the file. This is all using the `cpptraj` "remap" feature. The other thing this script does is take all the atoms named in the `restraints=` file and remaps them to the top of the list, regardless of the density they are in.  This way they are least likely to interfere with stripping the very worst waters out of the system. This requies re-naming all the restrained atoms, and that  list is provided in the output file `remapped_restraints.pdb`. Give this file to `dehydrate_amber_runme.com` as `notthese=remapped_restraints.pdb` along with the `remapped.rst7` AMBER restart file provided by this script.


#### Grafting a PDB into AMBER without LEAP
This is called by a few of the above scripts. In general, if you have made some changes to your system in PDB format (such as adding or subtracting waters), but you don't want to re-run LEAP. And, in fact, you want to keep the velocities from the previous run and just keep going. Then this is the script for you.
```
graft_atoms_runme.com edited.pdb previous.rst7 outprefix=grafted 
```
The protein atoms in `edited.pdb` must come in the same order as the ones in `previous.rst7`, but they can have different XYZ positions, and you can also have a differnet number of waters. If waters are missing from the end of `edited.pdb` then you will get a new system in `grafted.rst7` and `grafted.prmtop` that has the same number of waters. If you have extra waters in `edited.pdb`, and there is also a `padded.parm7` file available, then the output system in `grafted.rst7` and `grafted.prmtop` will have those new waters (albeit with zero velocity).  Using this new system, the `edited.pdb` file will effectively become the new `orignames.pdb`. However, as long as you don't mind the old names, the `rst2pdb_runme.com` will work very happily with an `orignames.pdb` that has more entries than it needs.
<br>
As a final check, this script will perform an energy minimization and display the total energy of the system before and after the teleport. Unless there is a bug, this step should be unexciting. You can disable these relatively time-consuming stages with `minimize=0` and `energycheck=0`<br>

#### converting and re-numbering pdb files
This is an awk jiffy program for performing various re-formatting functions on a PDB file. It is designed primarily to deal with common re-formatting issues between the default representations in AMBER vs those of crystallographic refinement programs like the Phenix and CCP4 suites.  With no command-line options, it tries to convert an AMBER formatted pdb file into something `refmac` or `phenix.refine` can read. For example, it converts `WAT` into `HOH`, and the N-terminal `H1` atom on each chain will be changed to just `H`, which `refmac` expects. It also changes the `HID/HIP/HIE` histidine residues into just plain `HIS`, and also changes `CYX` to `CYS` and `ASH` and `GLH` to `ASP` and `GLU`.  Run it like this:
```
% convert_pdb.awk -v output=amber refmacout.pdb > tleapme.pdb
```
This will do things like change the N-terminal `H` atom into the `H1` that tleap expects. It will also add `TER` records at the end of each protein, ligand and water molecule. It will further change `HOH` to `WAT`, and `NH4` to `AMM`, and `ACY` to `ACT`. Information about the protonation state can be provided with keywords prepended to the PDB file, such as:
```
echo "PROTON A324" | cat - refmacout.pdb | convert_pdb.awk -v output=amber > tleapme.pdb
```
You can also begin a line with `HID `, `HIE ` or `HIP ` to specify residue IDs (chain letter and residue number) to be given a non-default protonation state.<br>
You can also set other variables for different behaviors:<br>
* `only=protein` will print out only residues types known to be non-exotic amino acids
* `skip=water,H,EP` will skip all `HOH` residues, as well as leave out all hydrogen and extra-point atoms
* `fixEe=1` will justify 2-letter residue types, and element symbols as they are output by refmac.
* `fixEe=1` will justify 2-letter residue types, and element symbols as they are output by refmac.
* `append=ordresnum` will append the oridnal residue number at the end of each line
* `append=origid` will append the original, unmodified atom,type,chain,resnum string to the end of each line
* `renumber=ordinal` re-start residue numbering from the beginning
* `renumber=terify` add `TER` records
* `renumber=w8` use all eight spaces available for residue numbers
* `renumber=w4` try to make residue number fit into four digits, increment chain if neccesary
* `renumber=ordinal` start counting
* `renumber=watS` discard all chain information and re-chain, starting protein at chain A, start water at chain S.
* `renumber=chain` restart residue counter for each new chain encountered
* `renumber=dedupe` only check for and try to eliminate duplicate residue names
* `renumber=ordinal,watS,w4,chainrestart` re-start all residue numbers. Make protein start at chain A. Non-protein at chain L, water start at chain S. Use 4-character numbers and increment chain ID to prevent overflows.


#### compact residue numbers
This is an awk jiffy program for dealing with the problem of repesenting more than 9999 residues in a PDB format file. The output of programs like `cpptraj` curently just give you modulo 10000 values, which makes some downstream programs like `phenix.refine` or `refmac` unhappy with the duplicate residue names. One way to fix this is to use the extra space to the right of the 4-character residue number in the PDB file format, but few programs recognize this. Another, more widespread trick is called "hybrid 36" encoded values. This combines letters and numbers to form up to up to 2436111 residue numbers, and programs like gemmi, phenix and the CCP4 suite do work with these encoded values. So, I wrote a fast jiffy script for either sequentially re-numbering all residues as hy36, or for turning 8-character-wide residue numbers into hy36. Run it like this:
```
hy36_encode.awk -v ordinal=1 cpptraj_out.pdb > unique_resnums.pdb
```
This will re-number all residues starting at 1, ticking up the counter every time the 8-character string starting at column 23 changes from the previous line. The ordinal residue number is then packed into a 4-letter hy36 string in the usual residue number place. I also append the ordinal residue number in ordinary digits form at the end of each line.  If you do not set `ordinal=1` then the residue numbers are read by interpreting an 8-character wide string starting at column 23 as an integer.




#### rmsd
This is an awk jiffy program for computing the Root Mean Square Deviation between two PDB files. Unlike some other programs no attempt is made to align the models to minimize this difference. All this program does is look for pairs of atoms with the same name and compares the XYZ positions, and also any changes in B factor or occupancy. By "same name" I mean the 17-character string starting at column 12 on each line beginning with "ATOM" or "HETATM". The order does not matter, and it will complain if there is only one or more than two copies of a given atom.  Run it like this:
```
% rmsd original.pdb modified.pdb
2050 atom pairs found
RMSD(CA )= 0.00909198 (128 CA pairs)
RMSD(all)= 0.0117287 (2050 atom pairs)
RMSD(Bfac)= 0.0494906
MAXD(all)= 0.105669	for  HH12AARG A  18   
MAXD(Bfac)= -0.45	for   OE2AGLU A  32   
```
Because it is an awk program, it can also natively take the two PDB files as a stream on standard input. If you set the "xlog" variable the above summary will be listed on one line:
```
% rmsd -v xlog=1 original.pdb modified.pdb
0.00909198 0.01172872    0.00000   0.0495     0.1057    0.000   -0.450
```
And if you set the `debug` variable, each individual difference will be listed, allowing you to sort it.<br>
```
%egrep -vh "HOH|     H" original.pdb modified.pdb | rmsd -v debug=1 | sort -k1.25g | grep moved | tail
  CD AARG A  18    moved   0.0552 (XYZ)   0.00 (occ)   0.10 (B) at cen_x   -0.815   9.507 -13.204
  NE AARG A  18    moved   0.0556 (XYZ)   0.00 (occ)  -0.02 (B) at cen_x   -0.504   8.078 -13.274
  O  BLYS B  30    moved   0.0587 (XYZ)   0.00 (occ)   0.02 (B) at cen_x   -7.453  -1.368   4.529
  OE2AGLU A  32    moved   0.0591 (XYZ)   0.00 (occ)  -0.45 (B) at cen_x  -10.790   6.916   5.851
  O  ALYS A  50    moved   0.0628 (XYZ)   0.00 (occ)  -0.02 (B) at cen_x   -6.177   1.199   6.948
  NH2AARG A  18    moved   0.0665 (XYZ)   0.00 (occ)   0.04 (B) at cen_x   -0.888   5.829 -13.046
  O  BLYS B  50    moved   0.0685 (XYZ)   0.00 (occ)  -0.13 (B) at cen_x   -6.144   1.414   6.980
  CZ AARG A  18    moved   0.0695 (XYZ)   0.00 (occ)   0.02 (B) at cen_x   -1.303   7.084 -12.897
  OE1AGLU A  32    moved   0.0787 (XYZ)   0.00 (occ)  -0.29 (B) at cen_x   -8.762   7.280   5.083
  NH1AARG A  18    moved   0.0939 (XYZ)   0.00 (occ)   0.20 (B) at cen_x   -2.505   7.304 -12.378
```
There are indeed other ways to do all these things, but this program is needed in the `$PATH` so it can be called by the above scripts.


## Help

Let me know if you have any questions or find any bugs.  In general, debugging information is provided by adding the command-line option: `debug=1`, and possibly `tempfile=temp` to change the default temporary file prefix.
```
rst2pdb_runme.com amber.rst7 debug=1
```

