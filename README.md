# hydration

Amber-oriented scripts for adjusting water count on-the-fly in a molecular dynamics simulation

## Motivation

Getting the number of waters right in a constant-volume amber run is tricky. The pressure is
a traditional indicator, but cannot always be trusted if you are using restraints. Trouble is,
especially in the beginning of the run you usually want to have restraints to keep the molecule from 
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
* `AMBER` and its tools installed, particularly cpptraj and parmed
* an `*.rst7` AMBER restart file for your system, and the accompanying `*.parm7` topology file
* an `orignames.pdb` file that contains all the original non-xyz information from the PDB file you initially gave to LEAP.
* an optional `current_restraints.pdb` file, containing the names and xyz positions of your restraint reference points
* You will also need to have prepared a "padded" topology file, which you can produce by adding ~50k dummy waters at the end of the PDB file you provide to LEAP. You only need to do this once. 

### Installing

* git clone this repo
* copy all the files into somewhere in your shell `$PATH`, and make them executable:
  chmod u+x *.com *.awk
Yes, I know the extension says `*.com`, but these are not Windows executables. The use of `.com` to denote shell scripts pre-dates Windows.

### Executing program

#### file format conversion 
The first program you may want to run is:
```
rst2pdb_runme.com amber.rst7 orignames=orignames.pdb parmfile=xtal.prmtop
```
This will serve as a good test to see if everything is working. Success is indicated by the output of a new pdb formatted file called `amber.pdb`, which will have the same names as the atoms in `orignames.pdb`, but with the
xyz coordinates contained in the amber.rst7 file. The options shown for `orignames` and `parmfile` are the internally-stored defaults, so are not required if you already have those files with those names. I recommend using these names, but if you don't like those
names for some reason, I recommed editing the defaults at the top of the `rst2pdb_runme.com` script.<br>
Another internal default is `paddedparm=padded.parm7`, which is the topology file you will want to generate yourself
using LEAP, but with a few thousand extra dummy waters tacked on to the end. If for some reason the `xtal.prmtop` is missing or otherwise incompatible with `amber.rst7`, the `rst2pdb_runme.com` script will look for `padded.parm7` and strip out waters until it can be used to read `amber.rst7`. The new, smaller, topology file will be written out
as `resized.parm7` as well as a shortened version of `orignames.pdb`. However, it is okay, and indeed recommended, for `orignames.pdb` to be longer than your active system. In fact, make it as big as `padded.parm7`.<br>

You can change the output file prefix using `outpreifx=whatever`, or by specifying a pdb file name on the command line. In general, for all these scripts, any variable you see at the beginning of the script can be changed on the command line using the `variable=value` syntax.

#### Abhoring the Vacuum
The next thing you probably want to do is check to see if your current simulation has vacuum bubbles:
```
bubble_check_runme.com amber.pdb Vwater=96 minvoid=5
```
This will use `gemmi` to look for spaces between atoms and measure their volumes. In bulk water, there is about 96 A^3 of volume for each water molecule, which is the default for `Vwater`, so you don't have to specify it to take the default. The `minvoid` variable is the number of water molecules that a vacuum bubble can hold before it is considered a "void". How big are vacuum bubbles in normal bulk water at STP? I have no idea. But 5 strikes me as a bit much, so it is the default. If you have a vacuum bubble that can hold thousands of waters, then you might do well by filling it in.  This script will print out a reocmmended fill-in command, such as:<br>

#### Hydrate
The next thing you probably want to do is check to see if your current simulation has vacuum bubbles:
```
hydrarte_runme.com amber.pdb nadd=1000 outrst=wetter.rst7 outtop=xtal.prmtop
```
The output will be the result of running the AmberTools AddToBox program.  The last molecule in amber.pdb is taken as the "water" to add, thus preserving any extra points. In the end, you will have a new restart file, including velocities, and topology file to continue your simulation, but now with a few more waters. The filename values in the above command are the defaults. The key to adding in new waters to an AMBER parm file is the unmentioned `padded.parm7` file. The script will complain and exit if you don't have `padded.parm7`. It works by not adding extra waters to your previous `*.parm7` file, but rather by stripping an appropriate number of waters out of `padded.parm7` so that it is compatible with the new number of waters.<br> 
All new waters are given zero velocity, but will quickly heat up when the simulation restarts. I have not seen any need for an explicit re-heating protocol in my simulations so far.<br>
Also note: you can provide an `*.rst7` file on the command line instead of a PDB file and `hydrate_runme.com` will call `rst2pdb_runme.com` internally to generate the `*.pdb` file needed by AddToBox.<br>
If you provide a `restraints=current_restraints.pdb` file on the command line, this script will also update the reference points file called `ref.crd` so that it is also the same size as the new `wetter.rst7`, using the XYZ coordinates in `current_restraints.pdb`.

#### De-hydrate
Removing waters from an AMBER simulation is much easier than adding them. The existing and popular `cpptraj` and `parmed` programs can do this with the `strip` command. A good question, however, is: which ones?
```
dehydrarte_amber_runme.com toowet.rst7 maxreject=10 outprefix=drier.rst7
```


## Help

Let me know if you have any questions or find any bugs.  In general, debugging information is provided by adding the command-line option: `debug=1`, and possibly `tempfile=temp` to change the default temporary file prefix.
```
rst2pdb_runme.com amber.rst7 debug=1
```

