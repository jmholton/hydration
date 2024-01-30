# hydration

Amber-oriented scripts for adjusting water count on-the-fly.

## Motivation

Getting the number of waters right in a constant-volume amber run is tricky. The pressure is
a traditional indicator, but cannot always be trusted if you are using restraints. Trouble is,
especially in the beginning of the run you want to have restraints to keep the molecule from 
straying too far from it starting point before everything settles down. Unfortunately, the 
number of waters is one of those things you need to decide at the very beginning, and then you 
are married to it. The purpose of this project is to allow modification of the number of waters
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
Fo-Fc difference map as a guide. The trick is to avoid any
"transporter accidents" where two molecules try to occupy the same space, and exploding the run.
Also, if you are using restraints on some of the waters (like I do), you will want to be sure that
these special waters dont get stripped out.

## Getting Started

### Dependencies

* These are linux c-shell scripts, so you will need a working <pre>/bin/tcsh</pre>
* gemmi to measure vacuum bubbles, which is distributed with the CCP4 suite
* the Phenix Suite to probe electron density map values from your mtz file
* AMBER and its tools installed, particularly cpptraj and parmed
* an *.rst7 AMBER restart file for your system, and the accompanying *.parm7 topology file
* an "orignames.pdb" file that contains all the original non-xyz information from the PDB file you initially gave to LEAP.
* an optional current_restraints.pdb file, containing the names and locations of your restraint reference points
* You will also need to have prepared a "padded" topology file, which you can produce by adding ~50k dummy waters at the end of the PDB file you provide to LEAP. You only need to do this once. 

### Installing

* git clone this repo
* copy all the files into somewhere in your shell $PATH

### Executing program

* How to run the program
* Step-by-step bullets
```
code blocks for commands
```

## Help

Any advise for common problems or issues.
```
command to run if program contains helper info
```

