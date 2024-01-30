# hydration

Amber-oriented scripts for adjusting water count on-the-fly.

## Description

Getting the number of waters right in a constant-volume amber run is tricky. These scripts implement
a method for taking the constant N out of a constant NVT simulation. That is, you can safely add,
 remove or teleport waters in a simulation that is alredy
flying. They work off of the *.rst7 restart file

## Getting Started

### Dependencies

* These are linux c-shell scripts, so you will need a working /bin/tcsh
* You will also need AMBER installed
* You will need gemmi to measure vacuum bubbles, which is distributed with the CCP4 suite
* You will need the Phenix Suite to probe electron density map values from your mtz file

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

