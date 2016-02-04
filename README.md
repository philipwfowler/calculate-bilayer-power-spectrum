# calculate-bilayer-power-spectrum.py

## Introduction

This simple python code accompanies the paper "Membrane Stiffness is Modified by Integral Membrane Proteins" which is currently under consideration at a scientific journal. By making this code public we hope that other researchers will be more easily be able to reproduce our results.

## How to use

### Requirements

You will need to install
- python 2.7.X
- numpy and scipy
- [MDAnalysis](http://www.mdanalysis.org). Version 0.8.1 used for development and testing.

On a Mac python, numpy and scipy can be easily installed using either [MacPorts](http://www.macports.org) or [Homebrew](http://brew.sh). The other python modules should already be present on your system (argparse and math). 

### Quick start

Although the code is setup with command line flags so you can analyse your own simulation data, coordinate and trajectory files of the last 20% (i.e. 4-5µs) of the 1500 POPC coarse-grained [MARTINI](http://cgmartini.nl) simulation is included in the `example-trajectory/` folder. For simplicity these only contain the phosphate (i.e. PO4) beads as the code uses these to construct the two surfaces of the bilayer. Defaults in the code are setup so that if you simply run the following in a terminal

    ./calculate-bilayer-power-spectrum.py 
 
 It will use the same options as in the paper and analyse the small 1500 CG POPC simulation. The data files are written to 
 
    ls output-files/
    output-thickness.dat  output-undulation.dat

In both cases plotting column 1 against column 2 on a log-log plot will give you the power spectra shown in the upper panel of Figure 2C of the paper. 

### Options

The code is also setup with command line options. To see these

    ./calculate-bilayer-power-spectrum.py --help

and you should get
```
usage: calculate-bilayer-power-spectrum.py [-h] [--pdb PDB] [--traj TRAJ]
                                           [--phosphate PHOSPHATE]
                                           [--output OUTPUT]
                                           [--discard DISCARD] [--step STEP]
                                           [--bins BINS]

optional arguments:
  -h, --help            show this help message and exit
  --pdb PDB             the name of the coordinate file to read in (e.g.
                        foo.pdb)
  --traj TRAJ           the name of the trajectory file to read in (e.g.
                        foo.xtc)
  --phosphate PHOSPHATE
                        the selection text for the phosphate group (default is
                        'name PO4' for MARTINI CG sims)
  --output OUTPUT       the stem of the output file (default is 'output-
                        files/output')
  --discard DISCARD     the proportion of the trajectory to discard in the
                        range 0.0-1.0 (defaults to analyse the whole
                        trajectory i.e. 0.0)
  --step STEP           the size of the grid in nm (default=0.5nm)
  --bins BINS           divide the trajectory into this many bins (default is
                        to analyse the whole trajectory)
```
