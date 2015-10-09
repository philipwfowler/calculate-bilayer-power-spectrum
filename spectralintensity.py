#! /usr/bin/env python

import os
import MDAnalysis
import agpy
import math
from MDAnalysis.analysis.leaflet import LeafletFinder

def parse_trajectory_name(trajName):
    # parse the trajectory name
    tmp = trajName.split("/")[2]
    moo=tmp.split("-")
    return (moo[0],moo[1],moo[2],moo[3],moo[4])

def create_directory(path):
    try:
        os.makedirs(path)
    except OSError:
        pass

def determine_leaflets(universe,phosphateSelection,halfway):

    leaflets = {}

    bilayerCentre = universe.selectAtoms(phosphateSelection).centerOfGeometry()[2]

    if universe.trajectory.numatoms<30000:
    
        phosphates = LeafletFinder(universe,phosphateSelection)

        # check to see where the first leaflet lies
        if phosphates.group(0).centroid()[2] > bilayerCentre:
            leaflets[1] = phosphates.group(0)
            leaflets[2] = phosphates.group(1)
        else:
            leaflets[2] = phosphates.group(0)
            leaflets[1] = phosphates.group(1)
    else:		

        phosphates1 = LeafletFinder(universe,phosphateSelection + ' and prop x <= ' + str(halfway))
        phosphates2 = LeafletFinder(universe,phosphateSelection + ' and prop x > ' + str(halfway))

        if phosphates1.group(0).centerOfGeometry()[2] > bilayerCentre:
            if phosphates2.group(0).centerOfGeometry()[2] > bilayerCentre:
                leaflets[1] = phosphates1.group(0) + phosphates2.group(0)
                leaflets[2] = phosphates1.group(1) + phosphates2.group(1)            
            else:
                leaflets[1] = phosphates1.group(0) + phosphates2.group(1)
                leaflets[2] = phosphates1.group(1) + phosphates2.group(0)            
        else:
            if phosphates2.group(0).centerOfGeometry()[2] > bilayerCentre:
                leaflets[1] = phosphates1.group(1) + phosphates2.group(0)
                leaflets[2] = phosphates1.group(0) + phosphates2.group(1)                        
            else:
                leaflets[1] = phosphates1.group(1) + phosphates2.group(1)
                leaflets[2] = phosphates1.group(0) + phosphates2.group(0)            

    return leaflets

def write_output_file(power2D,fileName, radialBin,radiiFactor):
    
    radii,profile = agpy.azimuthalAverage(power2D,returnradii=True,binsize=radialBin)
    
    # if 'q_radius' not in locals():
    radii *= radiiFactor

    OUTPUT = open(fileName,"w")
    
    # don't write out the zero frequency power    
    for i in range(1,radii.size):
        
        # if the bin size is small, then some bins will contain NaN. Filter these out.
        if not math.isnan(profile[i]):
            print >> OUTPUT, "%8.4e %8.4e %8.4e" % (radii[i], profile[i],(profile[i]*(radii[i]**4)))
    OUTPUT.close()
    