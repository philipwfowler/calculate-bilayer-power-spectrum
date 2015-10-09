#! /usr/bin/env python

import MDAnalysis
import numpy
from scipy.interpolate import griddata
import argparse
import math
import matplotlib.pyplot as plt
#import progressbar
import time
from spectralintensity import *


parser = argparse.ArgumentParser()
parser.add_argument("--traj",help="XTC/DCD to read in",default="popc-1500-01-5us-dt600ns-nosol.xtc")
parser.add_argument("--phosphate",help="the selection text for the phosphate group",default="name PO4")
parser.add_argument("--tail",help="the selection text for the tail group",default="name C4A or name C4B")
parser.add_argument("--head",help="the selection text for the head group",default="name PO4 or name GL1")
parser.add_argument("--output",help="the path of the output file",default="output.dat")
parser.add_argument("--discard",help="the proportion of the trajectory to discard",default=0.0, type=float)
parser.add_argument("--step", type=float,default=1.0,dest="step",help="the size of the grid in nm")
parser.add_argument("--radialbin", type=float,default=0.5,dest="radialbin",help="the size of the radial averaging bin width in nm")
parser.add_argument("--bins",type=int,default=1,help="divide the trajectory into this many bins (default is to analyse the whole trajectory)")
parser.add_argument("--plot-fft",dest="plotFFT",action="store_true",help="plot the spectral densities in 2D")
parser.set_defaults(plotFFT=False)

options = parser.parse_args()

# parse the trajectory name in the standard format
(systemName, systemSize, systemRun, systemLength, systemStep) = parse_trajectory_name(options.traj)

# make the output path, creating folders if necessary
path = "dat/"
for i in [systemName,systemSize,systemRun]:    
    path += i + "/"
    create_directory(path)
    
if options.bins!=1:
    path += "bin" +str(options.bins)+"discard"+str(options.discard)+"/"
    create_directory(path)

# hack for parallel processing
if systemRun =='c36':
    options.phosphate = "name P"

# define the output file stem        
stem=systemName+"-"+systemSize+"-"+systemRun

pdbFile = "../gro/" + stem + "-small.gro"

# path+=stem+"-"+systemLength+"-"+systemStep+"-s"+str(options.step)+"-r"+str(options.radialbin)
path+=stem+"-"+systemLength+"-"+systemStep

# create a simple progress bar that will be updated on the command line
# progressBar = progressbar.AnimatedProgressBar(end=100, width=75)
# progressBar.show_progress()
# startTime = time.time()

# read in the coordinate and trajectory files
u = MDAnalysis.Universe(pdbFile,options.traj)

# defaults to zero
startFrame=int(u.trajectory.numframes*options.discard)
    
binFrames = int((u.trajectory.numframes*(1.0-options.discard))/options.bins)

# progressBar+2
# progressBar.show_progress()
# progressPerFrame = 95.0/(u.trajectory.numframes)

# find out the box dimensions
halfway = u.trajectory.ts.dimensions[0]/2.0

# identify the upper and lower leaflets
leaflets = determine_leaflets(u,options.phosphate,halfway)
            
# progressBar+3
# progressBar.show_progress()

binCounter = 0        
firstPass=True
Power = {}

# loop through the trajectory
for ts in u.trajectory:
    
    # progressBar+progressPerFrame
    # progressBar.show_progress()
    
    # proceed only if the time is in the correct range
    if ts.frame > startFrame:

        # what are the dimensions of the structure file? MDAnalysis reports this in Angstroms
        (box_x,box_y) = (u.dimensions[0]/10.,u.dimensions[1]/10.) # convert to nm

        # create a pixel coordinate mesh grid
        P,Q = numpy.mgrid[0:box_x:options.step, 0:box_y:options.step]                

        # how many pixels to a side?
        numberPixels = P.shape

        # create the q-vector (note the 2*pi)
        q = (2*math.pi*numpy.fft.fftfreq(numberPixels[0],options.step))
        
        # if this is the first frame, define some arrays and variables
        if firstPass:                        

            # setup the power arrays
            for i in ['undulation','thickness']:
                Power[i] = numpy.zeros(numberPixels)
            
            # store the size 
            grid_size=q.size
            
            # calculate the normalisation factor for the final array
            radiiFactor = 2*math.pi/(numberPixels[0]*options.step)
            
            # reset the variables
            firstPass=False
            frameCounter=0
            spectralIntensityFrames=0
        
        frameCounter+=1

        # check that the grids have the same size, otherwise skip this frame    
        if q.size == grid_size:
    
        	# what frame are we analysing?
            bin="%04i" % binCounter            

            # initialise the dictionaries
            coordinates={}
            mesh={}

            selectionText = {"phosphate": options.phosphate}

            for leaflet in [1,2]:

                # store the coordinates of the different groups (in Angstroms) in the coordinates dictionary
                for i in selectionText:            
                    coord = []
                    for j in leaflets[leaflet].residues:
                        c = j.selectAtoms(selectionText[i]).centerOfGeometry()
                        coord.append(c/10.) # convert to nm immediately                        
                    coordinates[i,leaflet]=numpy.array(coord)
                    
                # split out the xy and z coordinates separately of the Monge surface
                xy = numpy.delete(coordinates['phosphate',leaflet],-1,1)        
                z = coordinates['phosphate',leaflet][0:,2] - numpy.average(coordinates['phosphate',leaflet][0:,2])
                                
                # wrap them to avoid edge defects
                xy_wrapped = numpy.vstack((xy,xy+[box_x,0],xy-[box_x,0],xy+[0,box_y],xy-[0,box_y],xy+[box_x,box_y],xy-[box_x,box_y],xy+[box_x,-box_y],xy-[box_x,-box_y]))
                z_wrapped = numpy.hstack((z,z,z,z,z,z,z,z,z))
                
                # cubic interpolate to create a mesh describing the surface                
                #  although xy_wrapped, z_wrapped cover a 3x3 area, the result only covers the central box
                mesh['z',leaflet] = griddata(xy_wrapped, z_wrapped, (P, Q), method="cubic")
                
            # calculate the height and thickness fields
            mesh['undulation'] = (mesh['z',1] + mesh['z',2])/2.0
            mesh['thickness'] = (mesh['z',1] - mesh['z',2])/2.0
    
            FFT={}
            for i in ['undulation','thickness']:
                FFT[i] = numpy.fft.fft2(mesh[i])
                FFT[i] /= len(FFT[i])
                FFT[i] *= options.step                
                Power[i] += numpy.abs(numpy.fft.fftshift(FFT[i]))**2                    

            # add two since are adding the intensity in x and y
            spectralIntensityFrames += 1
            
        if ts.frame % binFrames == 0 and ts.time > 0 and options.bins != 1:

            # normalise by the number of frames
            for i in ['undulation','thickness']:
                Power[i] /= spectralIntensityFrames
                write_output_file(Power[i],path+"-"+i+"-"+bin+".dat",options.radialbin,radiiFactor)
            
            # reset the variables
            binCounter+=1
            firstPass=True

# catch the case when the bin size hasn't been specified so write all the data to file
if options.bins == 1:

    # normalise by the number of frames
    for i in ['undulation','thickness']:
        Power[i] /= spectralIntensityFrames
        write_output_file(Power[i],path+"-"+i+".dat",options.radialbin,radiiFactor)
        if options.plotFFT:
            plt.imshow(Power[i],origin='lower',extent=[min(q),max(q),min(q),max(q)],cmap="Greys")
            plt.savefig(path+"-power-"+i+"-lin.pdf")
            plt.imshow(numpy.log(Power[i]),origin='lower',extent=[min(q),max(q),min(q),max(q)],cmap="Greys")
            plt.savefig(path+"-power-"+i+"-log.pdf")
        
# print "\n",spectralIntensityFrames,frameCounter    
# endTime = time.time()
# print "Program took %7.1f seconds to run" % (endTime-startTime)
