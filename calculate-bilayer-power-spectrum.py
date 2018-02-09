#! /usr/bin/env python

from MDAnalysis.analysis.leaflet import LeafletFinder
import argparse
import MDAnalysis
import numpy
import math
from scipy.interpolate import griddata

def largest_groups(atoms):
    '''
    From a list of sizes, find out the indices of the two largest groups. These should correspond to the two leaflets of the bilayer.

    Keyword arguments:
    atoms -- list of sizes of clusters indentified by LeafletFinder
    '''
    largest=0
    second_largest=0

    for i in atoms:
        if atoms[i]>largest:
            largest_index=i
            largest = atoms[i]

    for i in atoms:
        if atoms[i]>second_largest and i!=largest_index:
            second_largest_index=i
            second_largest = atoms[i]

    return (largest_index,second_largest_index)


def determine_leaflets(universe,phosphateSelection):
    '''
    From a selection of phosphates, determine which belong to the upper and lower leaflets.

    Keyword arguments:
    universe -- an MDAnalysis Universe object
    phosphateSelection -- a string specifying how to select the phosphate atoms (e.g. "name P1")
    '''
    leaflets = {}

    # calculate the z value of the phosphates defining the bilayer (assumes bilayer is in x and y..)
    po4=universe.atoms.select_atoms(phosphateSelection)
    bilayerCentre = po4.center_of_geometry()[2]

    # apply the MDAnalysis LeafletFinder graph-based method to determine the two largest groups which
    #  should correspond to the upper and lower leaflets
    phosphates = LeafletFinder(universe,phosphateSelection)


    # find the two largest groups - required as the first two returned by LeafletFinder, whilst usually are the largest, this is not always so
    (a,b) = largest_groups(phosphates.sizes())

    # check to see where the first leaflet lies
    if phosphates.group(a).centroid()[2] > bilayerCentre:
        leaflets[1] = phosphates.group(a)
        leaflets[2] = phosphates.group(b)
    else:
        leaflets[2] = phosphates.group(a)
        leaflets[1] = phosphates.group(b)

    return leaflets

def write_output_file(power2D,fileName,radiiFactor):
    '''
    Write out a 1D power spectrum to file

    Keyword arguments:
    power2D -- a 2D numpy array of floats (NxN) containing the specified power spectrum
    filename -- string of the file to write the data to (including the PATH)
    radialBin -- float specifying the width of the radial histogram
    radiiFactor --  float of the normalisation factor specified in the main code. Takes into account 2pi etc.
    '''

    # find out the index of (0,0); this automatically accounts for odd/even behaviour
    centre = int(power2D.shape[0]/2.)

    # generate an array of (i,j) relative to the centre
    indices = numpy.indices(power2D.shape)-centre

    # now compute the square distance
    distance2 = indices[0]**2 + indices[1]**2

    r=radiiFactor*numpy.sqrt(distance2)
    bins=numpy.arange(0.01,10,0.01)

    digitized = numpy.digitize(r, bins)

    # occasionally get an error reported here
    bc = numpy.array([r[digitized == i].mean() for i in range(1, len(bins))])
    bm = numpy.array([power2D[digitized == i].mean() for i in range(1, len(bins))])

    bin_centers=bc[numpy.isfinite(bm)]
    bin_means=bm[numpy.isfinite(bm)]

    OUTPUT = open(fileName,"w")

    for i,j in enumerate(bin_centers):
        print >> OUTPUT, "%8.4e %8.4e %8.4e" % (bin_centers[i]+0.005, bin_means[i], ((bin_centers[i]**4)*bin_means[i]))

    OUTPUT.close()


if __name__ == "__main__":

    # use argparse to read in the options from the shell command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb",default="example-trajectory/popc-1500-CG-phosphates.gro",help="the name of the coordinate file to read in (e.g. foo.pdb)")
    parser.add_argument("--traj",default="example-trajectory/popc-1500-CG-phosphates.xtc",help="the name of the trajectory file to read in (e.g. foo.xtc)")
    parser.add_argument("--phosphate",help="the selection text for the phosphate group (default is \'name PO4\' for MARTINI CG sims)",default="name PO4")
    parser.add_argument("--output",help="the stem of the output file (default is 'output-files/output')",default="output-files/output")
    parser.add_argument("--discard",help="the proportion of the trajectory to discard in the range 0.0-1.0 (defaults to analyse the whole trajectory i.e. 0.0)",default=0.5, type=float)
    parser.add_argument("--step", type=float,default=0.5,dest="step",help="the size of the grid in nm (default=1nm)")
    parser.add_argument("--bins",type=int,default=1,help="divide the trajectory into this many bins (default is to analyse the whole trajectory)")
    options = parser.parse_args()

    # check the numeric options
    assert 0.0 <= options.discard < 1.0, "the proportion of the trajectory to discard must be >= 0.0 and < 1"
    assert options.bins >= 1, "the number of bins must be at least 1"
    assert options.step > 0.0, "the size of the grid must be greater than zero"

    # read in the coordinate and trajectory files
    u = MDAnalysis.Universe(options.pdb,options.traj)

    # calculate the frame to start calculating (defaults to the first frame)
    startFrame=int(u.trajectory.n_frames*options.discard)

    # calculate the number of frames in each bin (defaults to just one bin of the whole trajectory)
    binFrames = int((u.trajectory.n_frames*(1.0-options.discard))/options.bins)

    # identify the upper and lower leaflets using the MDAnalysis LeafletFinder method
    leaflets = determine_leaflets(u,options.phosphate)

    # initialise some variables
    binCounter = 0
    firstPass=True
    Power = {}

    # loop through the trajectory
    for ts in u.trajectory:

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
                            c = j.atoms.select_atoms(selectionText[i]).center_of_geometry()
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
                    output_file=options.output+"-"+bin+"-"+i+".dat"
                    write_output_file(Power[i],output_file,radiiFactor)

                # reset the variables
                binCounter+=1
                firstPass=True

    # catch the case when the bin size hasn't been specified so write all the data to file
    if options.bins == 1:

        # normalise by the number of frames
        for i in ['undulation','thickness']:
            Power[i] /= spectralIntensityFrames
            write_output_file(Power[i],options.output+"-"+i+".dat",radiiFactor)
