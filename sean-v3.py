#!/usr/bin/env python3

import argparse, h5py
import numpy as np
from losoto.h5parm import h5parm as losoto_h5parm
from losoto.h5parm import openSoltab

def evaluate_solutions(h5parm, master): # subroutine 1
    ''' input:    h5parm with phase solutions
        function: evaluate solutions for each antenna; use xx-yy statistic
                  described in the hybrid mapping section of the google doc and
                  determine validity
        output:   Append to master file with h5parm_name, ra, dec, and one
                  column per antenna with boolean for validity
    '''

    # get (ra, dec) for the source from the h5parm
    H = h5py.File(h5parm, 'r')
    ra_dec = H['/sol000/source'][0][1] # [ra, dec] in degrees

    # get the phase solutions for each antenna from the h5parm
    h = losoto_h5parm(h5parm, readonly = True)
    phase = openSoltab(h5parm, 'sol000', 'phase000')

    stations = phase.ant[:]

    # loop over all stations and for each polarisation, I have a list with one
    # element per per time step, and each time step is a list of values with one
    # for each subband
    for station in range(len(stations)):
        xx, yy = [], []

        # [polarisation (xx = 0, yy  = 1), direction (0), station, frequency, time]
        for value in phase.val[0, 0, station, :, :]:
            xx.append(value)

        for value in phase.val[1, 0, station, :, :]:
            yy.append(value)

    a = xx[0]
    print(len(xx))
    print(a)
    print(a[0])
    print(len(a))

    # scatter in XX-YY phase (get all phase solutions on a station, subtract
    # XX  and yy, reduce difference from -180 to 180, take the mean of the
    # absolute difference value and multiply by 2/pi)

    evaluations = [1, 0, 1]

    # append to master file
    with open(master, 'a') as f:
        f.write('{}, {}, {}'.format(h5parm, ra_dec[0], ra_dec[1]))
        for e in evaluations:
            f.write(', {}'.format(e))
        f.write('\n')

def make_h5parm(): # subroutine 2
    ''' input:    direction of measurement set to be self-calibrated; master
                  file with list of h5parms
        function: find nearest directions; construct a new h5parm that is a
                  combination of best phase solutions from nearest direction,
                  done per antenna. e.g. if the nearest direction has valid
                  solutions for all but the uk antenna, find the uk solutions
                  from the next nearest direction
        output:   a new h5parm to be applied to the measurement set
    '''
    pass

def applyh5parm(): # subroutine 3
    ''' input:    the output of make_h5parm; the measurement set for self-
                  calibration
        function: apply the output of make_h5parm to the measurement set
        output:   the measurement set for self-calibration with corrected data
    '''
    pass

def updatelist(): # subroutie 4
    ''' input:    the initial h5parm (i.e. from make_h5parm); the final h5parm
                  from loop 3
        function: combine the phase solutions from the initial h5parm and the
                  final h5parm; calls evaluate_solutions to update the master
                  file
        output:   a new h5parm that is a combination of both these h5parms; a
                  new line in the master file
    '''

    # update the master file
    # evaluate_solutions()
    pass

def main():
    ''' starting point: 1st iteration will have produced phase solutions
        independently for all directions, with each being its own h5parm
    '''

    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--h5parm', required = True)
    parser.add_argument('-m', '--master', required = True)
    args = parser.parse_args()
    h5parm = args.h5parm
    master = args.master

    # subroutine 1: evaluate the phase solutions in the h5parm
    evaluate_solutions(h5parm, master)

    make_h5parm()

    applyh5parm()

    # loop 3

    updatelist()

if __name__ == '__main__':
    main()
