#!/usr/bin/env python

import argparse, h5py
import numpy as np
from losoto.h5parm import openSoltab

def evaluate_solutions(h5parm, mtf):
    ''' input:    h5parm with phase solutions
        function: evaluate solutions for each antenna; use xx-yy statistic
                  described in the hybrid mapping section of the google doc;
                  determine validity
        output:   append h5parm_name, ra, dec, and one column per antenna with
                  boolean for validity to the master text file
    '''

    # get (ra, dec) for the source from the h5parm
    h = h5py.File(h5parm, 'r')
    ra_dec = h['/sol000/source'][0][1] # radians
    ra_dec = np.degrees(np.array(ra_dec))
    h.close()

    # get the phase solutions for each station from the h5parm
    phase = openSoltab(h5parm, 'sol000', 'phase000')
    stations = phase.ant[:]
    evaluations = {}

    for station in range(len(stations)): # loop over all stations
        xx, yy = [], []

        # [polarisation (xx = 0, yy  = 1), direction, station, frequency, time]
        for value in phase.val[0, 0, station, 0, :]:
            xx.append(value)

        for value in phase.val[1, 0, station, 0, :]:
            yy.append(value)

        # calculate xx-yy statistic
        xx = np.degrees(np.array(xx))
        yy = np.degrees(np.array(yy))
        xx_yy = xx - yy
        mean_xx_yy = np.nanmean(np.abs(xx_yy)) * (2 / np.pi)
        evaluations[stations[station]] = mean_xx_yy

    # append to master file
    with open(mtf, 'a') as f:
        f.write('{}, {}, {}'.format(h5parm, ra_dec[0], ra_dec[1]))
        for value in evaluations.values():
            f.write(', {}'.format(value)) # NB ensure written in correct order
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
    parser.add_argument('-m', '--mastertextfile', required = True)
    args = parser.parse_args()
    h5parm = args.h5parm
    mtf = args.mastertextfile

    # evaluate the phase solutions in the h5parm
    evaluate_solutions(h5parm, mtf)

    make_h5parm()

    applyh5parm()

    # loop 3

    updatelist()

if __name__ == '__main__':
    main()
