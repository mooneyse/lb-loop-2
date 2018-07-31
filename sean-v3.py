#!/usr/bin/env python

import argparse, csv, h5py
import numpy as np
from losoto.h5parm import openSoltab

def evaluate_solutions(mtf, threshold = 0.5):
    ''' input:    master text file
        function: evaluate solutions for each antenna; use xx-yy statistic
                  described in the hybrid mapping section of the google doc;
                  determine validity
        output:   append ra, dec, and a boolean for validity per station to the
                  master text file
    '''

    # get the h5parm from the bottom of the master text file
    h5parms = np.genfromtxt(mtf, delimiter = ',', unpack = True, dtype = str, usecols = 0)
    h5parm = h5parms[len(h5parms) - 1]

    # get the direction from the h5parm
    h = h5py.File(h5parm, 'r')
    ra_dec = h['/sol000/source'][0][1] # radians
    ra_dec = np.degrees(np.array(ra_dec))
    h.close()

    # get the phase solutions for each station from the h5parm
    phase = openSoltab(h5parm, 'sol000', 'phase000')
    stations = phase.ant[:]
    evaluations = {} # dictionary to hold the statistics for each station

    # calculate xx-yy statistic
    for station in range(len(stations)):
        xx, yy = [], []
        # [polarisation (xx = 0, yy  = 1), direction, station, frequency, time]
        for value in phase.val[0, 0, station, 0, :]:
            xx.append(value)

        for value in phase.val[1, 0, station, 0, :]:
            yy.append(value)

        xx = np.array(xx)
        yy = np.array(yy)
        xx_yy = xx - yy
        mean_xx_yy = np.nanmean(np.abs(xx_yy)) * (1 / (2 * np.pi))
        evaluations[stations[station]] = mean_xx_yy # 0 = best, 1 = worst

    # append to master file
    data = open(mtf, 'r')
    mtf_stations = list(csv.reader(data))[0][3:] # get the stations from the mtf

    with open(mtf, 'a') as f:
        f.write(', {}, {}'.format(ra_dec[0], ra_dec[1]))
        for mtf_station in mtf_stations:
            # look up the statistic for a station and determine if it is good
            try:
                value = evaluations[mtf_station[1:]]
            except KeyError:
                value = float('nan')

            if value < threshold: # pass
                f.write(', {}'.format(int(True)))
            elif np.isnan(value):
                f.write(', {}'.format('nan'))
            else: # fail
                f.write(', {}'.format(int(False)))

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
    parser.add_argument('-m', '--mastertextfile', help = 'the master text file', required = True)
    parser.add_argument('-t', '--threshold', type = float, help = 'the threshold determining the xx-yy statistic goodness', default = 0.5)
    args = parser.parse_args()
    mtf = args.mastertextfile
    threshold = args.threshold

    # evaluate the phase solutions in the h5parm
    evaluate_solutions(mtf, threshold)

    make_h5parm()

    applyh5parm()

    # loop 3

    updatelist()

if __name__ == '__main__':
    main()
