#!/usr/bin/env python

import argparse, csv, h5py
import numpy as np
import pyrap.tables as pt
from astropy.coordinates import SkyCoord
from losoto.h5parm import openSoltab

def evaluate_solutions(mtf, threshold = 0.25):
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
    direction = h['/sol000/source'][0][1] # radians
    direction = np.degrees(np.array(direction))
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
        f.write(', {}, {}'.format(direction[0], direction[1]))
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

def make_h5parm(mtf, ms):
    ''' input:    direction of measurement set to be self-calibrated; master
                  file with list of h5parms
        function: find nearest directions; construct a new h5parm that is a
                  combination of best phase solutions from nearest direction,
                  done per antenna. e.g. if the nearest direction has valid
                  solutions for all but the uk antenna, find the uk solutions
                  from the next nearest direction
        output:   a new h5parm to be applied to the measurement set
    '''

    # get the direction from the measurement set
    t  = pt.table(ms, readonly = True, ack = False)
    field = pt.table(t.getkeyword('FIELD'), readonly = True, ack = False)
    ms_direction = field.getcell('PHASE_DIR', 0)[0] # radians
    ms_direction = SkyCoord(ms_direction[0], ms_direction[1], unit = 'rad')
    field.close()
    t.close()

    # get the direction from the master text file
    h5parms, ras, decs = np.genfromtxt(mtf, delimiter = ',', unpack = True, dtype = str, usecols = (0, 1, 2))
    mtf_directions = []

    for h5parm, ra, dec in zip(h5parms, ras, decs):
        mtf_direction = SkyCoord(float(ra), float(dec), unit = 'deg')
        separation = ms_direction.separation(mtf_direction)
        print(separation.arcminute)

    # find the nearest good h5parm direction to the measurement set direction
    # and do this for each station


    # write the phase solutions for the nearest good h5parms to a new h5parm


def applyh5parm():
    ''' input:    the output of make_h5parm; the measurement set for self-
                  calibration
        function: apply the output of make_h5parm to the measurement set
        output:   the measurement set for self-calibration with corrected data
    '''
    pass

def updatelist():
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
    parser.add_argument('-m', '--mastertextfile', help = 'master text file', required = True)
    parser.add_argument('-f', '--measurementset', help = 'measurement set', required = True)
    parser.add_argument('-t', '--threshold', type = float, help = 'threshold determining the xx-yy statistic goodness', default = 0.25)
    args = parser.parse_args()
    mtf = args.mastertextfile
    ms = args.measurementset
    threshold = args.threshold

    # evaluate_solutions(mtf, threshold) # evaluate phase solutions in a h5parm

    make_h5parm(mtf, ms)

    applyh5parm()

    # loop 3

    updatelist()

if __name__ == '__main__':
    main()
