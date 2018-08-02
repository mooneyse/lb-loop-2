#!/usr/bin/env python

import argparse, csv, h5py, logging, os
import numpy as np
import pyrap.tables as pt
import losoto.h5parm as lh5
from astropy.coordinates import SkyCoord

# TODO add logging statements

def evaluate_solutions(mtf, threshold = 0.25):
    '''
    Description:
    - read the name of the h5parm from the master text file
    - get the direction from the h5parm
    - evaluate the phase solutions in the h5parm for each station using the xx-yy statistic
    - determine the validity of each xx-yy statistic that was calculated statistic
    - append the right ascension, declination, and validity to the master text file

    Parameters:
    - mtf      (str)            : master text file
    - thresold (float, optional): threshold to determine the goodness of the xx-yy statistic

    Returns:
    - None
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
    phase = lh5.openSoltab(h5parm, 'sol000', 'phase000') # TODO h5parm only closes on exit
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

    with open(mtf) as f:
        mtf_stations = list(csv.reader(f))[0][3:] # get stations from the mtf

    # append to master file
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
    '''
    Description:
    - get the direction from the measurement set
    - get the directions of the h5parms from the master text file
    - calculate the separation between the measurement set direction and the h5parm directions
    - for each station, find the h5parm of smallest separation which has valid phase solutions
    - create a new h5parm
    - write these phase solutions to this new h5parm

    Parameters:
    - mtf (str): master text file with list of h5parms
    - ms  (str): measurement set to be self-calibrated

    Returns:
    - new_h5parm (str): the new h5parm to be applied to the measurement set
    '''

    # get the direction from the measurement set
    t  = pt.table(ms, readonly = True, ack = False)
    field = pt.table(t.getkeyword('FIELD'), readonly = True, ack = False)
    ms_direction = field.getcell('PHASE_DIR', 0)[0] # radians
    ms_direction = SkyCoord(ms_direction[0], ms_direction[1], unit = 'rad')
    field.close()
    t.close()

    # get the direction from the master text file
    # BUG genfromtxt gives empty string for h5parms when names = True is used; importing them separately as a work around
    data = np.genfromtxt(mtf, delimiter = ',', unpack = True, dtype = float, names = True)
    h5parms = np.genfromtxt(mtf, delimiter = ',', unpack = True, dtype = str, usecols = 0)
    mtf_directions = {}

    # calculate the distance betweeen the ms direction and the h5parm directions
    for h5parm, ra, dec in zip(h5parms, data['ra'], data['dec']):
        mtf_direction = SkyCoord(float(ra), float(dec), unit = 'deg')
        separation = ms_direction.separation(mtf_direction)
        mtf_directions[separation] = h5parm # distances from ms to each h5parm

    # read in the stations from the master text file
    with open(mtf) as f: # get stations from the mtf
        mtf_stations = list(csv.reader(f))[0][3:] # skipping h5parm, ra, and dec
        mtf_stations = [x.lstrip() for x in mtf_stations] # remove leading space

    # find the closest h5parm which has an acceptable solution for each station
    # NOTE pandas could probably do better than this
    # these print statements are for testing only
    print('for this direction in the ms, make a new h5parm consisting of...')
    print('Station\t\tSeparation\th5parm\t\tRow\tBoolean')
    successful_stations = []

    for mtf_station in mtf_stations: # for each station
        for key in sorted(mtf_directions.keys()): # starting with shortest separations
            h5parm = mtf_directions[key]
            row = list(h5parms).index(h5parm) # row in mtf
            value = data[mtf_station][row] # boolean value for h5parm and station
            if value == 1 and mtf_station not in successful_stations:
                if mtf_station == 'ST001':
                    print('{}\t\t{}\t{}\t{}\t{}'.format(mtf_station, round(key.deg, 6), h5parm, row, int(value)))
                else:
                    print('{}\t{}\t{}\t{}\t{}'.format(mtf_station, round(key.deg, 6), h5parm, row, int(value)))
                successful_stations.append(mtf_station)

    # create a new h5parm
    ms = os.path.splitext(os.path.normpath(ms))[0]
    new_h5parm = '{}_{}_{}.h5'.format(ms, ms_direction.ra.deg, ms_direction.dec.deg)

    # write these best phase solutions to the new h5parm
    h = lh5.h5parm(new_h5parm, readonly = False)
    h.makeSolset(addTables = False)
    solset = h.getSolset('sol000')
    a = [0, 1] # dummy data
    b = np.array([a, a])
    c = solset.makeSoltab('phase', axesNames = ['freq', 'time'], axesVals = [a, a], vals = b, weights = b)
    h.close()

    return new_h5parm

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
    # evaluate_solutions(mtf, threshold)
    pass

def main():
    ''' starting point: 1st iteration will have produced phase solutions
        independently for all directions, with each being its own h5parm
    '''

    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mastertextfile', required = True, help = 'master text file')
    parser.add_argument('-f', '--measurementset', required = True, help = 'measurement set')
    parser.add_argument('-t', '--threshold', type = float, default = 0.25, help = 'threshold determining the xx-yy statistic goodness')
    args = parser.parse_args()
    mtf = args.mastertextfile
    ms = args.measurementset
    threshold = args.threshold

    evaluate_solutions(mtf, threshold) # evaluate phase solutions in a h5parm

    new_h5parm = make_h5parm(mtf, ms) # create a new h5parm of the best solutions

    applyh5parm(new_h5parm, ms)

    # loop 3

    updatelist(new_h5parm)

if __name__ == '__main__':
    main()
