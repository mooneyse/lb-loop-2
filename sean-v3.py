#!/usr/bin/env python

'''
A collection of functions for modifying HDF5 files.
'''

import argparse, csv, datetime, h5py, logging, os, subprocess, sys, threading
import numpy as np
import pyrap.tables as pt
import losoto.h5parm as lh5
from astropy.coordinates import SkyCoord
from pathlib import Path

def does_it_exist(the_file, clobber = False, append = False):
    # convenience function to check if a file already exists
    if Path(the_file).is_file():
        if append:
            logging.warn('{} already exists and it will be appended to (append = {})'.format(the_file, append))
        elif clobber:
            logging.warn('{} already exists but it will be overwritten (clobber = {})'.format(the_file, clobber))
            os.remove(the_file)
        else:
            logging.error('{} already exists and it will not be overwritten (clobber = {})'.format(the_file, clobber))
            sys.exit()
        return True
    else:
        logging.info('{} does not exist'.format(the_file))
        return False

def ascii(x):
    # convenience function to encode in ascii
    return x.encode('ascii', 'ignore')

def loop3():
    '''
    description:
    - calls loop 3

    parameters:

    returns:
    - loop3_h5parm (str): h5parm resulting from loop 3
    '''

    logging.info('executing loop3()')
    logging.info('loop3() completed')
    return '/data/scratch/sean/h5parms/loop3_h5parm.h5'

def evaluate_solutions(h5parm, mtf, threshold = 0.25):
    '''
    description:
    - get the direction from the h5parm
    - evaluate the phase solutions in the h5parm for each station using the xx-yy statistic
    - determine the validity of each xx-yy statistic that was calculated
    - append the right ascension, declination, and validity to the master text file

    parameters:
    - h5parm   (str)            : lofar hdf5 parameter file
    - mtf      (str)            : master text file
    - thresold (float, optional): threshold to determine the goodness of the xx-yy statistic

    returns:
    - none
    '''

    logging.info('executing evaluate_solutions(h5parm = {}, mtf = {}, threshold = {})'.format(h5parm, mtf, threshold))

    # get the direction from the h5parm
    h = h5py.File(h5parm, 'r') # TODO should probably use losoto for this
    direction = h['/sol000/source'][0][1] # radians
    direction = np.degrees(np.array(direction))
    h.close()

    # get the phase solutions for each station from the h5parm
    logging.info('opening {}'.format(h5parm))
    # NOTE the convenience function openSoltab does not close the h5parm so it is not used here
    lo = lh5.h5parm(h5parm, readonly = False)
    phase = lo.getSolset('sol000').getSoltab('phase000')
    logging.info('got the phase solution tab (phase000) from {}'.format(h5parm))
    stations = phase.ant[:]
    logging.info('got the stations (i.e. antennas) from {}'.format(h5parm))
    evaluations = {} # dictionary to hold the statistics for each station

    # calculate xx-yy statistic
    logging.info('evaluating the xx-yy statistic for {}'.format(h5parm))
    logging.info('a value < {} is good (1), otherwise the value is bad (0)'.format(threshold))
    for station in range(len(stations)):
        xx, yy = [], []
        # structure: phase.val[polarisation (xx = 0, yy  = 1), direction, station, frequency, time]
        for value in phase.val[0, 0, station, 0, :]:
            xx.append(value)

        for value in phase.val[1, 0, station, 0, :]:
            yy.append(value)

        xx = np.array(xx)
        yy = np.array(yy)
        xx_yy = xx - yy
        mean_xx_yy = np.nanmean(np.abs(xx_yy)) * (1 / (2 * np.pi))
        evaluations[stations[station]] = mean_xx_yy # 0 = best, 1 = worst

    # append to master file if it exists, else write
    if not does_it_exist(mtf, append = True):
        lo.close()
        sys.exit()
    logging.info('writing the results to the master text file {}'.format(mtf))
    with open(mtf) as f:
        mtf_stations = list(csv.reader(f))[0][3:] # get stations from the mtf

    with open(mtf, 'a') as f:
        f.write('{}, {}, {}'.format(h5parm, direction[0], direction[1]))
        for mtf_station in mtf_stations:
            # look up the statistic for a station and determine if it is good
            try:
                value = evaluations[mtf_station[1:]]
            except KeyError:
                value = float('nan')

            if value < threshold: # success
                f.write(', {}'.format(int(True)))
            elif np.isnan(value):
                f.write(', {}'.format('nan'))
            else: # fail
                f.write(', {}'.format(int(False)))

        f.write('\n')

    lo.close()
    logging.info('evaluate_solutions(h5parm = {}, mtf = {}, threshold = {}) completed'.format(h5parm, mtf, threshold))

def make_h5parm(mtf, ms, clobber = False):
    '''
    description:
    - get the direction from the measurement set
    - get the directions of the h5parms from the master text file
    - calculate the separation between the measurement set direction and the h5parm directions
    - for each station, find the h5parm of smallest separation which has valid phase solutions
    - create a new h5parm
    - write these phase solutions to this new h5parm

    parameters:
    - mtf (str): master text file with list of h5parms
    - ms  (str): measurement set to be self-calibrated

    returns:
    - new_h5parm (str): the new h5parm to be applied to the measurement set
    '''

    logging.info('executing make_h5parm(mtf = {}, ms = {}, clobber = {})'.format(mtf, ms, clobber))

    # get the direction from the measurement set
    t  = pt.table(ms, readonly = True, ack = False)
    field = pt.table(t.getkeyword('FIELD'), readonly = True, ack = False)
    ms_direction = field.getcell('PHASE_DIR', 0)[0] # radians
    ms_direction = SkyCoord(ms_direction[0], ms_direction[1], unit = 'rad')
    field.close()
    t.close()

    # get the direction from the master text file
    # HACK genfromtxt gives empty string for h5parms when names = True is used; importing them separately as a work around
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

    logging.info('for this direction in the ms, make a new h5parm consisting of the following:')
    logging.info('\tstation \tseparation\tboolean\trow\t5parm')
    successful_stations = []

    for mtf_station in mtf_stations: # for each station
        for key in sorted(mtf_directions.keys()): # starting with shortest separations
            h5parm = mtf_directions[key]
            row = list(h5parms).index(h5parm) # row in mtf
            value = data[mtf_station][row] # boolean value for h5parm and station
            if value == 1 and mtf_station not in successful_stations:
                logging.info('\t{}\t{}\t{}\t{}\t{}'.format(mtf_station.ljust(8), round(key.deg, 6), int(value), row, h5parm))
                successful_stations.append(mtf_station)

    # create a new h5parm
    ms = os.path.splitext(os.path.normpath(ms))[0]
    new_h5parm = '{}_{}_{}.h5'.format(ms, ms_direction.ra.deg, ms_direction.dec.deg)
    logging.info('making a new h5parm {}'.format(new_h5parm))
    does_it_exist(new_h5parm, clobber = clobber) # check if the h5parm exists

    # write these best phase solutions to the new h5parm
    h = lh5.h5parm(new_h5parm, readonly = False)
    h.makeSolset(addTables = False) # creates sol000
    # FIXME using 'addTables = False' because the default 'addTables = True' gives
    #       'NotImplementedError: structured arrays with columns with type description ``<U16`` are not supported yet, sorry'
    solset = h.getSolset('sol000')

    # get a h5parm with a result I am going to copy across
    # TODO extend this to loop over all h5parms to be copied
    my_h5parm = mtf_directions[sorted(mtf_directions.keys())[0]]

    # get the station for which the result is valid
    # TODO right now, just taking one station by way of example
    my_station = mtf_stations[0]
    logging.info('copying {} data from {} to {}'.format(my_station, my_h5parm, new_h5parm))

# NB THIS IS WHERE I AM AT THE MOMENT!
    # TODO use the h5parm and the station to get the relevant data
    lo = lh5.h5parm(my_h5parm, readonly = False)
    phase = lo.getSolset('sol000').getSoltab('phase000')
    for s in phase.ant[:]: # stations
        if s == my_station:
            logging.info('{} has {} dimensions, (pol, dir, ant, freq, time): {}'.format(my_h5parm, phase.val.ndim, phase.val.shape))
            polVals = phase.val[:,0,0,0,0] #  polarisation
            dirVals = phase.val[0,:,0,0,0] # direction
            antVals = phase.val[0,0,:,0,0] # station
            freqVals = phase.val[0,0,0,:,0] # frequency
            timeVals = phase.val[0,0,0,0,:] # time
            print(polVals, len(polVals))
            print(dirVals, len(dirVals))
            print(antVals, len(antVals))
            print(freqVals, len(freqVals))
            print(timeVals, len(timeVals))
    lo.close()

    # TODO copy this data into the new h5parm
    # TODO make sure this new h5parm has the same format as the standard lofar h5parms

    # dummy data
    # NOTE having a string here gives 'TypeError: Array objects cannot currently deal with void, unicode or object arrays'
    #      so encoding as ascii
    pol = [ascii('XX'), ascii('YY')]
    dir = [ascii('pointing')]
    ant = [ascii(mtf_station) for mtf_station in mtf_stations]
    freq = [1.3300628662109375E8]
    time = list(range(1686))
    vals = np.zeros((len(pol), len(dir), len(ant), len(freq), len(time)))
    weights = vals
    c = solset.makeSoltab('phase',
                          axesNames = ['pol', 'dir', 'ant', 'freq', 'time'],
                          axesVals = [pol, dir, ant, freq, time],
                          vals = vals,
                          weights = weights) # creates phase000
    h.close()


    logging.info('finished making the h5parm {}'.format(new_h5parm))
    logging.info('make_h5parm(mtf = {}, ms = {}, clobber = {}) completed'.format(mtf, ms, clobber))
    return new_h5parm

def applyh5parm(new_h5parm, ms, clobber = False):
    '''
    description:
    - create ndppp parset
    - apply the output of make_h5parm to the measurement set

    parameters:
    - new_h5parm (str): the output of make_h5parm
    - ms         (str): the measurement set for self-calibration

    returns:
    - ms (str): measurement set for self-calibration with corrected data
    '''

    logging.info('executing applyh5parm(new_h5parm = {}, ms = {}, clobber = {})'.format(new_h5parm, ms, clobber))
    # parset is saved in same directory as the h5parm
    parset = os.path.dirname(new_h5parm) + '/applyh5parm.parset'
    column_in = 'DATA'
    column_out = 'CORRECTED_DATA'

    does_it_exist(parset, clobber = clobber) # if parset already exists, warn user

    # create the parset
    with open(parset, 'w') as f:
        f.write('# applyh5parm function created this parset at {}\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        f.write('msin                = {}\n'.format(ms))
        f.write('msin.datacolumn     = {}\n'.format(column_in))
        f.write('msout               = .\n')
        f.write('msout.datacolumn    = %s\n' % column_out)
        f.write('steps               = [applycal]\n')
        f.write('applycal.type       = applycal\n')
        f.write('applycal.parmdb     = %s\n' % new_h5parm)
        f.write('applycal.correction = phase000\n')
    f.close()

    # apply the h5parm
    logging.info('apply the {} to {} with NDPPP'.format(new_h5parm, ms))
    ndppp_output = subprocess.check_output(['NDPPP', '--help']) # NOTE set to parset

    # format and log the output
    # NOTE might not be the best way of handling this given the text output could be large
    ndppp_output = ndppp_output.decode('utf-8')
    ndppp_output = ndppp_output.split('\n')
    for line in ndppp_output:
        if line: # do not print blank line
            logging.info(line)

    logging.info('finished applying {} to {}'.format(new_h5parm, ms))
    logging.info('applyh5parm(new_h5parm = {}, ms = {}, clobber = {}) completed'.format(new_h5parm, ms, clobber))
    return ms

def updatelist(new_h5parm, loop3_h5parm, mtf, clobber = False):
    '''
    description:
    - combine the phase solutions from the initial h5parm and the final h5parm
    - calls evaluate_solutions to update the master file with a new line appended

    parameters:
    - new_h5parm   (str): the initial h5parm (i.e. from make_h5parm)
    - loop3_h5parm (str): the final h5parm from loop 3
    - mtf          (str): master text file

    returns:
    - combined_h5parm (str): a new h5parm that is a combination of new_h5parm and loop3_h5par
    '''

    logging.info('executing updatelist(new_h5parm = {}, loop3_h5parm = {}, mtf = {}, clobber = {})'.format(new_h5parm, loop3_h5parm, mtf, clobber))
    # create new h5parm
    logging.info('combining phase solutions from {} and {}'.format(new_h5parm, loop3_h5parm))
    combined_h5parm = '{}-{}'.format(os.path.splitext(new_h5parm)[0], os.path.basename(loop3_h5parm))
    logging.info('new combined h5parm is called {}'.format(combined_h5parm))

    # combine two h5parms
    does_it_exist(combined_h5parm, clobber = clobber) # if h5parm already exists, then exit

    # write these best phase solutions to the new h5parm
    h = lh5.h5parm(combined_h5parm, readonly = False)
    h.makeSolset(addTables = False) # creates sol000
    # FIXME using 'addTables = False' because the default 'addTables = True' gives
    #       'NotImplementedError: structured arrays with columns with type description ``<U16`` are not supported yet, sorry'
    solset = h.getSolset('sol000')
    h.close()

    # evaluate the solutions and update the master file
    logging.info('updating {} with the {} solutions'.format(mtf, combined_h5parm))

    logging.info('evaluating the {} solutions'.format(new_h5parm))
    # evaluate_solutions(h5parm, mtf, threshold)
    logging.info('finished updating {}'.format(mtf))
    logging.info('updatelist(new_h5parm = {}, loop3_h5parm = {}, mtf = {}, clobber = {}) completed'.format(new_h5parm, loop3_h5parm, mtf, clobber))

    # One can then merge the two h5parms in a single file (this is needed if you want to interpolate/rescale/copy solutions in time/freq from the cal to the target):
    # https://support.astron.nl/LOFARImagingCookbook/losoto.html
    # H5parm_merge.py -v cal.h5:sol000 tgt.h5:cal000
    # One can create a h5parm also from a single SB
    # H5parm_importer.py -v LXXXXX_3c295.h5 LXXXXX_3c295.MS

    return combined_h5parm

def main():
    logging.basicConfig(format = '\033[1m%(asctime)s \033[31m%(levelname)s \033[00m%(message)s', datefmt = '%Y/%m/%d %H:%M:%S', level = logging.INFO)

    parser = argparse.ArgumentParser(description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-m', '--mtf', required = True, help = 'master text file')
    parser.add_argument('-p', '--h5parm', required = True, help = 'hdf5 file')
    parser.add_argument('-f', '--ms', required = True, help = 'measurement set')
    parser.add_argument('-t', '--threshold', type = float, default = 0.25, help = 'threshold determining the xx-yy statistic goodness')
    parser.add_argument('-c', '--clobber', help = 'overwrite the new h5parm if it exists', action = 'store_true')
    args = parser.parse_args()
    mtf = args.mtf
    h5parm = args.h5parm
    ms = args.ms
    threshold = args.threshold
    clobber = args.clobber

    logging.info('executing main()')
    # loop3() # run loop 3 to generate h5parm
    evaluate_solutions(h5parm, mtf, threshold) # evaluate phase solutions in a h5parm, append to mtf
    new_h5parm = make_h5parm(mtf, ms, clobber = clobber) # create a new h5parm of the best solutions
    applyh5parm(new_h5parm, ms, clobber = clobber) # apply h5parm to ms
    loop3_h5parm = loop3() # run loop 3, returning h5parm
    updatelist(new_h5parm, loop3_h5parm, mtf, clobber = clobber) # combine h5parms and update mtf
    logging.info('main() completed')

if __name__ == '__main__':
    main()
