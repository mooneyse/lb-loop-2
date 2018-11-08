#!/usr/bin/env python2.7

'''A collection of functions for modifying HDF5 files.'''

from __future__ import print_function
from functools import partial
from multiprocessing import Pool
from pathlib import Path  # pip install --user pathlib on CEP3
from astropy.coordinates import SkyCoord
import astropy.units as u
import pyrap.tables as pt
import losoto.h5parm as lh5  # module load losoto/2.0 on CEP3
import numpy as np
import argparse
import csv
import datetime
import logging
import multiprocessing
import os
import subprocess
import sys
import threading
import loop3A


def interpolate_nan(X):
    '''Interpolate NaN values using this answer from Stack Overflow:
    https://stackoverflow.com/a/6520696/6386612.'''
    nans, x = np.isnan(X), lambda z: z.nonzero()[0]
    X[nans] = np.interp(x(nans), x(~nans), X[~nans])
    return X


def coherence_metric(xx, yy):
    '''Calculates the coherence metric by comparing the XX and YY phases.'''
    xx, yy = interpolate_nan(xx), interpolate_nan(yy)
    return np.nanmean(np.gradient(abs(np.unwrap(xx - yy))) ** 2)


def evaluate_solutions(h5parm, mtf, threshold=0.25):
    '''Get the direction from the h5parm. Evaluate the phase solutions in the
    h5parm for each station using the coherence metric. Determine the validity
    of each xx-yy statistic that was calculated. Append the right ascension,
    declination, and validity to the master text file.

    Args:
    h5parm (str): LOFAR HDF5 parameter file.
    mtf (str): Master text file.
    threshold (float, optional): threshold to determine the goodness of the
        coherence metric.

    Returns:
    The coherence metric for each station. (dict)'''

    h = lh5.h5parm(h5parm)
    solname = h.getSolsetNames()[-1]  # only using the last solution set
    solset = h.getSolset(solname)
    soltabnames = solset.getSoltabNames()
    phase = solset.getSoltab('phase000')
    stations = phase.ant
    source = solset.getSou()  # dictionary
    direction = np.degrees(np.array(source[list(source.keys())[0]]))  # degrees
    generator = phase.getValuesIter(returnAxes=['freq', 'time'])
    evaluations, temporary = {}, {}  # evaluations holds the coherence metrics

    for g in generator:
        temporary[g[1]['pol'] + '_' + g[1]['ant']] = np.squeeze(g[0])

    for station in stations:
        xx = temporary['XX_' + station]
        yy = temporary['YY_' + station]
        evaluations[station] = coherence_metric(xx, yy)  # 0 = best

    with open(mtf) as f:
        mtf_stations = list(csv.reader(f))[0][3:]  # get stations from the mtf

    with open(mtf, 'a') as f:
        f.write('{}, {}, {}'.format(h5parm, direction[0], direction[1]))

        for mtf_station in mtf_stations:
            # look up the statistic for a station and determine if it is good
            try:
                value = evaluations[mtf_station[1:]]
            except KeyError:
                value = float('nan')

            if np.isnan(value):
                f.write(', {}'.format('nan'))
            elif value < threshold:  # success
                f.write(', {}'.format(int(True)))
            else:  # fail
                f.write(', {}'.format(int(False)))

        f.write('\n')
    h.close()
    return evaluations

def make_h5parm_multiprocessing(args):
    '''
    description:
    - wrapper for make_h5parm
    - runs make_h5parm on multiple cores in parallel

    parameters:
    - args (list): all arguments to be passed to the make_h5parm function with
                   each list entry being a tuple of parameters for one run

    returns:
    - make_h5parm (function): runs the make_h5parm in parallel
    '''

    # probably over-complicated and not very pythonic way of handling the
    # the different arguments passed
    if len(args) == 4:
        mtf, ms, clobber, directions = args

    elif len(args) == 3 and type(args[2]) == list:
        mtf, ms, directions = args
        clobber = False

    elif len(args) == 3 and type(args[2]) == bool:
        mtf, ms, clobber = args
        directions = []

    elif len(args) == 2:
        mtf, ms = args
        clobber, directions = False, []

    return make_h5parm(mtf, ms = ms, clobber = clobber, directions = directions)

def make_h5parm(mtf, ms = '', clobber = False, directions = []):
    '''
    description:
    - get the direction from the measurement set or list provided
    - get the directions of the h5parms from the master text file
    - calculate the separation between the measurement set direction and the h5parm directions
    - for each station, find the h5parm of smallest separation which has valid phase solutions
    - create a new h5parm
    - write these phase solutions to this new h5parm

    parameters:
    - mtf        (str) : master text file with list of h5parms
    - ms         (str) : measurement set to be self-calibrated
    - clobber    (bool): (default = false) overwrite new_h5parm if it exists
    - directions (list): (default = []) ra, dec of source (radians)

    returns:
    - new_h5parm (str): the new h5parm to be applied to the measurement set
    '''

    logging.info('executing make_h5parm(mtf = {}, ms = {}, clobber = {}, directions = {})'.format(mtf, ms, clobber, directions))

    # get the direction from the measurement set if source position is not given
    if not directions:
        if not ms:
            logging.info('no directions or ms supplied so exiting')
            sys.exit()
        else:
            t  = pt.table(ms, readonly = True, ack = False)
            field = pt.table(t.getkeyword('FIELD'), readonly = True, ack = False)
            directions = field.getcell('PHASE_DIR', 0)[0] # radians
            field.close()
            t.close()
            logging.info('no source positions given, using {} phase center'.format(ms))

    elif len(directions) != 2:
        logging.error('ra, dec not passed correctly')
        os._exit() # https://stackoverflow.com/a/23937527/6386612
        # not ideal as buffers are not flushed on exit
        # sys.exit()

    directions = SkyCoord(directions[0], directions[1], unit = 'rad')

    # get the direction from the master text file
    # HACK genfromtxt gives empty string for h5parms when names = True is used; importing them separately as a work around
    data = np.genfromtxt(mtf, delimiter = ',', unpack = True, dtype = float, names = True)
    h5parms = np.genfromtxt(mtf, delimiter = ',', unpack = True, dtype = str, usecols = 0)
    mtf_directions = {}

    # calculate the distance betweeen the ms direction and the h5parm directions
    # there is one entry in mtf_directions for each unique line in the mtf
    for h5parm, ra, dec in zip(h5parms, data['ra'], data['dec']):
        mtf_direction = SkyCoord(float(ra), float(dec), unit = 'deg')
        separation = directions.separation(mtf_direction)
        mtf_directions[separation] = h5parm # distances from ms to each h5parm

    # read in the stations from the master text file
    with open(mtf) as f:
        mtf_stations = list(csv.reader(f))[0][3:] # skipping h5parm, ra, and dec
        mtf_stations = [x.lstrip() for x in mtf_stations] # remove leading space

    # find the closest h5parm which has an acceptable solution for each station
    working_file = '{}/make_h5parm_{}_{}.txt'.format(os.path.dirname(os.path.dirname(ms)), directions.ra.deg, directions.dec.deg)
    logging.info('creating working file {}'.format(working_file))
    f = open(working_file, 'w')

    logging.info('for the ({}, {}) direction, make a new h5parm consisting of the following:'.format(directions.ra.deg, directions.dec.deg))
    logging.info('\tstation \tseparation\tboolean\trow\t5parm')
    successful_stations = []

    for mtf_station in mtf_stations: # for each station
        for key in sorted(mtf_directions.keys()): # starting with shortest separations
            h5parm = mtf_directions[key]
            row = list(h5parms).index(h5parm) # row in mtf
            value = data[mtf_station][row] # boolean value for h5parm and station
            if value == 1 and mtf_station not in successful_stations:
                working_information = '{}\t{}\t{}\t{}\t{}'.format(mtf_station.ljust(8), round(key.deg, 6), int(value), row, h5parm)
                # print this information and write it to a working file
                logging.info('\t{}'.format(working_information))
                f.write('{}\n'.format(working_information))
                successful_stations.append(mtf_station)
    f.close()

    # create a new h5parm
    ms = os.path.splitext(os.path.normpath(ms))[0]
    new_h5parm = '{}_{}_{}.h5'.format(ms, directions.ra.deg, directions.dec.deg)
    logging.info('creating {}'.format(new_h5parm))
    does_it_exist(new_h5parm, clobber = clobber) # check if the h5parm exists
    h = lh5.h5parm(new_h5parm, readonly = False)

    try:
        table = h.makeSolset() # creates sol000
    except:
        logging.error('could not make antenna and source tables')
        sys.exit()
        # we want the default 'addTables = True' but on my machine that gives
        # 'NotImplementedError: structured arrays with columns with type description ``<U16`` are not supported yet, sorry'

    solset = h.getSolset('sol000') # on the new h5parm

    # get data to be copied from the working file to which the information was written
    working_data = np.genfromtxt(working_file, delimiter = '\t', dtype = str)
    working_data = sorted(working_data.tolist()) # stations in h5parm are alphabetised
    val, weight = [], []

    for my_line in range(len(working_data)): # one line per station
        my_station = working_data[my_line][0]
        my_h5parm = working_data[my_line][len(working_data[my_line]) - 1]
        logging.info('copying {} data from {} to {}'.format(my_station, my_h5parm, new_h5parm))

        # use the station to get the relevant data to be copied from the h5parm
        lo = lh5.h5parm(my_h5parm, readonly = False) # NB try change this to True as it should be True
        phase = lo.getSolset('sol000').getSoltab('phase000')
        logging.info('{} has {} dimensions, (pol, dir, ant, freq, time): {}'.format(my_h5parm, phase.val.ndim, phase.val.shape))

        for s in range(len(phase.ant[:])): # stations
            if phase.ant[s] == my_station.strip():
                # copy values and weights
                v = phase.val[:, :, s, :, :]
                w = phase.weight[:, :, s, :, :]
                v_expanded = np.expand_dims(v, axis = 2)
                w_expanded = np.expand_dims(w, axis = 2)
                val.append(v_expanded)
                weight.append(w_expanded)

        # WARNING pol, dir, ant, time, freq should be the same in all h5parms so
        # using the last one in the loop for that information (could be a source of error in future)
        # also getting the antenna and source table from this last h5parm
        if my_line == len(working_data) - 1:
            soltab = lo.getSolset('sol000')
            antenna_soltab = soltab.getAnt() # dictionary
            source_soltab = soltab.getSou() # dictionary
            pol = phase.pol[:]
            dir = phase.dir[:]
            ant = phase.ant[:]
            time = phase.time[:]
            freq = phase.freq[:]

        lo.close()

    vals = np.concatenate(val, axis = 2)
    weights = np.concatenate(weight, axis = 2)

    # write these best phase solutions to the new h5parm
    c = solset.makeSoltab('phase',
                          axesNames = ['pol', 'dir', 'ant', 'freq', 'time'],
                          axesVals = [pol, dir, ant, freq, time],
                          vals = vals,
                          weights = weights) # creates phase000

    # copy source and antenna tables into the new h5parm
    source_table = table.obj._f_get_child('source')
    source_table.append(source_soltab.items()) # from dictionary to list
    antenna_table = table.obj._f_get_child('antenna')
    antenna_table.append(antenna_soltab.items()) # from dictionary to list

    h.close() # close the new h5parm

    # the end, tidying up
    logging.info('finished making the h5parm {}'.format(new_h5parm))
    logging.info('removing {}'.format(working_file))
    os.remove(working_file)
    logging.info('make_h5parm(mtf = {}, ms = {}, clobber = {}) completed'.format(mtf, ms, clobber))
    return new_h5parm

def apply_h5parm(h5parm, ms, column_out='DATA'):
    '''Creates an NDPPP parset. Applies the output of make_h5parm to the
    measurement set.

    Args:
    new_h5parm (str): The output of make_h5parm.
    ms (str): The measurement set for self-calibration.
    column_out (str, default = 'DATA'): The column for NDPPP to write to.

    Returns:
    The path to the measurement set for self-calibration with corrected data.
    (str)'''

    # parset is saved in same directory as the h5parm
    parset = os.path.dirname(h5parm) + '/applyh5parm.parset'
    column_in = 'DATA'
    now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    with open(parset, 'w') as f:  # create the parset
        f.write('# applyh5parm function created this parset at {}\n'.format(now)
        f.write('msin                = {}\n'.format(ms))
        f.write('msin.datacolumn     = {}\n'.format(column_in))
        f.write('msout               = {}\n'.format(ms))
        f.write('msout.datacolumn    = {}\n'.format(column_out)
        f.write('steps               = [applycal]\n')
        f.write('applycal.type       = applycal\n')
        f.write('applycal.parmdb     = {}\n'.format(h5parm))
        f.write('applycal.correction = phase000\n')
    f.close()

    ndppp_output = subprocess.check_output(['NDPPP', '--help']) # NOTE update
    return ms

def update_list(new_h5parm, loop3_h5parm, mtf, clobber = False, threshold = 0.25):
    '''
    description:
    - combine the phase solutions from the initial h5parm and the final h5parm
    - calls evaluate_solutions to update the master file with a new line appended

    parameters:
    - new_h5parm   (str): the initial h5parm (i.e. from make_h5parm)
    - loop3_h5parm (str): the final h5parm from loop 3
    - mtf          (str): master text file

    returns:
    - combined_h5parm (str): a new h5parm that is a combination of new_h5parm and loop3_h5parm
    '''

    logging.info('executing updatelist(new_h5parm = {}, loop3_h5parm = {}, mtf = {}, clobber = {})'.format(new_h5parm, loop3_h5parm, mtf, clobber))

    # get solutions from new_h5parm and loop3_h5parm

    # from new_h5parm
    h = lh5.h5parm(new_h5parm)
    phase = h.getSolset('sol000').getSoltab('phase000')
    pol = phase.pol[:]
    dir = phase.dir[:]
    ant = phase.ant[:]
    time = phase.time[:]
    freq = phase.freq[:]

    val_new_h5parm = phase.val[:]
    weight_new_h5parm = phase.weight[:]
    h.close()

    # from loop3_h5parm
    h = lh5.h5parm(loop3_h5parm)
    soltab = h.getSolset('sol000') # NB remove hard-coded names
    phase = soltab.getSoltab('phase000')
    antenna_soltab = soltab.getAnt().items() # dictionary to list
    source_soltab = soltab.getSou().items() # dictionary to list

    pol = phase.pol[:]
    dir = phase.dir[:]
    ant = phase.ant[:]
    time = phase.time[:]
    freq = phase.freq[:]

    val_loop3_h5parm = phase.val[:]
    weight_loop3_h5parm = phase.weight[:]
    h.close()

    # for comined_h5parm
    vals = val_new_h5parm + val_loop3_h5parm
    weights = weight_new_h5parm + weight_loop3_h5parm

    # create new h5parm
    logging.info('combining phase solutions from {} and {}'.format(new_h5parm, loop3_h5parm))
    combined_h5parm = '{}-{}'.format(os.path.splitext(new_h5parm)[0], os.path.basename(loop3_h5parm))
    logging.info('new combined h5parm is called {}'.format(combined_h5parm))

    # combine two h5parms
    does_it_exist(combined_h5parm, clobber = clobber) # if h5parm already exists, then exit

    # write these best phase solutions to the combined_h5parm
    h = lh5.h5parm(combined_h5parm, readonly = False)
    try:
        table = h.makeSolset() # creates sol000
    except NotImplementedError:
        logging.error('could not make antenna and source tables')
        sys.exit()
        # on my machine the default 'addTables = True' gives
        # 'NotImplementedError: structured arrays with columns with type description ``<U16`` are not supported yet, sorry'

    solset = h.getSolset('sol000')
    c = solset.makeSoltab('phase',
                          axesNames = ['pol', 'dir', 'ant', 'freq', 'time'],
                          axesVals = [pol, dir, ant, freq, time],
                          vals = vals,
                          weights = weights) # creates phase000

    # copy source and antenna tables into the new h5parm
    source_table = table.obj._f_get_child('source')
    source_table.append(source_soltab)
    antenna_table = table.obj._f_get_child('antenna')
    antenna_table.append(antenna_soltab) # from dictionary to list

    h.close()

    # evaluate the solutions and update the master file
    logging.info('updating {} with the {} solutions'.format(mtf, combined_h5parm))
    evaluate_solutions(combined_h5parm, mtf, threshold = threshold)
    logging.info('updatelist(new_h5parm = {}, loop3_h5parm = {}, mtf = {}, clobber = {}) completed'.format(new_h5parm, loop3_h5parm, mtf, clobber))

    return combined_h5parm

def main():
    '''Evaluates the h5parm phase solutions, makes a new h5parm of the best
    solutions, applies them to the measurement set, and updates the master text
    file with the best solutions after loop 3 is called.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-m',
                        '--mtf',
                        required=True,
                        type=str,
                        help='master text file')

    parser.add_argument('-p',
                        '--h5parm',
                        required=True,
                        type=str,
                        help='hdf5 file')

    parser.add_argument('-f',
                        '--ms',
                        required=True,
                        type=str,
                        help='measurement set')

    parser.add_argument('-t',
                        '--threshold',
                        required=False,
                        type=float,
                        default=0.25,
                        help='threshold for the xx-yy statistic goodness')

    parser.add_argument('-n',
                        '--cores',
                        required=False,
                        type=int,
                        default=4,
                        help='number of cores to use')

    parser.add_argument('-d',
                        '--directions',
                        type=float,
                        default=0,
                        nargs='+',
                        help='source positions (radians; ra1 dec1 ra2 dec2...)')

    args = parser.parse_args()
    mtf = args.mtf
    h5parm = args.h5parm
    ms = args.ms
    threshold = args.threshold
    cores = args.cores
    directions = args.directions

    evaluate_solutions(h5parm=h5parm, mtf=mtf, threshold=threshold)

    # if directions: # if a direction is given
    #     if len(directions) % 2 != 0: # should be even
    #         logging.error('uneven number of ra, dec given for source positions')
    #         sys.exit()
    #
    #     # if multiple ra, dec are given then do multiprocessing
    #     # first, some book-keeping to get things in the right place
    #     mtf_list, ms_list, clobber_list = [], [], []
    #     for i in range(int(len(directions) / 2)):
    #         mtf_list.append(mtf)
    #         ms_list.append(ms)
    #         clobber_list.append(clobber)
    #
    #     directions_paired = list(zip(directions[::2], directions[1::2])) # every second item is ra, dec
    #     multiprocessing = list(zip(mtf_list, ms_list, clobber_list, directions_paired))
    #
    #     pool = Pool(cores) # specify cores
    #     new_h5parms = pool.map(make_h5parm_multiprocessing, multiprocessing)
    #
    # else: # if directions are not given, use the ms phase centre
    #     new_h5parm = make_h5parm(mtf, ms = ms, clobber = clobber, directions = directions)
    #
    # for new_h5parm in new_h5parms:
    #     applyh5parm(new_h5parm, ms, clobber = clobber) # apply h5parm to ms
    #
    apply_h5parm(h5parm=h5parm, ms=ms)  # h5parm should be new_h5parm
    # loop3_h5parm = new_h5parm # for testing
    # updatelist(new_h5parm, loop3_h5parm, mtf, clobber = clobber, threshold = threshold) # combine h5parms and update mtf
    # logging.info('main() completed')

if __name__ == '__main__':
    main()
