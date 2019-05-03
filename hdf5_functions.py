#!/usr/bin/env python2.7

'''A collection of functions for modifying HDF5 files. These functions form
   loop 2 of the LOFAR long-baseline pipeline, which can be found at
   https://github.com/lmorabit/long_baseline_pipeline.'''

from __future__ import print_function
from functools import partial
from multiprocessing import Pool
from pathlib import Path  # on CEP3, "pip install --user pathlib"
from astropy.coordinates import SkyCoord
from losoto.lib_operations import reorderAxes
import losoto.h5parm as lh5  # on CEP3, "module load losoto/2.0"
import astropy.units as u
import pyrap.tables as pt
import numpy as np
import argparse
import csv
import datetime
import os
import subprocess

__author__ = 'Sean Mooney'
__date__ = '01 November 2018'

def make_blank_mtf(mtf):
    '''Create an empty master text file containing all of the LOFAR remote and
    international stations, and ST001.'''
    mtf_header = ('# h5parm, ra, dec, solutions, ST001, RS106HBA, RS205HBA, '
                  'RS208HBA, RS210HBA, RS305HBA, RS306HBA, RS307HBA, RS310HBA,'
                  ' RS404HBA, RS406HBA, RS407HBA, RS409HBA, RS410HBA, '
                  'RS503HBA, RS508HBA, RS509HBA, DE601HBA, DE602HBA, DE603HBA,'
                  ' DE604HBA, DE605HBA, FR606HBA, SE607HBA, UK608HBA, '
                  'DE609HBA, PL610HBA, PL611HBA, PL612HBA, IE613HBA\n')
    if not os.path.isfile(mtf):  # if it does not already exist
        with open(mtf, 'w+') as the_file:
            the_file.write(mtf_header)
    return mtf


def interpolate_nan(x_):
    '''Interpolate NaN values using this answer from Stack Overflow:
    https://stackoverflow.com/a/6520696/6386612. This works even if the first
    or last value is nan or if there are multiple nans. It raises an error if
    all values are nan.'''
    x_ = np.array(x_)
    if np.isnan(x_).all():  # if all values are nan
        raise ValueError('All values in the array are nan, so interpolation is'
                         ' not possible.')
    nans, x = np.isnan(x_), lambda z: z.nonzero()[0]
    x_[nans] = np.interp(x(nans), x(~nans), x_[~nans])
    return x_


def coherence_metric(xx, yy):
    '''Calculates the coherence metric by comparing the XX and YY phases.'''
    xx, yy = interpolate_nan(xx), interpolate_nan(yy)
    return np.nanmean(np.gradient(abs(np.unwrap(xx - yy))) ** 2)


def evaluate_solutions_wrapper(h5parm, mtf, solution_tables='phase',
                               threshold=0.25):
    '''Executes the evaluate_solutions function for each solution table. If a
    list of solution tables are given, a corresponding list of thresholds must
    also be given.'''

    if type(solution_tables) is str:
        evaluate_solutions(h5parm, mtf, solution_tables, threshold=threshold)

    elif type(solution_tables) is list:
        for solution_table, threshold in zip(solutions_tables, threshold):
            evaluate_solutions(h5parm, mtf, solution_table, threshold=threshold)


def evaluate_solutions(h5parm, mtf, solution_table='phase', threshold=0.25):
    '''Get the direction from the h5parm. Evaluate the phase solutions in the
    h5parm for each station using the coherence metric. Determine the validity
    of each coherence metric that was calculated. Append the right ascension,
    declination, and validity to the master text file.

    Args:
    h5parm (str): LOFAR HDF5 parameter file.
    mtf (str): Master text file.
    threshold (float, optional): threshold to determine the goodness of the
        coherence metric.

    Returns:
    The coherence metric for each station. (dict)'''

    h = lh5.h5parm(h5parm)
    solname = h.getSolsetNames()[0]  # set to -1 to use only the last solset
    solset = h.getSolset(solname)
    soltabnames = solset.getSoltabNames()
    soltab = solset.getSoltab(solution_table + '000')
    stations = soltab.ant
    source = solset.getSou()  # dictionary
    direction = np.degrees(np.array(source[list(source.keys())[0]]))  # degrees
    generator = soltab.getValuesIter(returnAxes=['freq', 'time'])
    evaluations, temporary = {}, {}  # evaluations holds the coherence metrics

    for g in generator:
        temporary[g[1]['pol'] + '_' + g[1]['ant']] = np.squeeze(g[0])

    for station in stations:
        xx = temporary['XX_' + station]
        yy = temporary['YY_' + station]
        evaluations[station] = coherence_metric(xx, yy)  # 0 = best

    with open(mtf) as f:
        mtf_stations = list(csv.reader(f))[0][4:]  # get stations from the mtf

    with open(mtf, 'a') as f:
        f.write('{}, {}, {}, {}'.format(h5parm, direction[0], direction[1], solution_table))

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


def dir2phasesol_wrapper(mtf, ms, directions=[], cores=4):
    '''Book-keeping to get the multiprocessing set up and running.'''

    mtf_list, ms_list = [], []
    for i in range(int(len(directions) / 2)):
        mtf_list.append(mtf)
        ms_list.append(ms)

    directions_paired = list(zip(directions[::2], directions[1::2]))
    multiprocessing = list(zip(mtf_list, ms_list, directions_paired))
    pool = Pool(cores)  # specify cores
    new_h5parms = pool.map(dir2phasesol_multiprocessing, multiprocessing)
    return new_h5parms


def dir2phasesol_multiprocessing(args):
    '''Wrapper to parallelize make_h5parm.'''
    mtf, ms, directions = args
    return dir2phasesol(mtf=mtf, ms=ms, directions=directions)


def dir2phasesol(mtf, ms='', directions=[]):
    '''Get the directions of the h5parms from the master text file. Calculate
    the separation between a list of given directions and the h5parm
    directions. For each station, find the h5parm of smallest separation which
    has valid phase solutions. Create a new h5parm. Write these phase solutions
    to this new h5parm.

    Args:
    mtf (str): Master text file with list of h5parms.
    ms (str): Measurement set to be self-calibrated.
    directions (list, default = []): RA, Dec of source (radians).

    Returns:
    The new h5parm to be applied to the measurement set. (str)'''

    # get the direction from the master text file
    # HACK genfromtxt gives empty string for h5parms when names=True is used
    # importing them separately as a work around
    data = np.genfromtxt(mtf, delimiter=',', unpack=True, dtype=float,
                         names=True)
    h5parms = np.genfromtxt(mtf, delimiter=',', unpack=True, dtype=str,
                            usecols=0)

    # calculate the distance betweeen the ms and the h5parm directions
    # there is one entry in mtf_directions for each unique line in the mtf
    directions = SkyCoord(directions[0], directions[1], unit='rad')
    mtf_directions = {}

    if h5parms.size == 1:
        # to handle mtf files with one row which cannot be iterated over
        mtf_direction = SkyCoord(float(data['ra']), float(data['dec']),
                                 unit='deg')
        separation = directions.separation(mtf_direction)
        mtf_directions[separation] = h5parms

    else:
        for h5parm, ra, dec in zip(h5parms, data['ra'], data['dec']):
            mtf_direction = SkyCoord(float(ra), float(dec), unit='deg')
            separation = directions.separation(mtf_direction)
            mtf_directions[separation] = h5parm  # distances from ms to h5parms

    # read in the stations from the master text file
    with open(mtf) as f:
        mtf_stations = list(csv.reader(f))[0][4:]  # skip h5parm, ra, and dec
        mtf_stations = [x.lstrip() for x in mtf_stations]  # remove first space

    # find the closest h5parm which has an acceptable solution for each station
    # a forward slash is added to the ms name in case it does not end in one
    parts = {'prefix': os.path.dirname(os.path.dirname(ms + '/')),
             'ra': directions.ra.deg,
             'dec': directions.dec.deg}

    working_file = '{prefix}/make_h5parm_{ra}_{dec}.txt'.format(**parts)
    f = open(working_file, 'w')
    successful_stations = []

    for mtf_station in mtf_stations:  # for each station
        for key in sorted(mtf_directions.keys()):  # shortest separation first
            h5parm = mtf_directions[key]

            # this try/except block is necessary because otherwise this function
            # crashes when the master text file only has one h5parm in it
            try:
                row = list(h5parms).index(h5parm)  # row in mtf
                value = data[mtf_station][row]  # boolean for h5parm and station

            except:
                row = 0
                value = data[mtf_station]

            if value == 1 and mtf_station not in successful_stations:
                w = '{}\t{}\t{}\t{}\t{}'.format(mtf_station.ljust(8),
                round(key.deg, 6), int(value),
                row, h5parm)
                f.write('{}\n'.format(w))
                successful_stations.append(mtf_station)

    f.close()

    # create a new h5parm
    ms = os.path.splitext(os.path.normpath(ms))[0]
    new_h5parm = '{}_{}_{}.h5'.format(ms, directions.ra.deg,
                                      directions.dec.deg)
    h = lh5.h5parm(new_h5parm, readonly=False)
    table = h.makeSolset()  # creates sol000
    solset = h.getSolset('sol000')  # on the new h5parm

    # get data to be copied from the working file
    working_data = np.genfromtxt(working_file, delimiter='\t', dtype=str)
    working_data = sorted(working_data.tolist())  # stations are alphabetised

    # NOTE working_data is the list of nearest stations with good solutions; if
    #      for a station there is no good solution in any h5parm the new h5parm
    #      will exclude that station

    val, weight = [], []
    time_check, freq_check, pol_check, ant_check, dir_check = [], [], [], [], []

    for my_line in range(len(working_data)):  # one line per station
        my_station = working_data[my_line][0]
        my_h5parm = working_data[my_line][len(working_data[my_line]) - 1]

        # use the station to get the relevant data to be copied from the h5parm
        lo = lh5.h5parm(my_h5parm, readonly=False)  # NB change this to True
        phase = lo.getSolset('sol000').getSoltab('phase000')

        axes_names = phase.getAxesNames()
        values = phase.val
        weights = phase.weight

        if 'dir' not in axes_names:  # add the direction dimension
            axes_names = ['dir'] + axes_names
            values = np.expand_dims(phase.val, 0)
            weights = np.expand_dims(phase.weight, 0)

        reordered_values = reorderAxes(values, axes_names,
                                       ['time', 'freq', 'ant', 'pol', 'dir'])
        reordered_weights = reorderAxes(weights, axes_names,
                                        ['time', 'freq', 'ant', 'pol', 'dir'])

        for s in range(len(phase.ant[:])):  # stations
            if phase.ant[s] == my_station.strip():
                # copy values and weights
                v = reordered_values[:, :, s, :, :]  # time, freq, ant, pol, dir
                w = reordered_weights[:, :, s, :, :]  # same order as v
                v_expanded = np.expand_dims(v, axis=2)
                w_expanded = np.expand_dims(w, axis=2)
                val.append(v_expanded)
                weight.append(w_expanded)


        time = phase.time[:]
        freq = phase.freq[:]
        pol = phase.pol[:]
        ant = phase.ant[:]

        time_check.append(time)
        freq_check.append(freq)
        pol_check.append(pol)
        ant_check.append(ant)
        lo.close()

    # check that every entry in the *_check lists are identical
    # TODO 1 better behaviour is needed here; we want to make sure the
    #      frequency, polarisations, and possibly directions are the same, but
    #      for time, it should use the minimum and maximum times and with the
    #      smallest interval from all the HDF5 files; for the antennas, by
    #      definition it should have a value for them all; then remove the
    #      NotImplementedError
    for my_list in [time_check, freq_check, pol_check, ant_check]:
        check = all(list(_) == list(my_list[0]) for _ in my_list)
        if not check:
            raise NotImplementedError('A new HDF5 file cannot be made from a '
                                      'few other HDF5 files because at least '
                                      'one of the time, frequency, polarisation'
                                      ', antenna, or direction axes do not '
                                      'match.')

    # direction of new h5parm
    dir = [str(directions.ra.rad) + ', ' + str(directions.dec.rad)]
    vals = np.concatenate(val, axis=2)
    weights = np.concatenate(weight, axis=2)

    # write these best phase solutions to the new h5parm
    c = solset.makeSoltab('phase',
                          axesNames=['time', 'freq', 'ant', 'pol', 'dir'],
                          axesVals=[time, freq, ant, pol, dir],
                          vals=vals,
                          weights=weights)  # creates phase000

    # copy source and antenna tables into the new h5parm
    source_soltab = {'POINTING':
                     np.array([directions.ra.rad, directions.dec.rad],
                              dtype='float32')}
    antenna_soltab = {'ST001': np.array([0, 0, 0], dtype='float32'),
                      'RS106HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS205HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS208HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS210HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS305HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS306HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS307HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS310HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS404HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS406HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS407HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS409HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS410HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS503HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS508HBA': np.array([0, 0, 0], dtype='float32'),
                      'RS509HBA': np.array([0, 0, 0], dtype='float32'),
                      'DE601HBA': np.array([0, 0, 0], dtype='float32'),
                      'DE602HBA': np.array([0, 0, 0], dtype='float32'),
                      'DE603HBA': np.array([0, 0, 0], dtype='float32'),
                      'DE604HBA': np.array([0, 0, 0], dtype='float32'),
                      'DE605HBA': np.array([0, 0, 0], dtype='float32'),
                      'FR606HBA': np.array([0, 0, 0], dtype='float32'),
                      'SE607HBA': np.array([0, 0, 0], dtype='float32'),
                      'UK608HBA': np.array([0, 0, 0], dtype='float32'),
                      'DE609HBA': np.array([0, 0, 0], dtype='float32'),
                      'PL610HBA': np.array([0, 0, 0], dtype='float32'),
                      'PL611HBA': np.array([0, 0, 0], dtype='float32'),
                      'PL612HBA': np.array([0, 0, 0], dtype='float32'),
                      'IE613HBA': np.array([0, 0, 0], dtype='float32')}

    source_table = table.obj._f_get_child('source')
    source_table.append(source_soltab.items())  # from dictionary to list
    antenna_table = table.obj._f_get_child('antenna')
    antenna_table.append(antenna_soltab.items())  # from dictionary to list
    h.close()  # close the new h5parm
    os.remove(working_file)
    return new_h5parm


def apply_h5parm(h5parm, ms, column_out='CORRECTED_DATA'):
    '''Creates an NDPPP parset. Applies the output of make_h5parm to the
    measurement set.

    Args:
    new_h5parm (str): The output of make_h5parm.
    ms (str): The measurement set for self-calibration.
    column_out (str, default='DATA'): The column for NDPPP to write to.

    Returns:
    None.'''

    # parset is saved in same directory as the h5parm
    parset = os.path.dirname(h5parm) + '/applyh5parm.parset'
    column_in = 'DATA'
    now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    with open(parset, 'w') as f:  # create the parset
        f.write('# created by apply_h5parm at {}\n'.format(now))
        f.write('msin                = {}\n'.format(ms))
        f.write('msin.datacolumn     = {}\n'.format(column_in))
        f.write('msout               = {}\n'.format(ms))
        f.write('msout.datacolumn    = {}\n'.format(column_out))
        f.write('steps               = [applycal]\n')
        f.write('applycal.type       = applycal\n')
        f.write('applycal.parmdb     = {}\n'.format(h5parm))
        f.write('applycal.correction = phase000\n')
    f.close()

    ndppp_output = subprocess.check_output(['NDPPP', '--help'])  # NOTE update
    os.remove(parset)


def update_list(new_h5parm, loop3_h5parm, mtf, soltab, threshold=0.25):
    '''Combine the phase solutions from the initial h5parm and the final
    h5parm. Calls evaluate_solutions to update the master file with a new line
    appended.

    Args:
    new_h5parm (str): The initial h5parm (i.e. from make_h5parm).
    loop3_h5parm (str): The final h5parm from loop 3.
    mtf (str): Master text file.
    soltab (str): Solution table (e.g. phase or TEC).
    threshold (float, default=0.25): Threshold determining goodness passed to
                                     evaluate_solutions.

    Returns:
    A new h5parm that is a combination of new_h5parm and loop3_h5parm (str).'''

    # get solutions from new_h5parm and loop3_h5parm
    h = lh5.h5parm(new_h5parm)  # from new_h5parm
    phase = h.getSolset('sol000').getSoltab(soltab + '000')
    pol = phase.pol[:]
    try:  #  may not contain a direction dimension
        dir = phase.dir[:]
    except:
        dir = ['0']  # NB this should be the position for the h5parm
    ant = phase.ant[:]
    time = phase.time[:]
    freq = phase.freq[:]

    val_new_h5parm = phase.val[:]
    weight_new_h5parm = phase.weight[:]
    h.close()

    h = lh5.h5parm(loop3_h5parm)  # from loop3_h5parm
    sol000 = h.getSolset('sol000')  # change to take highest solset?
    phase = sol000.getSoltab(soltab + '000')
    antenna_soltab = sol000.getAnt().items()  # dictionary to list
    source_soltab = sol000.getSou().items()  # dictionary to list

    pol = phase.pol[:]
    try:  #  may not contain a direction dimension
        dir = phase.dir[:]
    except:
        dir = ['0']  # NB this should be the position for the h5parm
    ant = phase.ant[:]
    time = phase.time[:]
    freq = phase.freq[:]

    val_loop3_h5parm = phase.val[:]
    weight_loop3_h5parm = phase.weight[:]
    h.close()

    # for comined_h5parm
    # TODO something more complicated needed here
    vals = val_new_h5parm + val_loop3_h5parm

    # complex_number = 1 + 1j
    weights = weight_new_h5parm + weight_loop3_h5parm
    combined_h5parm = (os.path.splitext(new_h5parm)[0] + '-' +
                       os.path.basename(loop3_h5parm))

    # write these best phase solutions to the combined_h5parm
    h = lh5.h5parm(combined_h5parm, readonly=False)

    table = h.makeSolset()  # creates sol000

    solset = h.getSolset('sol000')
    c = solset.makeSoltab('phase',
                          axesNames=['time', 'freq', 'ant', 'pol', 'dir'],
                          axesVals=[time, freq, ant, pol, dir],
                          vals=vals,
                          weights=weights)  # creates phase000

    # copy source and antenna tables into the new h5parm
    source_table = table.obj._f_get_child('source')
    source_table.append(source_soltab)
    antenna_table = table.obj._f_get_child('antenna')
    antenna_table.append(antenna_soltab)  # from dictionary to list
    h.close()

    # evaluate the solutions and update the master file
    evaluate_solutions(h5parm=combined_h5parm, mtf=mtf, threshold=threshold,
                       solution_table=soltab)
    return combined_h5parm


def main():
    '''First, evaluate the h5parm phase solutions. Then for a given direction,
    make a new h5parm of acceptable solutions from the nearest direction for
    each station. Apply the solutions to the measurement set. Run loop 3 to
    image the measurement set in the given direction. Updates the master text
    file with the new best solutions after loop 3 is called.'''

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-m',
                        '--mtf',
                        required=False,
                        type=str,
                        default='/data020/scratch/sean/letsgetloopy/mtf.txt',
                        help='master text file')

    parser.add_argument('-p',
                        '--h5parm',
                        required=False,
                        type=str,
                        default='/data020/scratch/sean/letsgetloopy/M1344+5503.ms_02_c0.h5',
                        help='hdf5 file')

    parser.add_argument('-f',
                        '--ms',
                        required=False,
                        type=str,
                        default='/data020/scratch/sean/letsgetloopy/L693725_SB256_uv.ndppp_prep_target',
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
                        default=[0.226893, 0.9512044, 0.244346, 0.9686577],
                        nargs='+',
                        help='source positions (radians; RA DEC RA DEC...)')

    parser.add_argument('-s',
                        '--soltabs',
                        required=False,
                        type=str,
                        default='phase',
                        help='measurement set')

    args = parser.parse_args()
    mtf = args.mtf
    h5parm = args.h5parm
    ms = args.ms
    threshold = args.threshold
    cores = args.cores
    directions = args.directions
    soltabs = args.soltabs

    make_blank_mtf(mtf=mtf)

    evaluate_solutions_wrapper(h5parm=h5parm,
                               mtf=mtf,
                               solution_tables=soltabs)

    new_h5parms = dir2phasesol_wrapper(mtf=mtf,
                                       ms=ms,
                                       directions=directions,
                                       cores=cores)

    for new_h5parm in new_h5parms:
        apply_h5parm(h5parm=new_h5parm, ms=ms)  # outputs a ms per direction

    # loop 3 goes here

    update_list(new_h5parm=new_h5parms[0], loop3_h5parm=new_h5parms[1],
                mtf=mtf, soltab='phase', threshold=threshold)  # new_h5parms used as a test


if __name__ == '__main__':
    main()
