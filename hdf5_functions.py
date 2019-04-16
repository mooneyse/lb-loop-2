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
import losoto.h5parm as lh5  # on CEP3, "module load losoto/2.0 on CEP3"
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

def interpolate_nan(x_):
    '''Interpolate NaN values using this answer from Stack Overflow:
    https://stackoverflow.com/a/6520696/6386612.'''
    nans, x = np.isnan(x_), lambda z: z.nonzero()[0]
    x_[nans] = np.interp(x(nans), x(~nans), x_[~nans])
    return x_


def coherence_metric(xx, yy):
    '''Calculates the coherence metric by comparing the XX and YY phases.'''
    xx, yy = interpolate_nan(xx), interpolate_nan(yy)
    return np.nanmean(np.gradient(abs(np.unwrap(xx - yy))) ** 2)


def evaluate_solutions_wrapper(h5parm, mtf, solution_tables, threshold):
    '''Executes the evaluate_solutions function for each solution table. If a
    list of solution tables are given, a corresponding list of thresholds must
    also be given.'''

    if type(solution_tables) is str:
        evaluate_solutions(h5parm, mtf, solution_tables, threshold=threshold)

    elif type(solution_tables) is list:
        for solution_table, threshold in zip(solutions_tables, threshold):
            evaluate_solutions(h5parm, mtf, solution_table, threshold=threshold)


def evaluate_solutions(h5parm, mtf, solution_table, threshold=0.25):
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

    print('MTF:', mtf)
    print('MS:', ms)
    print('DIR:', directions)

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

    for h5parm, ra, dec in zip(h5parms, data['ra'], data['dec']):
        mtf_direction = SkyCoord(float(ra), float(dec), unit='deg')
        separation = directions.separation(mtf_direction)
        mtf_directions[separation] = h5parm  # distances from ms to each h5parm

    # read in the stations from the master text file
    with open(mtf) as f:
        mtf_stations = list(csv.reader(f))[0][4:]  # skip h5parm, ra, and dec
        mtf_stations = [x.lstrip() for x in mtf_stations]  # remove first space

    # find the closest h5parm which has an acceptable solution for each station
    parts = {'prefix': os.path.dirname(os.path.dirname(ms)),
             'ra': directions.ra.deg, 'dec': directions.dec.deg}

    working_file = '{prefix}/make_h5parm_{ra}_{dec}.txt'.format(**parts)
    f = open(working_file, 'w')
    successful_stations = []

    for mtf_station in mtf_stations:  # for each station
        for key in sorted(mtf_directions.keys()):  # shortest separation first
            h5parm = mtf_directions[key]
            row = list(h5parms).index(h5parm)  # row in mtf
            value = data[mtf_station][row]  # boolean for h5parm and station

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
    val, weight = [], []

    for my_line in range(len(working_data)):  # one line per station
        my_station = working_data[my_line][0]
        my_h5parm = working_data[my_line][len(working_data[my_line]) - 1]

        # use the station to get the relevant data to be copied from the h5parm
        lo = lh5.h5parm(my_h5parm, readonly=False)  # NB change this to True
        phase = lo.getSolset('sol000').getSoltab('phase000')

        axes_names = phase.getAxesNames()
        if 'dir' not in axes_names:
            axes_names.append('dir')  # add the direction dimension
        print(axes_names)

        for s in range(len(phase.ant[:])):  # stations
            if phase.ant[s] == my_station.strip():
                # copy values and weights
                v = phase.val[:, :, s, :]
                w = phase.weight[:, :, s, :]
                v_expanded = np.expand_dims(v, axis=2)
                w_expanded = np.expand_dims(w, axis=2)
                val.append(v_expanded)
                weight.append(w_expanded)

        # WARNING pol, dir, ant, time, freq should be the same in all h5parms
        # so using the last one in the loop for that information (could be a
        # source of error in future)
        # also getting the antenna and source table from this last h5parm
        if my_line == len(working_data) - 1:
            soltab = lo.getSolset('sol000')
            antenna_soltab = soltab.getAnt()  # dictionary
            source_soltab = soltab.getSou()  # dictionary
            pol = phase.pol[:]
            print('POL:', pol)
            # dir = phase.dir[:]
            ant = phase.ant[:]
            time = phase.time[:]
            freq = phase.freq[:]

        lo.close()

    vals = np.concatenate(val, axis=2)
    weights = np.concatenate(weight, axis=2)

    # write these best phase solutions to the new h5parm
    print('list(vals.shape):', list(vals.shape))  # WARNING this leads to an assert dim == list(vals.shape) AssertionError if the order changes
    c = solset.makeSoltab('phase',
                          axesNames=['time', 'freq', 'ant', 'pol'],  # 'dir',
                          axesVals=[time, freq, ant, pol],  # dir,
                          vals=vals,
                          weights=weights)  # creates phase000

    # copy source and antenna tables into the new h5parm
    source_table = table.obj._f_get_child('source')
    source_table.append(source_soltab.items())  # from dictionary to list
    antenna_table = table.obj._f_get_child('antenna')
    antenna_table.append(antenna_soltab.items())  # from dictionary to list
    h.close()  # close the new h5parm
    return new_h5parm


def apply_h5parm(h5parm, ms, column_out='DATA'):
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


def update_list(new_h5parm, loop3_h5parm, mtf, threshold=0.25):
    '''Combine the phase solutions from the initial h5parm and the final
    h5parm. Calls evaluate_solutions to update the master file with a new line
    appended.

    Args:
    new_h5parm (str): The initial h5parm (i.e. from make_h5parm).
    loop3_h5parm (str): The final h5parm from loop 3.
    mtf (str): Master text file.

    Returns:
    A new h5parm that is a combination of new_h5parm and loop3_h5parm (str).'''

    # get solutions from new_h5parm and loop3_h5parm
    h = lh5.h5parm(new_h5parm)  # from new_h5parm
    phase = h.getSolset('sol000').getSoltab('phase000')
    pol = phase.pol[:]
    dir = phase.dir[:]
    ant = phase.ant[:]
    time = phase.time[:]
    freq = phase.freq[:]

    val_new_h5parm = phase.val[:]
    weight_new_h5parm = phase.weight[:]
    h.close()

    h = lh5.h5parm(loop3_h5parm)  # from loop3_h5parm
    soltab = h.getSolset('sol000')  # NB change to take highest solset
    phase = soltab.getSoltab('phase000')
    antenna_soltab = soltab.getAnt().items()  # dictionary to list
    source_soltab = soltab.getSou().items()  # dictionary to list

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
    combined_h5parm = (os.path.splitext(new_h5parm)[0] + '-' +
                       os.path.basename(loop3_h5parm))

    # write these best phase solutions to the combined_h5parm
    h = lh5.h5parm(combined_h5parm, readonly=False)

    table = h.makeSolset()  # creates sol000

    solset = h.getSolset('sol000')
    c = solset.makeSoltab('phase',
                          axesNames=['pol', 'dir', 'ant', 'freq', 'time'],
                          axesVals=[pol, dir, ant, freq, time],
                          vals=vals,
                          weights=weights)  # creates phase000

    # copy source and antenna tables into the new h5parm
    source_table = table.obj._f_get_child('source')
    source_table.append(source_soltab)
    antenna_table = table.obj._f_get_child('antenna')
    antenna_table.append(antenna_soltab)  # from dictionary to list
    h.close()

    # evaluate the solutions and update the master file
    evaluate_solutions(h5parm=combined_h5parm, mtf=mtf, threshold=threshold)
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
                        help='master text file')

    parser.add_argument('-p',
                        '--h5parm',
                        required=False,
                        type=str,
                        help='hdf5 file')

    parser.add_argument('-f',
                        '--ms',
                        required=False,
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
                        help='source positions (radians; RA DEC RA DEC...)')

    args = parser.parse_args()
    mtf = args.mtf
    h5parm = args.h5parm
    ms = args.ms
    threshold = args.threshold
    cores = args.cores
    directions = args.directions

    evaluate_solutions_wrapper(h5parm='/data020/scratch/sean/letsgetloopy/M1344+5503.ms_02_c0.h5',
                               mtf='/data020/scratch/sean/letsgetloopy/mtf.txt',
                               solution_tables='phase',
                               threshold=0.25)

    # new_h5parms = dir2phasesol_wrapper(mtf='/data020/scratch/sean/letsgetloopy/mtf.txt',
    #                                    ms='/data020/scratch/sean/letsgetloopy/L693725_SB256_uv.ndppp_prep_target',
    #                                    directions=[0.226893, 0.9512044, 0.244346, 0.9686577],
    #                                    cores=cores)
    new_h5parms = dir2phasesol(mtf='/data020/scratch/sean/letsgetloopy/mtf.txt',
                               ms='/data020/scratch/sean/letsgetloopy/L693725_SB256_uv.ndppp_prep_target',
                               directions=[0.226893, 0.9512044])
    new_h5parms = dir2phasesol(mtf='/data020/scratch/sean/letsgetloopy/mtf.txt',
                               ms='/data020/scratch/sean/letsgetloopy/L693725_SB256_uv.ndppp_prep_target',
                               directions=[0.244346, 0.9686577])
    # apply_h5parm(h5parm=new_h5parms[0], ms=ms)  # new_h5parms[0] used as a test

    # loop 3 goes here

    # update_list(new_h5parm=new_h5parms[0], loop3_h5parm=new_h5parms[1],
    #             mtf=mtf, threshold=threshold)  # new_h5parms used as a test

if __name__ == '__main__':
    main()
