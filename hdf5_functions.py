#!/usr/bin/env python2.7

'''A collection of functions for modifying HDF5 files. These functions form
   loop 2 of the LOFAR long-baseline pipeline, which can be found at
   https://github.com/lmorabit/long_baseline_pipeline.'''

from __future__ import print_function
from functools import partial
from multiprocessing import Pool
from pathlib import Path  # on CEP3, "pip install --user pathlib"
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord
from losoto.lib_operations import reorderAxes
import losoto.h5parm as lh5  # on CEP3, "module load losoto"
import astropy.units as u
import pyrap.tables as pt
import numpy as np
import argparse
import csv
import datetime
import os
import subprocess

__author__ = 'Sean Mooney'
__date__ = '01 May 2019'

def make_blank_mtf(mtf):
    '''Create an empty master text file containing all of the LOFAR remote and
    international stations, and ST001.

    Args:
    mtf (str): The master text file to be created.

    Returns:
    The name of the master text file (str).'''

    mtf_header = ('# h5parm, ra, dec, ST001, RS106HBA, RS205HBA, RS208HBA, '
                  'RS210HBA, RS305HBA, RS306HBA, RS307HBA, RS310HBA, RS404HBA,'
                  ' RS406HBA, RS407HBA, RS409HBA, RS410HBA, RS503HBA, '
                  'RS508HBA, RS509HBA, DE601HBA, DE602HBA, DE603HBA, DE604HBA,'
                  ' DE605HBA, FR606HBA, SE607HBA, UK608HBA, DE609HBA, '
                  'PL610HBA, PL611HBA, PL612HBA, IE613HBA\n')
    if not os.path.isfile(mtf):  # if it does not already exist
        with open(mtf, 'w+') as the_file:
            the_file.write(mtf_header)

    return mtf


def interpolate_nan(x_):
    '''Interpolate NaN values using this answer from Stack Overflow:
    https://stackoverflow.com/a/6520696/6386612. This works even if the first
    or last value is nan or if there are multiple nans. It raises an error if
    all values are nan.

    Args:
    x_ (list or NumPy array): Values to interpolate.

    Returns:
    The interpolated values. (NumPy array)'''

    x_ = np.array(x_)
    if np.isnan(x_).all():  # if all values are nan
        raise ValueError('All values in the array are nan, so interpolation is'
                         ' not possible.')
    nans, x = np.isnan(x_), lambda z: z.nonzero()[0]
    x_[nans] = np.interp(x(nans), x(~nans), x_[~nans])

    return x_


def coherence_metric(xx, yy):
    '''Calculates the coherence metric by comparing the XX and YY phases.

    Args:
    xx (list or NumPy array): One set of phase solutions.
    yy (list or NumPy array): The other set of phase solutions.

    Returns:
    The coherence metric. (float)'''

    xx, yy = interpolate_nan(xx), interpolate_nan(yy)

    return np.nanmean(np.gradient(abs(np.unwrap(xx - yy))) ** 2)


def evaluate_solutions(h5parm, mtf, threshold=0.25):
    '''Get the direction from the h5parm. Evaluate the phase solutions in the
    h5parm for each station using the coherence metric. Determine the validity
    of each coherence metric that was calculated. Append the right ascension,
    declination, and validity to the master text file.

    Args:
    h5parm (str): LOFAR HDF5 parameter file.
    mtf (str): Master text file.
    threshold (float; default = 0.25): threshold to determine the goodness of
        the coherence metric.

    Returns:
    The coherence metric for each station. (dict)'''

    h = lh5.h5parm(h5parm)
    solname = h.getSolsetNames()[0]  # set to -1 to use only the last solset
    solset = h.getSolset(solname)
    soltabnames = solset.getSoltabNames()
    soltab = solset.getSoltab('phase000')
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


def dir2phasesol_wrapper(mtf, ms, directions=[], cores=4):
    '''Book-keeping to get the multiprocessing set up and running.

    Args:
    mtf (str): The master text file.
    ms (str): The measurement set.
    directions (list; default = []): Directions in radians, of the form RA1,
        Dec1, RA2, Dec2, etc.
    cores (float; default = 4): Number of cores to use.

    Returns:
    The names of the newly created h5parms in the directions specified. (list)
    '''

    mtf_list, ms_list = [], []
    for i in range(int(len(directions) / 2)):
        mtf_list.append(mtf)
        ms_list.append(ms)

    directions_paired = list(zip(directions[::2], directions[1::2]))
    multiprocessing = list(zip(mtf_list, ms_list, directions_paired))
    pool = Pool(cores)  # specify cores
    new_h5parms = pool.map(dir2phasesol_multiprocessing, multiprocessing)

    return new_h5parms


def interpolate_time(the_array, the_times, new_times):
    '''Given a h5parm array, it will interpolate the values in the time axis
    from whatever it is to a given value.

    Args:
    the_array (NumPy array): The array of values or weights from the h5parm.
    the_times (NumPy array): The 1D array of values along the time axis.
    new_times (NumPy array): The 1D time axis that the values will be mapped
        to.

    Returns:
    The array of values or weights for a h5parm expanded to fit the new time
    axis. (NumPy array)'''

    # get the original data
    time, freq, ant, pol, dir_ = the_array.shape  # axes were reordered earlier

    # make the new array
    interpolated_array = np.ones(shape=(len(new_times), freq, ant, pol, dir_))

    for a in range(ant):  # for one antenna only
        old_x_values = the_array[:, 0, a, 0, 0]  # xx
        old_y_values = the_array[:, 0, a, 1, 0]  # yy

        # calculate the interpolated values
        x1 = interp1d(the_times, old_x_values, kind='nearest', bounds_error=False)
        y1 = interp1d(the_times, old_y_values, kind='nearest', bounds_error=False)

        new_x_values = x1(new_times)
        new_y_values = y1(new_times)

        # assign the interpolated values to the new array
        interpolated_array[:, 0, a, 0, 0] = new_x_values  # new x values
        interpolated_array[:, 0, a, 1, 0] = new_y_values  # new y values

    return interpolated_array


def dir2phasesol_multiprocessing(args):
    '''Wrapper to parallelise make_h5parm.

    Args:
    args (list or tuple): Parameters to be passed to the dir2phasesol
        function.

    Returns:
    The output of the dir2phasesol function, which is the name of a new
        h5parm file. (str)'''

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
    directions (list; default = []): RA, Dec of one source in radians.

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
        mtf_stations = list(csv.reader(f))[0][3:]  # skip h5parm, ra, and dec
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

    # working_data is the list of nearest stations with good solutions; if for
    # a station there is no good solution in any h5parm the new h5parm will
    # exclude that station
    val, weight = [], []
    time_mins, time_maxs, time_intervals = [], [], []
    frequencies = []

    # looping through the h5parms that will be used in the new h5parm to find
    # the shortest time interval of all h5parms being copied
    for my_line in range(len(working_data)):  # one line per station
        my_station = working_data[my_line][0]
        my_h5parm = working_data[my_line][len(working_data[my_line]) - 1]

        # use the station to get the relevant data to be copied from the h5parm
        lo = lh5.h5parm(my_h5parm, readonly=False)  # NB change this to True
        phase = lo.getSolset('sol000').getSoltab('phase000')
        time = phase.time[:]
        time_mins.append(np.min(time))
        time_maxs.append(np.max(time))
        time_intervals.append((np.max(time) - np.min(time)) / (len(time) - 1))
        frequencies.append(phase.freq[:])
        lo.close()

    # properties of the new h5parm
    # the time ranges from the lowest to the highest on the smallest interval
    num_of_steps = 1 + ((np.max(time_maxs) - np.min(time_mins)) /
                        np.min(time_intervals))
    new_time = np.linspace(np.min(time_mins), np.max(time_maxs), num_of_steps)

    # looping through the h5parms to get the solutions for the good stations
    # needed to build the new h5parm
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
                # TODO interpolate
                v_interpolated = interpolate_time(the_array=v_expanded,
                                                  the_times=phase.time[:],
                                                  new_times=new_time)
                w_interpolated = interpolate_time(the_array=w_expanded,
                                                  the_times=phase.time[:],
                                                  new_times=new_time)
                val.append(v_interpolated)
                weight.append(w_interpolated)

        frequencies.append(phase.freq[:])
        lo.close()

    # properties of the new h5parm
    # the time ranges from the lowest to the highest on the smallest interval
    freq = [np.average(frequencies)]  # all items in the list should be equal
    ant = successful_stations  # antennas that will be in the new h5parm
    pol = ['XX', 'YY']  # as standard
    dir = [str(directions.ra.rad) + ', ' + str(directions.dec.rad)]  # given

    vals = np.concatenate(val, axis=2)
    weights = np.concatenate(weight, axis=2)

    # write these best phase solutions to the new h5parm
    c = solset.makeSoltab('phase',
                          axesNames=['time', 'freq', 'ant', 'pol', 'dir'],
                          axesVals=[new_time, freq, ant, pol, dir],
                          vals=vals,
                          weights=weights)  # creates phase000

    # copy source and antenna tables into the new h5parm
    source_soltab = {'POINTING':
                     np.array([directions.ra.rad, directions.dec.rad],
                              dtype='float32')}
    # the X, Y, Z coordinates of the stations sould be in these arrays
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
    new_h5parm (str): The output of dir2phasesol.
    ms (str): The measurement set for self-calibration.
    column_out (str; default = 'CORRECTED_DATA'): The column NDPPP writes to.

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


def add_amplitude_and_phase_solutions(ampltides, amplitude_phases, phases):
    '''Convert amplitude and phase solutions into complex numbers, add them,
    and return the amplitude and phase components of the result. The solutions
    must be on the same time axis.

    Args:
    amplitudes (list or NumPy array): Amplitude solutions.
    amplitude_phases (list or NumPy array): Phase solutions from the amplitude
        solve.
    phases (list or NumPy array): Phase solutions.

    Returns:
    Amplitude solutions (NumPy array), phase solutions (NumPy array)'''

    amplitude_final, phase_final = [], []

    for A, theta_A, theta in zip(ampltides, amplitude_phases, phases):
        complex_amplitude = A * complex(np.cos(theta_A), np.sin(theta_A))
        complex_phase = complex(np.cos(theta), np.sin(theta))  # A is unity
        complex_ = complex_amplitude + complex_phase

        amplitude_final.append(abs(complex_))
        phase_final.append(np.arctan2(complex_.imag, complex_.real))

    return np.array(amplitude_final), np.array(phase_final)


def make_new_times(time1, time2):
    '''Make a new time axis from two others, going from the minimum to the
    maximum with the smallest time step.

    Args:
    time1 (list or NumPy array): Times.
    time2 (list or NumPy array): Other times.

    Returns:
    New time axis (list).'''

    times = [time1, time2]
    time_intervals = []
    for time in times:
        time_intervals.append((np.max(time) - np.min(time)) / (len(time) - 1))

    max_time = np.max([np.max(time1), np.max(time2)])
    min_time = np.min([np.min(time1), np.min(time2)])
    num_of_steps = 1 + (max_time - min_time) / np.min(time_intervals)
    new_time = np.linspace(min_time, max_time, num_of_steps)

    return new_time


def sort_axes(soltab):
    '''Add a direction axis if there is none and sort the axes
    into a predefined order.

    Args:
    soltab (Losoto object): Solution table.

    Returns:
    Values ordered, with a direction axis included (NumPy array.)
    Weights ordered, with a direction axis included (NumPy array.)'''

    axes_names = soltab.getAxesNames()
    if 'dir' not in axes_names:  # add the direction dimension
            axes_names = ['dir'] + axes_names
            values = np.expand_dims(soltab.val, 0)
            weights = np.expand_dims(soltab.weight, 0)

    reordered_values = reorderAxes(values, axes_names,
                                   ['time', 'freq', 'ant','pol', 'dir'])
    reordered_weights = reorderAxes(weights, axes_names,
                                    ['time', 'freq', 'ant','pol', 'dir'])

    return reordered_values, reordered_weights


def update_list(initial_h5parm, incremental_h5parm, mtf, threshold=0.25):
    '''Combine the phase solutions from the initial h5parm and the final
    h5parm. The initial h5parm contains the initial solutions and the final
    h5parm contains the incremental solutions so they need to be added to form
    the final solutions. Calls evaluate_solutions to update the master file
    with a new line appended.

    Args:
    new_h5parm (str): The initial h5parm (i.e. from dir2phasesol).
    loop3_h5parm (str): The final h5parm from loop 3.
    mtf (str): Master text file.
    threshold (float; default = 0.25): Threshold determining goodness passed to
        evaluate_solutions.

    Returns:
    A new h5parm that is a combination of new_h5parm and loop3_h5parm (str).'''

    # get solutions from new_h5parm and loop3_h5parm
    f = lh5.h5parm(initial_h5parm)  # from new_h5parm
    initial_phase = f.getSolset('sol000').getSoltab('phase000')
    try:  # h5parms from dir2phasesol have a direction, but in case not
        initial_dir = initial_phase.dir[:]
    except:
        initial_dir = ['0']  # if it is missing

    initial_time = initial_phase.time[:]
    initial_freq = initial_phase.freq[:]
    initial_val = initial_phase.val[:]
    initial_weight = initial_phase.weight[:]

    g = lh5.h5parm(incremental_h5parm)  # from loop3_h5parm
    sol000 = g.getSolset('sol000')  # change to take highest solset?
    incremental_phase = g.getSolset('sol000').getSoltab('phase000')
    antenna_soltab = g.getSolset('sol000').getAnt().items()  # dict to list
    source_soltab = g.getSolset('sol000').getSou().items()  # dict to list

    try:  #  may not contain a direction dimension
        dir = incremental_phase.dir[:]
    except:
        dir = initial_dir  # if none, take it from the other h5
    ant = incremental_phase.ant[:]
    incremental_time = incremental_phase.time[:]
    incremental_freq = incremental_phase.freq[:]
    incremental_val = incremental_phase.val[:]
    incremental_weight = incremental_phase.weight[:]

    # for comined_h5parm
    # make val_initial and val_incremental on the same time axis
    # first, build the new time axis and order the array
    new_times = make_new_times(initial_time, incremental_time)
    initial_sorted_val, initial_sorted_weight = sort_axes(initial_phase)
    incremental_sorted_val, incremental_sorted_weight = sort_axes(incremental_phase)

    # interpolate the solutions from both h5parms onto this new time axis
    initial_val_new = interpolate_time(initial_sorted_val, initial_time, new_times)
    initial_weight_new = interpolate_time(initial_sorted_weight, initial_time, new_times)
    incremental_val_new = interpolate_time(incremental_sorted_val, incremental_time, new_times)
    incremental_weight_new = interpolate_time(incremental_sorted_weight, incremental_time, new_times)

    # BEFORE ADDING WE HAVE TO MAKE SURE OF THE FOLLOWING. THE ANTENNAS ARE IN THE SAME ORDER,
    # OR MAKE THE NEW ARRAY OF VALUES HAVE ALL ANTENNAS AND WRITE NAN FOR ANTENNAS WHERE THERE IS NO
    # SOLUTION IN EITHER H5PARM.
    # FIRST, CHECK IF THE H5PARMS HAVE AMPLITUDE SOLUTION TABLES AND IF SO, GET THE AMPLITUDE AND
    # PHASES AND MAKE AN ARRAY OF VALUES WHICH ARE THE COMPLEX NUMBERS FOR EACH ARRAY, AND THEN
    # ADD THESE TOGETHER

    # add_amplitude_and_phase_solutions(ampltides, amplitude_phases, phases)
    vals = val_initial + val_incremental
    weights = val_initial + weight_incremental
    combined_h5parm = (os.path.splitext(initial_h5parm)[0] + '-' +
                       os.path.basename(incremental_h5parm))

    freq = np.array([np.mean([initial_freq, incremental_freq])])
    pol = np.array(['XX', 'YY'])

    # write these best phase solutions to the combined_h5parm
    h = lh5.h5parm(combined_h5parm, readonly=False)

    table = h.makeSolset()  # creates sol000

    solset = h.getSolset('sol000')
    c = solset.makeSoltab('phase',
                          axesNames=['time', 'freq', 'ant', 'pol', 'dir'],
                          axesVals=[new_times, freq, ant, pol, dir],
                          vals=vals,
                          weights=weights)  # creates phase000

    # copy source and antenna tables into the new h5parm
    source_table = table.obj._f_get_child('source')
    source_table.append(source_soltab)
    antenna_table = table.obj._f_get_child('antenna')
    antenna_table.append(antenna_soltab)  # from dictionary to list

    f.close()
    g.close()
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
                        default='/data020/scratch/sean/letsgetloopy/mtf.txt',
                        help='master text file')

    parser.add_argument('-p',
                        '--h5parm',
                        required=False,
                        type=str,
                        default='/data020/scratch/sean/letsgetloopy/phases.h5',
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

    args = parser.parse_args()
    mtf = args.mtf
    h5parm = args.h5parm
    ms = args.ms
    threshold = args.threshold
    cores = args.cores
    directions = args.directions

    # make_blank_mtf(mtf=mtf)
    #
    # evaluate_solutions(h5parm=h5parm, mtf=mtf)
    #
    # new_h5parms = dir2phasesol_wrapper(mtf=mtf,
    #                                    ms=ms,
    #                                    directions=directions,
    #                                    cores=cores)
    #
    # for new_h5parm in new_h5parms:
    #     apply_h5parm(h5parm=new_h5parm, ms=ms)  # outputs a ms per direction

    # loop 3 goes here

    # new_h5parms used as a test
    update_list(initial_h5parm='/data020/scratch/sean/letsgetloopy/phases.h5', incremental_h5parm='/data020/scratch/sean/letsgetloopy/amplitudes.h5',
                mtf=mtf, threshold=threshold)


if __name__ == '__main__':
    main()
