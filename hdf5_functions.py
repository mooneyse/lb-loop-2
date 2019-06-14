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
import uuid

__author__ = 'Sean Mooney'
__date__ = '01 May 2019'

def make_blank_mtf(mtf):
    '''Create an empty master text file containing all of the LOFAR core,
    remote, and international stations, and ST001.

    Args:
    mtf (str): The master text file to be created.

    Returns:
    The name of the master text file (str).'''

    mtf_header = ('# h5parm, ra, dec, ST001, RS106HBA, RS205HBA, RS208HBA, '
                  'RS210HBA, RS305HBA, RS306HBA, RS307HBA, RS310HBA, RS404HBA,'
                  ' RS406HBA, RS407HBA, RS409HBA, RS410HBA, RS503HBA, '
                  'RS508HBA, RS509HBA, DE601HBA, DE602HBA, DE603HBA, DE604HBA,'
                  ' DE605HBA, FR606HBA, SE607HBA, UK608HBA, DE609HBA, '
                  'PL610HBA, PL611HBA, PL612HBA, IE613HBA, '
                  'CS001HBA0, CS001HBA1, CS002HBA0, CS002HBA1, '
                  'CS003HBA0, CS003HBA1, CS004HBA0, CS004HBA1, '
                  'CS005HBA0, CS005HBA1, CS006HBA0, CS006HBA1, '
                  'CS007HBA0, CS007HBA1, CS011HBA0, CS011HBA1, '
                  'CS013HBA0, CS013HBA1, CS017HBA0, CS017HBA1, '
                  'CS021HBA0, CS021HBA1, CS024HBA0, CS024HBA1, '
                  'CS026HBA0, CS026HBA1, CS028HBA0, CS028HBA1, '
                  'CS030HBA0, CS030HBA1, CS031HBA0, CS031HBA1, '
                  'CS032HBA0, CS032HBA1, CS101HBA0, CS101HBA1, '
                  'CS103HBA0, CS103HBA1, CS201HBA0, CS201HBA1, '
                  'CS301HBA0, CS301HBA1, CS302HBA0, CS302HBA1, '
                  'CS401HBA0, CS401HBA1, CS501HBA0, CS501HBA1\n')
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

    try:
        xx, yy = interpolate_nan(xx), interpolate_nan(yy)
    except:
        return np.nan  # if the values are all nan, they cannot be interpolated
                       # so return a coherence value of nan also

    return np.nanmean(np.gradient(abs(np.unwrap(xx - yy))) ** 2)


def evaluate_solutions(h5parm, mtf, threshold=0.25):
    '''Get the direction from the h5parm. Evaluate the phase solutions in the
    h5parm for each station using the coherence metric. Determine the validity
    of each coherence metric that was calculated. Append the right ascension,
    declination, and validity to the master text file.

    Args:
    h5parm (str): LOFAR HDF5 parameter file.
    mtf (str): Master text file.
    threshold (float; default=0.25): threshold to determine the goodness of
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
        # TODO the first error arises here, because the h5parm has three
        #      axes and it was assumed that it would have only one, so for now
        #      as a HACK to get things moving, we'll take just the first axis,
        #      without which we get the following traceback:
        #      Traceback (most recent call last):
        #      File "./hdf5_functions.py", line 1057, in <module>
        #      main()
        #      File "./hdf5_functions.py", line 1034, in main
        #      evaluate_solutions(h5parm=h5parm0, mtf=mtf)
        #      File "./hdf5_functions.py", line 136, in evaluate_solutions
        #      evaluations[station] = coherence_metric(xx, yy)  # 0 = best
        #      File "./hdf5_functions.py", line 100, in coherence_metric
        #      return np.nanmean(np.gradient(abs(np.unwrap(xx - yy))) ** 2)
        #      TypeError: unsupported operand type(s) for ** or pow(): 'list' and 'int'
        evaluations[station] = coherence_metric(xx[:, 0], yy[:, 0])  # 0 = best

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


def build_soltab(soltab, working_data):
    '''Creates a solution table from many h5parms using data from the temporary
    working file.

    Args:
    soltab (str): The name of the solution table to copy solutions from.
    working_data (NumPy array): Data providing the list of good and bad
        stations, which was taken from the temporary working file, and the
        goodness relates to the coherence metric on the phase solutions.

    Returns:
    Values to populate the solution table (NumPy array).
    Weights to populate the solution table (NumPy array).
    Time axis to populate the solution table (NumPy array).
    Frequency axis to populate the solution table (NumPy array).
    Antenna axis to populate the solution table (NumPy array).'''

    for my_line in range(len(working_data)):  # one line per station
        my_station = working_data[my_line][0]
        my_h5parm = working_data[my_line][len(working_data[my_line]) - 1]
        lo = lh5.h5parm(my_h5parm, readonly=False)
        tab = lo.getSolset('sol000').getSoltab(soltab + '000')
        time_mins.append(np.min(tab.time[:]))
        time_maxs.append(np.max(tab.time[:]))
        time_intervals.append((np.max(tab.time[:]) - np.min(tab.time[:])) / (len(tab.time[:]) - 1))
        frequencies.append(tab.freq[:])
        lo.close()

    # the time ranges from the lowest to the highest on the smallest interval
    num_of_steps = 1 + ((np.max(time_maxs) - np.min(time_mins)) / np.min(time_intervals))
    new_time = np.linspace(np.min(time_mins), np.max(time_maxs), num_of_steps)

    # looping through the h5parms to get the solutions for the good stations
    for my_line in range(len(working_data)):  # one line per station
        my_station = working_data[my_line][0]
        my_h5parm = working_data[my_line][len(working_data[my_line]) - 1]
        lo = lh5.h5parm(my_h5parm, readonly=False)
        tab = lo.getSolset('sol000').getSoltab(soltab + '000')
        axes_names = tab.getAxesNames()
        values = tab.val
        weights = tab.weight

        if 'dir' not in axes_names:  # add the direction dimension
            axes_names = ['dir'] + axes_names
            values = np.expand_dims(tab.val, 0)
            weights = np.expand_dims(tab.weight, 0)

        reordered_values = reorderAxes(values, axes_names, ['time', 'freq', 'ant', 'pol', 'dir'])
        reordered_weights = reorderAxes(weights, axes_names, ['time', 'freq', 'ant', 'pol', 'dir'])

        for s in range(len(tab.ant[:])):  # stations
            if tab.ant[s] == my_station.strip():
                v = reordered_values[:, :, s, :, :]  # time, freq, ant, pol, dir
                w = reordered_weights[:, :, s, :, :]
                v_expanded = np.expand_dims(v, axis=2)
                w_expanded = np.expand_dims(w, axis=2)
                v_interpolated = interpolate_time(the_array=v_expanded, the_times=tab.time[:], new_times=new_time)
                w_interpolated = interpolate_time(the_array=w_expanded, the_times=tab.time[:], new_times=new_time)
                val.append(v_interpolated)
                weight.append(w_interpolated)
        lo.close()

    vals = np.concatenate(val, axis=2)
    weights = np.concatenate(weight, axis=2)

    return vals, weights, new_time, [np.average(frequencies)], successful_stations


def dir2phasesol(mtf, ms='', directions=[]):
    '''Get the directions of the h5parms from the master text file. Calculate
    the separation between a list of given directions and the h5parm
    directions. For each station, find the h5parm of smallest separation which
    has valid phase solutions. Create a new h5parm. Write these phase solutions
    to this new h5parm.

    Args:
    mtf (str): Master text file with list of h5parms.
    ms (str; default=''): Measurement set to be self-calibrated.
    directions (list; default=[]): Right ascension and declination of one
        source in radians.

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
            # TODO if there are multiple h5parms for one direction (which will
            #      be the case after the ms has been through loop 3 and the
            #      update_list function) then the best solutions will be at the
            #      bottom; however, here they will still both have the same
            #      direction and could be good solutions, getting a 1 in the
            #      mtf, but one solution table is better than another in reality
            #      and we need a way to distinguish this - it could involve
            #      changing the mtf to write the XX-YY statistic if the value is
            #      above the threshold, and zero otherwise, and checking to see
            #      which hdf5 has a higher value if it is the case that both
            #      directions are the same (or better, if the separation is
            #      equal, to cover this happening by chance with another hdf5
            #      with an equal separation but in a different direction)

            h5parm = mtf_directions[key]

            # this try/except block is necessary because otherwise this crashes
            # when the master text file only has one h5parm in it
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
    # the shortest time interval of all h5parms being copied, and the longest
    # time span
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

    # the time ranges from the lowest to the highest on the smallest interval
    num_of_steps = 1 + ((np.max(time_maxs) - np.min(time_mins)) /
                        np.min(time_intervals))
    new_time = np.linspace(np.min(time_mins), np.max(time_maxs), num_of_steps)
    stations_in_correct_order = []

    # looping through the h5parms again to get the solutions for the good
    # stations needed to build the new h5parm
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
                stations_in_correct_order.append(phase.ant[s])
                # copy values and weights
                v = reordered_values[:, :, s, :, :]  # time, freq, ant, pol, dir
                w = reordered_weights[:, :, s, :, :]  # same order as v
                v_expanded = np.expand_dims(v, axis=2)
                w_expanded = np.expand_dims(w, axis=2)
                v_interpolated = interpolate_time(the_array=v_expanded,
                                                  the_times=phase.time[:],
                                                  new_times=new_time)
                w_interpolated = interpolate_time(the_array=w_expanded,
                                                  the_times=phase.time[:],
                                                  new_times=new_time)
                val.append(v_interpolated)
                weight.append(w_interpolated)

        lo.close()

    # properties of the new h5parm
    freq = [np.average(frequencies)]  # all items in the list should be equal
    ant = stations_in_correct_order  # antennas that will be in the new h5parm
    pol = ['XX', 'YY']  # as standard
    dir_ = [str(directions.ra.rad) + ', ' + str(directions.dec.rad)]  # given

    vals = np.concatenate(val, axis=2)
    print('val shape', vals.shape, vals[0,0,23,0,0], ant[23])
    weights = np.concatenate(weight, axis=2)
    # TODO the HACK on the line below is necessary to get around the fact that
    #      there are three frequencies
    freq = [124900817.87109375, 131932067.87109375, 138963317.87109375]

    # write these best phase solutions to the new h5parm
    c = solset.makeSoltab('phase',
                          axesNames=['time', 'freq', 'ant', 'pol', 'dir'],
                          axesVals=[new_time, freq, ant, pol, dir_],
                          vals=vals,
                          weights=weights)  # creates phase000

    # WARNING the tec and amplitude soltab functionality has not been tested
    try:
        vals, weights, time, freq, ant = build_soltab(soltab='tec', working_data=working_data)
        c = solset.makeSoltab('tec',
                              axesNames=['time', 'freq', 'ant', 'pol', 'dir'],
                              axesVals=[time, freq, ant, pol, dir_],
                              vals=vals,
                              weights=weights)  # creates tec000
    except:
        pass

    try:
        vals, weights, time, freq = build_soltab(soltab='amplitude', working_data=working_data)
        c = solset.makeSoltab('amplitude',
                              axesNames=['time', 'freq', 'ant', 'pol', 'dir'],
                              axesVals=[time, freq, ant, pol, dir_],
                              vals=vals,
                              weights=weights)  # creates amplitude000
    except:
        pass

    # copy source and antenna tables into the new h5parm
    source_soltab = {'POINTING':
                     np.array([directions.ra.rad, directions.dec.rad],
                              dtype='float32')}
    # the x, y, z coordinates of the stations should be in these arrays
    tied = {'ST001': np.array([3826557.5, 461029.06, 5064908], dtype='float32')}

    core = {'CS001HBA0': np.array([3826896.235, 460979.455, 5064658.203], dtype='float32'),
            'CS001HBA1': np.array([3826979.384, 460897.597, 5064603.189], dtype='float32'),
            'CS002HBA0': np.array([3826600.961, 460953.402, 5064881.136], dtype='float32'),
            'CS002HBA1': np.array([3826565.594, 460958.110, 5064907.258], dtype='float32'),
            'CS003HBA0': np.array([3826471.348, 461000.138, 5064974.201], dtype='float32'),
            'CS003HBA1': np.array([3826517.812, 461035.258, 5064936.15], dtype='float32'),
            'CS004HBA0': np.array([3826585.626, 460865.844, 5064900.561], dtype='float32'),
            'CS004HBA1': np.array([3826579.486, 460917.48, 5064900.502], dtype='float32'),
            'CS005HBA0': np.array([3826701.16, 460989.25, 5064802.685], dtype='float32'),
            'CS005HBA1': np.array([3826631.194, 461021.815, 5064852.259], dtype='float32'),
            'CS006HBA0': np.array([3826653.783, 461136.440, 5064824.943], dtype='float32'),
            'CS006HBA1': np.array([3826612.499, 461080.298, 5064861.006], dtype='float32'),
            'CS007HBA0': np.array([3826478.715, 461083.720, 5064961.117], dtype='float32'),
            'CS007HBA1': np.array([3826538.021, 461169.731, 5064908.827], dtype='float32'),
            'CS011HBA0': np.array([3826637.421, 461227.345, 5064829.134], dtype='float32'),
            'CS011HBA1': np.array([3826648.961, 461354.241, 5064809.003], dtype='float32'),
            'CS013HBA0': np.array([3826318.954, 460856.125, 5065101.85], dtype='float32'),
            'CS013HBA1': np.array([3826402.103, 460774.267, 5065046.836], dtype='float32'),
            'CS017HBA0': np.array([3826405.095, 461507.460, 5064978.083], dtype='float32'),
            'CS017HBA1': np.array([3826499.783, 461552.498, 5064902.938], dtype='float32'),
            'CS021HBA0': np.array([3826463.502, 460533.094, 5065022.614], dtype='float32'),
            'CS021HBA1': np.array([3826368.813, 460488.057, 5065097.759], dtype='float32'),
            'CS024HBA0': np.array([3827218.193, 461403.898, 5064378.79], dtype='float32'),
            'CS024HBA1': np.array([3827123.504, 461358.861, 5064453.935], dtype='float32'),
            'CS026HBA0': np.array([3826418.227, 461805.837, 5064941.199], dtype='float32'),
            'CS026HBA1': np.array([3826335.078, 461887.696, 5064996.213], dtype='float32'),
            'CS028HBA0': np.array([3825573.134, 461324.607, 5065619.039], dtype='float32'),
            'CS028HBA1': np.array([3825656.283, 461242.749, 5065564.025], dtype='float32'),
            'CS030HBA0': np.array([3826041.577, 460323.374, 5065357.614], dtype='float32'),
            'CS030HBA1': np.array([3825958.428, 460405.233, 5065412.628], dtype='float32'),
            'CS031HBA0': np.array([3826383.037, 460279.343, 5065105.85], dtype='float32'),
            'CS031HBA1': np.array([3826477.725, 460324.381, 5065030.705], dtype='float32'),
            'CS032HBA0': np.array([3826864.262, 460451.924, 5064730.006], dtype='float32'),
            'CS032HBA1': np.array([3826947.411, 460370.066, 5064674.992], dtype='float32'),
            'CS101HBA0': np.array([3825899.977, 461698.906, 5065339.205], dtype='float32'),
            'CS101HBA1': np.array([3825805.288, 461653.869, 5065414.35], dtype='float32'),
            'CS103HBA0': np.array([3826331.59, 462759.074, 5064919.62], dtype='float32'),
            'CS103HBA1': np.array([3826248.441, 462840.933, 5064974.634], dtype='float32'),
            'CS201HBA0': np.array([3826679.281, 461855.243, 5064741.38], dtype='float32'),
            'CS201HBA1': np.array([3826690.821, 461982.139, 5064721.249], dtype='float32'),
            'CS301HBA0': np.array([3827442.564, 461050.814, 5064242.391], dtype='float32'),
            'CS301HBA1': np.array([3827431.025, 460923.919, 5064262.521], dtype='float32'),
            'CS302HBA0': np.array([3827973.226, 459728.624, 5063975.3], dtype='float32'),
            'CS302HBA1': np.array([3827890.077, 459810.483, 5064030.313], dtype='float32'),
            'CS401HBA0': np.array([3826795.752, 460158.894, 5064808.929], dtype='float32'),
            'CS401HBA1': np.array([3826784.211, 460031.993, 5064829.062], dtype='float32'),
            'CS501HBA0': np.array([3825568.82, 460647.62, 5065683.028], dtype='float32'),
            'CS501HBA1': np.array([3825663.508, 460692.658, 5065607.883], dtype='float32')}

    antenna_soltab = {'RS106HBA': np.array([3829205.598, 469142.533000, 5062181.002], dtype='float32'),
                      'RS205HBA': np.array([3831479.67, 463487.529000, 5060989.903], dtype='float32'),
                      'RS208HBA': np.array([3847753.31, 466962.809000, 5048397.244], dtype='float32'),
                      'RS210HBA': np.array([3877827.56186, 467536.604956, 5025445.584], dtype='float32'),
                      'RS305HBA': np.array([3828732.721, 454692.421000, 5063850.334], dtype='float32'),
                      'RS306HBA': np.array([3829771.249, 452761.702000, 5063243.181], dtype='float32'),
                      'RS307HBA': np.array([3837964.52, 449627.261000, 5057357.585], dtype='float32'),
                      'RS310HBA': np.array([3845376.29, 413616.564000, 5054796.341], dtype='float32'),
                      'RS404HBA': np.array([0.0, 0.0, 0.0], dtype='float32'),
                      'RS406HBA': np.array([3818424.939, 452020.269000, 5071817.644], dtype='float32'),
                      'RS407HBA': np.array([3811649.455, 453459.894000, 5076728.952], dtype='float32'),
                      'RS409HBA': np.array([3824812.621, 426130.330000, 5069251.754], dtype='float32'),
                      'RS410HBA': np.array([0.0, 0.0, 0.0], dtype='float32'),
                      'RS503HBA': np.array([3824138.566, 459476.972, 5066858.578], dtype='float32'),
                      'RS508HBA': np.array([3797136.484, 463114.447, 5086651.286], dtype='float32'),
                      'RS509HBA': np.array([3783537.525, 450130.064, 5097866.146], dtype='float32'),
                      'DE601HBA': np.array([4034101.522, 487012.757, 4900230.499], dtype='float32'),
                      'DE602HBA': np.array([4152568.006, 828789.153, 4754362.203], dtype='float32'),
                      'DE603HBA': np.array([3940295.706, 816722.865, 4932394.416], dtype='float32'),
                      'DE604HBA': np.array([3796379.823, 877614.13, 5032712.528], dtype='float32'),
                      'DE605HBA': np.array([4005681.02, 450968.643, 4926458.211], dtype='float32'),
                      'FR606HBA': np.array([4324016.708, 165545.525, 4670271.363], dtype='float32'),
                      'SE607HBA': np.array([3370271.657, 712125.881, 5349991.165], dtype='float32'),
                      'UK608HBA': np.array([4008461.941,-100376.609, 4943716.874], dtype='float32'),
                      'DE609HBA': np.array([3727217.673, 655109.175, 5117003.123], dtype='float32'),
                      'PL610HBA': np.array([3738462.416, 1148244.316, 5021710.658], dtype='float32'),
                      'PL611HBA': np.array([3850980.881, 1438994.879, 4860498.993], dtype='float32'),
                      'PL612HBA': np.array([3551481.817, 1334203.573, 5110157.41], dtype='float32'),
                      'IE613HBA': np.array([3801692.0, -528983.94, 5076958.0], dtype='float32')}

    # delete a key, value pair from the antenna table if it does not exist in
    # the antenna axis
    keys_to_remove = []
    for key in antenna_soltab:
        if key not in ant:
            keys_to_remove.append(key)

    for k in keys_to_remove:
        antenna_soltab.pop(k, None)

    for a in ant:
        if a[:2] == 'ST':
            antenna_soltab.update(tied)  # there will only be the tied station
        if a[:2] == 'CS':
            antenna_soltab.update(core)
            break  # only add the core stations to the antenna table once

    source_table = table.obj._f_get_child('source')
    source_table.append(source_soltab.items())  # from dictionary to list
    antenna_table = table.obj._f_get_child('antenna')
    antenna_table.append(antenna_soltab.items())  # from dictionary to list
    h.close()  # close the new h5parm
    os.remove(working_file)
    return new_h5parm


def apply_h5parm(h5parm, ms, column_out='DATA'):
    '''Creates an NDPPP parset. Applies the output of make_h5parm to the
    measurement set.

    Args:
    new_h5parm (str): The output of dir2phasesol.
    ms (str): The measurement set for self-calibration.
    column_out (str; default = 'DATA'): The column NDPPP writes to.

    Returns:
    None.'''

    # parset is saved in same directory as the h5parm
    parset = os.path.dirname(h5parm) + '/applyh5parm.parset'
    column_in = 'DATA'
    now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    msout =  ms + '-' + str(uuid.uuid4()) + '.MS'

    with open(parset, 'w') as f:  # create the parset
        f.write('# created by apply_h5parm at {}\n'.format(now))
        f.write('msin                = {}\n'.format(ms))
        f.write('msin.datacolumn     = {}\n'.format(column_in))
        f.write('msout               = {}\n'.format(msout))
        f.write('msout.datacolumn    = {}\n'.format(column_out))
        f.write('steps               = [applycal]\n')
        f.write('applycal.type       = applycal\n')
        f.write('applycal.parmdb     = {}\n'.format(h5parm))
        f.write('applycal.correction = phase000\n')
    f.close()

    ndppp_output = subprocess.check_output(['NDPPP', parset])
    os.remove(parset)
    return msout


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

    # convert nan to zero, otherwise nan + X = nan, not X
    ampltides = np.nan_to_num(ampltides)
    amplitude_phases = np.nan_to_num(amplitude_phases)
    phases = np.nan_to_num(phases)

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


def update_list(initial_h5parm, incremental_h5parm, mtf, threshold=0.25,
                amplitude_h5parm=''):
    '''Combine the phase solutions from the initial h5parm and the final
    h5parm. The initial h5parm contains the initial solutions and the final
    h5parm contains the incremental solutions so they need to be added to form
    the final solutions. Calls evaluate_solutions to update the master file
    with a new line appended.

    Args:
    new_h5parm (str): The initial h5parm (i.e. from dir2phasesol).
    loop3_h5parm (str): The final h5parm from loop 3.
    mtf (str): Master text file.
    threshold (float; default=0.25): Threshold determining goodness passed to
        evaluate_solutions.
    amplitude_h5parm (str): HDF5 file containing amplitude (and corresponding
        phase) solutions.

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
    initial_ant = initial_phase.ant[:]
    initial_val = initial_phase.val[:]
    initial_weight = initial_phase.weight[:]

    g = lh5.h5parm(incremental_h5parm)  # from loop3_h5parm
    sol000 = g.getSolset('sol000')  # change to take highest solset?
    incremental_phase = g.getSolset('sol000').getSoltab('phase000')
    antenna_soltab = g.getSolset('sol000').getAnt().items()  # dict to list
    source_soltab = g.getSolset('sol000').getSou().items()  # dict to list

    try:  #  may not contain a direction dimension
        dir_ = incremental_phase.dir[:]
    except:
        dir_ = initial_dir  # if none, take it from the other h5
    incremental_time = incremental_phase.time[:]
    incremental_freq = incremental_phase.freq[:]
    incremental_ant = incremental_phase.ant[:]
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

    # this protects against the antennas not being in the order in each h5parm
    all_antennas = sorted(list(set(initial_ant.tolist() + incremental_ant.tolist())))  # total unique list of antennas
    default_shape = (len(new_times), 1, 2, 1)
    summed_values, summed_weights = [], []

    for antenna in all_antennas:  # for each antenna in either h5parm
        # get values and weights from the first h5parm
        val1 = np.zeros(default_shape)
        wgt1 = np.zeros(default_shape)
        for ant1 in range(len(initial_ant)):
            if antenna == initial_ant[ant1]:
                val1 = initial_val_new[:, :, ant1, :, :]
                wgt1 = initial_weight_new[:, :, ant1, :, :]

        # get values and weights from the second h5parm
        val2 = np.zeros(default_shape)
        wgt2 = np.zeros(default_shape)
        for ant2 in range(len(incremental_ant)):
            if antenna == incremental_ant[ant2]:
                val2 = incremental_val_new[:, :, ant2, :, :]
                wgt2 = incremental_weight_new[:, :, ant2, :, :]

        # and add them, converting nan to zero
        val_new = np.expand_dims(np.nan_to_num(val1) + np.nan_to_num(val2), axis=2)
        wgt_new = np.expand_dims((np.nan_to_num(wgt1) + np.nan_to_num(wgt2)) / 2, axis=2)

        summed_values.append(val_new)
        summed_weights.append(wgt_new)

    vals = np.concatenate(summed_values, axis=2)
    weights = np.concatenate(summed_weights, axis=2)

    # if a h5parm is given with amplitude solutions, add this to our results
    if amplitude_h5parm != '':
        a = lh5.h5parm(amplitude_h5parm)
        amplitude = a.getSolset('sol000').getSoltab('amplitude000')
        amplitude_phases = a.getSolset('sol000').getSoltab('phase000')

        # get amplitude, amplitude_phases and phases onto a new time axis
        newest_times = make_new_times(new_times, amplitude.time[:])

        amp_val, amp_wgt = sort_axes(amplitude)  # adds dir and reorders
        amp_ph_val, amp_ph_wgt = sort_axes(amplitude_phases)

        ph_val_interp = interpolate_time(vals, new_times, newest_times)
        ph_wgt_interp = interpolate_time(vals, new_times, newest_times)

        amp_val_interp = interpolate_time(amp_val, amplitude.time[:], newest_times)
        amp_wgt_interp = interpolate_time(amp_wgt, amplitude.time[:], newest_times)

        amp_ph_val_interp = interpolate_time(amp_ph_val, amplitude.time[:], newest_times)
        amp_ph_wgt_interp = interpolate_time(amp_ph_wgt, amplitude.time[:], newest_times)

        # get list of antennas for the new array
        newest_ant = sorted(list(set(amplitude.ant.tolist() +
                                     amplitude_phases.ant.tolist() +
                                     list(all_antennas))))

        # add the amplitude/phases to the phases
        default_shape = np.zeros((len(newest_times), 1, 1, 1))  # time, freq, pol, dir
        empty_amp_val = np.zeros((len(newest_times), 1, len(newest_ant), 2, 1))  # time, freq, ant, pol, dir
        empty_amp_wgt = np.zeros((len(newest_times), 1, len(newest_ant), 2, 1))  # time, freq, ant, pol, dir
        empty_ph_val = np.zeros((len(newest_times), 1, len(newest_ant), 2, 1))  # time, freq, ant, pol, dir
        empty_ph_wgt = np.zeros((len(newest_times), 1, len(newest_ant), 2, 1))  # time, freq, ant, pol, dir

        summed_values, summed_weights = [], []

        for n in range(len(newest_ant)):  # for each antenna in either h5parm
            antenna = newest_ant[n]
            # set empty variables in case there is not data for all antennas
            amp_val_x, amp_val_y, amp_wgt_x, amp_wgt_y = default_shape, default_shape, default_shape, default_shape
            amp_ph_val_x, amp_ph_val_y, amp_ph_wgt_x, amp_ph_wgt_y = default_shape, default_shape, default_shape, default_shape
            ph_val_x, ph_val_y, ph_wgt_x, ph_wgt_y = default_shape, default_shape, default_shape, default_shape

            # get values and weights from the first h5parm
            for ant in range(len(amplitude.ant)):
                if antenna == amplitude.ant[ant]:
                    amp_val_x = amp_val_interp[:, 0, ant, 0, 0]
                    amp_val_y = amp_val_interp[:, 0, ant, 1, 0]
                    amp_wgt_x = amp_wgt_interp[:, 0, ant, 0, 0]
                    amp_wgt_y = amp_wgt_interp[:, 0, ant, 1, 0]

            for ant in range(len(amplitude_phases.ant)):
                if antenna == amplitude_phases.ant[ant]:
                    amp_ph_val_x = amp_ph_val_interp[:, 0, ant, 0, 0]
                    amp_ph_val_y = amp_ph_val_interp[:, 0, ant, 1, 0]
                    amp_ph_wgt_x = amp_ph_wgt_interp[:, 0, ant, 0, 0]
                    amp_ph_wgt_y = amp_ph_wgt_interp[:, 0, ant, 1, 0]

            # get values and weights from the second h5parm
            for ant in range(len(all_antennas)):
                if antenna == all_antennas[ant]:
                    ph_val_x = ph_val_interp[:, 0, ant, 0, 0]
                    ph_val_y = ph_val_interp[:, 0, ant, 1, 0]
                    ph_wgt_x = ph_wgt_interp[:, 0, ant, 0, 0]
                    ph_wgt_y = ph_wgt_interp[:, 0, ant, 1, 0]

            # and add them
            new_amp_val_x, new_ph_val_x = add_amplitude_and_phase_solutions(amp_val_x, amp_ph_val_x, ph_val_x)
            new_amp_val_y, new_ph_val_y = add_amplitude_and_phase_solutions(amp_val_y, amp_ph_val_y, ph_val_y)
            new_amp_wgt_x = (np.nan_to_num(amp_wgt_x) + np.nan_to_num(ph_wgt_x)) / 2
            new_ph_wgt_x = (np.nan_to_num(amp_ph_wgt_x) + np.nan_to_num(ph_wgt_x)) / 2
            new_amp_wgt_y = (np.nan_to_num(amp_wgt_y) + np.nan_to_num(ph_wgt_y)) / 2
            new_ph_wgt_y = (np.nan_to_num(amp_ph_wgt_y) + np.nan_to_num(ph_wgt_y)) / 2

            empty_amp_val[:, 0, ant, 0, 0] = new_amp_val_x
            empty_amp_val[:, 0, ant, 1, 0] = new_amp_val_y
            empty_amp_wgt[:, 0, ant, 0, 0] = new_amp_wgt_x
            empty_amp_wgt[:, 0, ant, 1, 0] = new_amp_wgt_y
            empty_ph_val[:, 0, ant, 0, 0] = new_ph_val_x
            empty_ph_val[:, 0, ant, 1, 0] = new_ph_val_y
            empty_ph_wgt[:, 0, ant, 0, 0] = new_ph_wgt_x
            empty_ph_wgt[:, 0, ant, 1, 0] = new_ph_wgt_y

        amp_vals = empty_amp_val
        amp_weights = empty_amp_wgt
        vals = empty_ph_val
        weights = empty_ph_wgt
        new_times = newest_times  # redefining these so the phase makeSoltab works correctly regardless
        all_antennas = newest_ant
        a.close()

    freq = np.array([np.mean([initial_freq, incremental_freq])])
    pol = np.array(['XX', 'YY'])

    combined_h5parm = (os.path.splitext(initial_h5parm)[0] + '-' +
                       os.path.basename(incremental_h5parm))

    # write these best phase solutions to the combined_h5parm
    h = lh5.h5parm(combined_h5parm, readonly=False)
    table = h.makeSolset()  # creates sol000
    solset = h.getSolset('sol000')
    c = solset.makeSoltab('phase',
                          axesNames=['time', 'freq', 'ant', 'pol', 'dir'],
                          axesVals=[new_times, freq, all_antennas, pol, dir_],
                          vals=vals,
                          weights=weights)  # creates phase000

    if amplitude_h5parm != '':
        d = solset.makeSoltab('amplitude',
                              axesNames=['time', 'freq', 'ant', 'pol', 'dir'],
                              axesVals=[new_times, freq, all_antennas, pol, dir_],
                              vals=amp_vals,
                              weights=amp_weights)  # creates amplitude000

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
                        '--h5parm0',
                        required=False,
                        type=str,
                        default='/data020/scratch/sean/letsgetloopy/SILTJ132737.15+550405.9_L693725_phasecal.apply_tec_02_c0.h5',
                        help='one hdf5 file')

    parser.add_argument('-P',
                        '--h5parm1',
                        required=False,
                        type=str,
                        default='/data020/scratch/sean/letsgetloopy/SILTJ133749.65+550102.6_L693725_phasecal.apply_tec_00_c0.h5',
                        help='another hdf5 file')

    parser.add_argument('-f',
                        '--ms',
                        required=False,
                        type=str,
                        default='/data020/scratch/sean/letsgetloopy/SILTJ135044.06+544752.7_L693725_phasecal.MS',
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
                        default=[-2.7043, 0.958154], # 0.226893, 0.9512044, 0.244346, 0.9686577
                        nargs='+',
                        help='source positions (radians; RA DEC RA DEC...)')

    args = parser.parse_args()
    mtf = args.mtf
    h5parm0 = args.h5parm0
    h5parm1 = args.h5parm1
    ms = args.ms
    threshold = args.threshold
    cores = args.cores
    directions = args.directions

    make_blank_mtf(mtf=mtf)

    # evaluate_solutions(h5parm=h5parm0, mtf=mtf)
    # evaluate_solutions(h5parm=h5parm1, mtf=mtf)

    # TODO the directions could be read from the ms in this case
    new_h5parms = dir2phasesol_wrapper(mtf=mtf,
                                       ms=ms,
                                       directions=directions,
                                       cores=cores)

    msouts = []
    # for new_h5parm in new_h5parms:
    #     msouts.append(apply_h5parm(h5parm=new_h5parm, ms=ms))  # outputs an ms per direction

    # TODO will this work if loop 2 is run from the directory with the ms?
    # for msout in msouts:  # loop 3
        # run_loop_3 = 'python /data020/scratch/sean/run1/git/long_baseline_pipeline/bin/loop3B_v1.py ' + msout
        # os.system(run_loop_3)

    # update_list(initial_h5parm=h5parm, incremental_h5parm=loop3_phases,
    #             mtf=mtf, threshold=threshold, amplitude_h5parm=loop3_amplitudes)


if __name__ == '__main__':
    main()
