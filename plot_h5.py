#!/usr/bin/env python2.7

'''Plot HDF5 solutions for a specified station from two HDF5 files.'''

from __future__ import print_function
import argparse
import numpy as np
import sys
import losoto.h5parm as lh5
from losoto.lib_operations import reorderAxes
import matplotlib.pyplot as plt

__author__ = 'Sean Mooney'
__date__ = '01 June 2019'

def get_values(h5, station='ST001', polarisation='XX'):
    station_exists=False
    h = lh5.h5parm(h5, readonly=True)
    phase = h.getSolset('sol000').getSoltab('phase000')
    time = phase.time
    for ant in range(len(phase.ant)):
        if phase.ant[ant] == station:
            station_exists = True
            break

    if not station_exists:
        print('{} does not exist in {}.'.format(station, h5))
        h.close()
        sys.exit()

    try:
        reordered_values = reorderAxes(phase.val, phase.getAxesNames(), ['time', 'freq', 'ant', 'pol', 'dir'])
        my_values_xx = reordered_values[:, 0, ant, 0, 0]
        my_values_yy = reordered_values[:, 0, ant, 1, 0]
    except:
        reordered_values = reorderAxes(phase.val, phase.getAxesNames(), ['time', 'freq', 'ant', 'pol'])
        my_values_xx = reordered_values[:, 0, ant, 0]
        my_values_yy = reordered_values[:, 0, ant, 1]

    h.close()
    if polarisation == 'XX':
        values = my_values_xx
    elif polarisation == 'YY':
        values = my_values_yy

    return values, time



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('h5a', type=str)
    parser.add_argument('h5b', type=str)
    parser.add_argument('h5c', type=str)
    parser.add_argument('-p', '--polarisation', required=False, type=str, default='XX')

    args = parser.parse_args()
    h5a = args.h5a
    h5b = args.h5b
    h5c = args.h5c
    polarisation = args.polarisation

    stations = ['RS106HBA', 'RS205HBA', 'RS208HBA', 'RS210HBA', 'RS305HBA', 'RS306HBA', 'RS307HBA',
                'RS310HBA', 'RS406HBA', 'RS407HBA', 'RS409HBA', 'RS503HBA', 'RS508HBA', 'RS509HBA',
                'DE601HBA', 'DE602HBA', 'DE603HBA', 'DE604HBA', 'DE605HBA', 'FR606HBA', 'SE607HBA',
                'UK608HBA', 'DE609HBA', 'PL610HBA', 'PL611HBA', 'PL612HBA', 'IE613HBA', 'ST001']

    fig = plt.figure(figsize=(17, 8.5))
    for i in range(len(stations)):
        values1, times1 = get_values(h5a, station=stations[i], polarisation=polarisation)
        values2, times2 = get_values(h5b, station=stations[i], polarisation=polarisation)
        values3, times3 = get_values(h5c, station=stations[i], polarisation=polarisation)
        plt.subplot(4, 7, i + 1)
        plt.plot(times1, values1, 'k-', lw=1, alpha=0.33)
        plt.plot(times2, values2, 'r-', lw=1, alpha=0.33)
        plt.plot(times3, values3, 'b-', lw=1, alpha=0.33)
        plt.xticks([])
        plt.yticks([])
        plt.xlabel('Time')
        plt.ylabel('Phase')
        plt.xlim(min([min(times1), min(times2), min(times3)]), max([max(times1), max(times2), max(times3)]))
        plt.ylim(-np.pi, np.pi)
        plt.title(stations[i])
    fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
