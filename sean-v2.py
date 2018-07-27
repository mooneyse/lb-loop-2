#!/usr/bin/env python

''' credit: Carole Roskowinski'''

import argparse, os
import numpy as np
import pyrap.tables as pt
import losoto.h5parm as losotoh5
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u

class ReadMs:
    ''' reading in the parameters of the target data with pyrap and putting them
        into directories for further use
    '''

    def __init__(self, ms):
        self.timepara = {'start':0, 'end':0, 'step':0, 'cent':0}
        self.freqpara = {'start':0, 'end':0, 'step':0, 'cent':0}
        self.msname = ms
        if not os.path.isdir(ms): sys.exit('Input does not exist!')

        t  = pt.table(ms, readonly = True, ack = False)

        # get time parameters
        t1 = t.sort('unique desc TIME')
        self.timepara['step'] = t.getcell('EXPOSURE', 0)
        self.timepara['start'] = np.min(t.getcol('TIME')) - self.timepara['step'] / 2.
        self.timepara['end'] = np.max(t.getcol('TIME')) + self.timepara['step'] / 2.
        self.timepara['cent'] = self.timepara['start'] + (self.timepara['end'] - self.timepara['start']) / 2.
        self.mstimevalues = t1.getcol('TIME')[::-1]
        t1.close()

        # get frequency parameters
        freq = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly = True, ack = False)
        self.fullband = freq.getcell('TOTAL_BANDWIDTH', 0)
        self.freqpara['cent'] = freq.getcell('REF_FREQUENCY', 0)
        self.freqpara['step'] = freq.getcell('CHAN_WIDTH', 0)[0]
        self.msfreqvalues = freq.getcell('CHAN_FREQ', 0)
        self.freqpara['start'] = self.msfreqvalues[0] - self.freqpara['step'] / 2.
        self.freqpara['end'] = self.msfreqvalues[-1] + self.freqpara['step'] / 2.
        freq.close()

        # get station names
        antennas = pt.table(t.getkeyword('ANTENNA'), readonly = True, ack = False)
        self.stations = antennas.getcol('NAME')
        self.positions = antennas.getcol('POSITION')
        antennas.close()

        # get pointing information
        pointing = pt.table(t.getkeyword('FIELD'), readonly = True, ack = False)
        self.direction = pointing.getcell('PHASE_DIR', 0) # in radians
        pointing.close()

        t.close()

    def GetTimepara(self, p = ''):
        if p != '': return self.timepara[p]
        else: return self.timepara

    def GetFreqpara(self, p = ''):
        if p != '': return self.freqpara[p]
        else: return self.freqpara

    def GetMSNamepara(self): return self.msname

def main(msname):
        ''' This is called by loop 1 (target selection) and calls loop 3 (self-
            calibration). It interfaces with a global h5parm to create a local
            h5parm for each target.
        '''

        # 1: housekeeping

        # get information on the MS, like pointing centre and frequency range
        tgtname = msname.split('.ms')[0].split('_')[0]
        msinfo = ReadMs(msname)
        target_direction = msinfo.direction[0]
        minfreq = np.min(msinfo.msfreqvalues)
        maxfreq = np.max(msinfo.msfreqvalues)
        bandwidth = (maxfreq - minfreq) / 1e6 # in MHz
        if bandwidth > freq_range: # if the bandwidth is too large, split into channels for WSClean
        	# split into channels for imaging
    	    nchan = int(np.ceil(bandwidth / freq_range))

        print(tgtname, target_direction, bandwidth)
        # 2: find the appropriate solutions and apply them
        # find the best h5parm solutions using the pointing center

if __name__ == '__main__':
    main('/data/scratch/sean/loop-2/L569711_SB051.ms/')
