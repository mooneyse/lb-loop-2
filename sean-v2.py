#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

def find_h5_solutions(target_coo, all_solutions_names):
    ''' target_coo:           - in radians
                              - from ReadMs
        all_soltutions_names: - list of the solutions already available
                              - [ra, dec], h5parm_name
                              - ra, dec in radians
                              - nomenclature for h5parm_name: dir_ra_dec.h5
                                with ra, dec in degrees with maximum precision
                              - from higher tiers
    '''

    # do phase only (i.e. no matching) if this is the first calibrator
    print('all solution names: %s' % all_solutions_names)
    if all_solutions_names: # we are not on the brightest calibrators anymore if this list is not empty
        # cone search
            # NB ReadMs gives coo in radian but skycoo can work with all if good unit given?
        print('target 0, target 1: %s, %s' % (target_coo[0], target_coo[1]))
        coo_tar = SkyCoord(target_coo[0], target_coo[1], unit = 'rad') # TODO check the frame - default: ircs
        # seems to return them in degrees by default

    #     coo_sol_min = SkyCoord(all_solutions_names[0][0][0], all_solutions_names[0][0][1], unit = 'deg')
    #     sep_min = coo_tar.separation(coo_sol_min)
    #     name_min = all_solutions_names[0][1]
    #
    #     for sol in all_solutions_names[1:]:
    #         coo = SkyCoord(sol[0][0], sol[0][1], unit = 'deg')
    #         sep = coo_tar.separation(coo)
    #
    #         if sep.degree < sep_min.degree: # TODO already take care now if pls exactly the same d?
    #             sep_min = sep
    #             coo_sol_min = coo
    #             name_min = sol[-1]
    #
    #     # print 'final'
    #     # print sep_min.degree
    #     # print coo_sol_min
    #     # print name_min
    #
    # # Continue the processing: common for all tiers
    # # calc, calc calc
    # # from Martin too
    #
    # # TODO rest of the processing common for all tiers?
    #     # REP wi
    #     # + cf how to combine / or at worst; put 2 times tests
    #
    #     # cone search if fainter sources
    #     # fg or check for 1st loop to avoid looking for previous solution
    #         # what is a best? A global fg or check if file exist?
    #         # other way for bookeeping; list of coo and if list empty; <=> q ​​on 1st loop
    #         # ms; need to know this list uk lvl sup ie loop 1
    #         # => have overall direction then check if empty; will be signal for 1st level of calibration
    #         # RQ: if h5parm ac N already; no need for direction
    #         # what should be normal pm to the level of the pipeline gene
    #
    #          # pb of multithreading though; if all start together okish
    #          # but if one finish before one start; not valid...
    #          # TODO deal ac fg
    #
    #          # salso deal with gl layer;
    #          # pw happen q list / dir contains close soil that is being solved or at best
    #          # => force keep solutions uk level sup
    #          # pw regler pb fg
    #
    #     # phase sol: TODO: even pt?
    #     # NB ds new fct
    #
    #     # merging
    #
    # # cone search on name
    #     # => list of solutions
    # # solving and merging
    #     # TODO solving? applying
    #
    # # return(0)

def main(msname, all_sol_names, freq_range = 10):
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
        if bandwidth > freq_range: # if the bandwidth is too large, split into channels for imaging in WSClean
    	    nchan = int(np.ceil(bandwidth / freq_range))

        print('Target name: %s' % tgtname)
        print('Target RA, dec: %s' % target_direction)
        print('Bandwidth: %s', bandwidth)

        # 2: find the appropriate solutions and apply them

        # find the best h5parm solutions using the pointing centre
        my_hy5parm = find_h5_solutions(target_direction, all_sol_names)


if __name__ == '__main__':
    sol_avl = [[[277.3825, 48.74611111], 'C10'], [[212.836625, 52.20219444], 'D9'], [[212.835375, 52.20297222], 'B3'], [[24.42208333, 33.15972222], 'A1']] # coo should be tuple and not list
    main('/data/scratch/sean/loop-2/L569711_SB051.ms/', sol_avl)
