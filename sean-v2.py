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
        coo_tar = SkyCoord(target_coo[0], target_coo[1], unit = 'rad') # degrees by default; ReadMs gives coordinates in radians
        # TODO check the frame - default: ircs

        coo_sol_min = SkyCoord(all_solutions_names[0][0][0], all_solutions_names[0][0][1], unit = 'deg')
        sep_min = coo_tar.separation(coo_sol_min)
        name_min = all_solutions_names[0][1]

        for sol in all_solutions_names[1:]:
            coo = SkyCoord(sol[0][0], sol[0][1], unit = 'deg')
            sep = coo_tar.separation(coo)

            if sep.degree < sep_min.degree: # TODO what if it is exactly the same distance?
                sep_min = sep
                coo_sol_min = coo
                name_min = sol[-1]

        print('minimum separation: %s' % sep_min.degree)
        print('minimum ra, dec:    %s' % coo_sol_min)
        print('minimum name:       %s' % name_min)

    # TODO continue the processing that is common to all tiers

def writeApplyH5parmParset(h5parmName, parsetname = 'ndppp_applyH5.parset', incol = 'DATA', outcol = 'CORRECTED_DATA', outms = '.'):
    with open(parsetname, 'w') as f: # overwrites previous files of the same name
        f.write('msin.datacolumn = %s\n' % incol)
        f.write('msout = %s\n' % outms)
        f.write('msout.datacolumn = %s\n' % outcol)
        f.write('steps = [applycal]\n')
        f.write('applycal.type = applycal\n')
        f.write('applycal.parmdb = %s\n' % h5parmName)
        f.write('applycal.correction = phase000\n')
    f.close()

def geographic_from_xyz(xyz_m): # https://github.com/brentjens/lofar-antenna-positions
    ''' Compute lonitude, latitude, and height.
    '''

    wgs84_a = 6378137.
    wgs84_f = 1. / 298.257223563
    wgs84_e2 = wgs84_f * (2. - wgs84_f)

    x_m, y_m, z_m = xyz_m
    lon_rad = np.arctan2(y_m, x_m)
    r_m = np.sqrt(x_m ** 2. + y_m ** 2.)
    # iterate to latitude solution
    phi_previous = 1e4
    phi = np.arctan2(z_m, r_m)
    while abs(phi - phi_previous) > 1.6e-12:
        phi_previous = phi
        phi = np.arctan2(z_m + wgs84_e2 * wgs84_a * normalized_earth_radius(phi) * np.sin(phi), r_m)
    lat_rad = phi
    height_m = r_m * np.cos(lat_rad) + z_m * np.sin(lat_rad) - wgs84_a * np.sqrt(1. - wgs84_e2 * np.sin(lat_rad) ** 2.)
    return {'lon_rad': lon_rad, 'lat_rad': lat_rad, 'height_m': height_m}

def closure(vis, tel, lastv = -1, plotfile = 'clplot.png'):
    # find which number is which antenna
    itels = np.sort(get_idx_tels (vis, tel))
    if itels == []:
        return -1

    # make 3 reference MSs with pointers
    print('itels: %s' % itels)
    d1, ut1, uvw = dget_t(vis, itels[0], itels[1])
    d2, ut2, uvw = dget_t(vis, itels[1], itels[2])
    d3, ut3, uvw = dget_t(vis, itels[0], itels[2])
    a1, p1 = getap(d1[:lastv])
    a2, p2 = getap(d2[:lastv])
    a3, p3 = getap(d3[:lastv])
    clph = p1 + p2 - p3
    np.putmask(clph, clph > np.pi, clph - 2 * np.pi)
    np.putmask(clph, clph < -np.pi, clph + 2 * np.pi)
    # return a statistic; 1.64 is random closure phase, < 1.64 is more coherent
    if len(plotfile):
        plt.plot(clph, 'b+')
        plt.savefig(plotfile)
    return np.nanmean(np.gradient(np.unwrap(clph)) ** 2.)

def gcirc_distance(dict1, dict2): # haversine formula
    lat_midpoint = 0.5 * (dict1['lat_rad'] - dict2['lat_rad'])
    lon_midpoint = 0.5 * (dict1['lon_rad'] - dict2['lon_rad'])
    a = np.sin(lat_midpoint) ** 2. + np.cos(dict1['lat_rad']) * np.cos(dict2['lat_rad']) * np.sin(lon_midpoint) ** 2.
    return 6378137. * normalized_earth_radius(dict1['lat_rad']) * 2. * np.arctan2(np.sqrt(a), np.sqrt(1. - a))

def main(msname, all_sol_names, freq_range = 10, mosaic_rad_arcsec = 0):
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

    print('Target name:     %s' % tgtname)
    print('Target RA, dec:  %s' % target_direction)
    print('Bandwidth (MHz): %s' % bandwidth)

    # 2: find the appropriate solutions and apply them

    # find the best h5parm solutions using the pointing centre
    my_hy5parm = find_h5_solutions(target_direction, all_sol_names)

    # apply these solutions
    applyParset = tgtname + '_ndppp_apply.parset'
    my_h5parm = '/data/scratch/sean/loop-2/ILTJ132737.2+550406.2.h5' # placeholder
    writeApplyH5parmParset(my_h5parm, parsetname = applyParset)
    s = '#NDPPP %s msin=%s' % (applyParset, msname)
    os.system(s)

    # 3. set imaging parameters

    # determine longest coherent baseline using closure phase
    stations = msinfo.stations # get a list of stations
    print(stations)
    ctel1 = [s for s in stations if 'ST' in s][0]
    rs_tels = [s for s in stations if 'RS' in s]
    ctel2 = rs_tels[0] # taking the 1st RS TODO update to find the best one
    intl_tels = [s for s in stations if 'RS' not in s and 'ST' not in s]

    # get station positions for use in calculating baseline lengths
    positions = msinfo.positions
    st_index = [i for i in np.arange(0, len(stations)) if 'ST' in stations[i]]
    st_pos = positions[st_index[0]]
    st_dict = geographic_from_xyz(st_pos)
    intl_index = [i for i, val in enumerate(stations) if val in intl_tels]
    intl_pos = positions[intl_index]

    # find the closure phase scatter
    cp_scatter = np.zeros(0)
    bl_lengths = np.zeros(0)
    for i, intl_tel in enumerate(intl_tels):
        closure_tels = [ctel1, ctel2, intl_tel]
        closure_phase_scatter = closure(msname, closure_tels, plotfile = '')
        cp_scatter = np.append(cp_scatter, closure_phase_scatter)
        intl_dict = geographic_from_xyz(intl_pos[i]) # get station position dictionary
        bl_len_m = gcirc_distance(st_dict, intl_dict)
        bl_lengths = np.append(bl_lengths, bl_len_m)

    # plot the scatter in closure phase versus baseline length
    resolution_element = (2.99e8 / 144e6) / np.max(baseline_lengths) * 206265. # lambda / max(distance) * (conversion to arcseconds)
    cell_size = resolution_element / 5.
    l_max_meters = np.max(baseline_lengths)

    # use LoTSS catalogue (this will always be avaliable, and more applicable than FIRST)
    # to get the largest size of the source (i.e. DC_Maj) and then set the image size, with (max, min) = (2", 20')
    # or to a standard size if in mosaic mode
    if mosaic_rad_arcsec > 0:
        max_imsize = 2 * mosaic_rad_arcsec
        ncells = np.ceiling(max_imsize / cell_size)

    # 4. get initial model

if __name__ == '__main__':
    sol_avl = [[[277.3825, 48.74611111], 'C10'], [[212.836625, 52.20219444], 'D9'], [[212.835375, 52.20297222], 'B3'], [[24.42208333, 33.15972222], 'A1']] # coo should be tuple and not list
    main('/data/scratch/sean/loop-2/L569711_SB051.ms/', sol_avl)
