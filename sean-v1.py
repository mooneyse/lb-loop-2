#!/usr/bin/env python

import os
# import argparse
import pyrap.tables as pt
import numpy as np
import losoto.h5parm as losotoh5
import matplotlib.pyplot as plt
from astropy import units as u # for cone search
from astropy.coordinates import SkyCoord # for cone search
# from astropy.coordinates import match_coordinates_sky

# reading in the the parameters of target data with pyrap and putting them into directories for further use

class ReadMs:
    def __init__(self, ms):
        self.timepara = {'start':0, 'end':0, 'step':0, 'cent':0}
        self.freqpara = {'start':0, 'end':0, 'step':0, 'cent':0}
        self.msname = ms
        if not os.path.isdir(ms): sys.exit('Input does not exist!')

        t  = pt.table(ms, readonly = True, ack = False)

        # getting time parameters first
        t1 = t.sort('unique desc TIME')
        self.timepara['step'] = t.getcell('EXPOSURE', 0)
        self.timepara['start'] = np.min(t.getcol('TIME')) - self.timepara['step'] / 2.
        self.timepara['end'] = np.max(t.getcol('TIME')) + self.timepara['step'] / 2.
        self.timepara['cent'] = self.timepara['start'] + (self.timepara['end'] - self.timepara['start']) / 2.
        self.mstimevalues = t1.getcol('TIME')[::-1]
        t1.close()

        # getting frequency parameters
        freq = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly = True, ack = False)
        self.fullband = freq.getcell('TOTAL_BANDWIDTH', 0)
        self.freqpara['cent'] = freq.getcell('REF_FREQUENCY', 0)
        self.freqpara['step'] = freq.getcell('CHAN_WIDTH', 0)[0]
        self.msfreqvalues = freq.getcell('CHAN_FREQ', 0)
        self.freqpara['start'] = self.msfreqvalues[0] - self.freqpara['step'] / 2.
        self.freqpara['end'] = self.msfreqvalues[-1] + self.freqpara['step'] / 2.
        freq.close()

        # getting station names
        antennas = pt.table(t.getkeyword('ANTENNA'), readonly = True, ack = False)
        self.stations = antennas.getcol('NAME')
        self.positions = antennas.getcol('POSITION')
        antennas.close()

        # getting pointing information
	    pointing = pt.table(t.getkeyword('FIELD'), readonly = True, ack = False)
        self.direction = pointing.getcell('PHASE_DIR', 0) # this returns radians
	    pointing.close()

        t.close()

    def GetTimepara(self, p = ''):
        if p != '': return self.timepara[p]
        else: return self.timepara

    def GetFreqpara(self, p = ''):
        if p != '': return self.freqpara[p]
        else: return self.freqpara

    def GetMSNamepara(self): return self.msname

def ReadH5(h5parm):
    # open the h5parm
    myh5 = losotoh5.h5parm(h5parm)
    # get the solution set
    mysolsets = myh5.getSolsetNames()

def make_h5_solutions(direction, h5parm):
    # where direction = [ra, dec] (in degrees?)
    print 'Plug in script from Carole'
    return(0)

def find_h5_solutions(target_coo, all_solutions_names):
    # nomenclature for the name of the h5parm: dir_ra_dec.h5
    # where RA and dec are the coordinates in degrees (with max precision)
    # RQ for management of the target_coo: no need to get lugged obj sky coo tt long; just do conv just for match

    ''' target_coo: provided in rad, rad (from MSread)
        all_soltutions_names: list of the solutions already available (from higher tiers - [[ra, dec], h5parm_name], ra & dec in deg])
    '''

    # if first calibrator then do phase only: no matching
    print all_solutions_names
    if all_solutions_names: # if this list is not empty, we are not on the brightest calibrators anymore
        # cone search
            # NB: readMS gives coo in radian but skycoo can work with all if good unit given?
        # print target_coo[0], target_coo[1]
        coo_tar = SkyCoord(target_coo[0], target_coo[1], unit = 'rad') # TODO check the frame! default: ircs
            # seems to return them in degrees by default

        coo_sol_min = SkyCoord(all_solutions_names[0][0][0], all_solutions_names[0][0][1], unit = 'deg')
        sep_min = coo_tar.separation(coo_sol_min)
        name_min = all_solutions_names[0][1]

        for sol in all_solutions_names[1:]:
            coo = SkyCoord(sol[0][0], sol[0][1], unit = 'deg')
            sep = coo_tar.separation(coo)

            if sep.degree < sep_min.degree: # TODO already take care now if pls exactly the same d?
                sep_min = sep
                coo_sol_min = coo
                name_min = sol[-1]

        # print 'final'
        # print sep_min.degree
        # print coo_sol_min
        # print name_min

    # Continue the processing: common for all tiers
    # calc, calc calc
    # from Martin too

    # TODO rest of the processing common for all tiers?
        # REP wi
        # + cf how to combine / or at worst; put 2 times tests

        # cone search if fainter sources
        # fg or check for 1st loop to avoid looking for previous solution
            # what is a best? A global fg or check if file exist?
            # other way for bookeeping; list of coo and if list empty; <=> q ​​on 1st loop
            # ms; need to know this list uk lvl sup ie loop 1
            # => have overall direction then check if empty; will be signal for 1st level of calibration
            # RQ: if h5parm ac N already; no need for direction
            # what should be normal pm to the level of the pipeline gene

             # pb of multithreading though; if all start together okish
             # but if one finish before one start; not valid...
             # TODO deal ac fg

             # salso deal with gl layer;
             # pw happen q list / dir contains close soil that is being solved or at best
             # => force keep solutions uk level sup
             # pw regler pb fg

        # phase sol: TODO: even pt?
        # NB ds new fct

        # merging

    # cone search on name
        # => list of solutions
    # solving and merging
        # TODO solving? applying

    # return(0)

def input2strlist_nomapfile(invar):
    ''' from bin/download_IONEX.py
        give the list of MSs from the list provided as a string
    '''

    str_list = None
    if type(invar) is str:
        if invar.startswith('[') and invar.endswith(']'):
            str_list = [f.strip(' \'\"') for f in invar.strip('[]').split(',')]
        else:
            str_list = [invar.strip(' \'\"')]
    elif type(invar) is list:
        str_list = [str(f).strip(' \'\"') for f in invar]
    else:
        raise TypeError('input2strlist: Type '+str(type(invar))+' unknown!')
    return str_list

def writeApplyH5parmParset(h5parmName, parsetname = 'ndppp_applyH5.parset', incol = 'DATA', outcol = 'CORRECTED_DATA', outms = '.'):
    # open write-only parset; this will overwrite previous files of the same name
    with open(parsetname, 'w') as f:
        f.write('msin.datacolumn=%s\n' % incol)
        f.write('msout=%s\n' % outms)
        f.write('msout.datacolumn=%s\n' % outcol)
        f.write('steps=[applycal]\n')
        f.write('applycal.type=applycal\n')
        f.write('applycal.parmdb=%s\n' % h5parmName)
        f.write('applycal.correction=phase000\n')
    f.close()

# borrowed directly from https://github.com/brentjens/lofar-antenna-positions
def normalized_earth_radius(latitude_rad):
    wgs84_f = 1. / 298.257223563
    return 1. / np.sqrt(np.cos(latitude_rad) ** 2. + ((1. - wgs84_f) ** 2.) * (np.sin(latitude_rad) ** 2.))

# borrowed directly from https://github.com/brentjens/lofar-antenna-positions
def geographic_from_xyz(xyz_m):
    ''' Compute lon, lat, and height
    '''
    wgs84_a = 6378137.
    wgs84_f = 1. / 298.257223563
    wgs84_e2 = wgs84_f * (2. - wgs84_f)

    x_m, y_m, z_m = xyz_m
    lon_rad = np.arctan2(y_m, x_m)
    r_m = np.sqrt(x_m ** 2 + y_m ** 2)
    # iterate to latitude solution
    phi_previous = 1e4
    phi = np.arctan2(z_m, r_m)
    while abs(phi -phi_previous) > 1.6e-12:
        phi_previous = phi
        phi = np.arctan2(z_m + wgs84_e2 * wgs84_a * normalized_earth_radius(phi) * np.sin(phi), r_m)
    lat_rad = phi
    height_m = r_m * np.cos(lat_rad) + z_m * np.sin(lat_rad) - wgs84_a * np.sqrt(1. - wgs84_e2 * np.sin(lat_rad) ** 2)
    return {'lon_rad': lon_rad, 'lat_rad': lat_rad, 'height_m': height_m}

def gcirc_distance(dict1, dict2):
    # haversine formula
    lat_midpoint = 0.5 * (dict1['lat_rad'] - dict2['lat_rad'])
    lon_midpoint = 0.5 * (dict1['lon_rad'] - dict2['lon_rad'])
    a = np.sin(lat_midpoint) ** 2. + np.cos(dict1['lat_rad']) * np.cos(dict2['lat_rad']) * np.sin(lon_midpoint) ** 2.
    return 6378137. * normalized_earth_radius(dict1['lat_rad']) * 2. * np.arctan2(np.sqrt(a), np.sqrt(1. - a))

# closure

def closure(vis, tel, lastv = -1, plotfile = 'clplot.png'):
    # find which number is which antenna
    itels = np.sort(get_idx_tels (vis, tel))
    if itels == []:
        return -1

    # make three reference MSs with pointers
    print 'itels', itels
    d1, ut1, uvw = dget_t(vis, itels[0], itels[1])
    d2, ut2, uvw = dget_t(vis, itels[1], itels[2])
    d3, ut3, uvw = dget_t(vis, itels[0], itels[2])
    a1, p1 = getap(d1[:lastv])
    a2, p2 = getap(d2[:lastv])
    a3, p3 = getap(d3[:lastv])
    clph = p1 + p2 - p3
    np.putmask(clph, clph > np.pi, clph - 2 * np.pi)
    np.putmask(clph, clph < -np.pi, clph + 2 * np.pi)
    # return a statistic - is 1.64 for random closure phase, less for coherent
    if len(plotfile):
        plt.plot(clph, 'b+')
        plt.savefig(plotfile)
    return np.nanmean(np.gradient(np.unwrap(clph)) ** 2.)

def get_idx_tels(data, tel):
    os.system('taql \'select NAME from %s/ANTENNA\' >closure_which' % data)
    f = open('closure_which')
    for line in f:
        if 'select' in line or not len(line.rstrip('\n')):
            continue
        try:
            a = int(line) # it's just a number, use STATION column instead
            f.close()
            os.system('taql \'select STATION from %s/ANTENNA\' >closure_which' % data)
            break
        except:
            f.close()
            break
    idx_tels, iline = [-1] * len(tel), 0
    f = open('closure_which')
    for line in f:
        if not 'select' in line:
            for i in range(len(tel)):
                if tel[i] == line[:len(tel[i])]:
                    idx_tels[i] = iline
            iline += 1
    f.close()
    if -1 in idx_tels:
        os.system('cat closure_which')
        print 'I did not find one or more of the telescopes.'
        print 'The telescopes present are those in the list above.'
        return []
    os.system('rm closure_which')
    return idx_tels

def dget_t(vis, tel1, tel2):
    os.system('taql \'select from %s where ANTENNA1==%d and ANTENNA2==%d giving %s\'' % (vis, tel1, tel2, 'cl_temp.ms'))
    t = pt.table('cl_temp.ms')
    ut = np.ravel(np.asarray([tuple(each.values()) for each in t.select('TIME')]))
    spw = np.ravel(np.asarray([tuple(each.values()) for each in t.select('DATA_DESC_ID')]))
    dc = t.select('DATA')
    d = np.asarray([tuple(each.values()) for each in dc])[:, 0, :, :]
    d = np.swapaxes(d, 0, 2) # makes the order pol - chan - time as in casa
    for u in t.select('UVW'):
        try:
            uvw = np.vstack ((uvw, u['UVW']))
        except:
            uvw = np.copy(u['UVW'])
    if spw.sum():
        for i in np.unique(spw):
            new = np.take(d, np.argwhere(spw == i), axis = 2)[:, :, :, 0]
            try:
                d_out = np.concatenate((d_out, new), axis = 1)
            except:
                d_out = np.copy(new)
        d = d_out
    return d, ut, uvw

def getap(d, pol = 0):
    ph = np.sum(d[pol], axis = 0)
    return np.sum(abs(d[pol]), axis = 0) / d.shape[1], np.arctan2(ph.imag, ph.real)

def main(msname, all_sol_names, freq_range = 10, data_col = 'DATA', mosaic_rad_arcsec = 0):

    ''' This script is called by loop 1 (target selection) and calls loop 3
        (self-calibration loop). It assumes that loop 1 has appropriately
        prepared the data and that calibration and imaging will be handled by
        loop 3. It interfaces with a global h5parm to create a local h5parm for
        each target it is run on.
    '''

    # section 1: housekeeping
    # send in part to the __main__
    # (due to issues with the parsing of a list and the input2strlist_nomapfile function)

    # get information on the measurement set, like pointing centre and frequency range
    tgtname = msname.split('.ms')[0].split('_')[0]
    msinfo = ReadMs(msname)
    target_direction = msinfo.direction[0]
    minfreq = np.min(msinfo.msfreqvalues)
    maxfreq = np.max(msinfo.msfreqvalues)
    bandwidth = (maxfreq - minfreq) / 1e6 # convert to MHz
    # if the bandwidth is too large, split into channels for wsclean
    if bandwidth > freq_range:
    	# split into channels for imaging
	    nchan = int(np.ceil(bandwidth / freq_range))

    # section 2: find the appropriate solutions and apply them
    # using the pointing center, find the best h5parm solutions
    # my_h5parm = make_h5_solutions(target_direction, h5parm)

    my_hy5parm = find_h5_solutions(target_direction, all_sol_names) # NB: to be rechfer pm; pw need the MS set

    # apply these solutions
        # pb: is this ok ac rest loop? must div ds fct?
        # NB: Each h5 sslt ac solution that is applied to MS
        # => h5 step should be ok for then?

    applyParset = 'ndppp_apply_' + tgtname + '.parset'
    writeApplyH5parmParset(my_h5parm, parsetname = applyParset)
    ss = 'NDPPP %s msin=%s' % (applyParset, msname)
    os.system(ss)

    # section 3: set imaging parameters
    # determine longest coherent baseline using closure phase
    # get a list of stations
    stations = msinfo.stations
    ctel1 = [s for s in stations if 'ST' in s][0]
    rs_tels = [s for s in stations if 'RS' in s]
    # just take the first RS - later can update to find best one
    ctel2 = rs_tels[0]
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
        closure_phase_scatter = closure(msname, closure_tels,plotfile = '')
        cp_scatter = np.append(cp_scatter, closure_phase_scatter)
        # get station position dictionary
        intl_dict = geographic_from_xyz(intl_pos[i])
        bl_len_m = gcirc_distance(st_dict, intl_dict)
        bl_lengths = np.append(bl_lengths, bl_len_m)

    # plot the scatter in closure phase versus baseline length

    # lambda / max_distance * conversion_to_arcsec
    resolution_element = (2.99e8 / 144e6) / np.max(baseline_lengths) * 206265.
    cell_size = resolution_element / 5.
    l_max_meters = np.max(baseline_lengths)

    # use LoTSS catalogue (this will always be avaliable, and more applicable than FIRST)
    # to get the largest size of the source (i.e. DC_Maj) to set the image size, with max = 2 arcmin, min = 20 arcsec
    # or a standard size if in 'mosaic' mode
    if mosaic_rad_arcsec > 0:
        max_imsize = 2 * mosaic_rad_arcsec
        ncells = np.ceiling(max_imsize / cell_size)

    # section 4: get initial model

    # get wsclean version - there are differences between 2.3 and 2.4 that cause it to fail
    ss = 'wsclean --version > wsclean.version'
    os.system(ss)
    with open('wsclean.version', 'r') as f:
        lines = f.readlines()
    f.close()
    ss = 'rm wsclean.version'
    os.system(ss)

    wsclean_version = lines[1].split()[2]
    if wsclean_version == '2.3':
	rms_bkg = '-rms-background-window 25'
        save_cc = '-save-component-list'
    else:
	rms_bkg = '-local-rms'
	save_cc = '-save-source-list'

    # iterate over loop 3 until convergence
    # combined solution from Martin

# if __name__ == '__main__':
    # import argparse

    # parser = argparse.ArgumentParser(description = 'LOFAR long-baseline pipeline: loop 2')
    # parser.add_argument( 'msname', type = str, help = 'Measurement set name')
    # parser.add_argument( 'h5parm', type = str, help = 'H5parm containing phase solutions') # TODO rather coo or just initialization of the N higher for remain as one wants / level global and pass the - or load the - the level infinite

    # parser.add_argument('all_solutions_names', type = str, nargs = '+', help = 'List of available solutions h5parm files ([[ra, dec], name])') TODO cf if or global
        # might remind an issue by calling it from command line, especially when no solution yet available

    # parser.add_argument('-f', '--freq_range', dest = 'freq_range', default = 10, help = 'Frequency range in MHz before breaking into channels to image (default 10)')
    # parser.add_argument('-d', dest = 'data_col', type = str, default = 'DATA', help = 'Which data column to start with the imaging')

    # parse the arguments # TODO to complete
    # args = parser.parse_args()
    # msname = args.msname
    # h5parm = args.h5parm
    # all_sol_names = args.all_solutions_names
    # freq_range = args.freq_range
    # data_col = args.data_col

    # main(msname, input2strlist_nomapfile(all_sol_names), freq_range, data_col) # TODO check type

sol_avl = [[[277.3825, 48.74611111], 'C10'], [[212.836625, 52.20219444], 'D9'], [[212.835375, 52.20297222], 'B3'], [[24.42208333, 33.15972222], 'A1']] # coo should be tuple and not list
     # coo deg
main('/data009/scratch/carosko/chess/data_sean/L569711_SB051_uv.MS/', sol_avl)
# TODO recheck syntax for parse ms qd list
# > python loop2_c_c_debug_a.py /data009/scratch/carosko/chess/data_sean/L569711_SB051_uv.MS/ 'bana', ' nana' -f 30 -d 'CORRECTED' ?
