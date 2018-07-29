#!/usr/bin/env python3

import argparse
import numpy as np

def evaluate_solutions(h5parm):
    ''' input:    h5parm with phase solutions
        function: evaluate solutions for each antenna; use XX-YY statistic
                  described in the Hybrid mapping section of the Google Doc and
                  determine validity
        output:   Append to master file with h5parm_name, RA, dec, and one
                  column per antenna with boolean for validity
    '''

    print(h5parm)

def make_h5parm():
    ''' input:    direction of measurement set to be self-calibrated; master
                  file with list of h5parms
        function: find nearest directions; construct a new h5parm that is a
                  combination of best phase solutions from nearest direction,
                  done per antenna. e.g. if the nearest direction has valid
                  solutions for all but the UK antenna, find the UK solutions
                  from the next nearest direction
        output:   a new h5parm to be applied to the measurement set
    '''
    pass

def apply_h5parm():
    ''' input:    the output of make_h5parm; the measurement set for self-
                  calibration
        function: apply the output of make_h5parm to the measurement set
        output:   the measurement set for self-calibration with corrected data
    '''
    pass

def comibine_h5parm():
    ''' input:    the initial h5parm (i.e. from make_h5parm); the final h5parm
                  from loop 3
        function: combine the phase solutions from the initial h5parm and the
                  final h5parm; calls evaluate_solutions to update the master
                  file
        output:   a new h5parm that is a combination of both these h5parms; a
                  new line in the master file
    '''
    pass

def main():
    ''' Starting point: 1st iteration will have produced phase solutions
        independently for all directions, with eachbeing its own h5parm
    '''

    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--h5parm', required = True)
    args = parser.parse_args()

    evaluate_solutions(args.h5parm)

    make_h5parm()

    apply_h5parm()

    # loop 3

    comibine_h5parm()

if __name__ == '__main__':
    main()
