# Script feeding into the LOFAR long-baseline pipeline

## How to use loop 2

* Log in to CEP3.
* Make a clean directory (`mkdir loop-2-tests` then `cd loop-2-tests`).
* Move the HDF5 files here (`mv [file.h5] .`).
* Clone the repository (`git clone https://github.com/mooneyse/lb-loop-2.git`).
* Load the LOFAR software and LoSoTo (`module load lofar losoto`).
* Loop 2 is the `hdf5_fuctions.py` script. It runs with Python 2.
* View the arguments with `lb-loop-2/hdf5_functions.py --help`.
* Do something like this:
```
$ lb-loop-2/hdf5_functions.py \
  -h5parm0 /data020/scratch/sean/letsgetloopy/SILTJ132737.15+550405.9_L693725_phasecal.apply_tec_02_c0.h5 \
  -h5parm1 /data020/scratch/sean/letsgetloopy/SILTJ133749.65+550102.6_L693725_phasecal.apply_tec_00_c0.h5 \
  --ms /data020/scratch/sean/letsgetloopy/SILTJ135044.06+544752.7_L693725_phasecal.MS \
  -mtf /data020/scratch/sean/letsgetloopy/mtf.txt \
  --direction -2.7043 0.958154  # RA, DEC in radians
```
* This will make the master text file if it does not exist, evaluate the solutions in the HDF5 files, make a new HDF5 file in the direction supplied, and apply this HDF5 file to the MS, outputting a new MS.
* Then run loop 3 on this new MS (`python2 lb-loop-2/loop3B_v1.py [file.ms]`).
* When loop 3 is done, edit `main()` in `hdf5_fuctions.py` to run the `update_list` function, passing it the HDF5 file originally applied to the MS, the HDF5 file outputted by loop 3, and the master text file. This function combines these incremental solutions, evaluates them, and writes the output to the master text file.
* For a list of issues with loop 2, see the [Issues on GitHub](https://github.com/mooneyse/lb-loop-2/issues).
