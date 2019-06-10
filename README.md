# Script feeding into the LOFAR long-baseline pipeline

## Issues that have arisen while testing loop 2

* There are three frequency axes in the h5parms from loop 3 now.
* When running loop 3, I encounter the following error:
```
    File "/data020/scratch/sean/run1/git/long_baseline_pipeline/bin/loop3B_v1.py", line 736, in main
    spec_info = casatb.table( vis + '::SPECTRAL_WINDOW')
    RuntimeError: Table SILTJ135044.06+544752.7_L693725_phasecal.MS-eef114bd-d760-48ed-b28e-442b57a55f04.MS::SPECTRAL_WINDOW does not exist
```
But from what I can see in `casabrowser`, it does exist. This error was arising because the command was not executed in the same directory where the MS is -- even though the full filepath was given.
