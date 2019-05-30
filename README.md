Run GMTB test cases using WRF
=============================


## Configuring experiments

Edit the WRF configuration in 'wrf_in' as desired. Files from here will be linked in to the run directory.

Experiment specific values such as the start date and initial sounding will be set up by the run script

## Running experiments

On Raijin, load the conda environment first:

    module use /g/data3/hh5/public/modules
    module load conda

Run all test cases:

    ./run_all.sh --gmtb-repo /path/to/gmtb-scm-release --wrf-main /path/to/WRFV3/main

Run a single test case

    ./run_test.py twpice --gmtb-repo /path/to/gmtb-scm-release --wrf-main /path/to/WRFV3/main 

Test cases can be any name from
[gtmb-scm-release/scm/etc/case_config](https://github.com/NCAR/gmtb-scm-release/tree/master/scm/etc/case_config)
(without the '.nml' extension)

## Experiment output

Output files from a run can be found in 'output/{testcase}-{timestamp}', e.g.
'output/twpice-20190529T013725'. The timestamp is intended to provide a
historical record of previous runs