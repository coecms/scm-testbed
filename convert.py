#!/usr/bin/env python
#
# Copyright 2019 Scott Wales
#
# Author: Scott Wales <scott.wales@unimelb.edu.au>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import xarray
import f90nml
import pandas
from pathlib import Path

# Constants from GMTB forcing_file_common.py
g = 9.80665
R_dry = 287.0
c_p = 1004.0
p0 = 100000.0
L_v = 2.5E6
L_s = 2.834E6

def wrf_initial_sounding(init_in, force_in, outputdir):
    T = force_in.T_nudge.isel(time=0)

    sounding = xarray.Dataset({
        'z': init_in.height,
        'u': init_in.u,
        'v': init_in.v,
        'theta': init_in.thetail/(1-(L_v/c_p*init_in.ql+L_s/c_p*init_in.qi)/T),
        'qv': init_in.qt - init_in.ql - init_in.qi,
        })

    Tsfc = force_in.T_surf.isel(time=0).data
    psfc = force_in.p_surf.isel(time=0).data

    # Grab the first level of the input sounding as the surface values
    zsfc = 0
    usfc = sounding.u.isel(levels=0).data
    vsfc = sounding.v.isel(levels=0).data
    qvsfc = sounding.qv.isel(levels=0).data

    with open(outputdir / 'input_sounding', 'w') as f:
        print(zsfc, usfc, vsfc, qvsfc, Tsfc, qvsfc, psfc, file=f)
        sounding.to_dataframe().to_csv(f, sep=' ', index=False, header=False)

def wrf_setup_namelist(casefile, lat, lon, wrfin, outputdir):
    # Translate namelists
    start = pandas.Timestamp(
            year=casefile['case_config']['year'],
            month=casefile['case_config']['month'],
            day=casefile['case_config']['day'],
            hour=casefile['case_config']['hour'])

    runtime = pandas.Timedelta(casefile['case_config']['runtime'],'seconds')

    end = start + runtime

    wrfin['time_control']['start_year'] = start.year
    wrfin['time_control']['start_month'] = start.month
    wrfin['time_control']['start_day'] = start.day
    wrfin['time_control']['start_hour'] = start.hour
    wrfin['time_control']['start_minute'] = start.minute
    wrfin['time_control']['start_second'] = start.second

    wrfin['time_control']['run_days'] = runtime.days
    wrfin['time_control']['run_hours'] = 0
    wrfin['time_control']['run_minutes'] = 0
    wrfin['time_control']['run_seconds'] = runtime.seconds

    wrfin['time_control']['end_year'] = end.year
    wrfin['time_control']['end_month'] = end.month
    wrfin['time_control']['end_day'] = end.day
    wrfin['time_control']['end_hour'] = end.hour
    wrfin['time_control']['end_minute'] = end.minute
    wrfin['time_control']['end_second'] = end.second

    wrfin['scm']['scm_lat'] = lat
    wrfin['scm']['scm_lat'] = lon

    with open(outputdir / 'namelist.input', 'w') as f:
        wrfin.write(f)

def main():
    casefile = 'gmtb-scm-release/scm/etc/case_config/twpice.nml'
    datadir = Path('gmtb-scm-release/scm/data/processed_case_input')

    wrfnml = 'WRF/test/em_scm_xy/namelist.input'
    outputdir = Path('wrf-scm')

    config_in = f90nml.read(casefile)
    init_in = xarray.open_dataset(datadir / (config_in['case_config']['case_name']+'.nc'), 'initial')
    force_in = xarray.open_dataset(datadir / (config_in['case_config']['case_name']+'.nc'), 'forcing')

    wrf_in = f90nml.read(wrfnml)
    
    wrf_initial_sounding(init_in, force_in, outputdir)
    wrf_setup_namelist(config_in, force_in.lat.data[0], force_in.lon.data[0], wrf_in, outputdir)

if __name__ == '__main__':
    main()
