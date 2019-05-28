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

from pathlib import Path
import f90nml
import xarray
import pandas
import subprocess
import dask
import numpy

# Constants from GMTB forcing_file_common.py
g = 9.80665
R_dry = 287.0
c_p = 1004.0
p0 = 100000.0
L_v = 2.5E6
L_s = 2.834E6

class GTMBWRFTest():
    def __init__(self, testcase, gtmbdir, wrfcfg, workdir, wrfmain):
        self.testcase = testcase
        self.gtmbdir = Path(gtmbdir)
        self.wrfcfg = Path(wrfcfg)
        self.workdir = Path(workdir)
        self.wrfmain = Path(wrfmain)

    def setup_workdir(self):
        self.workdir.mkdir()
        for c in self.wrfcfg.iterdir():
            if c.name != 'namelist.input':
                (self.workdir / c.name).symlink_to(c.resolve())

    def setup_config(self):
        casefile = self.gtmbdir / f'scm/etc/case_config/{self.testcase}.nml'
        casenml = f90nml.read(casefile)

        start = pandas.Timestamp(
                year=casenml['case_config']['year'],
                month=casenml['case_config']['month'],
                day=casenml['case_config']['day'],
                hour=casenml['case_config']['hour'])

        datadir = self.gtmbdir / 'scm/etc' / casenml['case_config']['case_data_dir']
        datafile = datadir / f'{casenml["case_config"]["case_name"]}.nc'

        base = xarray.open_dataset(datafile, decode_cf=False)
        base.time.attrs['units'] = f'seconds since {start.isoformat()}'
        base = xarray.decode_cf(base)

        initial = xarray.open_dataset(datafile, 'initial')
        forcing = xarray.open_dataset(datafile, 'forcing')
        
        initial = initial.update(base)
        forcing = forcing.update(base)

        self.initial_sounding(initial, forcing)

        self.setup_namelist(casenml, start, forcing.lat.data[0], forcing.lon.data[0])

        self.setup_forcing(forcing)

    def initial_sounding(self, initial, forcing):
        T = forcing.T_nudge.isel(time=0)

        sounding = xarray.Dataset({
            'z': initial.height,
            'u': initial.u,
            'v': initial.v,
            'theta': initial.thetail/(1-(L_v/c_p*initial.ql+L_s/c_p*initial.qi)/T),
            'qv': initial.qt - initial.ql - initial.qi,
            })

        Tsfc = forcing.T_surf.isel(time=0).data
        psfc = forcing.p_surf.isel(time=0).data

        # Grab the first level of the input sounding as the surface values
        zsfc = 0
        usfc = sounding.u.isel(levels=0).data
        vsfc = sounding.v.isel(levels=0).data
        qvsfc = sounding.qv.isel(levels=0).data

        with open(self.workdir / 'input_sounding', 'w') as f:
            print(zsfc, usfc, vsfc, Tsfc, qvsfc, psfc, file=f)
            sounding.to_dataframe().to_csv(f, columns=['z','u','v','theta','qv'], sep=' ', index=False, header=False)

    def setup_namelist(self, casenml, start, lat, lon):
        wrfnml = f90nml.read(self.wrfcfg / 'namelist.input')
        # Translate namelists
        runtime = pandas.Timedelta(casenml['case_config']['runtime'],'seconds')

        end = start + runtime

        wrfnml['time_control']['start_year'] = start.year
        wrfnml['time_control']['start_month'] = start.month
        wrfnml['time_control']['start_day'] = start.day
        wrfnml['time_control']['start_hour'] = start.hour
        wrfnml['time_control']['start_minute'] = start.minute
        wrfnml['time_control']['start_second'] = start.second

        wrfnml['time_control']['run_days'] = runtime.days
        wrfnml['time_control']['run_hours'] = 0
        wrfnml['time_control']['run_minutes'] = 0
        wrfnml['time_control']['run_seconds'] = runtime.seconds

        wrfnml['time_control']['end_year'] = end.year
        wrfnml['time_control']['end_month'] = end.month
        wrfnml['time_control']['end_day'] = end.day
        wrfnml['time_control']['end_hour'] = end.hour
        wrfnml['time_control']['end_minute'] = end.minute
        wrfnml['time_control']['end_second'] = end.second

        wrfnml['scm']['scm_lat'] = lat
        wrfnml['scm']['scm_lat'] = lon

        with open(self.workdir / 'namelist.input', 'w') as f:
            wrfnml.write(f)

    def setup_forcing(self, forcing):
        zero = (['time','levels'], dask.array.zeros((forcing.time.size, forcing.levels.size), dtype='f4'))

        Times = numpy.datetime_as_string(forcing.time, unit='s')
        Times = [t.replace('T','_') for t in Times]
        print(Times)

        wrf_forcing = xarray.Dataset({
                'Times': (['time'], Times),
                'Z_FORCE': zero, # TODO: Calculate heights
                'U_G': forcing.u_g.T,
                'V_G': forcing.v_g.T,
                'W_SUBS': forcing.w_ls.T,
                'TH_UPSTREAM_X': zero,
                'TH_UPSTREAM_Y': zero,
                'QV_UPSTREAM_X': zero,
                'QV_UPSTREAM_Y': zero,
                'U_UPSTREAM_X': zero,
                'U_UPSTREAM_Y': zero,
                'V_UPSTREAM_X': zero,
                'V_UPSTREAM_Y': zero,
                'Z_FORCE_TEND': zero,
                'U_G_TEND': zero,
                'V_G_TEND': zero,
                'W_SUBS_TEND': zero,
                'TH_UPSTREAM_X_TEND': zero,
                'TH_UPSTREAM_Y_TEND': zero,
                'QV_UPSTREAM_X_TEND': zero,
                'QV_UPSTREAM_Y_TEND': zero,
                'U_UPSTREAM_X_TEND': zero,
                'U_UPSTREAM_Y_TEND': zero,
                'V_UPSTREAM_X_TEND': zero,
                'V_UPSTREAM_Y_TEND': zero,
                'TAU_X': zero,
                'TAU_X_TEND': zero,
                'TAU_Y': zero,
                'TAU_Y_TEND': zero,
            })
        wrf_forcing.Times.encoding['dtype'] = 'S1'

        wrf_ideal = xarray.open_dataset(self.workdir / 'force_ideal.nc')
        #for k, v in wrf_forcing.items():
        #    if k in wrf_ideal:
        #        v.attrs = wrf_ideal[k].attrs
        wrf_forcing.attrs = wrf_ideal.attrs

        wrf_forcing.to_netcdf(self.workdir / 'force_test.nc')

        wrfnml = f90nml.read(self.workdir / 'namelist.input')
        wrfnml['time_control']['auxinput3_inname'] = 'force_test.nc'
        wrfnml['time_control']['auxinput3_interval_h'] = 3
        wrfnml.write(self.workdir / 'namelist.input', force=True)


    def run_testcase(self):
        self.setup_workdir()
        self.setup_config()
        self.run_ideal()
        self.run_wrf()

    def run_ideal(self):
        subprocess.run([self.wrfmain.resolve() / 'ideal.exe'], cwd=self.workdir)

    def run_wrf(self):
        subprocess.run([self.wrfmain.resolve() / 'wrf.exe'], cwd=self.workdir)


if __name__ == '__main__':
    ts = pandas.Timestamp.utcnow().strftime('%Y%m%dT%H%M%S')
    test = GTMBWRFTest(
            testcase = 'twpice',
            gtmbdir = '../gmtb-scm-release',
            wrfcfg = 'wrf_in',
            workdir = f'test/twpice-{ts}',
            wrfmain = '../WRF4/WRFV3/main',
            )
    test.run_testcase()
