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
import argparse

# Constants from GMTB forcing_file_common.py
g = 9.80665
R_dry = 287.0
c_p = 1004.0
p0 = 100000.0
L_v = 2.5E6
L_s = 2.834E6


def pressure_to_z(zsfc, T, p):
    """
    Convert GMTB pressure levels to heights

    Args:
        zsfc: Surface height
        T: Temperature
        p: GMTB pressure levels
    """
    # Following GMTB 'forcing_file_common.py'
    z = numpy.zeros(T.shape)
    z[:,0] = zsfc
    for k in range(p.size - 1):
        dphi = -0.5*R_dry*((T[:,k]+T[:,k+1])/(0.5*(p[k]+p[k+1])))*(p[k+1]-p[k])
        z[:,k+1] = z[:,k]+dphi/g
    return z


class gmtbWRFTest():
    def __init__(self, testcase, gmtbdir, wrfcfg, workdir, wrfmain):
        """
        Initialise paths needed by the test run

        Sets up a new WRF configuration based on the files in 'wrfcfg'

        Args:
            gmtbdir: Path to gmtb repository, used to find test case inputs
            wrfcfg: Path to WRF configuration
            workdir: Directory to run WRF in
            wrfmain: Parent path of ideal.exe and wrf.exe
        """
        self.testcase = testcase
        self.gmtbdir = Path(gmtbdir)
        self.wrfcfg = Path(wrfcfg)
        self.workdir = Path(workdir)
        self.wrfmain = Path(wrfmain)


    def setup_workdir(self):
        """
        Create the run directory and link files from the config directory

        The namelist is handled by 'setup_namelist()'
        """
        self.workdir.mkdir()
        for c in self.wrfcfg.iterdir():
            if c.name != 'namelist.input':
                (self.workdir / c.name).symlink_to(c.resolve())


    def setup_config(self):
        """
        Top level setup for a test case

        Reads the GMTB test configuration and input files, then sets up the WRF
        run appropriately
        """
        casefile = self.gmtbdir / f'scm/etc/case_config/{self.testcase}.nml'
        casenml = f90nml.read(casefile)

        start = pandas.Timestamp(
                year=casenml['case_config']['year'],
                month=casenml['case_config']['month'],
                day=casenml['case_config']['day'],
                hour=casenml['case_config']['hour'])

        datadir = self.gmtbdir / 'scm/etc' / casenml['case_config']['case_data_dir']
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

        zsfc = initial.height[0]
        self.setup_forcing(zsfc, forcing)


    def initial_sounding(self, initial, forcing):
        """
        Setup the WRF initial sounding

        The first height level in the forcing dataset is the surface height

        Theta and qv are found following equations 5.1 and 5.2 of the GMTB tech
        guide

        Args:
            inital: GMTB initial sounding
            forcing: GMTB forcing
        """
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
        zsfc = initial.height.isel(levels=0).data
        usfc = sounding.u.isel(levels=0).data
        vsfc = sounding.v.isel(levels=0).data
        qvsfc = sounding.qv.isel(levels=0).data

        with open(self.workdir / 'input_sounding', 'w') as f:
            print(zsfc, usfc, vsfc, Tsfc, qvsfc, psfc, file=f)
            sounding.isel(levels=slice(1,None)).to_dataframe().to_csv(f, columns=['z','u','v','theta','qv'], sep=' ', index=False, header=False)


    def setup_namelist(self, casenml, start, lat, lon):
        """
        Setup the WRF namelist to perform a GMTB run

        Gets the start date and duration from the GMTB namelist
        lat,lon location must be read from the forcing file

        Args:
            casenml (f90nml): GMTB namelist
            start (pandas.Timestamp): Start time
            lat, lon: Model location
        """
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

        wrfnml['time_control']['run_days'] = runtime.components.days
        wrfnml['time_control']['run_hours'] = runtime.components.hours
        wrfnml['time_control']['run_minutes'] = runtime.components.minutes
        wrfnml['time_control']['run_seconds'] = runtime.components.seconds

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


    def setup_forcing(self, zsfc, forcing):
        """
        Setup SCM forcing file by converting GMTB values

        Args:
            zsfc: Surface height (used to convert GMTB pressure levels to
                    height)
            forcing: GMTB forcing
        """
        # Convert pressure levels to heights
        z = pressure_to_z(zsfc, forcing.T_nudge.T, forcing.levels)

        # Setup WRF times
        Times = numpy.datetime_as_string(forcing.time, unit='s')
        Times = [t.replace('T','_') for t in Times]

        # Empty array
        zero = (['time','levels'], dask.array.zeros((forcing.time.size, forcing.levels.size), dtype='f4'))

        wrf_forcing = xarray.Dataset({
                'Times': (['time'], Times),
                'Z_FORCE': (['time', 'levels'], z),
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

        # Copy metadata from sample forcing file
        wrf_ideal = xarray.open_dataset(self.workdir / 'force_ideal.nc')
        for k, v in wrf_forcing.items():
            if k in wrf_ideal:
                v.attrs = wrf_ideal[k].attrs
        wrf_forcing.attrs = wrf_ideal.attrs

        # Save to file
        wrf_forcing.to_netcdf(self.workdir / 'force_test.nc')

        # Update namelist with forcing details
        interval = pandas.Timedelta((forcing.time[1] - forcing.time[0]).values).components

        wrfnml = f90nml.read(self.workdir / 'namelist.input')
        wrfnml['time_control']['auxinput3_inname'] = 'force_test.nc'
        wrfnml['time_control']['auxinput3_interval_d'] = interval.days
        wrfnml['time_control']['auxinput3_interval_h'] = interval.hours
        wrfnml['time_control']['auxinput3_interval_m'] = interval.minutes
        wrfnml['time_control']['auxinput3_interval_s'] = interval.seconds
        wrfnml.write(self.workdir / 'namelist.input', force=True)


    def run_testcase(self):
        """
        Setup and run the testcase
        """
        self.setup_workdir()
        self.setup_config()

        subprocess.run([self.wrfmain.resolve() / 'ideal.exe'], cwd=self.workdir)
        subprocess.run([self.wrfmain.resolve() / 'wrf.exe'], cwd=self.workdir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gmtb-repo', help='Path to local GMTB repository', required=True)
    parser.add_argument('--wrf-main', help='Path containing ideal.exe and wrf.exe', required=True)
    parser.add_argument('testcase', help='Testcase to run')
    args = parser.parse_args()

    ts = pandas.Timestamp.utcnow().strftime('%Y%m%dT%H%M%S')
    
    test = gmtbWRFTest(
            testcase = args.testcase,
            gmtbdir = args.gmtb_repo,
            wrfcfg = 'wrf_in',
            workdir = f'test/{args.testcase}-{ts}',
            wrfmain = args.wrf_main,
            )
    test.run_testcase()
