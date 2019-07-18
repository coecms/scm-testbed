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
import pandas
import numpy
import matplotlib.pyplot as plt
import stratify
from wrf.constants import Constants as wrf_constants

def compare_variable(gmtb, wrf, ax):

    vmax = max(gmtb.max(), wrf.max())
    vmin = min(gmtb.min(), wrf.min())
    if vmin < 0:
        vmax = max(-vmin, vmax)
        vmin = None

    #fig, ax = plt.subplots(1,3, sharey=True)
    gmtb.plot.pcolormesh('time','levels', ax=ax[0], vmin=vmin, vmax=vmax)

    _, t = numpy.meshgrid(range(wrf.shape[1]), wrf.XTIME)
    #ax[1].pcolormesh(t, wrf.PRES, wrf)
    wrf.coords['time2'] = (wrf.dims, t)
    wrf.plot.pcolormesh('time2', 'PRES', ax=ax[1], vmin=vmin, vmax=vmax)

    wrf_on_gmtb = stratify.interpolate(gmtb.levels, wrf.PRES, wrf, axis=1).T[:,::3]
    (gmtb - wrf_on_gmtb).plot.pcolormesh('time','levels', ax=ax[2])

    ax[0].set_ylim([100000,0])



def compare_models(gmtb, wrf):
    variables = {
        'u': ['u_nudge', 'U'],
        'v': ['v_nudge', 'V'],
        'qv': ['qt_nudge', 'QVAPOR'],
        'theta': ['thil_nudge', 'THETA'],
        'w': ['w_ls', 'W_G'],
    }

    fig, ax = plt.subplots(len(variables), 3, sharey='row', sharex=True)

    for i, (name, v) in enumerate(variables.items()):
        print(i)
        a = gmtb[v[0]]
        b = wrf[v[1]]

        compare_variable(a, b, ax[i,:])

    ax[0,0].set_title('GMTB Forcing')
    ax[0,1].set_title('WRF Output')
    ax[0,2].set_title('Difference')

    plt.show()

def main():
    wrf = xarray.open_dataset('output/twpice-latest/wrfout_d01_2006-01-19_03:00:00')
    wrf = wrf.isel(south_north=0, west_east=0, south_north_stag=0, west_east_stag=0)
    wrf.coords['PRES'] = wrf['P'] + wrf['PB']
    wrf['THETA'] = wrf['T'] + wrf_constants.T_BASE

    w_on_grid = stratify.interpolate(wrf.ZNU, wrf.ZNW, wrf.W, axis=1)
    wrf['W_G'] = (wrf.PRES.dims, w_on_grid)

    gmtb = xarray.open_dataset('../gmtb-scm-release/scm/data/processed_case_input/twpice.nc', 'forcing')
    gmtb.update(xarray.open_dataset('../gmtb-scm-release/scm/data/processed_case_input/twpice.nc'))
    gmtb.coords['time'] = wrf['XTIME'].values[0] + pandas.to_timedelta(gmtb['time'].values, 'seconds')

    compare_models(gmtb, wrf)

if __name__ == '__main__':
    main()
