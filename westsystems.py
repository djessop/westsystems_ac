#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# read_westsystems.py
# Created by: D E Jessop, DEC 2025
# Last modification: 2026-10-21
# TO DO:
#     - Read all metadata

from dataclasses import dataclass
from datetime import datetime as dt
from pathlib import Path

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import utm

M_CO2    = 44.009e-3    # molar_mass CO2/[kg/mol]
M_H2S    = 34.082e-3    # molar_mass H2S/[kg/mol]
M_SO2    = 64.066e-3    # molar_mass SO2/[kg/mol]
s_in_day = 3600 * 24    # seconds in a day/[s/day]

wsf_timestamp_fmt = '%d-%m-%Y %H:%M:%S'
columns = '''datetime,UTM ZONE,UTM LONGITUDE,UTM LATITUDE,ELEVATION,
PRESSURE [hPa],AIR TEMPERATURE [degC],AIR RELATIVE HUMIDITY [%],
ACCUMULATION CHAMBER,ACK,CO2_LLIMIT [sec],CO2_RLIMIT [sec],CO2_LCONC [ppm],
CO2_RCONC [ppm],CO2_SLOPE [ppm/s],CO2_SLOPE_ABS_ERR [ppm/s],
CO2_R^2,CO2_FLUX [mol/m2/day]'''.replace('\n', '').split(',')


# ---------------------------------------------------------------------------
# Result of linear fit over chosen time window
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class RangeFit:
    species       : str            # 'CO2' or 'H2S'
    llimit        : int            # window start [s]
    rlimit        : int            # window end [s]
    lconc         : float          # concentration at window start [ppm]
    rconc         : float          # concentration at window end [ppm]
    slope         : float          # [ppm/s]
    intercept     : float          # [ppm]
    slope_abs_err : float          # standard error of the slope [ppm/s]
    r2            : float
    flux          : float          # ACK * slope [mol/m2/day]
    n_points      : int

    def predict(self, x) -> np.ndarray:
        return self.intercept + self.slope * np.asarray(x, dtype=float)

    def to_row(self) -> dict:
        """Column names match the original output format."""
        s = self.species
        return {f'{s}_LLIMIT [sec]':          self.llimit,
                f'{s}_RLIMIT [sec]':          self.rlimit,
                f'{s}_LCONC [ppm]':           self.lconc,
                f'{s}_RCONC [ppm]':           self.rconc,
                f'{s}_SLOPE [ppm/s]':         self.slope,
                f'{s}_SLOPE_ABS_ERR [ppm/s]': self.slope_abs_err,
                f'{s}_R^2':                   self.r2,
                f'{s}_FLUX [mol/m2/day]':     self.flux}


# ---------------------------------------------------------------------------
# Pure functions: no plotting, no state, easy to test and batch-run
# ---------------------------------------------------------------------------

def clip_limits(data: pd.DataFrame, left: float, right: float):
    """Round to whole seconds, order them, and clip to the record length."""
    left, right = sorted((left, right))
    left  = max(0, int(round(left)))
    right = min(int(data['sec'].max()), int(round(right)))
    return left, right


def subset(data: pd.DataFrame, left: int, right: int) -> pd.DataFrame:
    return data[(data['sec'] >= left) & (data['sec'] <= right)]


def fit_range(data: pd.DataFrame, species: str,
              left: int, right: int, ack: float) -> RangeFit:
    """Ordinary least-squares line of concentration vs time on [left, right]."""
    sub = subset(data, left, right)
    n   = len(sub)
    if n < 3:
        raise ValueError(f'need at least 3 points in [{left}, {right}] s, '
                         f'got {n}')

    x = sub['sec'].to_numpy(dtype=float)
    y = sub[species].to_numpy(dtype=float)
    X = np.column_stack([np.ones_like(x), x])

    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid    = y - X @ beta
    ss_res   = float(resid @ resid)
    ss_tot   = float(((y - y.mean()) ** 2).sum())
    r2       = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    # standard error of the slope: sqrt(sigma^2 * [(X'X)^-1]_11)
    sigma2    = ss_res / (n - 2)
    slope_err = float(np.sqrt(sigma2 * np.linalg.inv(X.T @ X)[1, 1]))

    return RangeFit(species=species,
                    llimit=left, rlimit=right,
                    lconc=float(y[0]), rconc=float(y[-1]),
                    slope=float(beta[1]), intercept=float(beta[0]),
                    slope_abs_err=slope_err, r2=r2,
                    flux=float(ack * beta[1]),
                    n_points=n)


# ---------------------------------------------------------------------------
# File handling and interactive range selection
# ---------------------------------------------------------------------------

class WestsystemsFile:
    """Parse and analyse data from a Westsystems accumulation chamber.

    If neither `man_lims` nor `validate` is given, an interactive window opens
    so the fitting range can be picked with two mouse clicks.

    man_lims : (left, right) in seconds, used as the first guess (interactive)
               or as the final range (if validate=True).
    validate : if True, skip the interactive step.  Uses `man_lims`, or the
               whole record if `man_lims` is None.
    """

    def __init__(self, filename, gas_species='CO2',
                 ac_chamber=None, man_lims=None, validate=False):

        path = Path(filename)
        self.path      = path if path.is_absolute() else Path.cwd() / path
        self.pathname  = str(self.path.parent)
        self.filename  = self.path.name

        self.ac_chamber  = ac_chamber
        self.gas_species = gas_species
        self.man_lims    = man_lims
        self._validate   = validate
        self.fit         = None     # RangeFit, set once a range is chosen
        self.subdata     = None
        self.df          = None

        self._parse_file()

        # Keys are matched by prefix so a garbled unit such as "(°C)" cannot
        # break the lookup, whatever the file's encoding.
        self.elevation = self._meta_float('ELEVATION')
        self.longitude = self._meta_float('LONGITUDE')
        self.latitude  = self._meta_float('LATITUDE')
        self.T         = self._meta_float('TEMPERATURE') + 273.15   # K
        self.p         = self._meta_float('PRESSURE') * 100         # Pa
        self.RH        = None
        self.datetime  = dt.strptime(self.meta['TIME'], wsf_timestamp_fmt)

        site, fm, *_   = self.filename.split('_')
        self.sitename  = site
        self.fluxmeas  = fm

        self.ACK = calc_ack(self.p, self.T, self.ac_chamber)

        if man_lims is not None or validate:
            lims = man_lims if man_lims is not None \
                else (0, int(self.data['sec'].max()))
            self.set_range(*lims)
        if not validate:
            self.select_range()

    # -- parsing -------------------------------------------------------------

    def _meta_float(self, prefix: str) -> float:
        """Value of the metadata key whose name, ignoring any "(unit)" suffix,
        is `prefix`: "TEMPERATURE (°C)" -> "TEMPERATURE"."""
        try:
            for key, val in self.meta.items():
                if key.split('(')[0].strip().upper() == prefix:
                    return float(val)
        except ValueError:
            return np.nan
        # finally:
        #     raise KeyError(f'no metadata key named {prefix!r}; '
        #                f'found {list(self.meta)}')

    def _parse_file(self):
        meta = {}
        rows = []

        # Decode explicitly: the default encoding differs between Windows
        # (cp1252) and Linux/macOS (UTF-8).  Try UTF-8 first, then cp1252.
        try:
            text = self.path.read_text(encoding='utf-8-sig')
        except UnicodeDecodeError:
            text = self.path.read_text(encoding='cp1252')
        lines = iter(text.splitlines())

        meta['header'] = next(lines).strip()

        # metadata: "key: value" lines up to the first SENSOR# line
        for raw in lines:
            line = raw.strip()
            if 'SENSOR#' in line:
                break
            if line == '' or ':' not in line:
                continue
            key, val = line.split(':', 1)
            meta[key.strip()] = val.strip()

        # skip to the data block, then collect it
        read_data = False
        for raw in lines:
            read_data = read_data or ('FLUX RECORD TRACKS' in raw)
            if read_data:
                rows.append(raw.strip('\t\r\n').split('\t'))

        # rows[0] = "FLUX RECORD TRACKS" line, rows[1] = column names
        data = pd.DataFrame(rows[2:], columns=rows[1]).astype(float)
        data['sec'] = data['sec'].astype(int)

        self.data = data
        self.meta = meta

    # -- fitting -------------------------------------------------------------

    def set_range(self, left, right) -> RangeFit:
        """Fit the chosen window and refresh the output DataFrame."""
        left, right = clip_limits(self.data, left, right)
        self.subdata = subset(self.data, left, right)
        self.fit     = fit_range(self.data, self.gas_species,
                                 left, right, self.ACK)
        self._fmt_data()
        return self.fit

    def _fmt_data(self):
        if not np.isnan(self.latitude) or not np.isnan(self.longitude):
            utm_latlon = utm.from_latlon(self.latitude, self.longitude)
        else:
            utm_latlon = [np.nan, np.nan, np.nan]

        row = {'datetime':                   self.datetime,
               'UTM ZONE':                   utm_latlon[2],
               'UTM LONGITUDE':              utm_latlon[0],
               'UTM LATITUDE':               utm_latlon[1],
               'ELEVATION':                  self.elevation,
               'PRESSURE [hPa]':             self.p / 100.0,
               'AIR TEMPERATURE [degC]':     self.T - 273.15,
               'AIR RELATIVE HUMIDITY [%]':  self.RH,
               'ACCUMULATION CHAMBER':       self.ac_chamber,
               'ACK':                        self.ACK,
               **self.fit.to_row(),
               'FILENAME':                   str(self.path)}
        self.df = pd.Series(row).to_frame().T

    # -- interactive selection ----------------------------------------------

    def select_range(self):
        species = self.gas_species
        plt.close('all')
        fig, ax = plt.subplots()
        lims    = self.man_lims      # used for the first pass only
        validate = self._validate

        while not validate:
            ax.cla()
            ax.set_title(self.filename, fontsize='large')
            ax.scatter(self.data['sec'], self.data[species], c='C0')
            ax.set_xlabel('Time/[s]', fontsize='large')
            ax.set_ylabel(f'{species} concentration/[ppm]', fontsize='large')
            xlims, ylims = ax.get_xlim(), ax.get_ylim()

            if lims is None:
                pts = plt.ginput(2, timeout=0)
                if len(pts) < 2:
                    continue
                left, right = pts[0][0], pts[1][0]
            else:
                left, right = lims
                lims = None

            try:
                fit = self.set_range(left, right)
            except ValueError as err:
                print(err)
                continue

            ax.axvspan(fit.llimit, fit.rlimit, color='r', alpha=.2)
            ax.set_xlim(xlims)
            ax.set_ylim(ylims)

            ax.scatter(self.subdata['sec'], self.subdata[species], c='r')
            xs = self.subdata['sec'].to_numpy()
            ax.plot(xs, fit.predict(xs), '-k', lw=.75)

            text = (f'llimit: {fit.llimit:>3d} s\n'
                    f'rlimit: {fit.rlimit:>3d} s\n'
                    f'slope: {fit.slope:8.4f} ppm/s\n'
                    f'R2 : {fit.r2:11.4f}')
            ax.annotate(text, (.65, .05), xycoords='axes fraction',
                        fontsize='large',
                        bbox=dict(boxstyle='square,pad=0.3', fc='w', ec='k'))

            plt.pause(0.01)
            plt.show(block=False)

            inp      = input('Validate range (y/n)? ')
            validate = validate or not inp.strip().lower().startswith('n')
                #break

        self._validate = validate  #should be True
        self._fig, self._ax = fig, ax
        out = self.path.with_name(f'{self.path.stem}_{species}.png')
        fig.savefig(out, dpi=150)
        plt.close()


def calc_ack(p: float, T: float, ac_chamber: str='B') -> float:
    '''
    Calculates K, the constant of proportionality used to convert flux in
    ppm/s (CO2/H2S_SLOPE) to molar flux in mol/m2/day (CO2/H2S_FLUX)

    Output has been verified against WestSystems datasheet (see:
    https://www.westsystems.com/instruments/wp-content/uploads/2019/01/
    Handbook_Portable_9.1.pdf)
    
    Parameters
    ----------
    p: float or np.ndarray
        Ambient pressure, in Pa
    T: float or np.ndarray
        Air temperature, in K
    ac_chamber: str
        Type of accumulation chamber.  Must be one of 'A', 'B' or 'C'
    Note: if p and T are array_like, they must have the same dimension

    Returns
    -------
    ACK: floar or np.ndarray
        The K coefficient for conversion of fluxes from ppm/s to mol/m2/day
    '''
    from scipy.constants import R           #8.314 J/(mol K)

    # dict of accumulatation chamber data: inlet dia., height, vol. area.
    # Data from:
    # https://www.portablefluxmeter.com/#1499955177827-2d0b10ba-c11c
    ac_chamber_dict = {'A': (.2, .1, .002770, .0308),
                       'B': (.2, .2, .006186, .0317),
                       'C': (.3, .1, .006878, .0712)}

    if ac_chamber is None:
        return None
    # if ac_chamber != 'A' or ac_chamber != 'B' or ac_chamber != 'C':
    #     raise ValueError("'ac_chamber' must be one of 'A', 'B' or 'C'")
    try: 
        _, height, vol, area = ac_chamber_dict[ac_chamber]
    except KeyError:
        raise Warning('AC chamber type unknown.  Reverting to "B"')
        _, height, vol, area = ac_chamber_dict['B']
    secs_in_day = 24. * 3600
    return secs_in_day * p * vol / (R * T * area * 1_000_000)


## Run on all files in the given path
def batch_run(path: str | Path, gas_species: str='CO2', ac_chamber: str ='B',
              outfile: str='database.csv') -> pd.DataFrame:
    from glob import glob
    from numpy.linalg import LinAlgError
    from pathlib import Path

    import pandas as pd

    # Check that path is a pathlib.Path
    if isinstance(path, str):
        path = Path(path)

    if gas_species == 'H2S':
        for ind, col in enumerate(columns):
            columns[ind] = col.replace('CO2', 'H2S')

    filenames = sorted(path.glob('*.txt'))

    print('Found the following files:')
    print([filename.name for filename in filenames])

    df = pd.DataFrame(columns=columns)

    for filename in filenames:
        try:
            ws_data = WestsystemsFile(filename, gas_species, ac_chamber)
            df = pd.concat([df, ws_data.df])
            df.to_csv(outfile, index=False)
        except LinAlgError:
            pass

    return df.reset_index(drop=True)

if __name__ == '__main__':
    from pathlib import Path

    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description='Optional arguments for batch processing')
    parser.add_argument('path')
    parser.add_argument('-o', '--output', default='output.csv',
                        help=('Location of the output file'))
    parser.add_argument('-s', '--species', default='CO2',
                        help=('Gas species to process'))
    parser.add_argument('-c', '--chamber', default='B',
                        help=('Type of AC chamber used'))

    args = parser.parse_args()

    path        = Path(args.path)
    gas_species = args.species
    ac_chamber  = args.chamber
    output_name = args.output

    print(f'Looking for files in folder {path}...\nRunning for {gas_species} '
          + f'using AC chamber {ac_chamber}')
    path = Path(path)

    if path.is_dir():
        df = batch_run(path, gas_species, ac_chamber, output_name)
    else:
        raise FileNotFoundError(f"Path {path} not found")
