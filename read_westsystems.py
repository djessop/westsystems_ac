#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# read_westsystems.py
# Created by: D E Jessop, DEC 2025
# Last modification: 2026-01-23


from scipy.constants import R           #8.314 J/(mol K)
from datetime import datetime as dt

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


T0 = 273.15              # [K]
p0 = 101325.             # [bar]
molar_vol = R * T0 / p0  # vol occupied by 1 mol at STP ~ 0.0224/[m3/mol] 
M_CO2 = 44.009e-3        # molar_mass CO2/[kg/mol]
M_H2S = 34.082e-3        # molar_mass H2S/[kg/mol]
M_SO2 = 64.066e-3        # molar_mass SO2/[kg/mol]

# dict of accumulatation chamber data: inlet dia., height, vol. area.
# Data from:
# https://www.portablefluxmeter.com/#1499955177827-2d0b10ba-c11c
ac_chamber_dict = {'A': (.2, .1, .002770, .0308),
                   'B': (.2, .2, .006186, .0317),
                   'C': (.3, .1, .006878, .0712)}

wsf_timestamp_fmt = '%d-%m-%Y %H:%M:%S'

class WestsystemsFile:
    """A class providing methods for parsing and analysing data from 
    Westsystems accumulation chamber"""
    def __init__(self, filename, gas_species='CO2',
                 ac_chamber=None, man_lims=None, validate=False,
                 target_file='database.csv', source_file='database.csv'):
        *pathname, filename = filename.split('/')
        self.filename     = filename
        self.pathname     = '/'.join(pathname) + '/'
        self.ac_chamber   = ac_chamber
        self.gas_species  = gas_species
        self.target_file  = target_file
        self.source_file  = source_file
        self.CO2_SLOPE    = None
        self.H2S_SLOPE    = None
        self.CO2_FLUX     = None
        self.H2S_FLUX     = None
        self._beta        = None
        self._X           = None
        self._validate    = validate
        self.man_lims     = man_lims
        self._parse_file()
        self.elevation    = float(self.meta['ELEVATION (m)'])
        self.longitude    = float(self.meta['LONGITUDE'])
        self.latitude     = float(self.meta['LATITUDE'])
        self.T  	  = float(self.meta['TEMPERATURE (°C)']) + 273.15
        self.p  	  = float(self.meta['PRESSURE (HPa)']) * 100
        self.RH 	  = None,
        self.ACK          = self.calc_ack()
        self.datetime     = dt.strptime(self.meta['TIME'], wsf_timestamp_fmt)
        site, fm, *d_t    = filename.split('_')
        self.sitename     = site
        self.fluxmeas     = fm

        self._select_range()

    def _parse_file(self):
        fname = self.pathname + self.filename
        with open(fname, 'r') as f:
            data 	       = []
            meta 	       = {}
            meta['header'] = f.readline().strip('\n')
            linesread      = 1
            read_data      = False
            nest_depth     = 0
            ##--------------------READ METADATA--------------------##
            # Reading of metadata currently stops at the first line containing
            # "SENSOR#"
            try:
                while not read_data:
                    line = f.readline().strip()
                    linesread += 1
                    if line == '':
                        continue
                    elif 'SENSOR#' in line:
                        meta[line] = {}
                        nest_depth += 1  
                    else:
                        key, val  = line.split(':', 1)
                        meta[key] = val.strip()
                        
                    if nest_depth == 1:
                        key, val  = line.split(':', 1)
                        meta[head][key] = val.strip()
                        
                    read_data = read_data or ("FLUX RECORD TRACKS" in line) \
                        or ('SENSOR#' in line)
            except ValueError:
                pass
        
            # # Find current position in file
            # cur = f.tell()
            # f.seek(0, 2)
            # eof = f.tell()
            # f.seek(cur)
        
            # while f.tell() != eof:
            #     line = f.readline()
            for line in f:
                read_data = read_data or ("FLUX RECORD TRACKS" in line)
                if read_data:
                    data.append(line.strip('\t\n').split('\t'))
                    # else:
                    #     if line != '\n':
                    #         key, val = line.strip('\n').split(':')
                    #         meta[key] = val
        
            #     f.seek(0)
        
        data = pd.DataFrame(data[2:], columns=[data[1]]).astype(float)
        data['sec'] = data['sec'].astype(int)
        
        if isinstance(data.columns, pd.MultiIndex):
            data.columns = [x[0] for x in data.columns]

        ## NOT ALL WS DATA FILES HAVE THESE FIELDS
        # data['T'] = data['AirT'] + 273.15  # convert from °C to K
        # data['p'] = data['BarP'] * 100     # convert from hPa to Pa
    
        # return data, meta #data.drop('#', axis=1).set_index('sec'), meta
        self.data = data
        self.meta = meta
    
    def _fmt_data(self):
        import os
        import pandas as pd
        import utm

        if self.target_file == '':
            self.target_file = self.source_file

        # calculate molar flux
        self._ppm_to_molar_flux()

        utm_latlon  = utm.from_latlon(self.latitude, self.longitude)
        gas_species = self.gas_species

        ## There's probably a more intellegent way of doing this
        df_dict = {'datetime':                  self.datetime, 
                   'UTM ZONE':                  utm_latlon[2],
                   'UTM LONGITUDE':             utm_latlon[0],
                   'UTM LATITUDE':              utm_latlon[1],
                   'ELEVATION':                 self.elevation,
                   'PRESSURE [hPa]':            self.p,
                   'AIR TEMPERATURE [degC]':    self.T,
                   'AIR RELATIVE HUMIDITY [%]': None,
                   'ACCUMULATION CHAMBER':      self.ac_chamber,
                   'ACK':                       self.ACK,
                   'CO2_LLIMIT [sec]':  	self.CO2_LLIMIT,
                   'CO2_RLIMIT [sec]':  	self.CO2_RLIMIT,
                   'CO2_LCONC [ppm]':   	self.CO2_LCONC,
                   'CO2_RCONC [ppm]':   	self.CO2_RCONC,
                   'CO2_SLOPE [ppm/s]': 	self.CO2_SLOPE,
                   'CO2_SLOPE_ABS_ERR [ppm/s]': self.CO2_SLOPE_ABS_ERR,
                   'CO2_R^2':                   self.CO2_R2,
                   'CO2_FLUX [mol/m2/day]':     self.CO2_FLUX}

        if gas_species == 'H2S':
            #columns = columns.replace('CO2', 'H2S')

            df_dict = {'datetime': 		    self.datetime, 
                       'UTM ZONE': 		    utm_latlon[2],
                       'UTM LONGITUDE':             utm_latlon[0],
                       'UTM LATITUDE':              utm_latlon[1],
                       'ELEVATION':                 self.elevation,
                       'PRESSURE [hPa]':            self.p,
                       'AIR TEMPERATURE [degC]':    self.T,
                       'AIR RELATIVE HUMIDITY [%]': None,
                       'ACCUMULATION CHAMBER':      self.ac_chamber,
                       'ACK':                       self.ACK,
                       'H2S_LLIMIT [sec]':  	    self.H2S_LLIMIT,
                       'H2S_RLIMIT [sec]':  	    self.H2S_RLIMIT,
                       'H2S_LCONC [ppm]':   	    self.H2S_LCONC,
                       'H2S_RCONC [ppm]':   	    self.H2S_RCONC,
                       'H2S_SLOPE [ppm/s]': 	    self.H2S_SLOPE,
                       'H2S_SLOPE_ABS_ERR [ppm/s]': self.H2S_SLOPE_ABS_ERR,
                       'H2S_R^2':                   self.H2S_R2,
                       'H2S_FLUX [mol/m2/day]':     self.H2S_FLUX}

        # Check if data file exists.  If not, create one.
        if os.path.isfile(self.source_file):
            df = pd.read_csv(self.source_file)
        else:
            df = pd.DataFrame(columns=columns)
        
        self.df = pd.Series(df_dict).to_frame().T
        #df = pd.concat([df, pd.Series(df_dict).to_frame().T])
        #self.df = df
        #df.to_csv(self.target_file, index=False)

    def calc_ack(self):
        if self.ac_chamber is None:
            return None
        _, height, vol, area = ac_chamber_dict[self.ac_chamber]
        secs_in_day = 24. * 3600
        return secs_in_day * self.p * vol / (R * self.T * area * 1_000_000)

    def _ppm_to_molar_flux(self):
        if self.gas_species == 'CO2':
            self.CO2_FLUX = self.ACK * self.CO2_SLOPE
        if self.gas_species == 'H2S':
            self.H2S_FLUX = self.ACK * self.H2S_SLOPE

    def select_range(self):
        self._select_range()

    def _select_range(self):
        import matplotlib.pyplot as plt
        
        plt.close('all')
        fig, ax = plt.subplots()
        if self._validate:
            plt.close()

        while not self._validate: # or man_lims is not None:
            ax.cla()
            ax.set_title(self.filename, fontsize='large')
            gas_species = self.gas_species
    
            sca = plt.scatter(data=self.data, x='sec', y=gas_species, c='C0')
            ax.set_xlabel('Time/[s]', fontsize='large')
            ax.set_ylabel(f'{gas_species} concentration/[ppm]',
                          fontsize='large')
    
            xlims = ax.get_xlim()
            ylims = ax.get_ylim()

            if self.man_lims is None:
                left, right = plt.ginput(2, timeout=0)
            else:
                left  = [self.man_lims[0]]
                right = [self.man_lims[1]]
            left[0]  = max(0, np.round(left[0]))
            right[0] = min(self.data['sec'].max(), np.round(right[0]))
    
            spa = ax.axvspan(left[0], right[0], color='r', alpha=.2)

            ax.set_xlim(xlims)
            ax.set_ylim(ylims)

            self.subdata = self.data.query('@left[0] <= sec <= @right[0]')
            self._plot_model()
    
            rsc = plt.scatter(data=self.subdata, x='sec', y=gas_species, c='r')
    
            # text box with fitting data
            text = f'llimit: {left[0]:>3d} s\n'\
                +  f'rlimit: {right[0]:>3d} s\n'\
                +  f'slope: {eval(f"self.{gas_species}_SLOPE"):8.4f} ppm/s\n'\
                +  f'R2 : {eval(f"self.{gas_species}_R2"):11.4f}'
            
            ax.annotate(text, (.65, .05), xycoords='axes fraction',
                        fontsize='large',
                        bbox=dict(boxstyle="square,pad=0.3", fc="w", ec="k"))

            plt.pause(0.01)
            plt.show(block=False)
            
            self._validate = self._validate or (
                input('Validate range (y/n)? ').lower() == 'y')

            if gas_species == 'CO2':
                self.CO2_LLIMIT = left[0]
                self.CO2_RLIMIT = right[0]
                self.CO2_LCONC  = self.subdata['CO2'].iloc[0]
                self.CO2_RCONC  = self.subdata['CO2'].iloc[-1]
            if gas_species == 'H2S':
                self.H2S_LLIMIT = left[0]
                self.H2S_RLIMIT = right[0]
                self.H2S_LCONC  = self.subdata['H2S'].iloc[0]
                self.H2S_RCONC  = self.subdata['H2S'].iloc[-1]
                
            self._fmt_data()

        self._fig = fig
        self._ax  = ax
        
        fig.savefig(self.pathname + self.filename.split('.')[0] + '.png',
                    dpi=150)

    def _plot_model(self):
        self._lm()
        plt.plot(self._X.T[1], self._X@self._beta, '-k', lw=.75)

    def _lm(self):
        '''linear model'''
        from sklearn.metrics import r2_score

        data = self.subdata
        gas_species = self.gas_species
        x, y = data[['sec', gas_species]].values.T

        X = np.array([x**i for i in range(2)]).T  ## linear model
        self._beta, *foo = np.linalg.lstsq(X, y, rcond=-1)
        self._X = X
        self._intercept = self._beta[0]
        yhat = X@self._beta

        cov  = np.linalg.inv(X.T @ X)
        
        if gas_species == 'CO2':
            self.CO2_R2 = r2_score(y, yhat)
            self.CO2_SLOPE_ABS_ERR = np.sqrt(np.diag(cov))[1]
        else:
            self.H2S_R2 = r2_score(y, yhat)
            self.H2S_SLOPE_ABS_ERR = np.sqrt(np.diag(cov))[1]

        if gas_species == 'CO2':
            self.CO2_SLOPE = self._beta[1]
        if gas_species == 'H2S':
            self.H2S_SLOPE = self._beta[1]
        # print(f'{gas_species} flux: {eval(f"self.{gas_species}_SLOPE"):6.3f}'
        #       + ' ppm/s')
        

if __name__ == '__main__':
    import sys

    path = '/home/david/FieldWork/Soufriere-Guadeloupe/soil_degassing/data/' \
        + 'raw/'
    gas_species = 'CO2'
    if len(sys.argv) > 1:
        path        = sys.argv[1]
    if len(sys.argv) > 2:
        path        = sys.argv[1]
        gas_species = sys.argv[2]

    columns = '''datetime,UTM ZONE,UTM LONGITUDE,UTM LATITUDE,ELEVATION,
    PRESSURE [hPa],AIR TEMPERATURE [degC],AIR RELATIVE HUMIDITY [%],
    ACCUMULATION CHAMBER,ACK,CO2_LLIMIT [sec],CO2_RLIMIT [sec],
    CO2_LCONC [ppm],CO2_RCONC [ppm],CO2_SLOPE [ppm/s],
    CO2_SLOPE_ABS_ERR [ppm/s],CO2_r^2,
    CO2_FLUX [mol/m2/day]'''.replace('\n', '')

    if gas_species == 'H2S':
        columns = columns.replace('CO2', 'H2S')

    filenames = sorted(glob(path + '2025-*/**.txt'))

    print(filenames)

    df = pd.DataFrame(columns=columns)

    for filename in filenames:
        ws_data = WestsystemsFile(fname, gas_species, ac_chamber='B')
        df = pd.concat([df, ws_data.df])
