#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 11:30:51 2020

@author: kujawski
"""

import pandas as pd
import csv
import numpy as np
from os.path import join

calibpath = '.'
calibfile = '2021-05-05_calibration.csv'
calibdelimiter = ' '
serial_nbr = []
calib_value = []
calib_factor = []


def to_pa(level):
    return (10**(level/10))*(4e-10)


calibvalues = pd.read_csv(join(calibpath, calibfile),
                          header=None, delimiter=calibdelimiter)
calibfactors = to_pa(calibvalues[0])

# =============================================================================
# create file(s)
# =============================================================================

calib_id = "2021-05-05"
name = '2021-05-05_calibration'
filename = join(calibpath, name+'.xml')

with open(filename, 'w') as f:
    f.write(
        f'<?xml version="1.0" encoding="utf-8"?>\n<Calib name="{calib_id}">\n')
    for i in range(len(calibvalues)):
        channel_string = 'AnalogInput_'+"%03d" % (i+1)
        fac = calibfactors[i]
        f.write(f'	<pos Name="Point {channel_string}" factor="{fac}"/>\n')
    f.write('</Calib>')
