# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 15:23:10 2018

@author: Arron
"""

#!/usr/bin/env python

import sys
import os
from pandas import *
from pandasql import sqldf


def print_full(x):
    set_option('display.max_rows', len(x))
    print(x)
    reset_option('display.max_rows')

set_option('display.expand_frame_repr', False)

#snp = read_csv("C:/Users/arron/Documents/Phd/scripts/humvar3/humvar3_ss.csv",header=None)
#wild = read_csv("C:/Users/arron/Documents/Phd/scripts/humvarids/final.csv",header=None)

snp = read_csv("C:/Users/arron/Documents/Phd/scripts/epsnps/psipred.csv",header=None)
wild = read_csv("C:/Users/arron/Documents/Phd/scripts/epsnps/epidsmaster.csv",header=None)



snp.columns = ["rubbish","mut","ss","psi1","psi2","psi3","prot","pos","posNA"]
wild.columns = ["pos","wild","ss","psi1","psi2","psi3","prot","rubbish"]


merged = snp.merge(wild, how ='inner', left_on=['pos','prot'],right_on=['pos','prot'])
merged['diff1'] = merged.psi1_x - merged.psi1_y
merged['diff2'] = merged.psi2_x - merged.psi2_y
merged['diff3'] = merged.psi3_x - merged.psi3_y

merged.to_csv("C:/Users/arron/Documents/Phd/scripts/humvar3/merged.csv")