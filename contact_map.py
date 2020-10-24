# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 23:20:25 2020

@author: Max
"""

from main_library import contact_map, read_pdb
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", help="PDB input", type=str)
parser.add_argument("-o", help="Name of output", type=str)
parser.add_argument("-bw", help="Contact", action="store_false")
args = parser.parse_args()

df = read_pdb(args.f)

contact_map(df, args.bw, args.o)
    