# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 22:05:51 2020

@author: Max
"""

from main_library import find_salt_bridges, read_pdb
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", help="PDB input", type=str)
parser.add_argument("-o", help="Name of output", type=str)
parser.add_argument("-c", help="cutoff", type=int)
args = parser.parse_args()

out = args.o + ".txt"
df = read_pdb(args.f)

with open(out, "w+") as new:
    for i in find_salt_bridges(df,args.c):
        new.write(str(i)+"\n")