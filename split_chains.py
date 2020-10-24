# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 21:30:58 2020

@author: Max
"""

from main_library import split_chains
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", help="PDB input", type=str)
args = parser.parse_args()

split_chains(args.f)
print("Done")