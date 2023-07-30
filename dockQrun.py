import pathlib
import requests
import os
import pandas as pd
from subprocess import call

import subprocess

# PATH = "/DockQ/"
# '/Users/kpolonsky/PycharmProjects/rmsd/DockQ/7jmp_AGpred_ABpred.pdb' '/Users/kpolonsky/PycharmProjects/rmsd/DockQ/7jmp_native.pdb' -model_chain1 A B -model_chain2 C -native_chain1 H L -native_chain2 A -no_needle -perm1 -perm2

def rundockq(code, chainsAb, chainsAg, chainsAbNat, chainsAgNat):
    pred = code + "_AGpred_ABpred.pdb"
    # pred = code + "_AGpred_ABpred.pdb"
    native = code + "_native.pdb"
    # chainsAb = "A B "
    # chainsAg = "C "
    # chainsAbNat = "H L "
    # chainsAgNat = "P "
    args = "/DockQ.py" + " " + "'" + pred + "'" + " " + "'" + native + "'" + " " + "-model_chain1 " + chainsAb + "-model_chain2 " + chainsAg + "-native_chain1 " + chainsAbNat + "-native_chain2 " + chainsAgNat + " -no_needle -perm1 -perm2"

    cmd = "/usr/bin/python3 ." + args
    pr = open("/Users/kpolonsky/PycharmProjects/rmsd/DockQ/"+ code + "_dockQ.txt", "w")
    p = subprocess.Popen(cmd, stdout=pr, shell=True)
    # out, err = p.communicate()

rundockq("7sl5", "A B ", "C ", "A B ", "C ")