import os

from Bio import AlignIO
from crowelab_pyir import PyIR
import json
import sys
import subprocess
import pathlib
import run_TMalign
import rmsd_by_segment
import pandas as pd
import numpy as np
import shutil
from pathlib import Path
import re
import swalign
# import Bio.Align
# import Bio.AlignIO
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

# if __name__ == '__main__':

def copy_files(name):
    # path = "/Users/kseniapolonsky/Downloads/alphafold_results/AP2/" + name + "/"
    path = "/Users/kseniapolonsky/Downloads/RES/" + name + "/"
    # for name in list:
    cwd = os.getcwd()
    code = name.split("_")[0]
    chain = name.split("_")[1]
    target_path = os.path.join(cwd, code)
    # source_path = os.path.join(path, name)
    source_path = path
    fasta_name = "%s.fasta" % (name)
    pdb_name = "%s.pdb" % (code.lower())
    pdb_w_chain = "%s.pdb" % (code.lower() + "_" + chain)

    source_pdb = os.path.join(target_path, pdb_name)
    target_pdb = os.path.join(cwd, pdb_name)
    # my_file = Path(target_pdb)
    my_file = Path(target_pdb)
    if not my_file.exists():
        shutil.copyfile(source_pdb, target_pdb)

    # # source_fasta = os.path.join(source_path, fasta_name)
    # source_fasta = "/Users/kseniapolonsky/Downloads/alphafold_results/AP2/" + fasta_name
    source_fasta = os.path.join(target_path, fasta_name)
    target_fasta = os.path.join(cwd, fasta_name)
    my_file = Path(target_fasta)
    if not my_file.exists():
        shutil.copyfile(source_fasta, target_fasta)

    pred = "ranked_0.pdb"
    pred_full = name + "_" + pred

    source_pred = os.path.join(source_path, pred)
    target_pred = os.path.join(cwd, pred_full)
    my_file = Path(target_pred)
    if not my_file.exists():
        shutil.copyfile(source_pred, target_pred)

    return fasta_name, chain, pdb_name, pred_full, pdb_w_chain

def get_species(file):
    file_name = str(pathlib.Path.cwd()) + '/' + str(file)
    chain = file.split("_")[1]
    chain = chain.split(".")[0]
    f = open(file_name, 'r')
    species = "human"
    for line in f:
        if line.startswith('>'):
            # if (chain in line) and ("Mus musculus" in line):
            if (chain in line) and ("10090" in line):
                species = "mouse"
            # if (chain in line) and ("Homo sapiens" in line):
            if (chain in line) and ("9606" in line):
                species = "human"
            # if (chain in line) and ("Oryctolagus" in line):
            if (chain in line) and ("9986" in line):
                species = "rabbit"
            if (chain in line) and ("9544" in line):
            # if (chain in line) and ("Macaca mulatta" in line):
                species = "macaca"
            # if (chain in line) and ("Bos taurus" in line):
            if (chain in line) and ("9913" in line):
                species = "bovine"
            # if (chain in line) and ("Rattus norvegicus" in line):
            if (chain in line) and ("10116" in line):
                species = "rat"
    return species


