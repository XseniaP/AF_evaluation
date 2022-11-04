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
import re
import swalign
# import Bio.Align
# import Bio.AlignIO
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

def export(mode, df, PDB_PRED, SPECIES, igblast_entry, chain_type, rmsd_dict , tmalign_out, GDT_TS, GDT_HA, domains_dict ):
    cols = ['mode', 'pdb', 'chain', 'chain_type', 'light_type','species', 'rmsd_total', 'total_my_code', 'TM_score', 'coverage', 'cdr1','cdr2','cdr3','fr1','fr2','fr3','the_rest','sheet','helix','ssbond','evalue', 'bit-score','best_hit','gdt_ts','dgt_ha', 'variable_len', 'mismatches', 'cdr3_loc']
    df = pd.DataFrame(columns=cols)
    line = [0] * len(cols)
    temp = PDB_PRED.split('_')
    line[0] = mode
    line[1] = temp[0]
    line[2] = temp[1]
    line[5] = SPECIES
    if chain_type == 'lambda' or chain_type == 'kappa':
        line[3] = 'light'
        line[4] = chain_type
    else:
        line[3] = 'heavy'

    start = tmalign_out.find('RMSD=') + 5
    end = tmalign_out.find(', Seq_ID=', start)
    rmsd_total = tmalign_out[start:end].replace(" ","")
    line[6] = rmsd_total

    start = tmalign_out.find('TM-score=') + 9
    end = tmalign_out.find(' (if normalized', start)
    tm_score = tmalign_out[start:end].replace(" ","")
    line[8] = tm_score

    if 'total' in rmsd_dict:
        line[7] = rmsd_dict['total']
    if 'CDR1' in rmsd_dict:
        line[10] = rmsd_dict['CDR1']
    if 'CDR2' in rmsd_dict:
        line[11] = rmsd_dict['CDR2']
    if 'CDR3' in rmsd_dict:
        line[12] = rmsd_dict['CDR3']
    if 'FR1' in rmsd_dict:
        line[13] = rmsd_dict['FR1']
    if 'FR2' in rmsd_dict:
        line[14] = rmsd_dict['FR2']
    if 'FR3' in rmsd_dict:
        line[15] = rmsd_dict['FR3']
    if 'therest' in rmsd_dict:
        line[16] = rmsd_dict['therest']
    if 'SHEET' in rmsd_dict:
        line[17] = rmsd_dict['SHEET']
    if 'HELIX' in rmsd_dict:
        line[18] = rmsd_dict['HELIX']
    if 'SSBOND' in rmsd_dict:
        line[19] = rmsd_dict['SSBOND']

    line[20] = "{:e}".format(igblast_entry['Hits'][0]['e_value'])
    line[21] = igblast_entry['Hits'][0]['bit_score']
    line[22] = igblast_entry['Hits'][0]['gene']
    line[23] = float(GDT_TS)
    line[24] = float(GDT_HA)
    line[25] = igblast_entry['Total']['length']
    line[26] = igblast_entry['Total']['mismatches']
    if 'CDR3' in domains_dict:
        line[27] = domains_dict['CDR3']


    a_series = pd.Series(line, index=df.columns)
    df = df.append(a_series, ignore_index=True)

    df.to_csv('output.csv', mode='a', index=False, header=False)

    return line










