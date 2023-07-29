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
ABC = ["A","R","N","D","C", "E", "Q", "G", "H","I","L","K", "M","F","P","S","T","W","Y","V","X"]
# ABC_dict = {"ALA":"A", "ARG":"R","ASN":"N","ASP":"D","CYS":"C", "GLU":"E", "GLN":"Q", "GLY":"G", "HIS":"H","ILE":"I","LEU":"L","LYS":"K", "MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
ABC_dict = {"ALA":"A", "AALA":"A", "BALA":"A", "ARG":"R", "AARG":"R", "BARG":"R","ASN":"N", "AASN":"N", "BASN":"N","ASP":"D", "AASP":"D", "BASP":"D", "CYS":"C", "ACYS":"C", "BCYS":"C", "GLU":"E", "AGLU":"E", "BGLU":"E",  "GLN":"Q", "AGLN":"Q", "BGLN":"Q", "GLY":"G", "AGLY":"G", "BGLY":"G", "HIS":"H", "AHIS":"H", "BHIS":"H","ILE":"I", "AILE":"I", "BILE":"I","LEU":"L", "ALEU":"L", "BLEU":"L", "LYS":"K", "ALYS":"K", "BLYS":"K", "MET":"M", "AMET":"M", "BMET":"M","PHE":"F", "APHE":"F", "BPHE":"F","PRO":"P", "APRO":"P", "BPRO":"P","SER":"S", "ASER":"S", "BSER":"S", "THR":"T", "ATHR":"T", "BTHR":"T", "TRP":"W", "ATRP":"W", "BTRP":"W", "TYR":"Y", "ATYR":"Y", "BTYR":"Y", "VAL":"V", "AVAL":"V", "BVAL":"V"}

def getseq(df):
    seq1 = ""
    seq2 = ""
    residues_A_t = df.loc[df['Chain'] == 'A', 'Residue_type'].tolist()
    residues_B_t = df.loc[df['Chain'] == 'B', 'Residue_type'].tolist()

    for residue in residues_A_t:
        seq1 += ABC_dict[residue]
    for residue in residues_B_t:
        seq2 += ABC_dict[residue]
    return seq1, seq2

def sw_alignment(df):
    df_a = df.loc[df['Chain'] == 'A'].drop_duplicates(subset=['Residue_id'])
    df_b = df.loc[df['Chain'] == 'B'].drop_duplicates(subset=['Residue_id'])
    df = pd.concat([df_a, df_b], ignore_index=True, axis=0)
    residues_A = df.loc[df['Chain'] == 'A', 'Residue_id'].tolist()
    residues_B = df.loc[df['Chain'] == 'B', 'Residue_id'].tolist()
    seq1, seq2 = getseq(df)
    aligner = Align.PairwiseAligner()
    # aligner.mode = 'global'
    aligner.match_score = 3.0
    aligner.gap_score = -3.0
    aligner.extend_gap_score = -2.0
    aligner.mismatch_score = -1.0
    alignments = aligner.align(seq1, seq2)
    alignment = sorted(alignments)[0]
    # alignment = next(alignments)
    print("Score = %.1f:" % alignment.score)
    print(alignment)
    aligned = alignment.aligned
    aligned_A = alignment.aligned[0]
    aligned_B = alignment.aligned[1]
    # list1 = ['INS'] * len(seq1)
    # list2 = ['del'] * len(seq2)
    dict = {}
    dict2 = {}
    for tuple1, tuple2 in zip(aligned_A, aligned_B):
        temp = zip(residues_A[tuple1[0]:tuple1[1]], residues_B[tuple2[0]:tuple2[1]])
        dict.update(temp)
    print(dict)
    for tuple1, tuple2 in zip(aligned_B, aligned_A):
        temp = zip(residues_B[tuple1[0]:tuple1[1]], residues_A[tuple2[0]:tuple2[1]])
        dict2.update(temp)
    print(dict2)
    return dict, dict2, residues_A, residues_B



def parse(file):
    file_name = str(pathlib.Path.cwd()) + '/' + str(file)
    f = open(file_name, 'r')
    cols = ['At', 'SN', 'Atom', 'Residue_type', 'Chain', 'Residue_id', 'x', 'y', 'z']
    df = pd.DataFrame(columns = cols)
    dict = {}
    # for line in f:
    #     if line.startswith('REMARK Chain 1'):
    #         temp = line.split()
    #         size1 = int(temp[4])
    #     if line.startswith('REMARK Chain 2'):
    #         temp = line.split()
    #         size2 = int(temp[4])
    # size = min(size1, size2)
    # dict = {k : int(0) for k in range(1, size+1)}
    # f.seek(0)

    for line in f:
        if line.startswith('select') and 'select all' not in line and 'select *A' not in line and 'select *B' not in line:
            temp = line.split()
            a = int(temp[1].split(':')[0])
            b = int(temp[2].split(':')[0])
            dict[a] = b
        if line.startswith('ATOM'):
            line_set = line.split()
            # db2 = pd.DataFrame(line_set)
            df.loc[len(df)] = line_set
    # df['Residue_id'] = df['Residue_id'].astype(int)
    df['x'] = df['x'].astype(float)
    df['y'] = df['y'].astype(float)
    df['z'] = df['z'].astype(float)

    # while True:
    #     if dict[len(dict.keys())] == 0:
    #         del dict[len(dict.keys())]
    #     else:
    #         break
    return df, dict

def parse2(original, predicted):
    # cols = ['At', 'SN', 'Atom', 'Residue_type', 'Chain', 'Residue_id', 'x', 'y', 'z', 'na1', 'na2', 'na3']
    cols = ['At', 'SN', 'Atom', 'Residue_type', 'Chain', 'Residue_id', 'x', 'y', 'z']
    df1 = pd.DataFrame(columns=cols)
    df2 = pd.DataFrame(columns=cols)

    file_name1 = str(pathlib.Path.cwd()) + '/' + str(original)
    f = open(file_name1, 'r')

    for line in f:
        if line.startswith('ATOM'):
            line_set = line.split()
            df1.loc[len(df1)] = line_set[:9]
    df1['x'] = df1['x'].astype(float)
    df1['y'] = df1['y'].astype(float)
    df1['z'] = df1['z'].astype(float)
    df1['Chain'] = 'A'
    f.close()

    file_name2 = str(pathlib.Path.cwd()) + '/' + str(predicted)
    f = open(file_name2, 'r')

    for line in f:
        if line.startswith('ATOM'):
            line_set = line.split()
            df2.loc[len(df2)] = line_set[:9]
    df2['x'] = df2['x'].astype(float)
    df2['y'] = df2['y'].astype(float)
    df2['z'] = df2['z'].astype(float)
    df2['Chain'] = 'B'
    f.close()

    df = pd.concat([df1, df2], ignore_index=True, axis=0)
    df = df.loc[df['Atom'] == 'CA']

    return df



def get_alignment(out, fasta, df):
    # f = open(fasta, 'r')
    # for line in f:
    #     if not line.startswith(">"):
    #         sequence = str(line)

    residues_A = df.loc[df['Chain'] == 'A', 'Residue_id'].tolist()
    residues_B = df.loc[df['Chain'] == 'B', 'Residue_id'].tolist()

    true = 0
    out_by_line = out.split("\n")
    for line in out_by_line:
        if true == 2:
            true += 1
        elif true == 3:
            my_list2 = list(line)
            n = 0
            ins = 1
            for i in range(len(my_list2)):
                if my_list2[i] != '-':
                    my_list2[i] = residues_B[n]
                    n += 1
                else:
                    my_list2[i] = "gap" + str(ins)
                    ins += 1
            true = 0
            break
        elif true == 1:
            my_list = list(line)
            n = 0
            ins = 1
            for i in range(len(my_list)):
                if my_list[i] != '-':
                    my_list[i] = residues_A[n]
                    n += 1
                else:
                    my_list[i] = "INS" + str(ins)
                    ins += 1
            true += 1
        else:
            if "denotes residue pairs of d <  5.0 Angstrom" in line:
                true += 1
    if len(my_list) != len(my_list2):
        print("\n alignment is not accurate \n")
    zip_iterator = zip(my_list, my_list2)
    alignment = dict(zip_iterator)
    return alignment, residues_A, residues_B



def add_secondary(df, file, chain, residues_A):
    file_name = str(pathlib.Path.cwd()) + '/' + str(file)
    f = open(file_name, 'r')
    df['Secondary'] = np.nan
    df['SSbond'] = np.nan
    for line in f:
        if line.startswith('HELIX'):
            temp = line.split()
            x = temp[5]
            y = temp[8]
            if chain == temp[4] and chain == temp[7]:
                idx_x = df.index[(df['Residue_id'] == x)][0]
                idx_y = df.index[(df['Residue_id'] == y)][0]
                label = temp[0]
                df.loc[(df.index >= idx_x) & (df.index <= idx_y) & (df['Chain'] == 'A'), 'Secondary'] = label
                # df.loc[df['Residue_id'].astype(int).between(x, y), 'Secondary'] = label
        if line.startswith('SHEET'):
            temp = line.split()
            x = temp[6]
            y = temp[9]
            if chain == temp[5] and chain == temp[8]:
                idx_x = df.index[df['Residue_id'] == x][0]
                idx_y = df.index[df['Residue_id'] == y][0]
                label = temp[0]
                df.loc[(df.index >= idx_x) & (df.index <= idx_y) & (df['Chain'] == 'A'), 'Secondary'] = label
                # df.loc[df['Residue_id'].astype(int).between(x, y), 'Secondary'] = label
        if line.startswith('SSBOND'):
            temp = line.split()
            x = temp[4]
            y = temp[7]
            if chain == temp[3] and chain == temp[6]:
                idx_x = df.index[df['Residue_id'] == x][0]
                idx_y = df.index[df['Residue_id'] == y][0]
                label = temp[0]
                # df.loc[df['Residue_id'].astype(int).between(x, y), 'SSbond'] = label
                df.loc[(df.index >= idx_x) & (df.index <= idx_y) & (df['Chain'] == 'A'), 'SSbond'] = label
    return df

def color_code_results_domain(dict, new_dict):
    shutil.copyfile('TM_sup_atm.pml', 'TM_sup_atm_colorcoded.pml')
    full_name = str(pathlib.Path.cwd()) + '/' + str('TM_sup_atm_colorcoded.pml')
    f = open(full_name, 'a')
    f.write("color grey, chain A\n")
    f.write("color grey, chain B\n")
    low = float(1)
    med = float(1.6)
    high = float(2)
    n = 0
    for domain in dict:
        if domain in new_dict and float(new_dict[domain]) < low:
            text = "\nselect toBeColored" + str(n) + ", ((i. " + str(dict[domain][0]) + "-" + str(
                dict[domain][1]) + ") and c. A )\ncolor green, toBeColored" + str(n)
            f.write(text)
            n += 1
        if domain in new_dict and float(new_dict[domain]) > low and float(new_dict[domain]) < med:
            text = "\nselect toBeColored" + str(n) + ", ((i. " + str(dict[domain][0]) + "-" + str(
                dict[domain][1]) + ") and c. A )\ncolor blue, toBeColored" + str(n)
            f.write(text)
            n += 1
        if domain in new_dict and float(new_dict[domain]) > med and float(new_dict[domain]) < high:
            text = "\nselect toBeColored" + str(n) + ", ((i. " + str(dict[domain][0]) + "-" + str(
                dict[domain][1]) + ") and c. A )\ncolor orange, toBeColored" + str(n)
            f.write(text)
            n += 1
        if domain in new_dict and float(new_dict[domain]) > high:
            text = "\nselect toBeColored" + str(n) + ", ((i. " + str(dict[domain][0]) + "-" + str(
                dict[domain][1]) + ") and c. A )\ncolor red, toBeColored" + str(n)
            f.write(text)
            n += 1
    f.close()

def color_code_results(alignment, df):
    # shutil.copyfile('TM_sup_atm.pml', 'TM_sup_atm_colorcoded_ribbon.pml')
    shutil.copyfile('TM_sup_all_atm.pml', 'TM_sup_atm_colorcoded_ribbon.pml')
    full_name = str(pathlib.Path.cwd()) + '/' + str('TM_sup_atm_colorcoded_ribbon.pml')
    f = open(full_name, 'a')
    # f.write("hide all\n")
    # f.write("ribbon_width, 12\n")
    # f.write("show ribbon\n")
    f.write("color grey, chain A\n")
    f.write("color grey, chain B\n")

    red, blue, green, orange = [], [], [], []

    for residue in alignment:

        if 'INS' not in str(residue):

            if 'gap' not in str(alignment[residue]):

                x1 = df.loc[(df['Residue_id'] == residue) & (df['Chain'] == 'A'), 'x'].values[0]
                y1 = df.loc[
                    (df['Residue_id'] == residue) & (df['Chain'] == 'A'), 'y'].values[0]
                z1 = df.loc[
                    (df['Residue_id'] == residue) & (df['Chain'] == 'A'), 'z'].values[0]
                x2 = df.loc[
                    (df['Residue_id'] == alignment[residue]) & (df['Chain'] == 'B'), 'x'].values[0]
                y2 = df.loc[
                    (df['Residue_id'] == alignment[residue]) & (df['Chain'] == 'B'), 'y'].values[0]
                z2 = df.loc[
                    (df['Residue_id'] == alignment[residue]) & (df['Chain'] == 'B'), 'z'].values[0]
                rmsd = ((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)**(1/2)
                if rmsd <= 1.9:
                    green.append(residue)
                elif rmsd > 1.9 and rmsd <= 2.5:
                    blue.append(residue)
                elif rmsd > 2.5 and rmsd <= 3.0:
                    orange.append(residue)
                elif rmsd > 3.0:
                    red.append(residue)
            else:
        # elif alignment[residue] == 0:
                red.append(residue)
        else:
            red.append(residue)

    n = 0
    for color in [[green,"green"], [blue,"blue"], [orange,"orange"], [red,"red"]]:
        text = ""
        for residue in color[0]:
            if residue != color[0][-1]:
                text += str(residue) + "+"
            else:
                text += str(residue)
        text_f = "\nselect toBeColored" + str(n) + ", ((i. " + text + ") and c. A )\ncolor " + str(color[1]) +", toBeColored" + str(n)
        f.write(text_f)
        n += 1
    f.close()

def color_code_results_cartoon_GDT(alignment, df):
    # shutil.copyfile('TM_sup_atm.pml', 'TM_sup_atm_colorcoded_cartoon.pml')
    shutil.copyfile('TM_sup_all_atm.pml', 'TM_sup_atm_colorcoded_cartoon.pml')
    full_name = str(pathlib.Path.cwd()) + '/' + str('TM_sup_atm_colorcoded_cartoon.pml')
    f = open(full_name, 'a')
    f.write("hide all\n")
    f.write("show_as cartoon, chain A\n")
    f.write("color grey, chain A\n")

    red, blue, green, orange = [], [], [], []
    p05, p1, p2, p4, p8 = 0, 0, 0, 0, 0

    n = 0

    for residue in alignment:

        if 'INS' not in str(residue):

            if 'gap' not in str(alignment[residue]):

                x1 = df.loc[(df['Residue_id'] == residue) & (df['Chain'] == 'A'), 'x'].values[0]
                y1 = df.loc[
                    (df['Residue_id'] == residue) & (df['Chain'] == 'A'), 'y'].values[0]
                z1 = df.loc[
                    (df['Residue_id'] == residue) & (df['Chain'] == 'A'), 'z'].values[0]
                x2 = df.loc[
                    (df['Residue_id'] == alignment[residue]) & (df['Chain'] == 'B'), 'x'].values[0]
                y2 = df.loc[
                    (df['Residue_id'] == alignment[residue]) & (df['Chain'] == 'B'), 'y'].values[0]
                z2 = df.loc[
                    (df['Residue_id'] == alignment[residue]) & (df['Chain'] == 'B'), 'z'].values[0]
                rmsd = ((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) ** (1 / 2)
                n += 1

                if rmsd <= 0.5:
                    p05 += 1
                if rmsd <= 1:
                    p1 += 1
                if rmsd <= 2:
                    p2 += 1
                if rmsd <= 4:
                    p4 += 1
                if rmsd <= 8:
                    p8 += 1

                if rmsd <= 1.9:
                    green.append(residue)
                elif rmsd > 1.9 and rmsd <= 2.5:
                    blue.append(residue)
                elif rmsd > 2.5 and rmsd <= 3.0:
                    orange.append(residue)
                elif rmsd > 3.0:
                    red.append(residue)
            else:
                red.append(residue)
        else:
            red.append(residue)

    GDT_TS = (p1/n + p2/n + p4/n + p8/n)/4
    GDT_HA = (p05/n + p1/n + p2/n + p4/n) / 4

    n = 0
    for color in [[green, "green"], [blue, "blue"], [orange, "orange"], [red, "red"]]:
        text = ""
        for residue in color[0]:
            if residue != color[0][-1]:
                text += str(residue) + "+"
            else:
                text += str(residue)
        text_f = "\nselect toBeColored" + str(n) + ", ((i. " + text + ") and c. A )\ncolor " + str(
            color[1]) + ", toBeColored" + str(n)
        f.write(text_f)
        n += 1
    f.close()
    return GDT_TS, GDT_HA

# def convert_aligned_indices(pdb_2, alignment):
#     cols = ['At', 'SN', 'Atom', 'Residue_type', 'Chain', 'Residue_id', 'x', 'y', 'z', 'na1', 'na2', 'na3']
#     df = pd.read_csv(pdb_2, delim_whitespace=True, names=cols, index_col = False)
#     # df = pd.read_csv(pdb_2, delim_whitespace=True)
#     # df = pd.read_csv(pdb_2, sep=' ', delimiter=None, index_col = False, header = None)
#     rowEND = df.loc[(df['At'] == 'END')]
#     temp_at = df.loc[(df['At'] == 'TER'), 'At'].values[0]
#     temp_sn = df.loc[(df['At'] == 'TER'), 'SN'].values[0]
#     temp_rt = df.loc[(df['At'] == 'TER'), 'Atom'].values[0]
#     temp_ch = df.loc[(df['At'] == 'TER'), 'Residue_type'].values[0]
#     temp_rn = df.loc[(df['At'] == 'TER'), 'Chain'].values[0]
#     df = df.drop(df.index[df['At'] == 'TER'])
#     df = df.drop(df.index[df['At'] == 'END'])
#     df.loc[-1] = [temp_at, temp_sn, '', temp_rt, temp_ch, temp_rn, '', '', '', '', '', '']
#
#     df['SN'] = df['SN'].astype(int)
#     df['SN'] = df['SN'].astype(str)
#     df['Residue_id'] = df['Residue_id'].astype(int)
#     df['Residue_id'] = df['Residue_id'].astype(str)
#     for residue in alignment:
#         df.loc[df['Residue_id'] == alignment[residue], 'Residue_id'] = residue
#     df = pd.concat([df, rowEND], ignore_index=True, axis=0)
#     df.to_csv("test.pdb", header=None, index=None, sep=' ', mode='w')
def check_index(line_set):
    index1 = -1
    index2 = -1
    for i in range(len(line_set)):
        if line_set[i] != '' and index1 != -1:
            index2 = i
            temp = line_set[i]
            return temp, index1, index2
        if line_set[i] == 'A':
            index1 = i


def convert_aligned_indices(pdb_2, alignment):
    file_name = str(pathlib.Path.cwd()) + '/' + str(pdb_2)
    with open(file_name, 'r') as file:
        data = file.readlines()

    for i in range(len(data)):
        line_set = data[i].split()
        line_set2 = data[i].split(' ')

        if (data[i].startswith('TER')):
            temp, index1, index2 = check_index(line_set2)
            if temp in alignment:
                temp = alignment[temp]
            else:
                temp = 'N' * len(str(temp))

        elif data[i].startswith('END'):
            break
        else:
            temp, index1, index2 = check_index(line_set2)
            if temp in alignment:
                temp = alignment[temp]
            else:
                temp = 'N' * len(str(temp))
                # temp = str(line_set[5]) + 'X'

        n = len(temp) - len(line_set2[index2])
        line_set2[index2] = temp

        joined_string_1 = " ".join(line_set2[:index1])
        joined_string_2 = " ".join(line_set2[index1:])

        string_base = " " + str(temp) + " "

        # e.g. 99 to 99A : 1 space less on the right
        if (n == 1) & (not str(temp).isnumeric()):
            string = string_base + " " * n
            joined_string_2 = re.sub(string, string_base, joined_string_2)
        # e.g. 99 to 100A : 1 space less on the right and 1 less on the left
        elif (n == 2) & (not str(temp).isnumeric()):
            string = " " + string_base + " "
            joined_string_2 = re.sub(string, string_base, joined_string_2)
        # e.g. 100 to 99A : 1 more on left, 1 less on right
        elif (n == 0) & (not str(temp).isnumeric()):
            string = string_base + " "
            string_base = " " + string_base
            joined_string_2 = re.sub(string, string_base, joined_string_2)
        # e.g. 99 to 100 : n spaces less on the left
        elif (n >= 1) & (str(temp).isnumeric()):
            string = n * " " + string_base
            joined_string_2 = re.sub(string, string_base, joined_string_2)
        # e.g. 100 to 99 : n spaces more on the left
        elif (n <= -1) & (str(temp).isnumeric()):
            string = " " * (n * -1) + string_base
            joined_string_2 = re.sub(string_base, string, joined_string_2)

        joined_string = joined_string_1 + " " + joined_string_2
        data[i] = joined_string

    file.close()

    file_name = str(pathlib.Path.cwd()) + '/' + "x_" + str(pdb_2)
    # and write everything back
    with open(file_name, 'w') as file2:
        file2.writelines(data)
    file2.close()

def chain_pdb(pdb1, pdb2, CHAIN):
    file_name = str(pathlib.Path.cwd()) + '/' + str(pdb1)
    with open(file_name, 'r') as file:
        data = file.readlines()

    for i in range(len(data)):
        line_set = data[i].split()
        if (not data[i].startswith('ATOM')) or (not line_set[4] == CHAIN):
            if (not data[i].startswith('TER')) and (not data[i].startswith('END')):
                data[i] = ''
            elif (data[i].startswith('TER')) and (line_set[3] != CHAIN):
                data[i] = ''

    data = [x for x in data if x != '']

    file.close()

    file_name = str(pathlib.Path.cwd()) + '/' + str(pdb2)
    # and write everything back
    with open(file_name, 'w') as file2:
        file2.writelines(data)
    file2.close()

def check_chain_type(temp):
    lc_type = ''
    if "IGLV" in temp[1]:
        lc_type = "lambda"
    elif "IGKV" in temp[1]:
        lc_type = "kappa"
    elif "IGHV" in temp[1]:
        lc_type = "heavy"
    return lc_type


def cd3_location2(dict, chain_type, df, alignment):
    index = 0
    if 'CDR3' not in dict:
        dict['CDR3'] = [0, 0]
        cd3_start = 0
    elif 'CDR3' in dict:
        cd3_start = int(dict['CDR3'][0])
    seq1, seq2 = getseq(df)
    if 'FR3' in dict:
        index = int(dict['FR3'][1])

    seq = seq2[index - 3:]

    if chain_type == '':
        print("\n the chain is not defined \n")
    if chain_type == 'lambda':
        # start = ['TACTGT']
        start = ['YC']
        # end = ['GTATTC']
        end = ['VF','FGG']
    elif chain_type == 'kappa':
        # start = ['CAACAG']
        start = ['YC', 'QQ']
        # end = ['TTTCGG']
        end = ['FR', 'FGG']
    elif chain_type == 'heavy':
        # start = ['TGTGC','TGAG']
        start = ['YC', 'FC']
        # start = ['CA', 'CVK', 'CCR', 'CVR']
        # end = ['TGGGG']
        end = ['WG']

    # align sequence to the relevant start and and and see the location with the best score
    # start from the end of the known location of the FR3
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    # aligner.mode = 'global'
    aligner.match_score = 3.0
    aligner.gap_score = -3.0
    aligner.extend_gap_score = -2.0
    aligner.mismatch_score = -1.0

    cdr3_start_index = 0
    cdr3_end_index = 0

    if cd3_start == 0:
        score = 0
        for code in start:
            alignments = aligner.align(seq, code)
            alignment = sorted(alignments)[0]
            # alignment = next(alignments)
            # aligned = alignment.aligned
            for aligned_A in alignment.aligned:
                if aligned_A[0][1] < 7:
                    if alignment.score > score:
                        cdr3_start_index = aligned_A[0][1] + (index-3)
                        score = alignment.score
                    # print("Score = %.1f:" % alignment.score)
                    # print(alignment)
                    if score == aligner.match_score * len(code):
                        break
    else:
        cdr3_start_index = cd3_start

    score = 0
    for code in end:
        alignments = aligner.align(seq, code)
        alignment = sorted(alignments)[0]
    # alignment = next(alignments)
        for aligned_A in alignment.aligned:
            if aligned_A[0][0] < 30 and aligned_A[0][0] > 2:
                if alignment.score > score:
                    cdr3_end_index = aligned_A[0][0] + (index - 3)
                    score = alignment.score
                # print("Score = %.1f:" % alignment.score)
                # print(alignment)
                if score == aligner.match_score * len(code):
                    break

    dict['CDR3'][0] = cdr3_start_index
    dict['CDR3'][1] = cdr3_end_index



    # if 'FR3' in dict:
    #     print('\n' + str(dict['FR3']) + '\n')
    if 'CDR3' in dict:
        print('\n' + str(dict['CDR3']) + '\n')
        cd3seq = seq2[cdr3_start_index - 1: cdr3_end_index]
        print('CDR3 is: ' + cd3seq)

    if dict['CDR3'][0] == 0 or dict['CDR3'][1] == 0:
        del dict['CDR3']

    print(dict)
    return dict, cd3seq
