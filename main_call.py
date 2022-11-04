from crowelab_pyir import PyIR
import json
import sys
import subprocess
import pathlib
import run_TMalign
import rmsd_by_segment
import parse_superimposed
import numpy as np
import abangle
import shutil

from rmsd.export_to_spreadsheet import export
from predict_from_list import copy_files
from predict_from_list import get_species

# arguments = sys.argv
#
# # arg_dict = {}
# # for i in range(len(arguments)):
# #     if i%2 != 0:
# #         key = arguments[i]
# #     elif (i%2 == 0) & (i != 0):
# #         arg_dict[key] = arguments[i]
#
#
# # FILE = arg_dict['-fa']
# # CHAIN = arg_dict['-ch']
# # PDB_F = arg_dict['-fll']
# # PDB_PRED = arg_dict['-pr']
# # PDB_ORIG = arg_dict['-or']
#
# # name = "1NCW_L"
# name = arguments[1]
#
# FILE, CHAIN, PDB_F, PDB_PRED, PDB_ORIG = copy_files(name)
#
# # create a pdb for this specific chain *\n"
# parse_superimposed.chain_pdb(PDB_F, PDB_ORIG, CHAIN)
#
# # if '-sp' in arg_dict:
# #     SPECIES = arg_dict['-sp']
# # else:
# #     SPECIES = 'human'



# if __name__ == '__main__':
def main_call_function(name, FILE, CHAIN, PDB_F, PDB_PRED, PDB_ORIG):

    parse_superimposed.chain_pdb(PDB_F, PDB_ORIG, CHAIN)

    " *********************************************************************\n"
    " * get a pdb of a relevant chain only *\n"
    " * get the name and code of the species *\n"
    " * write into the log *\n"
    " *********************************************************************\n"

    SPECIES = get_species(FILE)
    print(FILE + " " + CHAIN + " " + PDB_F + " " + PDB_PRED + " " + PDB_ORIG + " " + SPECIES)

    fp = open('log.txt', 'a')
    fp.write(FILE + " " + CHAIN + " " + PDB_F + " " + PDB_PRED + " " + PDB_ORIG + " " + SPECIES + '\n')
    fp.close()

    dict = {}

    " *********************************************************************\n"
    " * create a dataframe with original chain as A and predicted chain as B *\n"
    " * pairwise align A and B, change residues in B to the aligned names from A *\n"
    " *********************************************************************\n"

    df1 = parse_superimposed.parse2(PDB_ORIG, PDB_PRED)
    alignment, alignment_reverse, residues_A, residues_B = parse_superimposed.sw_alignment(df1)
    parse_superimposed.convert_aligned_indices(PDB_PRED, alignment_reverse)


    " *********************************************************************\n"
    " * send the prediction from AlphaFord to TMalign *\n"
    " * parameter 0: sequence independent TMalign, 1: sequence dependent *\n"
    " *********************************************************************\n"
    # mode = 0
    mode = 1
    out = run_TMalign.run_TMalign(mode, name)
    # out = run_TMalign.run_TMalign(1)

    " *********************************************************************\n"
    " * send the original fasta sequence of the chain to IgBlast and *\n"
    " * get the dictionary with indices of each Ab segment *\n"
    " *********************************************************************\n"
    # pyirexample = PyIR(query=FILE, args=['--sequence_type', 'prot', '--legacy', '--gzip', 'False', '--germlineV',
    #                                      './db/human_igh_v'])

    if SPECIES == 'human':
    # homo sapiens
        pyirexample = PyIR(query=FILE, args=['--sequence_type', 'prot', '--legacy', '--gzip', 'False', '--germlineV',
                                         './db/human_igh_v', '--germlineD','./db/human_igh_d','--germlineJ','./db/human_igh_j'])
    if SPECIES == 'rabbit':
    # Oryctolagus cuniculus
        pyirexample = PyIR(query=FILE, args=['--sequence_type', 'prot', '--legacy', '--gzip', 'False', '--germlineV',
                                         './db/rabbit_igh_v', '--germlineD','./db/rabbit_igh_d','--germlineJ','./db/rabbit_igh_j'])
    if SPECIES == 'mouse':
    # mus muscus
        pyirexample = PyIR(query=FILE, args=['--sequence_type', 'prot', '--legacy', '--gzip', 'False', '--germlineV',
                                         './db/mouse_igh_v', '--germlineD','./db/mouse_igh_d','--germlineJ','./db/mouse_igh_j'])

    if SPECIES == 'bovine':
    # bovine
        pyirexample = PyIR(query=FILE, args=['--sequence_type', 'prot', '--legacy', '--gzip', 'False', '--germlineV',
                                         './db/bovine_igh_v', '--germlineD','./db/bovine_igh_d','--germlineJ','./db/bovine_igh_j'])

    if SPECIES == 'macaca':
    # Macaca mulatta
        pyirexample = PyIR(query=FILE, args=['--sequence_type', 'prot', '--legacy', '--gzip', 'False', '--germlineV',
                                         './db/macaca_igh_v', '--germlineD','./db/macaca_igh_d','--germlineJ','./db/macaca_igh_j'])

    if SPECIES == 'rat':
    # Rattus norvegicus
        pyirexample = PyIR(query=FILE, args=['--sequence_type', 'prot', '--legacy', '--gzip', 'False', '--germlineV',
                                         './db/rat_igh_v', '--germlineD','./db/rat_igh_d','--germlineJ','./db/rat_igh_j'])


    result = pyirexample.run()


    " *********************************************************************\n"
    " * define residues of each domain (CDR1/2/3, FR1/2/3/4) *\n"
    " *********************************************************************\n"

    chain_type = ''

    fp = open('log.txt', 'a')

    # define residues range for each of the domains
    f = open(result, 'r')
    for line in f:
        temp = line.split('Hits')
        print(temp[1] + "\n")

        fp.write(temp[1] + "\n")
        fp.close()

        chain_type = parse_superimposed.check_chain_type(temp)

        entry = json.loads(line)
        if 'CDR1' in entry:
            dict['CDR1'] = [entry['CDR1']['from'], entry['CDR1']['to']]
        if 'CDR2' in entry:
            dict['CDR2'] = [entry['CDR2']['from'], entry['CDR2']['to']]
        if 'CDR3' in entry:
            dict['CDR3'] = [entry['CDR3']['from'], entry['CDR3']['to']]
        if 'FR1' in entry:
            dict['FR1'] = [entry['FR1']['from'], entry['FR1']['to']]
        if 'FR2' in entry:
            dict['FR2'] = [entry['FR2']['from'], entry['FR2']['to']]
        if 'FR3' in entry:
            dict['FR3'] = [entry['FR3']['from'], entry['FR3']['to']]
        if 'FR4' in entry:
            dict['FR4'] = [entry['FR4']['from'], entry['FR4']['to']]
    f.close()

    " *********************************************************************\n"
    " * take the superimposed molecules from TM_sup_all  file produced *\n"
    " * by TMalign and calculate RMSD by segment *\n"
    " *********************************************************************\n"

    df, alignment = parse_superimposed.parse('TM_sup_all')

    # these are two options of the alignment: 1 based on the sequence, second one based on TM-align distance alignment
    # alignment, residues_A, residues_B = parse_superimposed.sw_alignment(df)
    alignment, residues_A, residues_B = parse_superimposed.get_alignment(out, FILE, df)

    # parse_superimposed.convert_aligned_indices(PDB_ORIG, alignment)

    # column = df["Residue_id"]
    # max_index = column.max()
    max_index = len(df.index)
    # df["Domain"] = np.nan
    df["Domain"] = 'therest'
    my_list = [None] * int(max_index)
    new_dict = {}
    for key in dict:
        x = int(dict[key][0])
        y = int(dict[key][1])
        # df.loc[(df['Residue_id'].astype(int).between(residues_A[x-1], residues_A[y-1])) & (df['Chain'] == 'A'), 'Domain'] = key
        idx_x = df.index[(df['Residue_id'] == residues_A[x-1]) & (df['Chain']=='A') ][0]
        idx_y = df.index[(df['Residue_id'] == residues_A[y-1]) & (df['Chain']=='A') ][0]
        df.loc[(df.index >= idx_x) & (df.index <= idx_y) & (df['Chain'] == 'A'), 'Domain'] = key


    # add CDR3 to the dict
    # parse_superimposed.cd3_location(dict, chain_type, df, alignment)

    classificator = 'Domain'
    rmsd = rmsd_by_segment.rmsd_for_domain(df, alignment, classificator, 'all')
    new_dict['total'] = rmsd


    dict = parse_superimposed.cd3_location2(dict, chain_type, df, alignment)



    for domain in dict:
        rmsd = rmsd_by_segment.rmsd_for_domain(df, alignment, classificator, domain)
        new_dict[domain] = rmsd
    rmsd = rmsd_by_segment.rmsd_for_domain(df, alignment, classificator,  'therest')
    new_dict['therest'] = rmsd
    # rmsd = rmsd_by_segment.rmsd_for_domain(df, alignment, classificator, 'all')
    # new_dict['total'] = rmsd


    # if 'CDR3' in new_dict:
    #     del new_dict['CDR3']



    # print(new_dict)
    # file = arguments[4]
    file = PDB_ORIG
    df = parse_superimposed.add_secondary(df, PDB_F, CHAIN, residues_A)
    classificator = 'Secondary'
    # for structure in ['SHEET','HELIX','SSBOND']:
    for structure in ['SHEET', 'HELIX']:
        rmsd = rmsd_by_segment.rmsd_for_domain(df, alignment, classificator, structure)
        new_dict[structure] = rmsd
    # print(new_dict)
    classificator = 'SSbond'
    for structure in ['SSBOND']:
        rmsd = rmsd_by_segment.rmsd_for_domain(df, alignment, classificator, structure)
        new_dict[structure] = rmsd
    if 'CDR3' not in new_dict:
        new_dict['CDR3'] = 0
    print("\n")
    print(new_dict)

    fp = open('log.txt', 'a')
    fp.write(str(new_dict))
    fp.write("\n")
    fp.close()

    parse_superimposed.color_code_results(alignment, df)
    GDT_TS, GDT_HA = parse_superimposed.color_code_results_cartoon_GDT(alignment, df)
    print("\nGDT Total Score is: " + str(GDT_TS))
    print("GDT High Accuracy Score is: " + str(GDT_HA) + "\n")
    export(mode, df, PDB_PRED, SPECIES, entry, chain_type, new_dict, out, GDT_TS, GDT_HA, dict)
    print("\n DONE \n")

    fp = open('log.txt', 'a')
    fp.write("\nGDT Total Score is: " + str(GDT_TS))
    fp.write("GDT High Accuracy Score is: " + str(GDT_HA) + "\n")
    fp.write("\n DONE \n")
    fp.close()