from crowelab_pyir import PyIR
import json
import sys
import subprocess
import pathlib
import run_TMalign
import rmsd_by_segment
import parse_superimposed
import numpy as np


" *********************************************************************\n"
" * Given two sets of n points v and w the RMSD is defined as follows *\n"
" *  \sqrt {{\frac {1}{n}}\sum _{i=1}^{n}((v_{ix}-w_{ix})^{2}+(v_{iy}-w_{iy})^{2}+(v_{iz}-w_{iz})^{2}}})\end{aligned}}} *\n"
" *********************************************************************"

def rmsd_for_domain(df, alignment, classificator, domain):
    sum = 0
    n = 0
    for residue in alignment:
        if 'INS' not in str(residue):
            x1 = df.loc[(df['Residue_id'] == residue) & (df['Chain'] == 'A'), 'x'].values[0]
            y1 = df.loc[(df['Residue_id'] == residue) & (df['Chain'] == 'A'), 'y'].values[0]
            z1 = df.loc[(df['Residue_id'] == residue) & (df['Chain'] == 'A'), 'z'].values[0]
        else:
            x1, y1, z1 = 0, 0, 0
        if ('gap' not in str(alignment[residue])) & ('N' not in str(alignment[residue])):
            x2 = df.loc[
                (df['Residue_id'] == alignment[residue]) & (df['Chain'] == 'B'), 'x'].values[0]
            y2 = df.loc[
                (df['Residue_id'] == alignment[residue]) & (df['Chain'] == 'B'), 'y'].values[0]
            z2 = df.loc[
                (df['Residue_id'] == alignment[residue]) & (df['Chain'] == 'B'), 'z'].values[0]
        elif ('gap' in alignment[residue]) or ('N' in str(alignment[residue])):
            x2, y2, z2 = 0, 0, 0
        if classificator == 'Domain':
            if ('INS' not in str(residue)) & ('gap' not in str(alignment[residue])) & ('N' not in str(alignment[residue])):
                if (df.loc[(df['Residue_id'] == residue) & (df['Chain'] == 'A'), classificator].values[
                        0] == domain) or (domain == 'all'):
                    sum += (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
                    n += 1
        if classificator == 'Secondary':
            if ('INS' not in str(residue)) & ('gap' not in str(alignment[residue])) & ('N' not in str(alignment[residue])):
                if (df.loc[(df['Residue_id'] == residue) & (df['Chain'] == 'A'), classificator].values[0] == domain):
                    sum += (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
                    n += 1
        if classificator == 'SSbond':
            if ('INS' not in str(residue)) & ('gap' not in str(alignment[residue])) & ('N' not in str(alignment[residue])):
                if (df.loc[(df['Residue_id'] == residue) & (df['Chain'] == 'A'), classificator].values[0] == domain):
                    sum += (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
                    n += 1
    if n != 0:
        rmsd = (sum/n)**(1/2)
    else:
        rmsd = 0
    return rmsd




# def rmsd_for_domain(df, alignment, classificator, domain):
#     sum = 0
#     # n = len(alignment)
#     n = 0
#     if domain == 'all':
#         domain = '*'
#     for residue in alignment:
#         x1 = df.loc[(df['Residue_id'] == int(residue)) & (df['Chain'] == 'A'), 'x'].values[0]
#         y1 = df.loc[
#             (df['Residue_id'].astype(int) == int(residue)) & (df['Chain'] == 'A'), 'y'].values[0]
#         z1 = df.loc[
#             (df['Residue_id'].astype(int) == int(residue)) & (df['Chain'] == 'A'), 'z'].values[0]
#         if alignment[residue] != 0:
#             x2 = df.loc[
#                 (df['Residue_id'] == int(alignment[residue])) & (df['Chain'] == 'B'), 'x'].values[0]
#             y2 = df.loc[
#                 (df['Residue_id'].astype(int) == int(alignment[residue])) & (df['Chain'] == 'B'), 'y'].values[0]
#             z2 = df.loc[
#                 (df['Residue_id'].astype(int) == int(alignment[residue])) & (df['Chain'] == 'B'), 'z'].values[0]
#         elif alignment[residue] == 0:
#             x1, y1, z2 = 0, 0, 0
#         if classificator == 'Domain':
#             if (df.loc[(df['Residue_id'] == int(residue)) & (df['Chain'] == 'A'), classificator].values[
#                     0] == domain) or (domain == '*'):
#                 sum += (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
#                 n += 1
#         if classificator == 'Secondary':
#             if (df.loc[(df['Residue_id'] == int(residue)) & (df['Chain'] == 'A'), classificator].values[0] == domain):
#                 sum += (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
#                 n += 1
#     if n != 0:
#         rmsd = (sum/n)**(1/2)
#     else:
#         rmsd = 0
#     return rmsd





# def rmsd_for_domain(df, alignment, classificator, domain):
#     length = len(alignment)
#     n = 0
#     sum = 0
#     for i in range(1, length+1):
#         x1 = df.loc[(df['Residue_id'] == int(i)) & (df['Chain'] == 'A'), 'x'].values[0]
#         y1 = df.loc[(df['Residue_id'].astype(int) == int(i)) & (df['Chain'] == 'A'), 'y'].values[0]
#         z1 = df.loc[(df['Residue_id'].astype(int) == int(i)) & (df['Chain'] == 'A'), 'z'].values[0]
#         x2 = df.loc[(df['Residue_id'] == int(i)) & (df['Chain'] == 'B'), 'x'].values[0]
#         y2 = df.loc[(df['Residue_id'].astype(int) == int(i)) & (df['Chain'] == 'B'), 'y'].values[0]
#         z2 = df.loc[(df['Residue_id'].astype(int) == int(i)) & (df['Chain'] == 'B'), 'z'].values[0]
#
#         if classificator == 'Domain':
#             if (df.loc[(df['Residue_id'] == int(i)) & (df['Chain'] == 'A'), classificator].values[0] == domain) or (domain == '*'):
#                     sum += (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
#                     n += 1
#             if classificator == 'Secondary':
#                 if (df.loc[(df['Residue_id'] == int(i)) & (df['Chain'] == 'A'), classificator].values[0] == domain):
#                     sum += (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
#                     n += 1
#     if n != 0:
#         rmsd = (sum/n)**(1/2)
#     else:
#         rmsd = 0
#     return rmsd