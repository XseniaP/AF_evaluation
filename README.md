# Evaluation of AlphaFold prediction accuracy

This script allows to evalute the prediction quality of the antibody heavy and light chains, to perform the calculations and combine all accuracy evaluation steps into a single pipeline. 

![Figure1_revised](https://github.com/XseniaP/AF_evaluation/assets/50076292/b61ee140-361e-4d7c-a37d-684aef8a440c)


The script is used to call the externally developed open-source products for alignment, superimposition, metrics calculation, and visualization, and then parses the results and moves the relevant information to the next step. In addition, the script provides on-demand calculations of the key metrics per domain and color coding, which are not covered by the externally developed packages. 

## Prerequisites

To run this script the following products and their prerequisites should be installed/compiled with the compilation files added to the main project folder:

* (1) **BioPython** - it is used to align the sequences
* (2) **TM-align** - the TM-align.cpp can be downloaded from the official website (https://zhanggroup.org/TM-align/),following the compilation the executable file should be saved in the root directory of the project under the name "TMalign"
* (3) **MM-align** - the MM-align.cpp can be downloaded from the official website (https://zhanggroup.org/MM-align/),following the compilation the executable file should be saved in the root directory of the project under the name "MMalign"
* (4) **PyIR** should be installed according to the guidelines found of the official page https://github.com/crowelab/PyIR, "from crowelab_pyir import PyIR" is called in the script. Please note that PyIR requires a set of BLAST germline databases. Please follow the official instructions to set it up
* (5) **DockQ** -   https://github.com/bjornwallner/DockQ
* (6) **ABangle** -  https://github.com/wjs20/ABangle/tree/main/abangle
* (7) **Others** - json, sys, subprocess, pathlib, pandas, numpy, shutil, re and other packages imported within the script

## Running the script
For each of the chains to be analyzed make sure that the following files are in the cwd: 
* (1) pdb file of the native chain structure with the filename in the format [name_chain.pdb] - e.g. 7mzi_H.pdb
* (2) fasta file with the sequence of the chain with the filename in the following format [name_chain.fasta] - e.g. 7MZI_H.fasta
* (3) pdb file with the prediction of the chain structure with the filename in the format [name_chain_ranked_0.pdb] - e.g. 7MZI_H_ranked_0.pdb

Run main_wrapper.py with the list of chains provided in a format [name_chain] as parameters separated by space, e.g.: 

python3 ./main_wrapper.py '6LDY_H 6LDY_L 6WAS_H 6WAS_L 6WKL_A 6WKL_B'

OR for a single chain

python3 ./main_wrapper.py 7n4j_L

You can expect the following output for a single chain: 

(1) files  
7n4j_L_TM_sup_all_m0 - TM-align output,   
7n4j_L.json - PyIR output,   
TM_sup_atm_colorcoded_cartoon.pml - color-coded PyMol file,   
output.csv - the results for each chain will be added as a separate line in the spreadsheet  

(2) console - below if the example of the output for a single chain (7n4j_L)

7n4j_L.fasta L 7n4j.pdb 7n4j_L_ranked_0.pdb 7n4j_L.pdb human

Score = 643.0:

target            0 QSVLMQPPSVSGAPGQRVTISCTGSSSNIGAGYDVHWYQQLPGTAPKLLIYGNNNRPSGV
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
query             0 QSVLMQPPSVSGAPGQRVTISCTGSSSNIGAGYDVHWYQQLPGTAPKLLIYGNNNRPSGV

target           60 PDRFSGSKSGTSASLAITGLQADDEADYYCQSYDSSLSGSKVFGGGTKLTVLGQPKAAPS
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
query            60 PDRFSGSKSGTSASLAITGLQADDEADYYCQSYDSSLSGSKVFGGGTKLTVLGQPKAAPS

target          120 VTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
query           120 VTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAA

target          180 SSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTE-- 216
                180 ||||||||||||||||||||||||||||||||||||-- 218
query           180 SSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS 218

 *********************************************************************
 * TM-align (Version 20210224): protein structure alignment          *
 * References: Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005) *
 * Please email comments and suggestions to yangzhanglab@umich.edu   *
 *********************************************************************
 
Name of Chain_1: 7n4j_L.pdb (to be superimposed onto Chain_2)

Name of Chain_2: 7n4j_L_ranked_0.pdb

Length of Chain_1: 216 residues

Length of Chain_2: 218 residues

Aligned length= 216, RMSD=   2.22, Seq_ID=n_identical/n_aligned= 1.000

TM-score= 0.88683 (if normalized by length of Chain_1, i.e., LN=216, d0=5.46)

TM-score= 0.87944 (if normalized by length of Chain_2, i.e., LN=218, d0=5.49)

(You should use TM-score normalized by length of the reference structure)

(":" denotes residue pairs of d <  5.0 Angstrom, "." denotes other aligned residues)

QSVLMQPPSVSGAPGQRVTISCTGSSSNIGAGYDVHWYQQLPGTAPKLLIYGNNNRPSGVPDRFSGSKSGTSASLAITGLQADDEADYYCQSYDSSLSGSKVFGGGTKLTVLGQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTE--
..::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::....::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  
QSVLMQPPSVSGAPGQRVTISCTGSSSNIGAGYDVHWYQQLPGTAPKLLIYGNNNRPSGVPDRFSGSKSGTSASLAITGLQADDEADYYCQSYDSSLSGSKVFGGGTKLTVLGQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS

Total CPU time is  0.19 seconds
Splitting input fasta file 7n4j_L.fasta
1 sequences successfully split into 1 pieces
Starting process pool using 8 processors
100%|██████████| 1/1 [00:04<00:00,  4.92s/seq]
1 sequences processed in 5.05 seconds, 0 sequences / s
Analysis complete, result file: 7n4j_L.json

": [{"gene": "IGLV1-40*01 unnamed protein product", "bit_score": 199.0, "e_value": 1e-68}, {"gene": "IGLV1-40*02 unnamed protein product", "bit_score": 198.0, "e_value": 3e-68}, {"gene": "IGLV1-40*03 unnamed protein product", "bit_score": 196.0, "e_value": 2e-67}], "FR1": {"from": 1.0, "to": 25.0, "length": 25.0, "matches": 24.0, "mismatches": 1.0, "gaps": 0.0, "percent identity": 96.0, "AA": "QSVLMQPPSVSGAPGQRVTISCTGS", "AA_Length": 25, "NT": "....T...................."}, "CDR1": {"from": 26.0, "to": 34.0, "length": 9.0, "matches": 9.0, "mismatches": 0.0, "gaps": 0.0, "percent identity": 100.0, "AA": "SSNIGAGYD", "AA_Length": 9, "NT": "........."}, "FR2": {"from": 35.0, "to": 51.0, "length": 17.0, "matches": 17.0, "mismatches": 0.0, "gaps": 0.0, "percent identity": 100.0, "AA": "VHWYQQLPGTAPKLLIY", "AA_Length": 17, "NT": "................."}, "CDR2": {"from": 52.0, "to": 54.0, "length": 3.0, "matches": 2.0, "mismatches": 1.0, "gaps": 0.0, "percent identity": 66.7, "AA": "GNN", "AA_Length": 3, "NT": "..S"}, "FR3": {"from": 55.0, "to": 90.0, "length": 36.0, "matches": 35.0, "mismatches": 1.0, "gaps": 0.0, "percent identity": 97.2, "AA": "NRPSGVPDRFSGSKSGTSASLAITGLQADDEADYYC", "AA_Length": 36, "NT": "............................E......."}, "CDR3": {"from": 91.0, "to": 99.0, "length": 9.0, "matches": 9.0, "mismatches": 0.0, "gaps": 0.0, "percent identity": 100.0}, "Total": {"from": "N/A", "to": "N/A", "length": 99.0, "matches": 96.0, "mismatches": 3.0, "gaps": 0.0, "percent identity": 97.0}, "Alignments": {"strings": ["<-------FR1-IMGT--------><CDR1-IM><---FR2-IMGT----><C><-------------FR3-IMGT------------->         ", "QSVLMQPPSVSGAPGQRVTISCTGSSSNIGAGYDVHWYQQLPGTAPKLLIYGNNNRPSGVPDRFSGSKSGTSASLAITGLQADDEADYYCQSYDSSLSG", "....T................................................S............................E................", "...VT................................................S............................E................", "...VT................................................S................A...........E................"], "keys": ["header", "Query_1", "IGLV1-40*01", "IGLV1-40*02", "IGLV1-40*03"]}, "NT-Trimmed": "QSVLMQPPSVSGAPGQRVTISCTGSSSNIGAGYDVHWYQQLPGTAPKLLIYGNNNRPSGVPDRFSGSKSGTSASLAITGLQADDEADYYCQSYDSSLSG"}
[91, 102]

CDR3 is: QSYDSSLSGSKV

{'CDR1': [26.0, 34.0], 'CDR2': [52.0, 54.0], 'CDR3': [91, 102], 'FR1': [1.0, 25.0], 'FR2': [35.0, 51.0], 'FR3': [55.0, 90.0]}

{'total': 2.246301958940681, 'CDR1': 2.1138627622856165, 'CDR2': 2.551029661398184, 'CDR3': 4.856094761568532, 'FR1': 3.1569509340501316, 'FR2': 2.754274666981576, 'FR3': 1.7425274029409117, 'therest': 1.5094966861884906, 'SHEET': 1.7471037957910018, 'HELIX': 1.752972760769545, 'SSBOND': 1.8750372874845358}

GDT Total Score is: 0.71875

GDT High Accuracy Score is: 0.474537037037037

DONE 
