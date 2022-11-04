import sys
import subprocess
import pathlib
import re
from subprocess import check_output

" *********************************************************************\n"
" * References: Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005) *\n"
" * TM-align arguments: original_model_pdb predicted_model_pdb [options] *\n"
" * take from the program arguments -or, -pr                             *\n"
" *********************************************************************"

arguments = sys.argv

# arg_dict = {}
# for i in range(len(arguments)):
#     if i%2 != 0:
#         key = arguments[i]
#     elif (i%2 == 0) & (i != 0):
#         arg_dict[key] = arguments[i]
#
#
# FILE = arg_dict['-fa']
# CHAIN = arg_dict['-ch']
# PDB_F = arg_dict['-fll']
# PDB_PRED = arg_dict['-pr']
# PDB_ORIG = arg_dict['-or']
# PDB_PRED_NEW = str(pathlib.Path.cwd()) + '/' + "x_" + PDB_PRED

# FILE = "%s.fasta" % (arguments[1])
# CHAIN = arguments[1].split("_")[1]
# PDB_F = "%s.pdb" % (arguments[1].split("_")[0].lower())
# PDB_PRED = arguments[1] + "_ranked_0.pdb"
# PDB_ORIG = "%s.pdb" % (arguments[1].split("_")[0].lower() + "_" + CHAIN)
# PDB_PRED_NEW = str(pathlib.Path.cwd()) + '/' + "x_" + PDB_PRED

def run_TMalign(mode, name):
    FILE = "%s.fasta" % (name)
    CHAIN = name.split("_")[1]
    PDB_F = "%s.pdb" % (name.split("_")[0].lower())
    PDB_PRED = name + "_ranked_0.pdb"
    PDB_ORIG = "%s.pdb" % (name.split("_")[0].lower() + "_" + CHAIN)
    PDB_PRED_NEW = str(pathlib.Path.cwd()) + '/' + "x_" + PDB_PRED

    if mode == 0:
        set_args = [None] * 6
        set_args[2] = PDB_PRED
    if mode == 1:
        set_args = [None] * 8
        set_args[6] = '-byresi'
        set_args[7] = '1'
        set_args[2] = PDB_PRED_NEW

    exec_name = str(pathlib.Path.cwd()) + "/TMalign"
    set_args[0] = exec_name
    set_args[1] = PDB_ORIG
    set_args[3] = '-o'
    set_args[4] = 'TM_sup'
    set_args[5] = '-v'


    # flag = 'false'
    # arguments = sys.argv
    # exec_name = str(pathlib.Path.cwd()) + "/TMalign"
    # set_args[0] = exec_name
    # for argument in arguments:
    #     if flag == 'true_original':
    #         set_args[1] = argument
    #         flag = 'false'
    #     if flag == 'true_prediction':
    #         set_args[2] = argument
    #         flag = 'false'
    #     if argument == '-or':
    #         flag = 'true_original'
    #     elif argument == '-pr':
    #         flag = 'true_prediction'
    #     else:
    #         flag = 'false'

    # set_args[3] = '-o'
    # set_args[4] = 'TM_sup'
    # set_args[5] = '-v'
    # set_args[6] ='-byresi'
    # set_args[7] = '1'


    res = subprocess.run(set_args, capture_output=True, text=True)
    p = subprocess.Popen(set_args, stdout=subprocess.PIPE)
    out, err = p.communicate()
    out = out.decode('utf-8')
    print(out)
    print("\n")

    fp = open('log.txt', 'a')
    fp.write(out)
    fp.write("\n")
    fp.close()

    # p.kill()
    p.terminate()



    return out


"Additional options:\n"
"    -dir     Perform all-against-all alignment among the list of PDB\n"
"             chains listed by 'chain_list' under 'chain_folder'. Note\n"
"             that the slash is necessary.\n"
"             $ TMalign -dir chain_folder/ chain_list\n"
"\n"
"    -dir1    Use chain2 to search a list of PDB chains listed by 'chain1_list'\n"
"             under 'chain1_folder'. Note that the slash is necessary.\n"
"             $ TMalign -dir1 chain1_folder/ chain1_list chain2\n"
"\n"
"    -dir2    Use chain1 to search a list of PDB chains listed by 'chain2_list'\n"
"             under 'chain2_folder'\n"
"             $ TMalign chain1 -dir2 chain2_folder/ chain2_list\n"
"\n"
"    -suffix  (Only when -dir1 and/or -dir2 are set, default is empty)\n"
"             add file name suffix to files listed by chain1_list or chain2_list\n"
"\n"
"    -atom    4-character atom name used to represent a residue.\n"
"             Default is \" CA \" for proteins\n"
"             (note the spaces before and after CA).\n"
"\n"
"    -ter     Strings to mark the end of a chain\n"
"             3: (default) TER, ENDMDL, END or different chain ID\n"
"             2: ENDMDL, END, or different chain ID\n"
"             1: ENDMDL or END\n"
"             0: (default in the first C++ TMalign) end of file\n"
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             0: (default) treat the whole structure as one single chain\n"
"             1: treat each MODEL as a separate chain (-ter should be 0)\n"
"             2: treat each chain as a separate chain (-ter should be <=1)\n"
"\n"
"    -outfmt  Output format\n"
"             0: (default) full output\n"
"             1: fasta format compact output\n"
"             2: tabular format very compact output\n"
"            -1: full output, but without version or citation information\n"
"\n"
"    -byresi  Whether to assume residue index correspondence between the\n"
"             two structures.\n"
"             0: (default) sequence independent alignment\n"
"             1: (same as TMscore program) sequence-dependent superposition,\n"
"                i.e. align by residue index\n"
"             2: (same as TMscore -c, should be used with -ter <=1)\n"
"                align by residue index and chain ID\n"
"             3: (similar to TMscore -c, should be used with -ter <=1)\n"
"                align by residue index and order of chain\n"
"\n"
"    -TMcut   -1: (default) do not consider TMcut\n"
"             Values in [0.5,1): Do not proceed with TM-align for this\n"
"                 structure pair if TM-score is unlikely to reach TMcut.\n"
"                 TMcut is normalized is set by -a option:\n"
"                 -2: normalized by longer structure length\n"
"                 -1: normalized by shorter structure length\n"
"                  0: (default, same as F) normalized by second structure\n"
"                  1: same as T, normalized by average structure length\n"
"\n"
"    -mirror  Whether to align the mirror image of input structure\n"
"             0: (default) do not align mirrored structure\n"
"             1: align mirror of chain1 to origin chain2\n"
"\n"
"    -het     Whether to align residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
"\n"
"    -infmt1  Input format for chain1\n"
"    -infmt2  Input format for chain2\n"
"            -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             1: SPICKER format\n"
"             2: xyz format\n"
"             3: PDBx/mmCIF format\n"
"\n"
"Usage: TMalign PDB1.pdb PDB2.pdb [Options]\n"
"\n"
"Options:\n"
"    -u    TM-score normalized by user assigned length (the same as -L)\n"
"          warning: it should be >= minimum length of the two structures\n"
"          otherwise, TM-score may be >1\n"
"\n"
"    -a    TM-score normalized by the average length of two structures\n"
"          T or F, (default F)\n"
"\n"
"    -i    Start with an alignment specified in fasta file 'align.txt'\n"
"\n"
"    -I    Stick to the alignment specified in 'align.txt'\n"
"\n"
"    -m    Output TM-align rotation matrix\n"
"\n"
"    -d    TM-score scaled by an assigned d0, e.g. 5 Angstroms\n"
"\n"
"    -o    Output the superposition to 'TM_sup*'\n"
"            $ TMalign PDB1.pdb PDB2.pdb -o TM_sup\n"
"          View superposed C-alpha traces of aligned regions by RasMol or PyMOL:\n"
"            $ rasmol -script TM_sup\n"
"            $ pymol -d @TM_sup.pml\n"
"          View superposed C-alpha traces of all regions:\n"
"            $ rasmol -script TM_sup_all\n"
"            $ pymol -d @TM_sup_all.pml\n"
"          View superposed full-atom structures of aligned regions:\n"
"            $ rasmol -script TM_sup_atm\n"
"            $ pymol -d @TM_sup_atm.pml\n"
"          View superposed full-atom structures of all regions:\n"
"            $ rasmol -script TM_sup_all_atm\n"
"            $ pymol -d @TM_sup_all_atm.pml\n"
"          View superposed full-atom structures and ligands of all regions\n"
"            $ rasmol -script TM_sup_all_atm_lig\n"
"            $ pymol -d @TM_sup_all_atm_lig.pml\n"
"\n"
" -fast    Fast but slightly inaccurate alignment by fTM-align algorithm\n"
"\n"
"   -cp    Alignment with circular permutation\n"
"\n"
"    -v    Print the version of TM-align\n"
"\n"
"    -h    Print the full help message, including additional options\n"
"\n"
"    (Options -u, -a, -d, -o will not change the final structure alignment)\n\n"
"Example usages:\n"
"    TMalign PDB1.pdb PDB2.pdb\n"
"    TMalign PDB1.pdb PDB2.pdb -u 100 -d 5.0\n"
"    TMalign PDB1.pdb PDB2.pdb -a T -o PDB1.sup\n"
"    TMalign PDB1.pdb PDB2.pdb -i align.txt\n"
"    TMalign PDB1.pdb PDB2.pdb -m matrix.txt\n"
"    TMalign PDB1.pdb PDB2.pdb -fast\n"
"    TMalign PDB1.pdb PDB2.pdb -cp\n"
