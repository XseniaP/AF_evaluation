import main_call
import run_TMalign
import rmsd_by_segment
import parse_superimposed
import predict_from_list
import sys
from predict_from_list import copy_files
from main_call import main_call_function



if __name__ == '__main__':

    fp = open('log.txt', 'w')
    fp.close()

    fp = open('exceptions.txt', 'w')
    fp.close()

    arguments = sys.argv
    my_list = arguments[1:]
    # my_list = sys.argv[1].replace('[', ' ').replace(']', ' ').replace(',', ' ').split()
    for name in my_list:
        try:
            FILE, CHAIN, PDB_F, PDB_PRED, PDB_ORIG = copy_files(name)
            parse_superimposed.chain_pdb(PDB_F, PDB_ORIG, CHAIN)
            main_call_function(name, FILE, CHAIN, PDB_F, PDB_PRED, PDB_ORIG)
        except Exception:
            fp = open('exceptions.txt', 'a')
            fp.write(str(name) + "\n")
            fp.close()
            continue
    sys.exit(0)



