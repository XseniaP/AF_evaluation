from crowelab_pyir import PyIR
import json

FILE = 'rcsb_pdb_2X7L.fasta'

if __name__ == '__main__':

    pyirexample = PyIR(query=FILE, args = ['--sequence_type', 'prot', '--legacy', '--gzip', 'False', '--germlineV', './db/human_igh_v'])
    result = pyirexample.run()


    # f = open(result, 'r')
    # for line in f:
    #     entry = json.loads(line)


    #Prints jason name
    print(result)

