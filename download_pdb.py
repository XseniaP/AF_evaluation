import os
import ssl
import sys
import urllib
from Bio import Entrez, SeqIO

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO

#
# $ wget -O - https://files.rcsb.org/download/1A80.pdb 2>/dev/null \
#    | python -c "import sys; from Bio import SeqIO; SeqIO.convert(sys.stdin, 'pdb-atom', sys.stdout, 'fasta')"

def fetch(id_list):

    Entrez.email = "kseniap@mail.tau.ac.il"
    path = "/Users/kseniapolonsky/Downloads"

    for my_id in id_list:
        name = my_id + ".pdb"
        full_path = path + "/" + my_id.upper()
        link = "http://files.rcsb.org/download/" + name
        os.mkdir(full_path)
        urllib.request.urlretrieve(link, full_path + "/" + name)

    return name

fetch(["1osp","5mo3","6plh","5sy8","4tpr","4ocr","4ocs","5jue","1q9o","5lgh","6obz","1mju","4qyo","3vfg","4qy8","1mj8","2ajx","2ajs","2aju","2ajv","3o41","5alc","3qpq","3qq9","6rco","5ob5","3tcl","1a6v","3sob","6p4c","4rbp","3se8","3mxv","3mxw","3nps","6azm"])
