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

        # try:
        #     handle = Entrez.efetch(db="protein", id=my_id, retmode="text")
        #     res = handle.read()
        #     print(res)
        #
        # except Exception:
        #     sys.stderr.write("Error! Cannot fetch: %s        \n" % my_id)

        # fasta = Entrez.efetch(db="protein", id=my_id, rettype="fasta")
        # fasta_record = SeqIO.read(fasta, "fasta")
        # name = my_id.upper() + ".fasta"
        # print(f'>{fasta_record.id}|{fasta_record.description}\n{fasta_record.seq}')
        # with open(name, 'w') as f:
        #     f.write(f'>{fasta_record.id}|{fasta_record.description}\n{fasta_record.seq}')
        #     f.close()
        # https: // www.rcsb.org / fasta / entry / 5MYO
        # link = "https://www.rcsb.org/fasta/entry/" + my_id.upper()
        # name = my_id.upper() + ".fasta"
        # urllib.request.urlretrieve(link, full_path + "/" + name)


    return name

# fetch(["4jfx","4jfz","5hi4","5hi5", "6k0y", "4hie", "6x1s","6x1u","6x1w","2fr4","6b14","3lex","3ley","3iu4","4n0y","5cp3","6dg2","5xli","5tpp","5l7x","3qcu","3sy0","5ewi","3cx5","6ell","6elj","6ele","5v6m","6pdr","6pds","6pdu","4jdv","4okv","5ik3","4r90","3v52"])
# fetch(["6k65","3mcl","4dgy","5n4j","4dgv","5k9j","5i30","7ean","7eam","5ye3","5ye4","5gz0","6fab","6lun","1q0x","6ca7","5drn","1rzf","1kiq","5ayu"])
# fetch(["5dr5","3qeg","5ur8","3t65","3cvi","4g6k","4g6m","6z2l","6pbv","5csz","1zea","6jp7","6xud","6xuk","3m8o","5w5z","1uac","1ua6","5i1l","5i1k"])
# fetch(["7km6","5i18","5i16","7kmi","7kmh","3fo9","6meh","6mee","6meg","6b5s","6b5r","6me1","6uud","6b5n","6b5m","6dkj","5en2","5obf","1n7m","5dtf","5xhg","5dt1","3kdm","3t4y","4xbe","7nx8","7nx2","4jha","3whx","4hgw"])
# fetch(["1mvu","2hvk","2w60","6wyt","6wyr","2fx7","5i76","6qcu","4k3d","5u6a","2vxt","2vxv","2vxq","3hc4","3hc0","3hc3","5vpg","5vpl","4jn1","4jn2","2xkn","3ntc","5jw5","1s3k","4r3s","8fab","6j9o","6ddm","6b9z","3s96","5u4r","4n8c","4x7s","7doh","2w9d","6tkd","6sf6"])
fetch(["1osp","5mo3","6plh","5sy8","4tpr","4ocr","4ocs","5jue","1q9o","5lgh","6obz","1mju","4qyo","3vfg","4qy8","1mj8","2ajx","2ajs","2aju","2ajv","3o41","5alc","3qpq","3qq9","6rco","5ob5","3tcl","1a6v","3sob","6p4c","4rbp","3se8","3mxv","3mxw","3nps","6azm"])