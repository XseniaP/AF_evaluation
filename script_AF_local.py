import urllib
# import pip
# pip.main(['install','pandas'])
import pandas as pd
import urllib.request
import subprocess
from Bio import Entrez, SeqIO, PDB


def get_seqeunce(code, chain_id):
    PDB_file_path = f'{code}.pdb'
    query_chain_id = f'{code}:{chain_id}'
    # PDB_file_path = "%s.pdb" % (code)
    # query_chain_id = "%s:%s" % (code, chain_id)
    chain = {record.id: record.seq for record in SeqIO.parse(PDB_file_path, 'pdb-seqres')}
    query_chain = chain[query_chain_id]
    print('done')
    return query_chain


def run_AF():
    filename = 'rmsd_results_130.xlsx'
    sheetname = 'Sheet1'
    df = pd.read_excel(filename, sheet_name=sheetname)
    pdb_code = df['pdb'].tolist()
    for code in pdb_code:
        parser = PDB.PDBParser()
        pdf_file = urllib.request.urlretrieve(f"http://files.rcsb.org/download/{code}.pdb", f"FILES/{code}.pdb")
        # pdf_file = urllib.request.urlretrieve("http://files.rcsb.org/download/%s.pdb" % (code), "%s.pdb" % (code))
        structure = parser.get_structure(code, f"{code}.pdb")
        # structure = parser.get_structure(code, "%s.pdb" % (code))
        # for model in structure:
        #     for chain in model:
        #         print(chain)
        #         fasta = get_fasta(code, chain)
        counter = 0
        for chain in structure.get_chains():
            ID = chain.id
            sequence = get_seqeunce(code, ID)
            if counter <= 1:
                fasta_name = f"FILES/{code}_{chain.id}.fasta"
                # fasta_name = "%s_%s.fasta" % (code, chain.id)

                with open(fasta_name, 'w') as out:
                    line1 = '>1'
                    out.write('{}\n{}\n'.format(line1, sequence))
                    out.close()
                counter += 1
                # rc = subprocess.call(f"qsub script.txt ~\AF2\{code} ~\{fasta_name}")
                # rc = subprocess.call("qsub script.txt ~\AF2\%s ~\%s" % (code, fasta_name))
            else:
                break

        # pdf_file = urllib.request.urlretrieve(f'http://files.rcsb.org/download/{code}.pdb', f'{code}.pdb')
        # fasta = urllib.request.urlretrieve(f'https://www.rcsb.org/fasta/entry/{code}', f'{code}.fasta')
        # fasta = get_fasta(code, chain_id)


if __name__ == '__main__':
    run_AF()
