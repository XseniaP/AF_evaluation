
# from google.cloud import storage
#
# bucket_name = 'alpha_f_eval'
# # Instantiates a client
# storage_client = storage.Client()
#
# # Get GCS bucket
# bucket = storage_client.get_bucket(bucket_name)
#
# # Get blobs in bucket (including all subdirectories)
# blobs_all = list(bucket.list_blobs())
# print(blobs_all)
#
# # Get blobs in specific subirectory
# # blobs_specific = list(bucket.list_blobs(prefix='path/to/subfolder/'))
import os
import pathlib
import re
from pandas import *

to_run = set()

# file_name1 = str(pathlib.Path.cwd()) + str("pocket_list.txt")
f = open('/Users/kseniapolonsky/rmsd/rmsd/pocket_list.txt', 'r')
my_list = list()


for line in f:
    if "fasta" in line:
        import re

        line_set = re.split('PocketDisk/|.fasta', line)
        # line_set = line.split('PocketDisk/')
        my_list.append(line_set[1])

# print(my_list)

# f = open('/Users/kseniapolonsky/Downloads/output_set2.csv', 'r')
data = read_csv('/Users/kseniapolonsky/Downloads/output_set2.csv', header=0)
# data['pdb_w_chain'] = data[['pdb', 'chain']].agg('_'.join, axis=1)
# data['pdb_w_chain'] = data['pdb'] + "_" + data['chain']
data['pdb_w_chain'] = data.iloc[:, 1].astype(str) + "_" + data.iloc[:, 2].astype(str)
pdb = set(data['pdb_w_chain'].to_list())
print(pdb)

for item in my_list:
    if item not in pdb:
        to_run.add(item)

my_list = [f.name for f in os.scandir('/Volumes/PocketDisk/NEW/AP2') if f.is_dir()]
# print(my_list)


for item in my_list:
    if item not in pdb:
        to_run.add(item)

my_list = [f.name for f in os.scandir('/Volumes/PocketDisk') if f.is_dir()]
# print(my_list)

for item in my_list:
    if item not in pdb:
        to_run.add(item)


my_list = [f.name for f in os.scandir('/Users/kseniapolonsky/Downloads/RES') if f.is_dir()]
# print(my_list)

for item in my_list:
    if item not in pdb:
        to_run.add(item)

print(to_run)







