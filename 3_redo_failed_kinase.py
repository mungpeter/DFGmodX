#!/usr/bin/env python3

# extract sequence from database, based on a list of kinases marked "failed"

import re,sys,os
from Bio import SeqIO
from x_pir_multi_align import CacheSeqDatabase

failed_kinase = sys.argv[1]
fasta_database = sys.argv[2]
fasta_out = sys.argv[3]

with open(failed_kinase, 'r') as fi:
  Failed = [x for x in fi if not None]

# Cache the fasta database
Database, db_order = CacheSeqDatabase( fasta_database )
Data_Keys = Database.keys()

Output = []
for name in Failed:
  for key in Data_Keys:
    if re.search(r'{0}'.format(key.split('/')[1]), name):
      print(key)
      Output.append(Database[key])

fo = open(fasta_out, 'w')
for fasta in Output:
  SeqIO.write(fasta, fo, 'fasta')

