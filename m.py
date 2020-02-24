#!/usr/bin/python

import sys,re,os
from Bio import SeqIO

data = list(SeqIO.parse('stdy_kinase.raw.clean.nogap.151102.fasta', 'fasta'))
sets = {}
for fasta in data:
  sets[fasta.id.split('|')[0]] = fasta
names = [fasta.id.split('|')[0] for fasta in data]
print(len(names))
print(names)

with open('3.list', 'r') as fi:
  pdb = [l.rsplit()[0] for l in fi]
print(len(pdb))

out = list(set(names).intersection(set(pdb)))
print(len(out))

with open('xxxxxx.fasta', 'w') as fo:
  for name in out:
    if sets[name]:
      SeqIO.write(sets[name], fo, "fasta")

