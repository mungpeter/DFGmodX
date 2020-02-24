#!/usr/bin/env python3

import re,sys,os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from x_pir_multi_align import CacheSeqDatabase

##########################################################################

msg = '''
\t> {0}
\t\t[ Typical Kinase Seq Alignment | fasta ]
\t\t[ Edited Sequence for replace  | fasta ]
\t\t[ Output filename ]\n'''.format(sys.argv[0])
if len(sys.argv) != 4: sys.exit(msg)

##########################################################################
#
# primarily used to append the manual editing to the activation loop segment
#

Kinases, k_order = CacheSeqDatabase(sys.argv[1]) 
ALoops,  l_order = CacheSeqDatabase(sys.argv[2])

print('# Kinase Seq: '+str(len(Kinases)))
print('# Replace:    '+str(len(ALoops)))

loop_keys = ALoops.keys()
kin_keys  = Kinases.keys()

l_key = {}
for k in kin_keys:
  k_head = k.split('/')[0]
  for l in loop_keys:
    if re.search(r'{0}'.format(k_head), l):
      l_key[k] = l
      break 

record = []  
for k_key in k_order:
  kinase = Kinases[k_key]
  loop   = ALoops[l_key[k_key]]
  
  newseq = kinase.seq[ :2140] + loop.seq + kinase.seq[2370: ]
  
  itm = SeqRecord( id=k_key, seq=newseq, description=kinase.description,
                   name=kinase.name )
  record.append(itm)

print(len(record))
with open(sys.argv[3], 'w') as fo:
  SeqIO.write(record, fo, 'fasta')
