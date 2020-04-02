#!/usr/bin/env python3

##########################################################################
#
#  Peter M.U. Ung @ MSSM/Yale
#
#  v1   20.04.01
#
#  Purpose: go through the list of model directories to look at the homology
#      model zDOPE scores and detect if any directories did not generate 
#      model. Usually models with zDOPE > -0.3 will need some checking.
#
#      e.g. KDR/PDGFR/FLTs have large insertion between C/D helices, either
#           need to go change the .pir to redo it, or ignore. Many Others
#           kinases will have high zDOPE since it is structurally different
#           from most other typical families.
#
##########################################################################

import sys
msg = '''\n  > {0}
      [ FASTA database with Protein Names used to generate Model ]
      [ Conformation (used in filename) ]
      [ Output filename ]
'''.format(sys.argv[0])
if len(sys.argv) != 4: sys.exit(msg)

import os,re
import pandas as pd

from Bio import SeqIO

##########################################################################

def main( fasta_file, conf, out_file ):

  prot_list = [f.id.split('|')[0] for f in SeqIO.parse(fasta_file, 'fasta')]

  failed = []
  finish = []
  ## Read in zDOPE file if it is there, collect the avg zDOPE for each protein
  for name in prot_list:
    zdope= []
    if os.path.isfile('{0}/1_result/{1}.{0}.zDOPE.txt'.format(name, conf)):
      with open('{0}/1_result/{1}.{0}.zDOPE.txt'.format(name, conf), 'r') as fi:
        # only collect zDOPE from first 5 rows of ranked zDOPE data
        for i, l in enumerate(fi):
          if i < 2 or i > 6:
            continue 
          else:
            zdope.append(float(l.split()[5]))
        
        if len(zdope) == 0:
          failed.append(name)
        else:
          finish.append([ name, sum(zdope)/len(zdope) ])
    else:
      failed.append(name)

  ## Sort protein with worst zDOPE scores first
  ordered = sorted(finish, key=lambda x: x[1])
  with open(out_file, 'w') as fo:
    for p in ordered:
      fo.write('{0:10s}{1:8.4f}\n'.format(p[0],p[1]))

  ## List those failed to model
  with open('failed.'+out_file, 'w') as fo:
    for f in failed:
      fo.write('{0}\n'.format(f))


##########################################################################
if __name__ == '__main__':
  main( sys.argv[1], sys.argv[2], sys.argv[3])
