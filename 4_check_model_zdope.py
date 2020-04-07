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

  > need to be in the home directory where all kinase subdirectories are
  > and <conf>.<name>.zDOPE.txt are in the /1_result folder in kinase directory
  e.g.>  ---- /home_directory
                    |-- /AKT1
                          |-- /1_result
                                  |-- cido.AKT1.zDOPE.txt
                    |-- /AKT2
                          |-- /1_result
                                  |-- cido.AKT2.zDOPE.txt
                    |-- /AKT3
                    |...\n  
'''.format(sys.argv[0])
if len(sys.argv) != 4: sys.exit(msg)

import os,re

from Bio import SeqIO

##########################################################################

def main( fasta_file, conf, out_file ):

  prot_list = [f.id.split('|')[0] for f in SeqIO.parse(fasta_file, 'fasta')]

  if conf == 'codi':
    num = ['.1', '.2','.3','.4']
    [ ParseFile(prot_list, conf, out_file, n) for n in num ]
  else:
    num = ''
    ParseFile(prot_list, conf, out_file, num)


##########################################################################
def ParseFile( prot_list, conf, out_file, num ):
  finish = []
  ## Read in zDOPE file if it is there, collect the avg zDOPE for each protein
  for name in prot_list:
    zdope= []
    zdope_file = '{0}/1_result/{1}.{0}{2}.zDOPE.txt'.format(name, conf, num)
    if os.path.isfile(zdope_file):
      with open(zdope_file, 'r') as fi:
        ## skip first 2 rows, take the 5th element, which is the zDOPE score
        for i, l in enumerate(fi):
          if i < 2:
            continue 
          else:
            zdope.append(float(l.split()[5]))

      ## sort the zDOPE from lowest (good) to highest (bad), take best 5
      top_5 = sorted(zdope)[:5]
      if len(zdope) == 0:
        finish.append([ name, 10. ])  # zDOPE > 2.0 means failed result
      else:
        finish.append([ name, sum(top_5)/len(top_5) ]) # average of best 5
    else:
      finish.append([ name, 10. ])

  ## Sort protein with lowest zDOPE scores first
  ordered = sorted(finish, key=lambda x: x[1])
  with open(out_file+num, 'w') as fo:
    for p in ordered:
      fo.write('{0:10s}{1:8.4f}\n'.format(p[0],p[1]))


##########################################################################
if __name__ == '__main__':
  main( sys.argv[1], sys.argv[2], sys.argv[3])
