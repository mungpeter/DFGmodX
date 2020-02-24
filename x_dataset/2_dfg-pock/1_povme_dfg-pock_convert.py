#!/usr/bin/python

##
##  16.02.24
##  convert POVME-generated DFG-pocket PDB into multiple-frame PDB by adding
##  MODEL flag

import sys,re

if len(sys.argv) != 3: sys.exit('\n   > x.py [input PDB] [output PDB]\n')
fo = open(sys.argv[2], 'wh')
with open (sys.argv[1], 'rh') as fi:
  for idx, l in enumerate(fi):
    if re.search(r'REMARK Frame', l):
      if idx > 1:  fo.write('ENDMDL\n')
      fo.write(re.sub(r'REMARK Frame', 'MODEL', l))
    else:
      fo.write(l)
fo.write('ENDMDL')
