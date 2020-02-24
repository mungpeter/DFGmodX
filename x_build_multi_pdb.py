#!/usr/bin/env python3

##########################################################################
#
#   Peter M.U. Ung @ MSSM
#
#   v1.0 - 14.05.09
#   v2.0 - 16.07.15 - minor change in 'print'
#   v3.0 - 17.07.14 - change into Object
#
#   Purpose: Read in a list of PDB files (.gz and .bz2 accepted) and write
#            out a multi-model PDB file.
#
##########################################################################

import sys,re
from CommonUtility import *

msg = '''\n    > {0}\n        [List of PDB: .list] [Output PDB]
'''.format(sys.argv[0])
#if len(sys.argv) != 3: sys.exit(msg)

def BuildMultiPDB( mdl_list, output_name ):
  with open(mdl_list, 'r') as fi, open(output_name, 'w') as fo:
    for idx, pdb in enumerate(fi):
      print(' -- Write structure [{0}] into {1} --'.format(idx+1, output_name))
      fo.write('MODEL {0}\n'.format(idx+1))
      handle = file_handle(pdb.rstrip())
      for line in handle: 
        fo.write(line)
      fo.write('ENDMDL\n')
