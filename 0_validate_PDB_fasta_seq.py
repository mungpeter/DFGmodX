#!/usr/bin/env python3

import sys
import os,re
import pandas as pd

from tqdm import tqdm
from x_pir_multi_align import FASTA_Gen
from x_pir_multi_align import RunClustalO
from x_pir_multi_align import CacheSeqDatabase

msg = '''\n  > {0}
      [ List of PDB file to be validated ]
      [ No-gap FASTA sequences downloaded from RCSB for the PDBs ]
'''.format(sys.argv[0])
if len(sys.argv) != 3: sys.exit(msg)

def main( infile, pdb_seq ):

  df = pd.read_csv(infile, header=None, comment='#')
  df.columns= ['pdb_file']
  df['pdb_id'] = df.pdb_file.apply(lambda x: x.split('.')[0])
#  df['pdb_seq'] = df.pdb_file.apply(lambda x: str(FASTA_Gen(x, x.split('.')[0])))
  df['pdb_seq'] = [str(FASTA_Gen(row, row.split('.')[0])) for row in tqdm(df.pdb_file.values)]

  fasta = FASTA(fasta=CacheSeqDatabase(pdb_seq)[0])
  df['fas_seq'] = df.pdb_id.apply(fasta)
  df['pdb_len'] = df.pdb_seq.apply(len)
  df['fas_len'] = df.fas_seq.apply(len)
  diff_df = df[ df.pdb_len - df.fas_len != 0 ]

  if len(diff_df) > 0:
    ## use ClustalO to align sequences with different length
    ## write out a temp file with PDB that need to fix the FASTA
    for idx, inp in tqdm(diff_df.iterrows()):
      with open('_tmp_.{0}.fasta'.format(inp.pdb_id),'w') as fo:
        fo.write('>fasta {0}\n{1}\n'.format(inp.fas_len, inp.fas_seq))
        fo.write('>pdb {0}\n{1}'.format(inp.pdb_len, inp.pdb_seq))
      RunClustalO('_tmp_.{0}.fasta'.format(inp.pdb_id),'_tmp_.{0}.rst'.format(inp.pdb_id))
      os.system('rm _tmp_.{0}.fasta'.format(inp.pdb_id))

    print('\033[36m>> Look for "_tmp_.<pdb_id>.fasta" for FASTA differences\n\033[0m')

  print('\033[34m>> Number of Seqs that are different in length: \033[0m{0}'.format(len(diff_df)))


########################################
class FASTA(object):
  def __init__(self, fasta=''):
    self.fasta = fasta
  def __call__(self, idx):
    return self.get_fasta(idx)
  def get_fasta(self, idx):
    if idx in self.fasta:
      return str(self.fasta[idx].seq)
    else:
      return ''

########################################
if __name__ == "__main__":
  main(sys.argv[1], sys.argv[2])

##########################################
#
#  Peter M.U. Ung @ MSSN/Yale
#
#  v1  20.03.28
#
#  PDB sequence downloaded from RCSB reflects the full-length
#  sequence used in crystallography, but actually resolved
#  PDB structure can have unresolved residues, resulting in
#  Xtal FASTA sequence that is shorter by a few residues.
#  Sometimes there are extra residues in Xtal structures that
#  are not included in the published FASTA sequences too.
#
#  This script compares the FASTA sequences published in RCSB and
#  the xtal-FASTA extracted from the actual PDB structures.
#  If the FASTA sequences have different lengths that indicates a
#  discrepancy, they and aligned with Clustalo to show where
#  the missing/inserted residues are. Use this detection result
#  to correct the sequence in the FASTA database manually.
#