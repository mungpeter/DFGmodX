#!/usr/bin/env python3

import io
import os
import re
import sys
import pandas as pd

from tqdm import tqdm
from pathos import multiprocessing

from Bio import PDB
from Bio.PDB import Select

from x_homolog_templ_check import BlastpPairwiseIdentity

#######################################################################################
#
#  v1.0  20.03.12
#
#
############################################################################################

def RunBlastDB( data ):
  kinase_db, pdb_id, chain_id, pdb_idx, seq = data

  fasta_file = '{0}_{1}.fasta'.format(pdb_id, chain_id)
  with open(fasta_file, 'w') as fo:
    fo.write('>'+pdb_idx+'\n'+seq+'\n')

  Identity = BlastpPairwiseIdentity( '.', fasta_file, kinase_db )
  if Identity is not None:
    imat_df = pd.DataFrame(Identity, columns=['kinase','length','ident','simi'])
  else:
    imat_df = None

  return [data, imat_df]


############################################################################################
## Check FASTA is kinase by comparing sequence identity to all canonical human kinases (kinome)
## output is a Dictionary of pdb_id to chain_id
def CheckKinaseSeqIdentity( uniq_f_df, kinase_db, len_cutoff, idt_cutoff, outpref ):

  # good:  length match > 175, seq ident > 40%
  # check: length match > 175, seq ident 30-40%
  # bad:   length match > 175, seq ident < 30%
  # no:    length match <= 175
  f_ok  = open(outpref+'.good_seq_ident.txt', 'w')
  f_ok.write('pdb_idx\tlength\tident\tsimi\n')
  f_chk = open(outpref+'.check_seq_ident.txt', 'w')
  f_chk.write('pdb_idx\tlength\tident\tsimi\n')
  f_bad = open(outpref+'.bad_seq_ident.txt', 'w')
  f_bad.write('pdb_idx\tlength\tident\tsimi\n')
  f_no  = open(outpref+'.no_seq_ident.txt', 'w')
  f_no.write('pdb_idx\tlength\tident\tsimi\n')

  new_pdb_ids = []    # collection of confirmed kinases
  non_kin_idx = []    # collection of non-kinase pdb_idx
  print('\n \033[34m## Comparing sequence identity of unique PDB FASTA ##\033[0m')

  ## Run Blastp in parallel
  Data = [[kinase_db, r.pdb_id, r.chain_id, r.pdb_idx, r.seq] for idx, r in uniq_f_df.iterrows()]
  mpi  = multiprocessing.Pool()
  df_l = [x for x in tqdm(mpi.imap(RunBlastDB, Data), total=len(Data))]
  mpi.close()
  mpi.join()

  for items in df_l:
    data, imat_df = items
    kinase_db, pdb_id, chain_id, pdb_idx, seq = data
    fasta_file = '{0}_{1}.fasta'.format(pdb_id, chain_id)

    if imat_df is None:
      f_no.write(fasta_file+'\t no output\n')
      continue

    ## Required matching to at least 200 residues as cutoff, if fewer than 200 res,
    ## unlikely to be a kinase catalytic domain. if seq ident < 40% of any 
    ## known human kinases, need to check and confirm if it is a bacterial
    ## or viral kinases; mammalian kinases are very similar, even to chicken/fish
    sele_df = imat_df[ imat_df.length > (len_cutoff - 40) ]
    info = '{0}_{1}\t{2}\t{3}\t{4}\n'.format(pdb_id, chain_id, 
                imat_df.length[0], imat_df.ident[0], imat_df.simi[0])

    ## catch a strange error with ident.iloc[0] key error
    try:
      test = sele_df.ident
    except KeyError:
      print('  \033[31mERROR:\033[0m "sele_df.ident[0]" - '+pdb_idx)
      continue

    if len(sele_df) == 0 or len(sele_df.ident) == 0:
      non_kin_idx.append('{0}_{1}'.format(pdb_id, chain_id))
      f_no.write(info)
      continue
    elif sele_df.ident.iloc[0] > idt_cutoff:
      f_ok.write(info)  
      new_pdb_ids.append( [pdb_id, chain_id] )
    elif sele_df.ident.iloc[0] > (idt_cutoff - 10.):
      f_chk.write(info)
    else:
      non_kin_idx.append('{0}_{1}'.format(pdb_id, chain_id))
      f_bad.write(info)

  f_bad.close()
  f_chk.close()
  f_ok.close()
  f_no.close()

  nonkin_seq_df = uniq_f_df[ uniq_f_df.pdb_idx.isin(non_kin_idx) ]

  return new_pdb_ids, nonkin_seq_df


########################################################################
