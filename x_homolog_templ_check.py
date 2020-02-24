#!/usr/bin/env python3

import re,sys,os
import os.path
import pandas as pd

from Bio import SeqIO

from x_pir_multi_align import FASTA_Gen
from x_pir_multi_align import CacheSeqDatabase

##########################################################################
# Calculate percent identity among the input sequence and available structure
# sequence to identify the most similar kinase. Use that kinase name/structure
# to proceed
def SearchKinaseStruct(
        script_directory, home_directory, work_directory, result_directory,
        pdb_directory, kinase_profile, ident_thres, template_pdb,
        mdl_prot_fasta, Settings ):

  Vars = [ 'ScriptDirectory', 'HomeDirectory', 'WorkingDirectory',
           'ResultDirectory', 'KinaseProfile', 'KinaseStructInput', 
           'ModelKinaseFasta' ]

  print('\n  ** KinaseStructInput is None - search for structure or homolog **')
  print('  -- Initiate Kinase Structure Search in Database --')
  for var in Vars:
    if Settings[var] is None:
      sys.exit('\n  > #2# ERROR: \'{0}\' is not specified: {0}: {1}'.format(
                      var, mdl_prot_fasta ))

  #################################################
  print(os.getcwd())
  os.chdir(work_directory)

  # Pairwise comparison of query sequence against a kinase Fasta database
  # Output Data lists sequence identity in descending order
  Data = FastaDatabaseSearch(result_directory, mdl_prot_fasta, kinase_profile)

  # If the input Fasta sequence has corresponding structure in database, use
  # that structure, otherwise, override the fasta ID and select the top 
  # scoring kinase and output the path to the structure
#  for name in list(zip(*Data))[0]:
#    if re.search(mdl_prot_fasta, name):
#      k_name = mdl_prot_fasta
#    else:
#      k_name = Data[0][0]

  # If the sequence length is < 200, the true percent identity for the full
  # kinase would be uncertain. If the best structure is length < 200, use the
  # the next best structure
#  Top = Data[0]
  for idx in range(0, len(Data)):
    if Data[idx][1] > 200:
      Top    = Data[idx]
      k_name = Data[idx][0]
    # If the input Fasta sequence has corresponding structure in database, 
    # e.g. 1ATP_E.fasta, use that structure, otherwise, override the fasta name
    # and select the top scoring kinase and output the path to the structure
    # To avoid rare cases of fasta filename matching part of PDBID, use strict
    # comparison instead of re.search()
      if re.search(r'.fasta', mdl_prot_fasta):
        for name in list(zip(*Data))[0]:
          if mdl_prot_fasta.split('/')[-1].split('.fasta')[0] == name:
            k_name = mdl_prot_fasta.split('/')[-1].split('.fasta')[0]
      break

  # Locate the chosen structure in the database and return the path
  Paths = pdb_directory.split(',')
  for path in Paths:
    if os.path.isfile('{0}/{1}.{2}'.format(path, k_name, template_pdb)):
      print('\n  > Use {0}.{1} as structure template\n  > Length: {2}\n  > Percent Identity: {3:4.1f}\tPercent Similarity: {4:4.1f}\n'.format(k_name, template_pdb, Top[1], Top[2], Top[3]))
    else:
      continue

  # Check percentage identity between input target sequence and input template 
  # structure. Pass if better than threshold % identity. Warning and quit if 
  # below.
  #  percent_ident = TCoffeePercentIdentity('_TEMP.homolog.fasta')
    if Data[0][2] > ident_thres:
      print('  SUCCESS: Protein structure input: {0}: {2}.{1}\n  {2} has enough sequence identity to the query sequence\n    Seq Identity: {3:4.1f} %\n'.format(k_name, template_pdb, mdl_prot_fasta, Top[2]))
      print('  ** Use the following Homolog Structure as Template **\n{0}/{1}.{2}'.format(path, k_name, template_pdb))
      return '{0}/{1}.{2}'.format(path, k_name, template_pdb), Data[0][2]
    else:
      print('\n  > #2# WARNING: Input sequence and template structure {0}.{1}: {4}\n  > #2# WARNING: have low sequence identity: {2:4.1f} %\n  > #2# WARNING: Choose another template structure that is a closer homolog to the target kinase\n  > #2# WARNING: and with higher sequence identity (> {3}%)'.format(
                  k_name, template_pdb, Data[0][2], ident_thres, mdl_prot_fasta ))

  #    sys.exit('\n  #2# WARNING: Input sequence and template structure {0}.{1}\n       have low sequence identity: {2:4.1f} %\n       Choose another template structure that is a closer homolog to the target kinase\n       and with higher sequence identity (> {3}%)'.format(k_name, template_pdb, Data[0][2], ident_thres))
      print('\n  > #2# WARNING: Input sequence and template structure {0}.{1}: {4}\n  > #2# WARNING: have low sequence identity: {2:4.1f} %\n  > #2# WARNING: Choose another template structure that is a closer homolog to the target kinase\n  > #2# WARNING: and with higher sequence identity (> {3}%)'.format(
                  k_name, template_pdb, Data[0][2], ident_thres, mdl_prot_fasta ))
      return '{0}/{1}.{2}'.format(path, k_name, template_pdb), Data[0][2]

  sys.exit('  > #2# ERROR: {4} Cannot find protein structure of seq: "{0}.{1}" in directories:\n  > #2# ERROR: "{2}"\n  > #2# ERROR: "{3}"'.format(
                k_name, template_pdb, Paths[0], Paths[1], mdl_prot_fasta ))



##########################################################################
# Use Blastp to generate pairwise percent identity between a query sequence
# and a database of sequences. Do not generate a matrix of pairwise identity
# to save on time
def FastaDatabaseSearch( result_directory, mdl_prot_fasta, kinase_profile ):

  print('\n  ** Search Fasta Database for similar kinases **')
  # Check if query fasta file exists. If not, find in fasta database and
  # print it out for blastp
  if re.search(r'.fasta', mdl_prot_fasta):
    if not os.path.isfile(mdl_prot_fasta):
      sys.exit('\n  > #2# ERROR: {0} (file) does not exist.'.format(mdl_prot_fasta))
  else:
    Database, db_order = CacheSeqDatabase(kinase_profile)
    print('  # Number of Fasta in Database: '+str(len(Database))+'\n')
    if Database[mdl_prot_fasta]:
      print('  Use FASTA: {0}'.format(mdl_prot_fasta))
      print('        Seq: {0}'.format(Database[mdl_prot_fasta].seq))
      with open(mdl_prot_fasta+'.fasta', 'w') as fo:
        SeqIO.write(Database[mdl_prot_fasta], fo, 'fasta')
    else:
      sys.exit('\n  > #2# ERROR: {0} does not exist in fasta database.'.format(mdl_prot_fasta))

  return BlastpPairwiseIdentity( result_directory, 
              mdl_prot_fasta, kinase_profile )


##########################################################################
# Calculate sequence identity and similarity of a query seq to a library of 
# sequence (or single seq) and output a list with the best one at the 1st row
def BlastpPairwiseIdentity( result_directory, mdl_prot_fasta, kinase_profile ):

  # If input Fasta is a file, reconfigure to only the fasta name  
  if os.path.isfile(mdl_prot_fasta):
    fasta_name = mdl_prot_fasta.split('.fasta')[0]
  else:
    fasta_name = mdl_prot_fasta

  print('\n  ** Calculate Sequence Identity between Query and Profile Sequences **')
  print('  Query Fasta:    '+fasta_name+'.fasta')
  print('  Kinase Profile: '+kinase_profile)
  # blastp to output: Name, AA_length, percent identity, percent positive
  # result in .csv format, omit other irrelevant data
  print('blastp -query "{0}" -subject "{1}" -max_target_seqs 5000 -out "{2}"/{3}.idmat.txt -outfmt "6 sseqid length pident ppos"'.format(fasta_name+'.fasta', kinase_profile, result_directory, fasta_name.split('/')[-1]))
  os.system('blastp -query "{0}" -subject "{1}" -max_target_seqs 5000 -out "{2}"/{3}.idmat.txt -outfmt "6 sseqid length pident ppos"'.format(fasta_name+'.fasta', kinase_profile, result_directory, fasta_name.split('/')[-1]))

  # Parse percent identity result generated by BlastP. Did not use clustalo or
  # t_coffee because they do redundent pairwise identity calculation for other 
  # kinases to create a true matrix and that is not needed; only need 1 set of
  # pairwise identity between query sequence and the database sequences
  if not os.path.isfile('{0}/{1}.idmat.txt'.format(
                        result_directory, fasta_name.split('/')[-1])):
    print('\n  > #2# Alignment Warning: Cannot find Blastp output. Seq identity to kinase too low? '+fasta_name)
    return None
  elif os.stat('{0}/{1}.idmat.txt'.format(
                result_directory, fasta_name.split('/')[-1])).st_size == 0:
    print('\n  > #2# Alignment Warning: Blastp failed. Seq Identity to kinase too low? '+fasta_name)
    return None


  ## Extract the identity information from Blastp result; sometimes a chain is
  ## broken into fragments, need to combine them according to residue ratios
  Ident = {}
  with open('{0}/{1}.idmat.txt'.format(result_directory, fasta_name.split('/')[-1]), 'rU') as fi:
    for line in fi:
      Items = line.split('\t')
      name, aa, identity, positive = (Items[0].split('|')[0], int(Items[1]), 
                                      float(Items[2]), float(Items[3]) )
      if name in Ident:
        Ident[name].append([name, aa, identity, positive])
      else:
        Ident[name] = [ [name, aa, identity, positive] ]

  # Convert dictionary into Tulip data. If a Fasta name has multiple lines, 
  # the alignment/identity calculation is broken down into pieces for 1 seq.
  # Summarize the pieces into 1 by adding up the ratio
  Data = []
  for name in Ident:
    length = sum(list(zip(*Ident[name]))[1])  # rearrange tulip groups 
    x, y = 0.0, 0.0
    nm   = name.split('_')
    if enumerate(nm) != 2:
      nm.append('A')
    for row in Ident[name]:
      x += row[1] * row[2]
      y += row[1] * row[3]

    Data.append( [nm[0], nm[1], length, (x/length), (y/length)] )

  # sort the dataset by percent identity or positive, then by available length,
  # then by filename to prefer A or B, etc
  pdata = pd.DataFrame(Data)  
  pdata.columns = ['pdb_id', 'chain', 'length', 'identity', 'similarity']  
  pdata['pdb_full'] = pdata['pdb_id']+'_'+pdata['chain']

  pdata = pdata.sort_values( by=['identity', 'length', 'chain'], 
                              ascending=[False, False, True] )
  pdata_temp = pdata.drop('pdb_id',1).drop('chain',1)
  col = pdata_temp.columns.tolist()
  col = col[-1:] + col[:-1]
  pdata_temp = pdata_temp[col]
  pdata=pdata_temp

  pdata.to_csv('{0}/{1}.idmat.sort.txt'.format(result_directory, fasta_name.split('/')[-1]), sep='\t', encoding='utf-8', float_format='%4.2f', index=False)

  Data = []
  for index, row in pdata.iterrows():
    Data.append( [row['pdb_full'], row['length'], row['identity'], 
                  row['similarity'] ])

  return Data


##########################################################################
##
##  Peter M.U. Ung @ MSSM
##  
##  v1.0    16.11.25 -
##  v2.0    16.11.27 -	added function to find kinase structure most related
##			to the input kinase sequence within the internal 
##			structure database
##  v3.0    16.12.05 -  remove CheckHomologTemplate as it superseded by
##                      SearchKinaseStruct
##  v3.1    17.03.28 -  fix bug with pdb_directory path
##  v4.0    17.03.28 -  fix bug with hierarchical sorting of ident, len, pdb_id
##                      by using Pandas dataframe
##  v5.0    17.07.14    fix bug when comparing fasta file name to existing PDB
##                      when fasta name is a part of PDB name
##  v6.0    18.03.28    update blastp to use precompiled fasta database

##  Originally thought this script would take in a target sequence and a 
##  crystal structure of a close homolog as template to generate a primary
##  model of the target sequence. But by using a primary model to gererate
##  a secondary DFG-out model of the target sequence, that would be a
##  homology model of another homology model, which would be introducing 
##  another layer of error/noise. 
##  Since the idea is to use a crystal structure of a close homolog as the 
##  primary template anyways, the template with high sequence identity, 
##  preferrably > 50%, should have sufficiently similar structure as the
##  target, hence the homolog template should be used directly in the DFG-out
##  homology modeling.
##  This script will now intend to check the sequence identity between the
##  input sequence and the input homolog template. The default threshold for 
##  similarity is set at 50%.

##  originally thought of using Biopython 1.61+ SeqIO (pdb-atom) function to 
##  read in a PDB file and generate the Fasta sequence, but this SeqIO actually
##  just read the HEADER and SEQRES handles to get the fasta. For PDB without
##  those handles, SeqIO will fail. So fall back onto just reading the protein,
##  force the fasta generation by getting each residue.
