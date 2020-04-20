#!/usr/bin/env python3

import sys,os
import re,glob
import subprocess
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

from aa_residue import AA

##########################################################################
# ## Use already-aligned sequences from database to put together a FASTA file
## For each crystal-FASTA, do clustalo to correct for the missing loops
def SequenceAlign(  Database, NoGapDB, kinome_database, pdb_directory, 
                    Tmpl_List, tget_pdb, mdl_prot_fasta, best_match_struc,
                    pc_ident, align_switch, mdl_pir_file, mdl_output_pref ):

  ali_prefix = mdl_pir_file.split('.pir')[0]

#######################

  ## Place (already aligned) fasta sequence of multiple templates from Database,
  ## all with gaps into Seq; if CIDI modeling, skip this step coz only 1 template
  Seq  = []
  if Tmpl_List is None:
    # for cut-site purpose, will include a hard-coded reference 1ATP for placeholder
    Seq.append( '>{0}\n{1}\n'.format('1ATP_E',Database['1ATP_E'].seq) )
  else:
    for tmpl_name in Tmpl_List:
      tmpl_id = tmpl_name.split('.')[0]
      print('  \033[34m** Template Protein:\033[0m '+tmpl_id)
      Seq.append( '>{0}\n{1}\n'.format(tmpl_id, Database[tmpl_id].seq) )

#####################
  ## Establish a gapped alignment of "Base" kinase for later structure building.
  ## Append Base (best-fit) structure for model building. If 'tget_pdb' is 
  ## manually supplied, use it. If found through searching, use 'best_pdb_id'

  pdb_id      = tget_pdb.split('/')[-1].split('.')[0]    # input PDB_ID from filename
  best_pdb_id = best_match_struc.split('/')[-1].split('.')[0] # Best match ID from filename
  print('\n  \033[34m** Target kinase:\033[0m       '+pdb_id)
  print('  \033[34m** Best matched kinase:\033[0m '+best_pdb_id)

  # Canonical kinase sequence from Database
  with open('_TEMP.best-gap.fasta', 'w') as db:
    SeqIO.write(Database[best_pdb_id], db, 'fasta')

  # Check if kinase structure exist in database; if not, do alignment to best-match seq
  if pdb_id in Database:
    base_seq = Database[pdb_id]
  else:
    with open('_TEMP.tget-pdb.fasta', 'w') as tg:
      tg.write('>{0}\n{1}'.format(pdb_id, str(FASTA_Gen(tget_pdb, pdb_id))))

    base_seq = MuscleProfileAlign('_TEMP.best-gap.fasta', '_TEMP.tget-pdb.fasta',
                                  '_TEMP.x2.fasta')
    if len(base_seq.seq) != len(Database[best_pdb_id].seq):
      sys.exit('  \033[34mERROR: MUSCLE alignment of input PDB failed with wrong length (_TEMP.x2.fasta): \033[0m'+tget_pdb)

  print(' \033[31m> base_seq:\033[0m '+base_seq.id)

  # Append the Matched Base kinase sequence to the end of list of template kinase
  Seq.append('\n>{0}\n{1}\n'.format(base_seq.id, base_seq.seq))

########################

  # print out Target kinase sequence and full-length of closest kinase
  MissingLoopCorrection(pdb_id, NoGapDB[best_pdb_id].seq, 
                        FASTA_Gen(tget_pdb, pdb_id))

  # Check if phospho- or unnatural amino acid is in the PDB
  CheckUnnaturalAA(tget_pdb, pdb_id)

#######################
  ## First, check if incoming PDB 'tget_pdb' is a known structure by inquiring
  ## Database. If exists, it should be the best-matching template. Do single-
  ## seq MUSCLE-profile alignment to map incoming sequence to this best template
  ## and generate templates-included alignment FASTA file

  # check if input fasta a file and if a known seq in the database
  if not os.path.isfile(mdl_prot_fasta):
    # if 'mdl_prot_fasta' is not Fasta file and is an ID in Database
    if mdl_prot_fasta in Database:
      print('\n \033[34m** FASTA input is known ID -- use FASTA in database:\033[0m '+mdl_prot_fasta)
      tget_seq = Database[mdl_prot_fasta]
    else:
    # If not input FASTA file, get fasta from input structure, then profile align
      print('\n \033[34m** FASTA input is "None" -- use Input Kinase Structure for FASTA:\033[0m '+tget_pdb)
      tget_seq = base_seq
  else:
    print('\n  \033[34m** FASTA input is a file:\033[0m \n'+mdl_prot_fasta+'\n')
    mdl_id = mdl_prot_fasta.split('/')[-1].split('.')[0]
    tget_seq = list(SeqIO.parse(mdl_prot_fasta, 'fasta'))[0]
    tget_seq.id = mdl_id
    print('\033[34m> input fasta seq length:\033[0m '+str(len(tget_seq.seq)))
    if len(tget_seq.seq) == len(base_seq.seq):
      print('  \033[31m> Input FASTA file appears to be pre-aligned to MD-kinome seq with equal length:\033[0m '+mdl_prot_fasta)
    else:
      # Do single-seq MUSCLE-profile alignment
      print('  \033[31m> Input FASTA file length differs from length in Database. Do profile alignment for:\033[0m '+mdl_prot_fasta)
      tget_seq = MuscleProfileAlign('_TEMP.best-gap.fasta', mdl_prot_fasta,
                                    '_TEMP.x3.fasta')

#####################

  # Append the target kinase sequence to a list of template sequence for later use
  Seq.append('\n>{0}\n{1}\n'.format(mdl_output_pref, tget_seq.seq))
  with open('_TEMP.{0}.y1.fasta'.format(mdl_output_pref), 'w') as fo:
    for fasta in Seq:  fo.write(fasta)

  # Check if the alignment of target seq is different from template
  # if checks out okay, convert the temp alignment file to full fasta file for use
  if len(base_seq.seq) != len(tget_seq.seq):
    x = '  base_seq.seq: {0}\n'.format(len(base_seq.seq))
    x = x+('  tget_pdb.seq: {0}\n'.format(len(tget_seq.seq)))
    sys.exit(x+'  \033[31m> #4# FATAL: Alignment of Target seq to Template seq has different length:\033[0m '+mdl_output_pref)
  else:
    print('    \033[34m- #1# Okay: len(base_seq.seq) == len(tget_seq):\033[0m '+str(len(tget_seq.seq)))
    # Generates a full FASTA file for use later
    TemplYCheck(pdb_id, ali_prefix, '_TEMP.{0}.y1.fasta'.format(mdl_output_pref), 
                pc_ident, align_switch, kinome_database)


##########################################################################
##########################################################################
# Read in the multiple sequence alignment file generated by Muscle/Tcoffee
# with the last PDB as the Model kinase
def ParseFastaSeqForPIR( mdl_pir_file ):
  ali_prefix = mdl_pir_file.split('.pir')[0]
  print('\n\033[34m## Converting \033[31mSingle-Piece\033[34m template alignment .fasta to .pir:\033[0m\n{0}'.format(ali_prefix))

  ## Read all templates + matched base + target Fasta and format them the same way 
  Templ = []
  for item in list(SeqIO.parse(ali_prefix+'.fasta', 'fasta')):
    temp = item.id.split()[0]    # Expresso fasta format contains extra element
    flnm = '>P1;chimera_'+temp
    pdb_id = flnm.split(';')[1].rstrip()

    # Extract RES_ID and CHAIN information from first line of template PDBs
    if os.path.isfile('{0}/{1}.pdb'.format(os.getcwd(),pdb_id)):
      print(' \033[33m# Extract 1st position and chain info from PDB: \033[32m{0}\033[0m'.format(pdb_id))
      line = os.popen('head -1 "{0}.pdb" '.format(pdb_id)).readline()
      resi = line[22:26]
      chan = line[21:22]
    else:  
      print(' \033[31m# Target PDB {0} may not exist but that is fine #\033[0m'.format(pdb_id))
      resi = '1   '
      chan = 'A'
    struct = '\nstructureX:{0}:{1}:{2}:LAST:{2}:{0}::-1.00:-1.00\n'.format(
                    pdb_id, resi, chan)
    Templ.append([flnm, struct, str(item.seq)])
  
  return Templ


##########################################################################
## Convert _TEMP.y.fasta into the final form before converting into .pir
def TemplYCheck( pdb_id, ali_prefix, fasta_to_pir, pc_ident, align_switch, kinome_database ):

  # If target seq ident is lower than MUSCLE-to-EXPRESSO switching threshold
  if pc_ident < align_switch:
    print('  \033[31m> #1# WARNING:\033[0m {0} Input FASTA has low Identity to best match < 50%, use T_Coffee: {1}'.format(
              ali_prefix, pdb_id))
    RunTCoffeeExpresso(fasta_to_pir, ali_prefix, kinome_database)
  else:
    RemoveFastaGapColumn(fasta_to_pir, ali_prefix+'.fasta')


##########################################################################
## When _TEMP.{x}.y1.fasta needs to be cleaned up by JalView, clean the sequence
## name, where JalView adds '/Start-End' tag to the name, which breaks the code
def CleanFASTAName( fasta_file, work_directory, mdl_pir_file ):

  name       = fasta_file.split('.y1.')[0]
  ali_prefix = mdl_pir_file.split('.pir')[0] 
  Data = []

  for fas in SeqIO.parse(fasta_file, 'fasta'):
    Data.append( SeqRecord( id=fas.id.split('/')[0], seq=fas.seq,
                      description=fas.description, name=fas.name ) )

  with open('{0}.y2.fasta'.format(name), 'w') as fo:
    SeqIO.write(Data, fo, 'fasta')

  RemoveFastaGapColumn( '{0}.y2.fasta'.format(name), ali_prefix+'.fasta' )


##########################################################################
## Align PDB-FASTA to full-seq FASTA to detect/correct missing loops in PDB
def MissingLoopCorrection( pdb_id, full_seq, xtal_seq ):

  missing = len(full_seq) - len(xtal_seq)

  print( '** There are \033[31m{0:3d}\033[0m residues missing in {1} **'.format(missing, pdb_id))
  with open('_TEMP.fasta', 'w') as w:
    w.write('>{0}|full-seq\n{1}\n\n>{0}\n{2}'.format(pdb_id, full_seq, xtal_seq))

  RunClustalO('_TEMP.fasta', '_TEMP.corrected')
  seq_record = list(SeqIO.parse('_TEMP.corrected.fasta', 'fasta'))
  for seq in seq_record:
    print(seq.format('fasta'))

  return seq_record[1].seq


##########################################################################
## Alignment using a pre-existing MSA profile. T-Coffee has issue with this. 
## MUSCLE and ClustalO came out about same time, 2004 and 2003, MUSCLE performs
## the best especially when doing single-seq profile alignment. Output the 
## profile-aligned fasta object
## Penalities for gap opening and extension are modified to enforce the new 
## alignment adopt the same gapping as the profile sequence. See examples:
## http://www.drive5.com/muscle/muscle_userguide3.8.html
## https://www.dnastar.com/manuals/MegAlignPro/15.3/en/topic/muscle-alignment-options
def MuscleProfileAlign( fasta_database, fasta_file, temp_file ):

  os.system('muscle -profile -in1 "{0}" -in2 "{1}" -out "{2}" -maxiters 500 -seqtype protein -gapopen -5.0 -gapextend -2.0 -center 0.0 -quiet'.format(
              fasta_database, fasta_file, temp_file ) )

  Tget_List = list(SeqIO.parse(temp_file, 'fasta'))

  return Tget_List[1]


#########################################################################
## Basic setup for running T-coffee Expresso / new use with 3d T-coffee
def RunTCoffeeExpresso( fasta_file, ali_prefix, kinome_database ):

  fasta_name = fasta_file.split('.fasta')[0]
#  os.system('t_coffee -in "{0}" -output=fasta,clustalw,html -mode expresso -max_n_proc 10 -method=mafft_msa,t_coffee_msa,dialigntx_msa,muscle_msa,kalign_msa -email {1}'.format(fasta_file, 'pmung@umich.edu'))
  os.system('t_coffee -in "{0}" -output=fasta,clustalw,html -template_file {1} -max_n_proc 10 -method=mafft_msa,t_coffee_msa,dialigntx_msa,muscle_msa,kalign_msa -email {2}'.format(fasta_file, kinome_database, 'pmung@umich.edu'))
  os.system('mv "{0}_aln" "{1}.fasta"'.format(fasta_file, ali_prefix))
  os.system('mv "{0}.html" "{1}.html"'.format(fasta_name, ali_prefix))
  os.system('mv "{0}.clustalw" "{1}.ali"'.format(fasta_name, ali_prefix))


##########################################################################
# Reformat and no alignment: Remove empty gap column from aligned FASTA file
# the '-action +rm_gap <% empty>' tag indicate which columns to remove, if
# column contains <% empty> seq that are empty in that column
# if not include '-action +rm_gap <>' flag, all '-' will be removed
def RemoveFastaGapColumn( fasta_input, fasta_output ):

  os.system('t_coffee -other_pg seq_reformat -in "{0}" -action +rm_gap 100 -output=fasta > {1}'.format(fasta_input, fasta_output))


##########################################################################
## Align with Clustal Omega
def RunClustalO( fasta_file, ali_prefix ):
  os.system('clustalo -i "{0}" -o "{1}.fasta" --full --force'.format(
              fasta_file, ali_prefix))


##########################################################################
## Align crystal-FASTA to full-seq FASTA to detect missing loops in PDB
def CacheSeqDatabase( fasta_database ):
  if type(fasta_database) is not list:
    db_list = [fasta_database]
  else:
    db_list = fasta_database
  print('\n  \033[34m## Caching sequence database:\033[0m\n',db_list)
  Database = {}
  Order    = []
  for db_name in db_list:
    for seq_record in SeqIO.parse(db_name, 'fasta'):
      pdb_id = seq_record.id.split('/')[0].split('|')[0].replace(':','_')
      seq_record.id = pdb_id
      seq_record.description = ''
      Database[pdb_id] = seq_record
      Order.append(pdb_id)
  return Database, Order


#########################################################################
## Use BioPython to generate FASTA sequence from PDB structure
def FASTA_Gen( pdb_name, pdb_id ):

  #print(pdb_name)
  try:
    m = PDBParser(PERMISSIVE=1).get_structure(pdb_id, pdb_name)
    peptides = PPBuilder().build_peptides(m)

    seq = ''
    for p in peptides:
      seq = seq + p.get_sequence()
  except FileNotFoundError:
    print('\033[31m ERROR: PDB File not found: \033[0m'+pdb_name)
    seq = None

  return seq


##########################################################################
## Build tethered template PDB file for Modeller alignment/model generation
def BuildMultiPieceTemplatePDB( pdb_directory, Tmpl_List, tget_pdb, chimera_tmpl_list ):

  print('\n\033[34m## Building Chimera PDB for homology modeling ##\033[0m')

  with open(chimera_tmpl_list, 'w') as w:
    for pdb_name in Tmpl_List:
      print('  \033[34m# Building chimera PDB based on:\033[0m '+pdb_name)
      pdb_id = pdb_name.split('.')[0]

      pdb_file = '{0}/{1}'.format(pdb_directory, pdb_name)
      w.write('chimera_{0}\n'.format(pdb_id))

      print(pdb_file)
      print(tget_pdb)
      os.system('cat "{0}" "{1}" "{0}" "{1}" | grep -v "HETATM" | grep -v "END" | grep -v "CONECT" | sed -E "s/^(.{2})./\{3}/" > {4}_{5}.pdb ; wait'.format(
                  pdb_file, tget_pdb, '{21}', '1A','chimera', pdb_id))


##########################################################################
## Check if phospho- or unnatural amino acid is in the PDB
def CheckUnnaturalAA( pdb_name, pdb_id ):
  m = PDBParser(PERMISSIVE=1).get_structure(pdb_id, pdb_name)
  UAA = []
  for residue in m.get_residues():
    if AA(residue.get_resname()) == '':
      UAA.append(residue.get_resname())
  print('\033[34m>> Phospho- or unnatural amino acid, and ligand in Target PDB\033[31m {0}\033[0m:'.format(pdb_id))
  print('{0}\n'.format(UAA))


##########################################################################
