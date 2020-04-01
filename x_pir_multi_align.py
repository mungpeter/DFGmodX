#!/usr/bin/env python3

import sys,os
import re,glob
import subprocess
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.PDB.PDBParser   import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

from aa_residue import AA
#from CommonUtility import *
from x_variables import per_line
from x_pir_edit  import CheckPIR
from x_pir_edit  import ModifyPIR

##########################################################################
#
#	Peter M.U. @ MSSM
#
#   v1.0 - 14.03.08
#   v2.0 - 14.04.02 -   added function to take in the full-length sequence
#                       of the target protein and compare it to crystal-
#			generated sequence to find out the missing loops,
#       		and incorporate the info into model generation.
#   v3.0 - 14.04.09 --  added function to take in mTOR/PI3K atypical kinase
#   v4.0 - 14.04.28 --  alert for unnatural amino acids or heteroatoms
#   v5.0 - 14.07.29 --  allow manually supplied alignment file
#   v6.0 - 14.07.31 --  change main alignment program to TCoffee/Expresso
#   v7.0 - 14.10.02 --  added the .pir modification function to parse Model
#                       Removed the 'atypical' function, 1) not compatible
#                       with the auto-parsing of Model, 2) not accurate
#   v8.0 - 16.12.05 --  mark out sys.exit; modify Biopython SeqIO
#   v9.0 - 17.03.30 --  fix a bug when using automated template search, rewrite
#                       pir writeout and tether-PDB writeout
#   v10. - 17.05.03 --  fix a bug with modeller cannot take PDB without uniform
#                       chain ID; fix a bug with generating .pir with correct
#                       number of chimera-template sequences; remove duplicate
#                       fasta seqs to avoid t_coffee problem
#   v11. - 17.06.28 --  use pre-aligned kinase profile as scaffod for alignment
#   v12. - 17.07.14 --  fix bug not generating de-HETATM pdb file in rare cases
#   v13. - 18.03.09 --  fix bug with wrong sequence counting, and adopt CODI
#                       multiple runs
#   v14. - 18.03.29     use Muscle to do alignment
#   v15  - 18.06.09     check change between '_TEMP.x2.fasta' n '_TEMP.x3.fasta'
#                       if any insertion due to missing residues in xtal struct
#   v16  - 20.02.20     edit MUSCLE gap scaling to enforce profile gap matching
#
#	Purpose:
#	This script generating the chimera PDB and aligned FASTA file input
#	for DFG-in --> DFG-out kinase modelling by MODELLER.
#
#	This script reads in the list of PDB files that will be used as
#	DFG-out template and generate the FASTA from the PDB.
#	A problem with BioPython FASTA generation with proteins with missing
#	loops is that the resultent crystal-FASTA does not have indication of
#	the missing loop; other programs looking at this crystal-FASTA may
#	interprete the FASTA as a single continuous chain. To make sure the 
#	missing loops in the PDB structure does not affect the model generated 
#	by MODELLER, the crystal-FASTA is aligned to the full-seq FASTA of
#	the corresponding protein. The resultent crystal-FASTA will be 
#	corrected for the missing loops with '-'. ClustalO will recognite the
#	'-' designation and align the sequences accordingly.
#	
###	TCoffee is called to align all templates to the kinase being modeled.
###	TCoffee mode 'Expresso' is used. This function retreives 3D structures
###	similar to the given sequences and uses the structural information to
###	generate final alignments.
#
#   # TCoffee is decipated for use and now will use Muscle to perform
#     single-sequence profile alignment - faster and mostly accurate
#
#       The MODELLER-alignment file (.pir) is parsed by the function DFGModify,
#       by recognizing several consistent sites for cut/paste (S/G next to 
#       hinge region, before DFG motif, and before APE motif).
#	(Previously, the MODELLER-alignment file was partially done and 
#       required additional modification to indicate the region of small lobe 
#       and DFG loop being modelled)
#
#	This script caternates the *superposed-PDB* files to create a tethered 
#	PDB used by MODELLER to model the regions in the kinase:
#		1) Small lobe (based on the knonw DFG-out template PDB)
#		2) large lobe (based on the DFG-in of the modelled PDB)
#	        3) DFG-loop   (based on the known DFG-out template PDB)
#	        4) large lobe (based on the DFG-in of the modelled PDB)
#		5) ligand (extracted from template PDB and cat to end of PDB)
#
#	Required:
#	- Superposed PDB files
#	- Remove non-ligands HETATM (e.g. solvent, salts)
#	- Post-processing of MODELLER alignment file
#	- Full-length sequence of the target kinase
#
##########################################################################

msg = '''
    > {0}
      [FASTA Database] [PDB Directory]
      [List of Template PDB Name] 
      [Target PDB]
      [Target FASTA sequence: .fasta | or | Target PDBID in FASTA Database]
      [Tethered PDB List Output File Name]
      [Alignment Output prefix]
      [Output prefix for DFG-out models]\n
      e.g. <python> Y_kinase.140129.filter.clusto.checked.nogap.fasta 
              2_ykinase/pdb ykinase_template.txt 2_ykinase/case_1/2X2L.pdb 
              '2X2L_A' chimera_2x2l_template.list 2x2l_align DFGmod.2x2l
      '''.format(sys.argv[0])
#if len(sys.argv) < 9 or len(sys.argv) > 11: sys.exit(msg)

##########################################################################

def ModellerMultiAlignGen(  pdb_directory, work_directory,
                            struct_database, struct_nogap, kinome_database, kinome_nogap,
                            template_list, tget_pdb, mdl_prot_fasta, 
                            best_match_struc, pc_ident, align_switch, correct_fasta,
                            chimera_tmpl_list, mdl_pir_file, mdl_output_pref ):

  # Build database of full-seq FASTA
  Database, db_order = CacheSeqDatabase([struct_database,kinome_database])
  NoGapDB,  ng_order = CacheSeqDatabase([struct_nogap,kinome_nogap])

  # Read in the list of PDBs to be used as templates
  Tmpl_List = pd.read_csv(template_list, comment='#', sep='\s+', header=None).iloc[:,0].to_numpy()
  print('\n\033[34m## The following PDBs will be used as templates:\033[0m\n'+template_list)
  print(Tmpl_List)
  print('\n')


  # Generate '_TEMP.{x}.y1.fasta' as intermediate, will check for percent
  # identity. This is also the insert point for corrected fasta if alignment
  # did not work out and need manual intervention '*y1.corr.fasta'. Since
  # it is corrected and need to maintain the exact format, no alignment will
  # be done to the corrected fasta
  # Convert the intermediate/corrected fasta to '_TEMP.{x}.y2.fasta'
  if re.search('None', correct_fasta, re.IGNORECASE):
    print('\033[34m## Running with Normal \033[31m"_TEMP.{0}.y1.fasta"\033[0m'.format(mdl_output_pref))
    AlignSequences( Database, NoGapDB, kinome_database, pdb_directory, Tmpl_List, 
                    tget_pdb, mdl_prot_fasta, best_match_struc,
                    pc_ident, align_switch, mdl_pir_file, mdl_output_pref )
  else:
    print('\033[31m## Running with CORRECTED Fasta file:\n\033[31m{0}\033[0m'.format(correct_fasta))
    CleanFASTAName( correct_fasta, work_directory, mdl_pir_file )

  BuildTemplatePDB(pdb_directory, Tmpl_List, tget_pdb, chimera_tmpl_list)

  # Generate the Modeller .pir file from alignment .fasta file
  GenerateModellerAlignmentFile(mdl_pir_file, mdl_output_pref)

  print('\n## Check the Target FASTA to make sure starting and ending\n##  residues, and missing loops are accounted for.\n##  Check the presence of phospho- or unnatural amino acid residue.\n')
  for fasta in SeqIO.parse('_TEMP.corrected.fasta', 'fasta'):
    print(fasta.format('fasta'))



###########################################################################
## Convert multi seq alignment into Modeller .pir alignment file
def GenerateModellerAlignmentFile( mdl_pir_file, mdl_output_pref ):

  # Read in the multiple sequence alignment file generated by TCoffee
  # with the last PDB as the Model kinase
  ali_prefix = mdl_pir_file.split('.pir')[0]
  print('\n\033[34m## Converting alignment .fasta to .pir:\033[0m\n{0}'.format(ali_prefix))


  ## Read all templates + matched base + target Fasta and format them the same way 
  Templ = []
  for item in list(SeqIO.parse(ali_prefix+'.fasta', 'fasta')):
    print('\033[31min generatemodelleralignment \033[0m',item.id)
    temp = item.id.split()[0]    # Expresso fasta format contains extra element
    flnm = '>P1;chimera_'+temp
    pdb_id = flnm.split(';')[1].rstrip()
    # Extract RES_ID and CHAIN information from first line of template PDBs
    if os.path.isfile('{0}.pdb'.format(pdb_id)):
      line = os.popen('head -1 "{0}.pdb" '.format(pdb_id)).readline()
      resi = line[22:26]
      chan = line[21:22]
    else:   # Target PDB may not exist but that's fine
      resi = '1   '
      chan = 'A'
    struct = '\nstructureX:{0}:{1}:{2}:LAST:{2}:{0}::-1.00:-1.00\n'.format(
                    pdb_id, resi, chan)
    Templ.append([flnm, struct, str(item.seq)])

  # Pop out Target Fasta ([-1]) and Base (best-fit [-2] from imported Fasta file
  Target = Templ.pop(-1)
  Base   = Templ.pop(-1)

  print(Base)
  ## Write out partially-prepared Modeller file (*.pir.prep)from alignment file
    # First write out the chimera-template Fasta. Alter the header lines for
    # the Target Fasta and append to last of the file
  with open('{0}.prep'.format(mdl_pir_file), 'w') as w:
    # Write out the chimera-template kinases first,
    for idx, Ali in enumerate(Templ):
      w.write('### {0}: {1} ###\n'.format(idx, Ali[0]))
      # Write the tethered template Fasta, 4 times
      for ln in Ali:          # Write template small lobe seq with header
        w.write(ln)
      w.write('\n')
      for ln in Base[2:]:	  # Write 1st large lobe part without header
        w.write(ln+'\n')
      for ln in Ali[2:]:      # Write DFG-loop without header
        w.write(ln+'\n')
      for ln in Base[2:]:	  # Write 2nd large lobe part without header
        w.write(ln+'\n')
      w.write("*\n\n")		    # Write the Modeller sequence Ending

    # Write the Target seq in tethered form, reformat the header lines
    w.write('### {0}: {1} ###\n'.format('final', Target[0]))
    for line in Target:
      if re.search(r'>P1;', line):
        old_name = line.rstrip().split(';')[1]
      line = re.sub(old_name, mdl_output_pref, line)
      line = re.sub(r'structureX', 'sequence', line)
      line = re.sub(r'-1\.00:-1\.00', ':', line)
      w.write(line)
    w.write('\n')
    for x in range(3):  # write the sequence 3 more times without header
      for ln in Target[2:]:
        w.write(ln+'\n')
    w.write("*")

  # Write out the partially-prepared Modeller .pir.prep file
  ModifyPIR(  CheckPIR( mdl_pir_file+'.prep', mdl_output_pref ), 
              per_line(), mdl_pir_file, 
              mdl_output_pref )


##########################################################################
## Call T-Coffee to perform multiple sequence alignment
## For each crystal-FASTA, do clustalo to correct for the missing loops
def AlignSequences( Database, NoGapDB, kinome_database, pdb_directory, Tmpl_List, 
                    tget_pdb, mdl_prot_fasta, best_match_struc,
                    pc_ident, align_switch, mdl_pir_file, mdl_output_pref ):

  ali_prefix = mdl_pir_file.split('.pir')[0]

#######################

  ## Place Template proteins in the beginning of the .pir file; all with gaps
  Seq  = []
  for tmpl_name in Tmpl_List:
    tmpl_id = tmpl_name.split('.')[0]
    print('  \033[34m** Template Protein:\033[0m '+tmpl_id)

    # Get the (already aligned) fasta sequence of kinase templates from Database
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
def BuildTemplatePDB( pdb_directory, Tmpl_List, tget_pdb, chimera_tmpl_list ):

  print('\n\033[34m## Building Chimera PDB for homology modeling ##\033[0m')

  with open(chimera_tmpl_list, 'w') as w:
    for pdb_name in Tmpl_List:
      print('  \033[34m# Building chimera PDB based on:\033[0m '+pdb_name)
      pdb_id = pdb_name.split('.')[0]

      pdb_file = '{0}/{1}'.format(pdb_directory, pdb_name)
      w.write('chimera_{0}\n'.format(pdb_id))

      print(pdb_file)
      print(tget_pdb)
      os.system('cat "{0}" "{1}" "{0}" "{1}" | grep -v "HETATM" | grep -v "END" | grep -v "CONECT" | sed -E "s/^(.{2})./\{3}/" > {4}_{5}.pdb'.format(pdb_file, tget_pdb, '{21}', '1A','chimera', pdb_id))


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