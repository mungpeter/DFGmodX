#!/usr/bin/env python3

import re,sys
from x_variables import DFGTempls
from x_variables import NCTermTempls

##########################################################################
##
##  Peter M.U. Ung @ MSSM
## 
##  v1.0    14.10.01
##  v2.0    17.05.04    add beginning and ending cut-sites to only model
##                      the essential kinase catalytic domain
##  v3.0    17.06.13    retire 1XKK (too short) and use 4YBK as template
##  v4.0    18.03.30    loosen the cutsite requirement for NC-termini
##  v5.0    18.04.06    extend rendered DFG-motif length to 5 after DFG
##                      include cutsites for additional conformations
##
##  Purpose:    Take in .pir file (model sequence in the last) and parse
##              the model sequence suitable for DFGmodel DFG-out Kinase
##              remodeling.
##
##  CutSite templates for DFGmodel
##  3KQ7:MGA|DL/CEL|KI/YRA|PE   (c-in,  d-out, s/t) p38a
##  3NPC:MDA|NL/CTL|KI/YRA|PE   (c-in,  d-out, s/t) jnk2
##  3UGC:PYG|SL/NRV|KI/GES|PI   (c-in,  d-out, y)   jak2
##  3ETA:AHG|DL/FTV|KI/GLL|PV   (c-in,  d-out, y)   igfr
##  4YBK:SKG|SL/LVC|KV/WTA|PE   (c-out, d-in,  y)   src
##  5FD2:EGS|SL/LTV|KI/WMA|PE   (c-out, d-in,  s/t) braf
##  3EFJ:KHG|DL/FTV|KV/WMA|LE   (c-out, d-out, y)   cmet
##  3NAX:KNG|EL/MHI|QI/YVS|PE   (c-out, d-out, s/t) pdk1
##  
##  CutSites for essential kinase catalytic domain,
##  remove residues before N-cutsite and after C-cutsite
##  3KQ7:VPE|RY/FAQ|YH          (c-in,  d-out, s/t) p38a
##  3NPC:VLK|RY/ITV|WY          (c-in,  d-out, s/t) jnk2
##  3UGC:EER|HL/QIR|DN          (c-in,  d-out, y)   jak2
##  3ETA:SRE|KI/DDL|HP          (c-in,  d-out, y)   igfr
##  4YBK:PRE|SL/DYF|TS          (c-out, d-in,  y)   src     EDY/FT
##  5FD2:PDG|QI/LLA|RS          (c-out, d-in,  s/t) braf    ELL/AR
##  3EFJ:GPS|SL/AIF|ST          (c-out, d-out, y)   cmet
##  3NAX:RPE|DF/FES|VT          (c-out, d-out, s/t) pdk1
##
##  Retired:
##  1XKK (C-ter too short for cutsite alignment)
##    DFG:      1XKK:PFG|CL/QHV|KI/WMA|LE   (c-out, d-in, y) egfr1
##    NC-ter:   1XKK:KET|EF/KMA|RD          (c-out, d-in, y) egfr1
##
##########################################################################

##########################################################################
## Read in the aligned .fasta file to determine the positions of the 'cuts'
## and generate the Modeller .pir with corrected target sequence
def ModifyMultiPiecePIR(  PIRSequences, per_line, 
                          out_filename, mdl_output_pref, cidi_model=False ):

  # Parse reference PDB fasta for columns in aligned FASTA file and use
  # that column to find corresponding seq in target kinases
  DFGCuts = CutSiteTemplates( DFGTempls(),    'DFG'    )
  NCCuts  = CutSiteTemplates( NCTermTempls(), 'NC-ter' )

  # Take in list of Sequences for editing
  # Sequences format: Item[0] (>P1;) 
  #                   Item[1] (structureX|sequence)
  #                   Item[2] (aligned fasta) 
  #                   Item[3] (fasta length)
  #                   Item[4] (protein name)
  RefSeqs = {}
  Model   = []
  print( '\n\033[34m## Read-in PIR Sequences: \033[31m{0}\033[0m'.format(len(PIRSequences)))
  for index, Item in enumerate(PIRSequences):
    print(Item[0], '\t', str(Item[3]))

    # Extract model protein sequence for editing
    if re.search(r'sequence:', Item[1]): 
      Model = Item
      del PIRSequences[index]

    # Extract reference sequences for referencing based on protein name
    for key in DFGCuts.keys():
      if re.search('{0}'.format(key), Item[0]):
        RefSeqs[key] = Item[2]

#############

  if len(Model) == 0:
    sys.exit('\n  > #2# ERROR: No Model Sequence was found in .pir file: '+mdl_output_pref)

  # Found out the positions of DFG cut-sites in reference sequences
  print('\n\033[34m## DFGCuts ##\033[0m\n', (DFGCuts.keys()))
  print('\n\033[34m## RefSeqs ##\033[0m\n', (RefSeqs.keys()))

  DFGCutSites = FindCutSites(DFGCuts, RefSeqs, 'DFG', mdl_output_pref)
  print('\n\033[34m## Finding DFG Cut-Site position in reference column:\033[0m\n', DFGCutSites)

  if len(DFGCutSites) != 2: 
    sys.exit('\n\033[31m  > #2# ERROR: Require 2 DFG-CutSite references: only {0} found in\033[0m\n{1}: {2}'.format(
                  len(DFGCutSites), RefSeqs.keys(), mdl_output_pref )) 

  NCTCutSites = FindCutSites(NCCuts, RefSeqs, 'NC-ter', mdl_output_pref)
  if len(NCTCutSites) != 2:
    sys.exit('\n\033[31m  > #2# ERROR: Require 2 NC-CutSite references: only {0} found in\033[30m\n{1}: {2}'.format(
                  len(NCTCutSites), RefSeqs.keys(), mdl_output_pref ))

  # Consistency check: compare 2 reference positions
  # then reformat/parse the model sequence to DFGmodel specification
  CutSiteCheck(DFGCutSites, RefSeqs, 'DFG',    mdl_output_pref)
  CutSiteCheck(NCTCutSites, RefSeqs, 'NC-ter', mdl_output_pref)
  ModifiedModel = CutAndReprintModel(DFGCutSites[0], NCTCutSites[0], Model)


  # Reprint the .pir file with modified(cut) model sequence in it
  # if cidi_model is True, only print out last 2 PIRSequences; first 2 are only placeholders
  ReprintPIRFile(PIRSequences, ModifiedModel, cidi_model, per_line, out_filename)


##########################################################################
## Parse CutSite template file for referencing
## Must be formatted. Require 3 CutSites, each is 5 residue long, 3 + 2.
## CutSite patterns are separated by '/', while CutSite is separated by '|'
## NAME:XXX|xx/YYY|yy/ZZZ|zz

def CutSiteTemplates( CutSites, txt ):

  Templs = {}
  for line in CutSites:
    Temp = line.split(':')
    Templs[Temp[0]] = list(Temp[1].split('/'))
  print('\n## \033[34m{0} {1}\033[0m CutSite references are found:'.format(len(Templs), txt))

  for templ in Templs: 
    print(' \033[33m* '+templ+' * \033[0m')
    print(Templs[templ])

  return Templs


##########################################################################
## From the template file of cut sites, read in reference aligned sequences
## and find the locations of those cut sites. 5-residue sites are used for
## recognition, in the format of 3-2
## Template with cut site is formated
## >PDBName:XXX|YY/AAA|BB/III|JJ

def FindCutSites( Templs, RefSeqs, site, mdl_output_pref ):
  # For each reference sequence, identify the position of residues for cut 
  Sites = []
  for key in RefSeqs.keys():
    Positions = []

    # Retrieve the template sequences for cutting
    Cuts = Templs[key]
    for tmpl_seq in Cuts:
      front, back = tmpl_seq.split('|')
      TmplSeq = list(front+back)

      # Check the reference sequences to locate the cut position
      RefSeq  = RefSeqs[key]
      for index, resi in enumerate(RefSeq):
        tmpl_posi = 0

        # Recognize resid #1
        if resi == TmplSeq[tmpl_posi]:
          posit = index
          while RefSeq[posit+1] is '-': posit += 1
          tmpl_posi += 1
          # Recongize resid #2
          if RefSeq[posit+1] == TmplSeq[tmpl_posi]:
            posit += 1
            while RefSeq[posit+1] is '-': posit += 1
            tmpl_posi += 1
            # Recognize resid #3
            if RefSeq[posit+1] == TmplSeq[tmpl_posi]:
              posit += 1
              while RefSeq[posit+1] is '-': posit += 1
              tmpl_posi += 1
              # Recognize resid #4
              if RefSeq[posit+1] == TmplSeq[tmpl_posi]:
                posit += 1
                while RefSeq[posit+1] is '-': posit += 1
                tmpl_posi += 1
                cut_posit  = posit  # Record the position for cutting
                # Recognize resid #5
                if RefSeq[posit+1] == TmplSeq[tmpl_posi]:
                  posit += 1
                  # Save cut positions for each reference set
                  Positions.append(cut_posit)
                  break # Force reading only 1st of 4 tethered sequences
                else: continue
              else: continue
            else: continue
          else: continue
        else: continue
    Sites.append(Positions)
    if len(Positions) < 2: 
      sys.exit('\n  > #2# ERROR: {0} has only {1} of {2} {3} Cutsites -- {4}: {5}'.format(
                      key, len(Positions), len(Cuts), site, Positions, mdl_output_pref ))

  return Sites


##########################################################################
## 
def CutAndReprintModel( DFGCutSites, NCTCutSites, Model ):
  ModelInfo = [Model[0], Model[1]]
  ModelSeq  = Model[2]
#  print(CutSites)
  # Basic length of individual component of an aligned sequence
  seq_length = int(len(ModelSeq)/4)

  # Cut sites position in the full, aligned, tethered sequence
  # Since the full sequence is a tether of 4 identical sequences, with iiii 
  # jjjj as the N-/C-termini needed to be trimmed, edit the seq
  # with '-' in between the required fragments can be visualized as:
  #   0|iiii|______||s1xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|jjjj|m
  # m+1|iiiiixxxxxx||ms1___||ms2xxxxxxxxxxxxxxxxxxxxxxxxx|jjjj|n
  # n+1|iiiiixxxxxxxxxxxxx||ns2___||ns3xxxxxxxxxxxxxxxxxx|jjjj|o
  # o+1|iiiiixxxxxxxxxxxxxxxxxxxxx||os3__________________|jjjj|p

  i   = 0                               # N-term on 1st tethered element
  s0  = NCTCutSites[0]                  # N-term on 1st tethered element
  s1  = DFGCutSites[0]                  # on 1st tethered element
  ms1 = DFGCutSites[0] + seq_length     # on 2nd tethered element
  ms2 = DFGCutSites[1] + seq_length     # on 2nd tethered element
  ns2 = DFGCutSites[1] + (seq_length*2) # on 3rd tethered element
  ns3 = DFGCutSites[2] + (seq_length*2) # on 3rd tethered element
  os3 = DFGCutSites[2] + (seq_length*3) # on 4th tethered element
  j   = NCTCutSites[1] + (seq_length*3) # C-term on 4th tethered element
  p   = len(ModelSeq)                   # C-term on 4th tethered element


  Positions = []
  # Append the 1st fragment: N-lobe
  Positions.extend('-' for m in range(i, s0))
  Positions.extend(ModelSeq[m] for m in range(s0, s1))
  Positions.extend('-' for m in range(0, (ms1 - s1)))

  # Append the 2nd fragment: C-lobe #1
  Positions.extend(ModelSeq[m] for m in range(ms1, ms2))
  Positions.extend('-' for m in range(0, (ns2 - ms2)))

  # Append the 3rd fragment: DFG-loop -- extend to 5 residues after DFG
  Positions.extend(ModelSeq[m] for m in range(ns2, ns2+12))
  Positions.extend('/')
  Positions.extend('-' for m in range(0, (os3 - ns2-12)-1))

  # Append the 4th fragment: C-lobe #2
  Positions.extend(ModelSeq[m] for m in range(os3, j))
  Positions.extend('-' for m in range(0, (p - j)))


  # Format the tethered sequence into 4 separate sequences
  NewPosit = []
  NewPosit.append([Positions[m] for m in range(0,           seq_length  )])
  NewPosit.append([Positions[m] for m in range(seq_length,  seq_length*2)])
  NewPosit.append([Positions[m] for m in range(seq_length*2,seq_length*3)])
  NewPosit.append([Positions[m] for m in range(seq_length*3,seq_length*4)])
  ModelInfo.append(NewPosit)

  return ModelInfo


##########################################################################
## Format and print the parsed .pir for Modeller remodeling
## if cidi_model is True, only print last 2 PIRSequences; first 2 are only placeholders
def ReprintPIRFile( PIRSequences, Model, cidi_model, per_line, out_filename ):
  # Print out the reference sequences
  Prints = []

  # if dealing with CIDI modelling, use only last 2 seq - best-match and target; 
  # otherwise, use all
  if cidi_model is True:
    PIRSequences = PIRSequences[-1:]

  for Item in PIRSequences:
    FormatSeqForPrint(FormatRefSeq(Item), Prints, per_line)

  # Print out the model sequence
  FormatSeqForPrint(Model, Prints, per_line)

  with open(out_filename, 'w') as fo:
    for line in Prints: 
      fo.write(line)


##########################################################################
## Separate the single chain into 4 equally-sized chains for Ref. Sequences
def FormatRefSeq( Seq ):
  length = int(len(Seq[2])/4)
  Ref    = [Seq[m] for m in range(0,2)]
  RefSeq = []

  RefSeq.append([Seq[2][m] for m in range(0,        length  )])
  RefSeq.append([Seq[2][m] for m in range(length,   length*2)])
  RefSeq.append([Seq[2][m] for m in range(length*2, length*3)])
  RefSeq.append([Seq[2][m] for m in range(length*3, length*4)])
  Ref.append(RefSeq)
  return Ref


##########################################################################
## Put the modified sequences into .pir format, each line contains X residues
def FormatSeqForPrint( Fasta, Prints, per_line ):
  Prints.append(Fasta[0]+'\n')   # P1; header line
  Prints.append(Fasta[1]+'\n')   # Xtal:: info line

  Seq = Fasta[2]
  for Residues in Seq:
    Line = []
    for index, residue in enumerate(Residues):
      if index == 0:
        Line.append(residue)
      elif index == (len(Residues)-1):
        Line.append(residue)
        Prints.append(''.join(Line)+'\n')
        Line = []
      elif (index+1) % per_line == 0:
        Line.append(residue)
        Prints.append(''.join(Line)+'\n')
        Line = []
      else:
        Line.append(residue)
  Prints.append('*\n\n')

  return Prints


##########################################################################
# Consistency check: compare 2 reference positions
# then reformate/parse the model sequence to DFGmodel specification
# For C-terminus, it is less important and some chains are shorter. Ignore C-ter
def CutSiteCheck( CutSites, RefSeqs, txt, mdl_output_pref ):

  print('\n## Using \033[34m{0} {1}\033[0m CuteSite references:'.format(len(CutSites), txt))
  for key in RefSeqs:
    print('\n\033[35m- CutSite locations:\033[0m')
    print(key, '=', CutSites)

  if txt == 'DFG':
    if CutSites[0] != CutSites[1]:
      sys.exit('\n  > #2# ERROR: Parsing reference DFG-sequences {0} and {1}\n  > #2# ERROR: do not give same CutSite profile: {2} -- {3} vs {4}: {5}'.format(
                      RefSeqs.keys()[0], RefSeqs.keys()[1], txt,
                      CutSites[0], CutSites[1], mdl_output_pref ) ) 
  else:
    if CutSites[0][0] != CutSites[1][0]:
      sys.exit('\n  > #2# ERROR: Parsing reference NC-sequences {0} and {1}\n  >#2# ERROR: do not give same CutSite profile: {2} -- {3} vs {4}: {5}'.format(
                      RefSeqs.keys()[0], RefSeqs.keys()[1], txt,
                      CutSites[0], CutSites[1], mdl_output_pref ) )


##########################################################################
## Check if input files are in correct .pir format, then parse and output 
## the sequences for further modification 
def CheckPIR( seq_file, mdl_output_pref ):
  # Parse the sequence (fasta) file
  Sets = []
  Temp = []
  new_set = False
  with open(seq_file, 'r') as fi:
    for l in fi:   # skip blank lines
      line = l.rstrip()
      if re.search('##', line):
        continue
      if line is False:
        continue
      if new_set:
        Temp.append(line)
      if re.search(r'\*', line):
        Sets.append(Temp)
        new_set = False
        Temp = []
      if re.search(r'>', line):
        new_set = True
        Temp.append(line)

  # Check the format of each entry and save the sequence dataset
  # Sequences format: Item[0] (>P1;)
  #                   Item[1] (structureX|sequence)
  #                   Item[2] (aligned fasta)
  #                   Item[3] (fasta length)
  #                   Item[4] (protein name)
  Pir = []
  for index, Prot in enumerate(Sets):
    Seq = []
    if re.search(r'>P1;', Prot[0]): 
      Seq.append(Prot[0])     # Seq[0]
    else:
      sys.exit('\n  > #2# ERROR: PIR file Missing "Header >P1;" for entry no.{0}: {1}'.format(
                      index+1, mdl_output_pref ) )

    if re.search('(sequence:|structureX:)', Prot[1]):
      Seq.append(Prot[1])     # Seq[1]
    else:
      sys.exit('\n  > #2# ERROR: PIR file Missing "Info sequence:" for entry no.{0}: {1}'.format(
                      index+1, mdl_output_pref ) )

    fasta = ''  
    for a in range(2, len(Prot)):
      fasta += re.sub('\*|\s+', '', Prot[a])    # remove extra '*' TER card
    Seq.append(fasta)                               # Seq[2]
    Seq.append(len(fasta))                          # Seq[3]
    Seq.append(re.split(r';|\||\s+', Prot[0])[1])   # Seq[4]
    Pir.append(Seq)

  # Extract the length of modeled protein
  seq_num = 0
  for Seq in Pir:
    if re.search(r'sequence', Seq[1]):
      seq_num = Seq[3]

  for Seq in Pir:
    if Seq[3] == seq_num: 
      # Aligned sequence should be divisible by 4 since it is a tether of 4
      if seq_num % 4 == 0:
        continue
      else:
        sys.exit('\n  > #2# ERROR: The sequence has {0} elements and is not dividable by 4: {1}'.format(
                        seq_num, mdl_output_pref ))
    else:
      print('{0} {1}'.format(Seq[3], seq_num))
      sys.exit('\n  > #2# ERROR: Misalignment betweeen model and target sequence {0}: {1}'.format(
                      Seq[0], mdl_output_pref ))

  PIRSequences = Pir
  return PIRSequences    # Output the dataset of entry
