#!/usr/bin/env python3

import re, os
import subprocess

from CommonUtility import *
from x_variables import DefaultVariables
from x_variables import TemplatePDBList as TL

##########################################################################
## Generate a template setup script with all keywords and default settings
def GenerateTemplSetupScript( setting_file, V=DefaultVariables() ):

  os.chdir(V['HomeDirectory'])
  f = open(setting_file, 'w')

  f.write('## Output prefix for all intermediate files\n')
  f.write('{0:20s} {1}\n'.format('OutputPrefix', 'XXX'))
  f.write('{0:20s} {1}\n'.format('NumberOfCPU', V['NumberOfCPU']))

  f.write('\n\n### Step 1 - Generate modified PIR\t# -pir flag ###\n')
  f.write('# Generate C-helix/DFG model:\n')
  f.write('# C-in/DFG-in: cidi; C-in/DFG-out: cido; C-out/DFG-in: codi; C-out/DFG-out: codo; Random: wcd\n')
  f.write('{0:20s} {1}\n'.format('CHelixDFGModel', V['CHelixDFGModel']))
  f.write('# Kinase family: (serine|st|tyrosine|y|not-specified|na)\n')
  f.write('{0:20s} {1}\n\n'.format('KinaseFamily', V['KinaseFamily']))

  f.write('# Kinase PDB to be modeled. If structure of the target kinase is not available \n# use a closely related homolog kinase (%ident > {0}%) as input structure\n'.format(V['SeqIdentThres']))
  f.write('# If left as "None", kinase with the closest match will be used\n')
  f.write('{0:20s} {1}\n'.format('KinaseStructInput', V['KinaseStructInput']))
  f.write('# FASTA file of model kinase, aligned to MD-kinome MSA (or just ID if already in Fasta Databases)\n')
  f.write('{0:20s} {1}\n\n'.format('ModelKinaseFasta', V['ModelKinaseFasta']))
  f.write('# Corrected FASTA alignment file, "None" unless "_TEMP.<OutputPrefix>.y1.fasta" has issues\n')
  f.write('{0:20s} {1}\n'.format('CorrectFASTAFile', V['CorrectFASTAFile']))

  f.write('\n\n### Step 2 - Generate DFGmodel\t# -mod flag ###\n')
  f.write('# Number of initial models to generate\n')
  f.write('{0:20s} {1}\n'.format('NumberOfModel', V['NumberOfModel']))

  f.write('\n\n# Step 3 - Select top modelsl\t# -vol flag\n')
  f.write('# Number of models with largest volume saved\n')
  f.write('{0:20s} {1}\n'.format('NumberOfTopModel', V['NumberOfTopModel']))


  f.write('\n\n### Basic Settings -- Most are default and do not need to change ###\n')
  f.write('# Directory of DFGmodx scripts\n')
  f.write('{0:20s} {1}\n'.format('ScriptDirectory', V['ScriptDirectory']))
  f.write('# Directory of internal data used by DFGmodx scripts\n')
  f.write('{0:20s} {1}\n'.format('DatasetDirectory', V['DatasetDirectory']))
  f.write('# Directory of internal superposed RCSB PDB structure databases\n')
  f.write('{0:20s} {1}\n'.format('PDBDirectory', V['PDBDirectory']))
  f.write('# PyMOL executable\n')
  f.write('{0:20s} {1}\n'.format('PymolExecutable', V['PymolExecutable']))
  f.write('# POVME2.1 executible\n')
  f.write('{0:20s} {1}\n\n'.format('POVMEExecutable', V['POVMEExecutable']))
  f.write('# Home Directory for work\n')
  f.write('{0:20s} {1}\n'.format('HomeDirectory', V['HomeDirectory']))
  f.write('# Result Directory for work\n')
  f.write('{0:20s} {1}\n'.format('ResultDirectory', V['ResultDirectory']))
  f.write('# Temporary Working Directory for work\n')
  f.write('{0:20s} {1}\n'.format('WorkingDirectory', V['WorkingDirectory']))
  f.write('# Temporary POVME Directory for work\n')
  f.write('{0:20s} {1}\n\n'.format('POVMEDirectory', V['POVMEDirectory']))

  f.write('# Modi/Dunbrack-aligned human kinome (Typical kinases, 495-seq) FASTA file\n')
  f.write('{0:20s} {1}\n'.format('KinomeDatabase', V['KinomeDatabase']))
  f.write('# Modi/Dunbrack-aligned human kinome (Typical kinases, 495-seq) FASTA file - No Gap\n')
  f.write('{0:20s} {1}\n'.format('KinomeNoGap', V['KinomeNoGap']))
  f.write('# Modi/Dunbrack-aligned RCSB Typical S/T-kinase structure FASTA file\n')
  f.write('{0:20s} {1}\n'.format('KinomeDatabase', V['KinomeDatabase']))
  f.write('# Modi/Dunbrack-aligned RCSB Typical S/T-kinase structure FASTA file - No Gap\n')
  f.write('{0:20s} {1}\n'.format('StructNoGap', V['StructNoGap']))
  f.write('# Reference PDB structure for reference superposition\n')
  f.write('{0:20s} {1}\n'.format('ReferencePDB', V['ReferencePDB']))
  f.write('# C-lobe catalytic domain beta-sheet on reference PDB for superposition\n')
  f.write('{0:20s} {1}'.format('SuperposeRefResi', V['SuperposeRefResi']))

  f.close()


##########################################################################
## Extract the variables in the setting file
def ParseInputVariables( setting_file ):

  Variables = DefaultVariables()
  Lines = remove_remark(file_handle(setting_file))

  for line in Lines:
    if re.search(r'ScriptDirectory',  line, re.IGNORECASE):
      Variables['ScriptDirectory']    = line.split()[1]
    if re.search(r'DatasetDirectory', line, re.IGNORECASE):
      Variables['DatasetDirectory']   = line.split()[1]
    if re.search(r'HomeDirectory',    line, re.IGNORECASE):
      Variables['HomeDirectory']      = line.split()[1]
    if re.search(r'WorkingDirectory', line, re.IGNORECASE):
      Variables['WorkingDirectory']   = line.split()[1]
    if re.search(r'ResultDirectory',  line, re.IGNORECASE):
      Variables['ResultDirectory']    = line.split()[1]

    if re.search(r'StructDatabase',   line, re.IGNORECASE):
      Variables['StructDatabase']     = line.split()[1]
    if re.search(r'StructNoGap',      line, re.IGNORECASE):
      Variables['StructNoGap']        = line.split()[1]
    if re.search(r'KinomeDatabase',   line, re.IGNORECASE):
      Variables['KinomeDatabase']     = line.split()[1]
    if re.search(r'KinomeNoGap',      line, re.IGNORECASE):
      Variables['KinomeNoGap']        = line.split()[1]

    if re.search(r'RefKinasePath',    line, re.IGNORECASE):
      Variables['RefKinasePath']      = line.split()[1]
    if re.search(r'SuperposeRefResi', line, re.IGNORECASE):
      Variables['SuperposeRefResi']   = line.split('RefResi')[1]  ## cannot use ' ' to split

    if re.search(r'KinaseStructInput',line, re.IGNORECASE):
      Variables['KinaseStructInput']  = line.split()[1]
    if re.search(r'ModelKinaseFasta', line, re.IGNORECASE):
      Variables['ModelKinaseFasta']   = line.split()[1]
    if re.search(r'ChimeraTemplList', line, re.IGNORECASE):
      Variables['ChimeraTemplList']   = line.split()[1]
    if re.search(r'CorrectFASTAFile', line, re.IGNORECASE):
      Variables['CorrectFASTAFile']   = line.split()[1]
    if re.search(r'SeqIdentity',      line, re.IGNORECASE):
      Variables['SeqIdentity']        = line.split()[1]
    if re.search(r'BestMatchStruc',   line, re.IGNORECASE):
      Variables['BestMatchStruc']     = line.split()[1]

    if re.search(r'ModifiedPIRFile',  line, re.IGNORECASE):
      Variables['ModifiedPIRFile']    = line.split()[1]
    if re.search(r'CHelixDFGModel',   line, re.IGNORECASE):
      Variables['CHelixDFGModel']     = line.split()[1]
    if re.search(r'KinaseFamily',     line, re.IGNORECASE):
      Variables['KinaseFamily']       = line.split()[1]

    if re.search(r'PymolExecutable',  line, re.IGNORECASE):
      Variables['PymolExecutable']    = line.split()[1]
    if re.search(r'POVMEExecutable',  line, re.IGNORECASE):
      Variables['POVMEExecutable']    = line.split()[1]

    if re.search(r'NumberOfModel',    line, re.IGNORECASE):
      Variables['NumberOfModel']      = int(line.split()[1])
    if re.search(r'NumberOfTopModel', line, re.IGNORECASE):
      Variables['NumberOfTopModel']   = int(line.split()[1])
    if re.search(r'NumberOfCPU',      line, re.IGNORECASE):
      Variables['NumberOfCPU']        = int(line.split()[1])
    if re.search(r'OutputPrefix',     line, re.IGNORECASE):
      Variables['OutputPrefix']       = line.split()[1]

  Variables = VariableSetup(Variables, '', '',
                            Variables['CHelixDFGModel'],
                            Variables['NumberOfCPU'],
                            Variables['NumberOfModel'],
                            Variables['NumberOfTopModel'],
                            Variables['OutputPrefix'] )

  # Print Variables on Standard IO for checking
  for key in Variables:
    print('* \033[34m{0:20s} \033[31m{1}\033[0m'.format(key ,Variables[key]))

  return Variables


#####################################################################
## generate specific DFGmod input setup file, with modifications
def VariableSetup( Vars, family, name, conf, cpu_num, 
                    num_mod, num_top, prefix ):

  if re.search(r'TK:', family):
    Vars['KinaseFamily'] = 'tyrosine'

  ## Generate prefix name with different options, use when doing for whole kinome
  if prefix is None:
    prefix = conf+'.'+name
    Vars['OutputPrefix']      = prefix
#    Vars['ModelKinaseFasta']  = str(os.getcwd())+'/'+name+'.fasta'
    Vars['ModelKinaseFasta']  = name
    Vars['KinaseStructInput'] = 'None'
    Vars['CHelixDFGModel']    = conf

  Vars['NumberOfModel']     = num_mod
  Vars['NumberOfTopModel']  = num_top
  Vars['NumberOfCPU']       = cpu_num
  Vars['ChimeraTemplList']  = prefix+'.chimera.list'
  Vars['ModifiedPIRFile']   = prefix+'.chimera.pir'
  Vars['AlignedModelList']  = prefix+'.aligned_pdb.list'
  Vars['PymolAlignPrefix']  = prefix+'.align'
  Vars['POVMEStructure']    = prefix+'.multi_frame.pdb'


  dataset_dir = Vars['DatasetDirectory']
  ## For each conformation type, select the corresponding template sets for modeling

  # If typical DFG-in/out (C-in), select either S/T or Y template set
  if re.search(r'cido', conf, re.IGNORECASE):
    if re.search(r'(st)|(serine)|(ser/ter)|(s/t)', Vars['KinaseFamily'], re.IGNORECASE):
      Vars['TemplateList'] = '{0}/{1}'.format(dataset_dir, TL('cido_s_templ'))
    if re.search(r'(y)|(tyr)|(tyrosine)',  Vars['KinaseFamily'], re.IGNORECASE):
      Vars['TemplateList'] = '{0}/{1}'.format(dataset_dir, TL('cido_y_templ'))

  # If C-out or wCD, S/T or Y kinase set does not matter since not enough samples anyways
  # but there are C-Out specific conformations
  if re.search(r'codi', conf, re.IGNORECASE):
    Vars['TemplateList'] = ['{0}/{1}'.format(dataset_dir,codi) for codi in TL('codi_templ').split(',')]

  if re.search(r'III', conf, re.IGNORECASE):
    Vars['TemplateList'] = '{0}/{1}'.format(dataset_dir, TL('III_templ'))
  if re.search(r'MET', conf, re.IGNORECASE):
    Vars['TemplateList'] = '{0}/{1}'.format(dataset_dir, TL('met_templ'))
  if re.search(r'EGFR', conf, re.IGNORECASE):
    Vars['TemplateList'] = '{0}/{1}'.format(dataset_dir, TL('egfr_templ'))

  if re.search(r'cidi', conf, re.IGNORECASE):
    Vars['TemplateList'] = '{0}/{1}'.format(dataset_dir, TL('cidi_templ'))

  if re.search(r'codo', conf, re.IGNORECASE):
    Vars['TemplateList'] = '{0}/{1}'.format(dataset_dir, TL('codo_templ'))

  if re.search(r'wcd', conf, re.IGNORECASE):
    Vars['TemplateList'] = ['{0}/{1}'.format(dataset_dir, wcd) for wcd in TL('wcd_templ').split(',')]


  return Vars


#####################################################################
#
#   v1.0    15.11.02
#   v2.0    16.06.07    add C-out parameters
#   v3.0    16.11.27    add auto structure search parameters
#   v4.0    16.12.05    remove homolog check
#   v5.0    17.02.03    use new sequence alignment 'resi 124-138+160-183'
#   v6.0    17.03.30    internalize PDBDirectory by type of stkinase/ykinase
#   v7.0    17.06.28    add pre-aligned seq profile
#   v8.0    18.03.06    allow modeling for multiple types of conformations,
#                       and multiple subconformations for CODI and wCD
#   v9.0    18.03.27    simplify printout and readin options
#   v10.0   18.03.28    new function to handle variable setup
#   v11.0   18.04.06    include fasta cutsites for kinase conformations
#
#   v1.0    18.04.08    separate from "x_variables.py"
#   v2.0    20.03.12    add kinome database
