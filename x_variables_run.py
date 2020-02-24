#!/usr/bin/env python3

from CommonUtility import *
from x_variables import DefaultVariables
from x_variables import TemplatePDBList as TL
import re, os
import subprocess

##########################################################################
## Generate a template setup script with all keywords and default settings
def GenerateTemplSetupScript( setting_file, V=DefaultVariables() ):

  os.chdir(V['HomeDirectory'])
  f = open(setting_file, 'w')

  f.write('\n\n### Step 1 - Generate modified PIR\t# -pir flag ###\n')
  f.write('\n# Generate C-helix/DFG model:\n')
  f.write('# C-in/DFG-in: cidi; C-in/DFG-out: cido; C-out/DFG-in: codi; C-out/DFG-out: codo; Random: wcd\n')
  f.write('{0:15s}\t{1}\n\n'.format('CHelixDFGModel', V['CHelixDFGModel']))
  f.write('# Kinase family: (serine|st|tyrosine|y|not-specified|na)\n\n')
  f.write('{0:15s}\t{1}\n'.format('KinaseFamily', V['KinaseFamily']))

  f.write('# Kinase PDB to be modeled. If structure of the target kinase is not available \n# use a closely related homolog kinase (%ident > '+V['SeqIdentThres']+'%) as input structure\n')
  f.write('# If left as "None", kinase with the closest match will be used\n')
  f.write('{0:15s}\t{1}\n'.format('KinaseStructInput', V['KinaseStructInput']))
  f.write('# FASTA file of model kinase (just PDB code if already in Fasta Database)\n')
  f.write('{0:15s}\t{1}\n'.format('ModelKinaseFasta', V['ModelKinaseFasta']))
  f.write('# Output prefix for all intermediate files\n')
  f.write('{0:15s}\t{1}\n'.format('OutputPrefix', 'XXX'))

  f.write('\n\n### Step 2 - Generate DFGmodel\t# -mod flag ###\n')
  f.write('# PyMOL executable for superposition\n')
  f.write('{0:15s}\t{1}\t'.format('PyMOL', V['PymolExecutable']))
  f.write('# Number of initial models to generate\n')
  f.write('{0:15s}\t{1}\n'.format('NumberOfModel', V['NumberOfModel']))

  f.write('\n\n# Step 3 - Select top modelsl\t# -vol flag\n')
  f.write('# Number of models with largest volume saved\n')
  f.write('{0:15s}\t{1}\n'.format('NumberOfTopModel', V['NumberOfTopModel']))
  f.write('{0:15s}\t{1}\n'.format('NumberOfCPU', V['NumberOfCPU']))
  f.close()


##########################################################################
## Extract the variables in the setting file
def ParseInputVariables( setting_file ):

  Variables = DefaultVariables()
  Lines = remove_remark(file_handle(setting_file))

  for line in Lines:
    if re.search(r'ScriptDirectory',  line, re.IGNORECASE):
      Variables['ScriptDirectory']    = line.split('\t')[1]
    if re.search(r'DatasetDirectory', line, re.IGNORECASE):
      Variables['DatasetDirectory']   = line.split('\t')[1]
    if re.search(r'HomeDirectory',    line, re.IGNORECASE):
      Variables['HomeDirectory']      = line.split('\t')[1]
    if re.search(r'WorkingDirectory', line, re.IGNORECASE):
      Variables['WorkingDirectory']   = line.split('\t')[1]
    if re.search(r'ResultDirectory',  line, re.IGNORECASE):
      Variables['ResultDirectory']    = line.split('\t')[1]

    if re.search(r'FastaDatabase',    line, re.IGNORECASE):
      Variables['FastaDatabase']      = line.split('\t')[1]
    if re.search(r'KinaseProfile',    line, re.IGNORECASE):
      Variables['KinaseProfile']      = line.split('\t')[1]

    if re.search(r'RefKinasePath',    line, re.IGNORECASE):
      Variables['RefKinasePath']      = line.split('\t')[1]
    if re.search(r'SuperposeRefResi', line, re.IGNORECASE):
      Variables['SuperposeRefResi']   = line.split('\t')[1]

    if re.search(r'KinaseStructInput',line, re.IGNORECASE):
      Variables['KinaseStructInput']  = line.split('\t')[1]
    if re.search(r'ModelKinaseFasta', line, re.IGNORECASE):
      Variables['ModelKinaseFasta']   = line.split('\t')[1]
    if re.search(r'ChimeraTemplList', line, re.IGNORECASE):
      Variables['ChimeraTemplList']   = line.split('\t')[1]
    if re.search(r'CorrectFASTAFile', line, re.IGNORECASE):
      Variables['CorrectFASTAFile']   = line.split('\t')[1]
    if re.search(r'SeqIdentity',      line, re.IGNORECASE):
      Variables['SeqIdentity']        = line.split('\t')[1]

    if re.search(r'ModifiedPIRFile',  line, re.IGNORECASE):
      Variables['ModifiedPIRFile']    = line.split('\t')[1]
    if re.search(r'CHelixDFGModel',   line, re.IGNORECASE):
      Variables['CHelixDFGModel']     = line.split('\t')[1]
    if re.search(r'KinaseFamily',     line, re.IGNORECASE):
      Variables['KinaseFamily']       = line.split('\t')[1]
    if re.search(r'UseHomologTempl',  line, re.IGNORECASE):
      Variables['UseHomologTempl']    = line.split('\t')[1]

    if re.search(r'PymolExecutable',  line, re.IGNORECASE):
      Variables['PymolExecutable']    = line.split('\t')[1]
    if re.search(r'NumberOfModel',    line, re.IGNORECASE):
      Variables['NumberOfModel']      = int(line.split()[1])
    if re.search(r'NumberOfTopModel', line, re.IGNORECASE):
      print(line.split())
      Variables['NumberOfTopModel']   = int(line.split()[1])
    if re.search(r'NumberOfCPU',      line, re.IGNORECASE):
      Variables['NumberOfCPU']        = int(line.split()[1])
    if re.search(r'OutputPrefix',     line, re.IGNORECASE):
      Variables['OutputPrefix']       = line.split('\t')[1]

  Variables = VariableSetup( Variables, '', '',
                             Variables['CHelixDFGModel'],
                             Variables['NumberOfCPU'],
                             Variables['NumberOfModel'],
                             Variables['OutputPrefix'],
                             Variables['CheckFASTAFile'],
                             None )

  # Print Variables on Standard IO
  for key in Variables:
    print('  * {0:20s}\t{1}'.format(key ,Variables[key]))

  return Variables


#####################################################################
## generate specific DFGmod input setup file, with modifications
def VariableSetup( Vars, family, name, conf, cpu_num, num_mod, 
                        num_top, corr, prefix ):

  if re.search(r'TK:', family):
    Vars['KinaseFamily'] = 'tyrosine'

  ## Generate prefix name for different outputs/options
  if prefix is None:
    prefix = conf+'.'+name
    Vars['OutputPrefix']      = prefix
    Vars['ModelKinaseFasta']  = str(os.getcwd())+'/'+name+'.fasta'
    Vars['KinaseStructInput'] = 'None'
    Vars['CHelixDFGModel']    = conf

  Vars['CorrectFASTAFile']  = corr
  Vars['NumberOfModel']     = num_mod
  Vars['NumberOfTopModel']  = num_top
  Vars['ChimeraTemplList']  = prefix+'.chimera.list'
  Vars['ModifiedPIRFile']   = prefix+'.chimera.pir'
  Vars['AlignedModelList']  = prefix+'.aligned_pdb.list'
  Vars['PymolAlignPrefix']  = prefix+'.align'
  Vars['POVMEStructure']    = prefix+'.multi_frame.pdb'

  dataset_dir = Vars['DatasetDirectory']
  ## For each conformation type, select the corresponding template sets for modeling
  # If typical DFG-in/out (C-in), select either S/T or Y template set
  if re.search(r'cidi', conf, re.IGNORECASE):
    if re.search(r'(st)|(serine)|(ser/ter)|(s/t)', Vars['KinaseFamily'], re.IGNORECASE):
      Vars['TemplateList'] = '{0}/{1}'.format(dataset_dir, TL('cidi_s_templ'))
    if re.search(r'(y)|(tyr)|(tyrosine)',  Vars['KinaseFamily'], re.IGNORECASE):
      Vars['TemplateList'] = '{0}/{1}'.format(dataset_dir, TL('cidi_y_templ'))

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
