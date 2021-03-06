#!/usr/bin/env python3

import sys
msg = '''
    > {0}
        [ Setup Script: formatted setup file | Pickle ]\n
        [ Mode: -full   Perform modified alignment, DFGmodel, and selection
                -pir    Run modified alignment only
                -mod    Run DFGmodel only
                -vol    Select top models only
                -set    Generate template setup script\n
                -none   * use new paths in "x_variables"; Workstation -> HPC
      **# Always use -pir then -mod, and ALWAYS check .pir for additional chain
          breaks. Only use -full if you have checked the kinase and it has no
          additional chain break or missing residues\n
        [ Opt:   -force   Forced restart of /dfgworking directory ]
        [ Opt:   -restart Forced restart of /1_result   directory ]
        [ Opt:   -pass    Pass along; always use this unless otherwise ]
        [ Opt:   -paths   * Update directory paths in setup files; use with "-none" ]\n
      e.g. > {0} setup.file -pir -pass\n
'''.format(sys.argv[0])
if len(sys.argv) < 2 or len(sys.argv) > 4: sys.exit(msg)

import os,re
import shutil
import pickle

from x_variables import DefaultVariables
from x_pir_gen import GenerateAndModifyPIR
from x_pir_multi_align import CacheSeqDatabase
from x_modl_gen import GenerateAndAlignModeller
from x_vol_gen import GeneratePOVMEAndSortModels
from x_variables_run import ParseInputVariables
from x_variables_run import GenerateTemplSetupScript
from x_homolog_templ_check import SearchKinaseStruct

##########################################################################
def main( setting_file, running, **kwargs ):

  # Settings to run which function of DFGmodel
  pir, mod, vol = False, False, False
  if running == '-full':
    pir, mod, vol = True, True, True
  elif running == '-pir':
    pir = True
  elif running == '-mod':
    mod = True
    vol = True
  elif running == '-vol':
    vol = True
  elif running == '-set':
    GenerateTemplSetupScript(setting_file)
    sys.exit('\n  Template setup script: {0}\n'.format(setting_file))
  elif running == '-none':
    print('')
  else:
    sys.exit('\n  Error: [Running Option] is not specified: {0}'.format(running))

  # Retreive Settings variables
  name = setting_file.split('/')[-1].split('.')[0]
  if re.search(r'.pkl', setting_file):
    Settings = pickle.load( open(setting_file, 'rb') )
  else:
    Settings = ParseInputVariables(setting_file)

  ## Use only when transfer from workstation to HPC; update directory paths
  if running == '-none':
    New = DefaultVariables()
    Settings['ScriptDirectory']  = New['ScriptDirectory']
    Settings['DatasetDirectory'] = New['DatasetDirectory']
    Settings['PymolExecutable']  = New['PymolExecutable']
    Settings['HomeDirectory']    = New['HomeDirectory']
    Settings['ResultDirectory']  = New['ResultDirectory']
    Settings['WorkingDirectory'] = New['WorkingDirectory']
    Settings['POVMEDirectory']   = New['POVMEDirectory']
    Settings['StructDatabase']   = New['StructDatabase']
    Settings['StructNoGap']      = New['StructNoGap']
    Settings['KinomeDatabase']   = New['KinomeDatabase']
    Settings['KinomeNoGap']      = New['KinomeNoGap']
    Settings['ConfClassify']     = New['ConfClassify']
    Settings['PDBDirectory']     = New['PDBDirectory']

    # This sets the use of one best matched xtal structure as C-lobe base
    best_struct = Settings['BestMatchStruc'].split('/')[-1]
    New['BestMatchStruc'] = Settings['PDBDirectory']+best_struct
    Settings['BestMatchStruc']   = New['BestMatchStruc']

    # This sets the use of one of multple template_list
    if type(Settings['TemplateList']) is not list:
      templ_list = Settings['TemplateList'].split('/')[-1]
      New['TemplateList'] = Settings['DatasetDirectory']+templ_list
      Settings['TemplateList'] = New['TemplateList']
    else:
      New['TemplateList'] = [Settings['DatasetDirectory']+i.split('/')[-1]    
                              for i in Settings['TemplateList']]
      Settings['TemplateList'] = New['TemplateList']  


  script_directory = Settings['ScriptDirectory']  # Directory with all the scripts
  dataset_dir      = Settings['DatasetDirectory'] # Directory with all the datasets
  home_directory   = Settings['HomeDirectory']    # Current directory
  work_directory   = Settings['WorkingDirectory'] # Working directory
  result_directory = Settings['ResultDirectory']  # Result directory

  struct_database  = Settings['StructDatabase']   # all-structure alignment fasta database
  struct_nogap     = Settings['StructNoGap']      # all-structure no-gap fasta database
  kinome_database  = Settings['KinomeDatabase']   # aligned canonical human kinome seq database
  kinome_nogap     = Settings['KinomeNoGap']      # no-gap canonical human kinome seq database
  conf_database    = Settings['ConfClassify']     # database of xtal conformation

  pdb_directory    = Settings['PDBDirectory']     # PDB structure directory
  template_list    = Settings['TemplateList']     # List of C-helix/DFG templates
  superpose_resi   = Settings['SuperposeRefResi'] # reference residue superpose
  seq_ident_thres  = Settings['SeqIdentThres']    # threshold to flag seq alignment
  align_switch     = Settings['AlignSwitchThres'] # switch from MUSCLE to EXPRESSO
  conformation     = Settings['CHelixDFGModel']   # conformation to generate

  povme_directory  = Settings['POVMEDirectory']   # Directory to temp folder
  pymol_exec       = Settings['PymolExecutable']  # pymol executable
  reference_pdb    = Settings['ReferencePDB']    # Reference PDB for superpose

  # Create .pir alignment file for running
  prot_struc_input = Settings['KinaseStructInput']# Name of model kinase
  mdl_prot_fasta   = Settings['ModelKinaseFasta'] # fasta of model kinase
  chimera_tmpl_list= Settings['ChimeraTemplList'] # Name of tether chimera file
  best_match_struc = Settings['BestMatchStruc']   # Best matching PDB from Blast
  correct_fasta    = Settings['CorrectFASTAFile'] # need to correct _TEMP.y.fasta
  pc_ident         = Settings['SeqIdentity']      # best match sequence identity

  # Running Modeller and superpose models
  mdl_pir_file     = Settings['ModifiedPIRFile']  # PIR file for Modeller
  number_of_model  = Settings['NumberOfModel']    # Number of model to generate
  number_of_cpu    = Settings['NumberOfCPU']      # Number of CPU to use
  mdl_output_pref  = Settings['OutputPrefix']     # Prefix of generated models

  # Ranking models based on POVME volume calculation
  top_model        = Settings['NumberOfTopModel'] # Number of selected model
  povme_exec       = Settings['POVMEExecutable']  # Directory of POVME software
  povme_pdb        = Settings['POVMEStructure']   # Name of multi-struct protein

#################################################

  # Check presence of critical files before running
#  CheckEssentFiles(script_directory, dataset_dir, Settings)

  # Check presence of pre-existing result folder
  if   option == '-restart':
    if os.path.isdir(result_directory):
      print( '\n  ** Found [ResultDirectory] **\n    {0}\n     Forced removal of old directory'.format(result_directory))
      shutil.rmtree(result_directory)   # Remove the existing directory
    if os.path.exists(work_directory):
      print( '\n  ** Found [WorkingDirectory] **\n   {0}\n     Forced removal of old directory'.format(work_directory))
      shutil.rmtree(work_directory)   # Remove the existing directory
  elif option == '-force':
    if os.path.isdir(result_directory):
      print( '\n  ** Found [ResultDirectory] **\n    {0}\n     Ignore and proceed'.format(result_directory))
    if os.path.isdir(work_directory):
      print( '\n  ** Found [WorkingDirectory] **\n   {0}\n     Forced removal of old directory'.format(work_directory))
      shutil.rmtree(work_directory)   # Remove the existing directory
  elif option == '-pass':
    if os.path.isdir(result_directory):
      print( '\n  ** Found [ResultDirectory] **\n    {0}\n     Ignore and proceed'.format(result_directory))
    if os.path.exists(work_directory):
      print( '\n  ** Found [WorkingDirectory] **\n    {0}\n     Ignore and proceed'.format(work_directory))
  elif option == '-paths':
      print( '\n\033[33m  ** Just redo the setup file with directory paths; no change to existing files **\033[0m')
      ## update the setup files
      pickle.dump( Settings, open( name+'.setup.pkl', 'wb' ) )
      with open( name+'.setup', 'w') as fo:
        for key in Settings.keys():
          fo.write('{0:20s} {1}\n'.format(key, Settings[key]))
      sys.exit()
  else:
      sys.exit('\n  \033[31mERROR: Invalid Option ( -restart | -force | -pass | -paths ):\033[0m '+option)


  if not os.path.isdir(work_directory):
    os.makedirs(work_directory)
  if not os.path.isdir(result_directory):
    os.makedirs(result_directory)


##########################################################################
  
  # Generate modified alignment file (.pir) for Modeller
  # with structure search
  if pir:

    # If "KinaseStructInput" is "None" (not supplied), trigger the search for
    # *Kinase Structure* with highest sequence identity within the internal struct
    # library, and use it in the subsequent runs
    if re.search(r'None', prot_struc_input, re.IGNORECASE) and re.search(r'None', mdl_prot_fasta, re.IGNORECASE):
      sys.exit('\n  > #2# ERROR: Both "KinaseStructInput" and "ModelKinaseFasta" are "None: '+mdl_output_pref)
    elif re.search('None', correct_fasta, re.IGNORECASE):
      # return matched PDB file, and its percent identity
      best_match_struc, pc_ident = SearchKinaseStruct(
          script_directory, home_directory, work_directory, result_directory,
          pdb_directory, struct_nogap, kinome_nogap, conf_database, conformation,
          seq_ident_thres, reference_pdb, mdl_prot_fasta, Settings )
      os.chdir(home_directory)  # come back to home directory
      Settings['BestMatchStruc'] = best_match_struc
      Settings['SeqIdentity']    = pc_ident

      ## update the setup files
      pickle.dump( Settings, open( '{0}.setup.pkl'.format(name), 'wb' ) )
      with open('{0}.setup'.format(name), 'w') as fo:
        for key in Settings.keys():
          fo.write('{0:20s} {1}\n'.format(key, Settings[key]))


    # if template_list is > 1 (CODI), run the subconformations in iterations
    if type(template_list) is list:     # CODI use [list of template_lists]
      for idx, template in enumerate(template_list):
        GenerateAndModifyPIR(
            script_directory, dataset_dir, home_directory,
            work_directory, result_directory,pdb_directory,
            struct_database, struct_nogap, kinome_database, kinome_nogap,
            template, reference_pdb, superpose_resi, best_match_struc,
            prot_struc_input, mdl_prot_fasta,
            pc_ident, align_switch, correct_fasta,
            chimera_tmpl_list+'.'+str(idx+1), mdl_pir_file+'.'+str(idx+1),
            mdl_output_pref+'.'+str(idx+1), pymol_exec, Settings )
    else:                               # CIDI/CIDO/CODO use 1 template_list
      GenerateAndModifyPIR(
          script_directory, dataset_dir, home_directory,
          work_directory, result_directory, pdb_directory,
          struct_database, struct_nogap, kinome_database, kinome_nogap,
          template_list, reference_pdb, superpose_resi, best_match_struc,
          prot_struc_input, mdl_prot_fasta,
          pc_ident, align_switch, correct_fasta,
          chimera_tmpl_list, mdl_pir_file,
          mdl_output_pref, pymol_exec, Settings )


##########################################################################
  # Build Models based on .pir and superposed structure files
  if mod:
    print('*** '+Settings['BestMatchStruc']+' ***')
    print(best_match_struc)
    if type(template_list) is list:      # CODI use [list of template_lists]
      for idx in range( len(template_list) ):
        GenerateAndAlignModeller(
            script_directory, home_directory, work_directory, result_directory,
            chimera_tmpl_list+'.'+str(idx+1), mdl_pir_file+'.'+str(idx+1),
            number_of_model, number_of_cpu,
            mdl_output_pref+'.'+str(idx+1),
            reference_pdb, best_match_struc, superpose_resi,
            pymol_exec, Settings )
    else:                                # CIDI/CIDO/CODO use 1 template_list
      GenerateAndAlignModeller(
          script_directory, home_directory, work_directory, result_directory,
          chimera_tmpl_list, mdl_pir_file,
          number_of_model, number_of_cpu,
          mdl_output_pref,
          reference_pdb, best_match_struc, superpose_resi,
          pymol_exec, Settings )


##########################################################################
  # Select models with largest binding pocket
  if vol:                         # CODI use [list of template_lists]
    if type(template_list) is not str:
      for idx in range( len(template_list) ):
        GeneratePOVMEAndSortModels(
            script_directory, home_directory, work_directory, result_directory,
            povme_exec, povme_directory, povme_pdb,
            mdl_output_pref+'.'+str(idx+1),
            number_of_cpu, conformation,
            int(top_model), Settings )
    else:                         # CIDI/CIDO/CODO use 1 template_list
      GeneratePOVMEAndSortModels(
          script_directory, home_directory, work_directory, result_directory,
          povme_exec, povme_directory, povme_pdb,
          mdl_output_pref,
          number_of_cpu, conformation,
          int(top_model), Settings )


##########################################################################
if __name__ == "__main__":
  if len(sys.argv) == 2:
    option = None
  if len(sys.argv) == 3:
    option = None
  if len(sys.argv) == 4:
    option = sys.argv[3]
  main(sys.argv[1], sys.argv[2], f=option)

#########################################################################
#
#   Peter M.U. Ung @ MSSM
#
#   v1.0    14.10.17
#   v2.0    15.07.22    bug fixes
#   v3.0    16.06.07    added C-out modeling option -depreciated
#   v4.0    16.08.10    superpose residues come from global variables
#   v5.0    16.11.25    check homolog template for kinase with only sequence
#   v6.0    16.11.27    search for closest structure if no supplied
#   v7.0    16.12.05    integrate CheckHomologTemplate to SearchKinaseStruct
#   v8.0    18.03.06    allow modeling for multiple types of conformations,
#                       and multiple subconformations for CODI and wCD
#   v9.0    18.03.28    change the 'setting_file' to accept pre-defined Dict
#   v10     20.03.21    adopt alignment to kinome MSA profile, detect if input
#                       fasta is already aligned and then skip MUSCLE alignment 
#                       add CIDI modeling using Konformation-classifed xtals
#
