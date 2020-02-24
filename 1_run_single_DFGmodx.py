#!/usr/bin/env python3

import sys,os,re
import shutil
import pickle

from x_pir_gen  import GenerateAndModifyPIR
from x_check_scripts import CheckEssentFiles
from x_modl_gen import GenerateAndAlignModeller
from x_vol_gen  import GeneratePOVMEAndSortModels
from x_variables_run import ParseInputVariables
from x_variables_run import GenerateTemplSetupScript
from x_homolog_templ_check import SearchKinaseStruct

msg = '''
    > {0}
        [ Setup Script: formatted setup file | Pickle ]\n
        [ Mode: -full   Perform modified alignment, DFGmodel, and selection
                -pir    Run modified alignment only
                -mod    Run DFGmodel only
                -vol    Select top models only
                -set    Generate template setup script\n
      **# Always use -pir then -mod, and ALWAYS check .pir for additional chain
          breaks. Only use -full if you have checked the kinase and it has no
          additional chain break or missing residues\n
        [ Opt:   -force   Forced restart of /dfgworking directory ]
        [ Opt:   -restart Forced restart of /1_result   directory ]
        [ Opt:   -pass    Pass along; always use this unless otherwise ]\n
      e.g. > {0} setup.file -pir -pass\n
'''.format(sys.argv[0])
if len(sys.argv) < 2 or len(sys.argv) > 4: sys.exit(msg)

##########################################################################
def main(setting_file, running, **kwargs):

  # Settings to run which function of DFGmodel
  pir, mod, vol = False, False, False
  if re.search(r'-full', running):
    pir, mod, vol = True, True, True
  elif re.search(r'-pir', running):
    pir = True
  elif re.search(r'-mod', running):
    mod = True
    vol = True
  elif re.search(r'-vol', running):
    vol = True
  elif re.search(r'-set', running):
    GenerateTemplSetupScript(setting_file)
    sys.exit('\n  Template setup script: {0}\n'.format(setting_file))
  else:
    sys.exit('\n  Error: [Running Option] is not specified: {0}'.format(running))

  # Variables of Settings
  if re.search(r'.pkl', setting_file):
    Settings = pickle.load( open(setting_file, 'rb') )
  else:
    Settings = ParseInputVariables(setting_file)

  script_directory = Settings['ScriptDirectory']# Directory with all the scripts
  dataset_dir      = Settings['DatasetDirectory']# Directory with all the datasets
  home_directory   = Settings['HomeDirectory']  # Current directory
  work_directory   = Settings['WorkingDirectory'] # Working directory
  result_directory = Settings['ResultDirectory']  # Result directory

  fasta_database   = Settings['FastaDatabase']  # all-structure alignment database
  kinase_profile   = Settings['KinaseProfile']  # kinase profile 99% redundancy
  pdb_directory    = Settings['PDBDirectory']   # PDB structure directory
  template_list    = Settings['TemplateList']   # List of C-helix/DFG templates
  superpose_resi   = Settings['SuperposeRefResi'] # reference residue superpose
  seq_ident_thres  = Settings['SeqIdentThres']  # threshold to flag seq alignment
  align_switch     = Settings['AlignSwitchThres'] # switch from MUSCLE to EXPRESSO
  conformation     = Settings['CHelixDFGModel'] # conformation to generate

  # Create .pir alignment file
  prot_struc_input = Settings['KinaseStructInput']# Name of model kinase
  mdl_prot_fasta   = Settings['ModelKinaseFasta'] # fasta of model kinase
  chimera_tmpl_list= Settings['ChimeraTemplList'] # Name of tether chimera file
  best_match_struc = Settings['BestMatchStruc']   # Best matching PDB from Blast
  correct_fasta    = Settings['CorrectFASTAFile'] # need to correct _TEMP.y.fasta
  pc_ident         = Settings['SeqIdentity']      # best match sequence identity

  # Running Modeller and superpose models
  mdl_pir_file     = Settings['ModifiedPIRFile'] # PIR file for Modeller
  number_of_model  = Settings['NumberOfModel']  # Number of model to generate
  number_of_cpu    = Settings['NumberOfCPU']    # Number of CPU to use
  mdl_output_pref  = Settings['OutputPrefix'] # Prefix of generated models
  pymol_exec       = Settings['PymolExecutable'] # pymol executable

  # Ranking models based on POVME volume calculation
  povme_location   = Settings['POVMELocation']  # Directory of POVME software
  top_model        = Settings['NumberOfTopModel'] # Number of selected model
  povme_pdb        = Settings['POVMEStructure'] # Name of multi-struct protein
  povme_directory  = Settings['POVMEDirectory'] # Directory to temp folder
  povme_template   = Settings['POVMETemplateFile'] # Template POVME setup file

  template_pdb     = Settings['TemplatePDB']  # Template PDB for superpose

#################################################

  # Check presence of critical files before running
  CheckEssentFiles(script_directory, dataset_dir, Settings)

  # Check presence of pre-existing result folder
  if   option == '-restart':
    if os.path.isdir(result_directory):
      print( '\n  ** Found [ResultDirectory] **\n    {0}\n     Forced removal of old directory **'.format(result_directory))
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
  else:
      sys.exit('\n  Warning: Invalid Option: -restart | -force | -pass')


  if not os.path.isdir(work_directory):
    os.makedirs(work_directory)
  if not os.path.isdir(result_directory):
    os.makedirs(result_directory)


##########################################################################
  # Generate modified alignment file (.pir) for Modeller
  # with structure search
  if pir:
    # If "KinaseStructInput" is "None" (not supplied), trigger the search for
    # kinase structure with highest sequence identity within the internal struct
    # library, and use it in the subsequent runs
    if re.search(r'None', prot_struc_input, re.IGNORECASE) and re.search(r'None', mdl_prot_fasta, re.IGNORECASE):
      sys.exit('\n  > #2# ERROR: Both "KinaseStructInput" and "ModelKinaseFasta" are "None: '+mdl_output_pref)
    elif not correct_fasta:  
      best_match_struc, pc_ident = SearchKinaseStruct(
          script_directory, home_directory, work_directory, result_directory,
          pdb_directory, kinase_profile, seq_ident_thres, template_pdb,
          mdl_prot_fasta, Settings )
      Settings['BestMatchStruc'] = best_match_struc
      Settings['SeqIdentity']    = pc_ident
      pickle.dump( Settings, open( "../SetupVars.pkl", "wb" ) )


    # if template_list is > 1 (CODI), run the subconformations in iterations
    if type(template_list) is not str:
      for idx, template in enumerate(template_list):
        GenerateAndModifyPIR(
            script_directory, dataset_dir, home_directory, 
            work_directory, result_directory,
            fasta_database, kinase_profile, pdb_directory, template,
            template_pdb, superpose_resi, best_match_struc,
            prot_struc_input, mdl_prot_fasta, 
            pc_ident, align_switch, correct_fasta,
            chimera_tmpl_list+'.'+str(idx+1), mdl_pir_file+'.'+str(idx+1), 
            mdl_output_pref+'.'+str(idx+1), pymol_exec, Settings )
    else:
      GenerateAndModifyPIR(
          script_directory, dataset_dir, home_directory, 
          work_directory, result_directory,
          fasta_database, kinase_profile, pdb_directory, template_list,
          template_pdb, superpose_resi, best_match_struc,
          prot_struc_input, mdl_prot_fasta, 
          pc_ident, align_switch, correct_fasta,
          chimera_tmpl_list, mdl_pir_file, 
          mdl_output_pref, pymol_exec, Settings )


##########################################################################
  # Build Models based on .pir and superposed structure files
  if mod:
    print('*** '+Settings['BestMatchStruc']+' ***')
    print(best_match_struc)
    if type(template_list) is not str:
      for idx in range( len(template_list) ):
        GenerateAndAlignModeller(
            script_directory, home_directory, work_directory, result_directory,
            chimera_tmpl_list+'.'+str(idx+1), mdl_pir_file+'.'+str(idx+1), 
            number_of_model,
            number_of_cpu, mdl_output_pref+'.'+str(idx+1),
            template_pdb, best_match_struc, superpose_resi,
            pymol_exec, Settings )
    else:
      GenerateAndAlignModeller(
          script_directory, home_directory, work_directory, result_directory,
          chimera_tmpl_list, mdl_pir_file, number_of_model,
          number_of_cpu, mdl_output_pref,
          template_pdb, best_match_struc, superpose_resi,
          pymol_exec, Settings )


##########################################################################
  # Select models with largest binding pocket
  if vol:
    if type(template_list) is not str:
      for idx in range( len(template_list) ):
        GeneratePOVMEAndSortModels(
            script_directory, home_directory, work_directory, result_directory,
            povme_location, povme_directory,
            povme_template, povme_pdb,
            mdl_output_pref+'.'+str(idx+1), number_of_cpu, conformation, 
            int(top_model), Settings )
    else:
      GeneratePOVMEAndSortModels(
          script_directory, home_directory, work_directory, result_directory,
          povme_location, povme_directory, povme_template, povme_pdb,
          mdl_output_pref, number_of_cpu, conformation, 
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
#   v6.0    16.11.27	search for closest structure if no supplied
#   v7.0    16.12.05    integrate CheckHomologTemplate to SearchKinaseStruct
#   v8.0    18.03.06    allow modeling for multiple types of conformations,
#                       and multiple subconformations for CODI and wCD
#   v9.0    18.03.28    change the 'setting_file' to accept pre-defined Dict
#
#
