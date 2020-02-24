#!/usr/bin/env python3

import re,os,sys
import tarfile
from x_pir_multi_align import ModellerMultiAlignGen
from x_pymol_alignment import PyMOLSuperpose

#####################################################################
  # Generate modified alignment file (.pir) for Modeller
def GenerateAndModifyPIR(
        script_directory, dataset_dir, home_directory, 
        work_directory, result_directory,
        fasta_database, kinase_profile, pdb_directory, template_list, 
	template_pdb, superpose_resi, best_match_struc,
        prot_struc_input, mdl_prot_fasta, 
        pc_ident, align_switch, correct_fasta,
        chimera_tmpl_list, mdl_pir_file, mdl_output_pref, pymol_exec,
        Settings ):

  Vars = ['FastaDatabase', 'KinaseProfile', 'PDBDirectory', 'TemplateList',
          'OutputPrefix', 'KinaseStructInput', 'ModelKinaseFasta',
          'ChimeraTemplList', 'ScriptDirectory', 'DatasetDirectory',
          'PyMOL']

  print('\n  -- Generating Alignment File --')
  for var in Vars:
    if var == 'ModelKinaseFasta':   # ModelKinaseFasta can be from fasta DB
      if re.search(r'.fasta', Settings[var]):
        if not os.path.exists(Settings[var]):
          sys.exit('\n  ERROR: {0} does not exist.'.format(Settings[var]))

    if type(Settings[var]) is None:
      sys.exit('\n  ERROR: \'{0}\' is not specified: {0}'.format(var))

#################################################

  os.chdir(work_directory)

  with open('_TEMP.pdb.list', 'w') as fo:
    fo.write(prot_struc_input)

  # align input structure to reference and redirect it
  # if input kinase structure is 'None', will use best matching PDB as input
  if not re.search(r'None', prot_struc_input, re.IGNORECASE):
    tget_pdb = SuperposeStruct( dataset_dir, template_pdb, best_match_struc,
                                superpose_resi, prot_struc_input,
                                '_TEMP.pdb.list', '_TEMP.struct_super',
                                pymol_exec )
  else:
    tget_pdb = prot_struc_input

  print('  ** Running multiple template-target structure superposition **')

  ModellerMultiAlignGen( fasta_database, kinase_profile, pdb_directory, 
                         work_directory, template_list, tget_pdb, 
                         mdl_prot_fasta, best_match_struc, 
                         pc_ident, align_switch, correct_fasta,
                         chimera_tmpl_list, mdl_pir_file, mdl_output_pref )

  print('  ** Finished template-target structure superposition **')

  # replace existing tar.bz2
  if os.path.exists('{0}/{1}.chim_pdb.tar.bz2'.format(result_directory, mdl_output_pref)):
    os.remove('{0}/{1}.chim_pdb.tar.bz2'.format(result_directory, mdl_output_pref))
  tar = tarfile.open('{0}/{1}.chim_pdb.tar.bz2'.format(result_directory, 
        mdl_output_pref), mode='w:bz2')

  if os.path.isfile(chimera_tmpl_list):
    with open(chimera_tmpl_list, 'r') as fi:
      for chimera in fi:
        tar.add(chimera.strip()+'.pdb')
    tar.add(chimera_tmpl_list)
    tar.close()
  else:
    print('### No such file: '+chimera_tmpl_list)

  os.system('cp {0} {1} {2}'.format(mdl_pir_file, chimera_tmpl_list,
                       result_directory))
  os.chdir(home_directory)


##########################################################################
## Superpose the input kinase PDB to reference 1atp to ensure consistency
## and redirect "tget_pdb" kinase structure to the aligned version in the 
## dfgmod working directory
def SuperposeStruct( dataset_dir, template_pdb, best_match_struc,
        superpose_resi, prot_struc_input, org_pdb_list, 
        pymol_align_pref, pymol_exec ):

  print('\n  ** Initiating Structure Superposition **')

  PyMOLSuperpose( pymol_exec, dataset_dir+'/'+template_pdb, 
        best_match_struc, superpose_resi, org_pdb_list, 
        pymol_align_pref, '_TEMP.list', 1, 'mod.pdb' )

  # rename file; x_pymol_alignment.py will name the aligned structure file
  # to xxx.1atp.mod.pdb. Rename to xxx.1atp.pdb
  name = prot_struc_input.split('/')[-1].split('.')[0]
  os.system('mv {0}.*mod.pdb {0}.{1}'.format(name, template_pdb))

  return str(os.getcwd())+'/'+name+'.'+template_pdb


##########################################################################
##
##	Peter M.U. Ung @ MSSM
##
##	v1.0 -- 14.06.10
##	v2.0 -- 16.11.26 - add function to align all input structure before
##			   any homology modeling to ensure superposition
##  v3.0 -- 17.07.12 - x_pymol_alignment is now used as Object
##  v4.0 -- 18.03.09 - adopt to handle CODI multiple runs
##  v5.0 -- 18.03.30 - change the handling of 'None' structure
##  v6.0 -- 19.12.28 - PyMOLSuperpose with new output PDB extension
