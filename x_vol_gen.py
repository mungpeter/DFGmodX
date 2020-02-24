import sys,os
import re,glob
import tarfile
import pandas as pd

from povme import povme
from x_variables import POVMESetup
from x_build_multi_pdb import BuildMultiPDB

##########################################################################
# 
def GeneratePOVMEAndSortModels(
        script_directory, home_directory, work_directory, result_directory,
        povme_location, povme_directory, povme_template, povme_pdb,
        mdl_output_pref, number_of_cpu, conformation, top_model, Settings ):

    Vars = ['POVMELocation', 'POVMEDirectory', 'OutputPrefix', 
            'POVMEStructure', 'NumberOfCPU', 'NumberOfTopModel']

    print('\n  -- Selection of Models --')
    for var in Vars:
      if Settings[var] is None:
        sys.exit('\n  Error: \'{0}\' is not specified: {0}'.format(var))

#################################################

    print(os.getcwd(), '\n')
    os.chdir(work_directory)

    #### Hard-coded variable
    aligned_mdl_list = mdl_output_pref+'.align_pdb.list'


    if not os.path.isdir(povme_directory):
      os.makedirs(povme_directory)

    ## Generate a multi-model PDB file of all model 
    mdl_list_name  = '{0}*.mod.pdb'.format(mdl_output_pref)
    povme_pdb_name = '{0}.multi_mod.pdb'.format(mdl_output_pref)
    GenerateMultiStruct( script_directory, work_directory, 
        mdl_list_name, povme_pdb_name, aligned_mdl_list )

    os.system('cp {0} "{1}"'.format(povme_pdb_name, povme_directory))
    os.chdir(povme_directory)

    ## Calculate the volume of the binding site
    povme_sorted_file = '{0}.volumes.sorted.txt'.format(mdl_output_pref)
    RunPOVMEVolumes( 
        script_directory, work_directory, result_directory,
        povme_location, povme_directory, povme_template, 
        povme_pdb_name, povme_sorted_file, mdl_output_pref, 
        number_of_cpu, conformation )

    ## Select top models based on binding site volume
    os.chdir(work_directory)
    Model_Volumes = pd.read_csv(povme_sorted_file, sep='\s+', comment='#', header=None).to_numpy()
    Model_List    = pd.read_csv(aligned_mdl_list, sep='\s+', comment='#', header=None).to_numpy()

    Top_Mdl_Names = []
    for idx in range(0, top_model):
      try:
        Top_Mdl_Names.append(Model_List[ int(Model_Volumes[idx])-1 ])
      except IndexError:
        break

    # Get the models with top 10 volume
    Top_mdl_Numbers = []
    with open('{0}.top{1}.list'.format(mdl_output_pref, top_model), 'w') as fo:
      for name in Top_Mdl_Names:
        fo.write(name+'\n')
        # Extract the Model number
        if re.search(r'.B0', name):
          tmp_name = name.split('{0}.B0'.format(mdl_output_pref))[1]
        else:
          tmp_name = name.split('{0}.B'.format(mdl_output_pref))[1]
        Top_mdl_Numbers.append(tmp_name.split('.')[0])
      print(' -- Wrote {0}.top{1}.list --'.format(mdl_output_pref, top_model))

    print('\n  *** Top number ***')
    print(Top_mdl_Numbers)

    ## Build multi-model PDB for binding site volumes (all and top)
    os.chdir(povme_directory)
    Top_Frms = [glob.glob('*frame_{0}.*'.format(x))[0] for x in Top_mdl_Numbers]
    povme_pdb_name = '{0}.multi_vol.pdb'.format(mdl_output_pref)

    GenerateMultiStruct(  script_directory, work_directory,
                          Top_Frms, povme_pdb_name, 'top.frame.list')
    if os.path.isfile('{0}/{1}.top_vol.tar.bz2'.format(result_directory, mdl_output_pref)):
      os.remove('{0}/{1}.top_vol.tar.bz2'.format(result_directory, mdl_output_pref))
    tar = tarfile.open('{0}/{1}.top_vol.tar.bz2'.format(result_directory, 
        mdl_output_pref), mode='w:bz2')
    for frm in Top_Frms:
      tar.add(frm)
    tar.close()
    os.system('cp {0} "{1}"'.format('top.frame.list', result_directory))


    ## Save the Top models to Result Directory
    os.chdir(work_directory)
    if os.path.isfile('{0}/{1}.top_pdb.tar.bz2'.format(result_directory, mdl_output_pref)):
      os.remove('{0}/{1}.top_pdb.tar.bz2'.format(result_directory, mdl_output_pref))
    tar = tarfile.open('{0}/{1}.top_pdb.tar.bz2'.format(result_directory,
          mdl_output_pref), mode='w:bz2')
    for mdl in Top_Mdl_Names:
      tar.add(mdl)
    tar.add('{0}.top{1}.list'.format(mdl_output_pref, top_model))
    tar.close()


##########################################################################
## Generate multi-PDB structure for POVME
def GenerateMultiStruct(
        script_directory, work_directory, 
        mdl_list_input, povme_pdb_name, aligned_mdl_list ):

  print('\n - Generate multi-frame model kinase PDB for POVME')
  if type(mdl_list_input) is list:    # option to input as list or *
    Model_List = mdl_list_input
  else:
    try:
      Model_List = sorted(glob.glob(mdl_list_input))
    except IndexError:
      sys.exit('\n  > #2# ERROR: Cannot find structure to aggregate for POVME: {0}'.format(mdl_list_input))
    print(' -- Models to be aggregated for POVME --\n{0}'.format(Model_List))

  ## Write a list of models to be aggregated.
  with open(aligned_mdl_list, 'w') as fo:
    for name in Model_List: 
      fo.write(name.strip()+'\n')
  ## Write a congregated PDB file
  BuildMultiPDB( aligned_mdl_list, povme_pdb_name )


##########################################################################
## Calculate binding site volume for top model selection
def RunPOVMEVolumes(
        script_directory, work_directory, result_directory,
        povme_location, povme_directory, povme_template, 
        povme_pdb_name, povme_sorted_file, mdl_output_pref, 
        number_of_cpu, conformation ):

  Vars = [povme_pdb_name]
  for var in Vars:
    if not os.path.isfile(var):
      sys.exit('\n  > #2# ERROR: Cannot find {0}: {1}'.format(var, mdl_output_pref))

  povme_setup_file = '{0}.povme_setup.in'.format(mdl_output_pref)
  GeneratePOVMEInputFile(
        povme_template, povme_directory, povme_pdb_name, script_directory,
        povme_setup_file, mdl_output_pref, number_of_cpu, conformation )

  print('\n ** Using Python3 and POVME 2.1 **\n')
#  os.system('python {0} {1}'.format(povme_location, povme_setup_file))
  povme.RunPOVME(povme_setup_file)

  # Sort the models by order of volume (descending order)
  Vols = pd.read_csv(mdl_output_pref+'.volumes.tabbed.txt', sep='\s+',
                      comment='#', header=None).sort_values(1, ascending=False)
  Vols.to_csv(povme_sorted_file, header=False, index=False, sep='\t')


  # Copy POVME calculation results to result directory
  os.system('cp {0} "{1}"'.format(povme_sorted_file, result_directory))
  os.system('cp {0} "{1}"'.format(povme_sorted_file, work_directory))


##########################################################################
## Setup and modify POVME input file
def GeneratePOVMEInputFile(
        povme_template, povme_directory, povme_pdb_name, script_directory, 
        povme_setup_file, mdl_output_pref, number_of_cpu, conformation ):

  lines  = '## POVME 2.1 input file for conformation: {0} ##\n\n'.format(conformation)
  lines += 'OutputFilenamePrefix\t{0}/{1}.\n'.format(povme_directory, mdl_output_pref)
  lines += 'PDBFileName\t\t{0}/{1}\n'.format(povme_directory, povme_pdb_name)
  lines += 'NumProcessors\t\t{0}\n\n'.format(number_of_cpu)
  lines += POVMESetup(conformation)
  
  with open(povme_setup_file, 'w') as fo:
    fo.write(lines)


##########################################################################
#
#   Peter M.U. Ung @ MSSM
#
#   v2.0    16.07.15 - select top volume PDB and copy to result directory
#   v3.0    17.07.14 - use x_build_multi_pdb as Object
#   v4.0    18.03.09 - update
#   v5.0    18.04.08 - ditch template POVME input file, "use x_variables"
#   v6.0    20.02.22 - update to POVME2.1 on python3, call it as function
#
